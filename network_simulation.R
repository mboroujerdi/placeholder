###############################################################################
# Network Extension Simulation Framework
# Second-Order Negative Feedback Model for Endocrine Systems
# Based on Boroujerdi & Schmidt (2013) framework
###############################################################################

# Load required packages
library(deSolve)    # For ODE solving
library(ggplot2)    # For plotting
library(gridExtra)  # For multi-panel plots
library(tidyr)      # For data manipulation
library(dplyr)      # For data manipulation

###############################################################################
# PART 1: CORE FUNCTIONS - Transfer Functions and Building Blocks
###############################################################################

#' Create a second-order negative feedback element
#'
#' @param omega_n Natural frequency (rad/time)
#' @param zeta Damping ratio (dimensionless)
#' @param type "G" for G-type (feed-forward) or "H" for H-type (feedback)
#' @return List containing transfer function numerator and denominator coefficients
create_element <- function(omega_n, zeta, type = "G") {
  # Compute rate constants
  k2 <- 2 * zeta * omega_n
  k3 <- omega_n
  k4 <- omega_n
  
  if (type == "G") {
    # G-type: numerator has (s + k2) term
    num <- c(k3, k3 * k2)  # k3*s + k3*k2
    den <- c(1, k2, k3 * k4)  # s^2 + k2*s + k3*k4
  } else if (type == "H") {
    # H-type: numerator is constant
    num <- c(k2 * k4)  # k2*k4
    den <- c(1, k2, k3 * k4)  # s^2 + k2*s + k3*k4
  } else {
    stop("Type must be 'G' or 'H'")
  }
  
  return(list(
    num = num,
    den = den,
    omega_n = omega_n,
    zeta = zeta,
    type = type,
    k2 = k2,
    k3 = k3,
    k4 = k4
  ))
}

#' Polynomial evaluation at a point
#'
#' @param p Polynomial coefficients (highest degree first)
#' @param x Point at which to evaluate
#' @return Polynomial value at x
polyval <- function(p, x) {
  n <- length(p)
  sum(p * x^((n-1):0))
}

#' Polynomial multiplication (convolution)
#'
#' @param a First polynomial coefficients
#' @param b Second polynomial coefficients
#' @return Product polynomial coefficients
polymult <- function(a, b) {
  convolve(a, rev(b), type = "open")
}

#' Create cascade transfer function from multiple elements
#'
#' @param elements List of elements created by create_element()
#' @return List with numerator and denominator of cascade
cascade_tf <- function(elements) {
  # Start with first element
  num <- elements[[1]]$num
  den <- elements[[1]]$den
  
  # Multiply subsequent elements
  if (length(elements) > 1) {
    for (i in 2:length(elements)) {
      num <- polymult(num, elements[[i]]$num)
      den <- polymult(den, elements[[i]]$den)
    }
  }
  
  return(list(num = num, den = den, elements = elements))
}

###############################################################################
# PART 2: TIME-DOMAIN ANALYTICAL SOLUTIONS
###############################################################################

#' Analytical step response for second-order system
#'
#' @param t Time vector
#' @param omega_n Natural frequency
#' @param zeta Damping ratio
#' @param amplitude Step amplitude (default 1)
#' @return Response vector
step_response_analytical <- function(t, omega_n, zeta, amplitude = 1) {
  y <- numeric(length(t))
  
  if (zeta < 1) {
    # Underdamped case (oscillatory)
    omega_d <- omega_n * sqrt(1 - zeta^2)  # Damped natural frequency
    phi <- atan(sqrt(1 - zeta^2) / zeta)   # Phase angle
    
    y <- amplitude * (1 - exp(-zeta * omega_n * t) / sqrt(1 - zeta^2) * 
                        sin(omega_d * t + phi))
    
  } else if (abs(zeta - 1) < 1e-10) {
    # Critically damped case
    y <- amplitude * (1 - exp(-omega_n * t) * (1 + omega_n * t))
    
  } else {
    # Overdamped case
    sqrt_term <- sqrt(zeta^2 - 1)
    s1 <- -omega_n * (zeta + sqrt_term)
    s2 <- -omega_n * (zeta - sqrt_term)
    
    y <- amplitude * (1 - ((zeta + sqrt_term) * exp(s1 * t) - 
                            (zeta - sqrt_term) * exp(s2 * t)) / 
                        (2 * sqrt_term))
  }
  
  return(y)
}

#' Compute performance metrics from step response
#'
#' @param t Time vector
#' @param y Response vector
#' @param final_value Expected final value (default 1)
#' @return List of performance metrics
compute_metrics <- function(t, y, final_value = 1) {
  # Rise time (10% to 90%)
  idx_10 <- which(y >= 0.1 * final_value)[1]
  idx_90 <- which(y >= 0.9 * final_value)[1]
  rise_time <- t[idx_90] - t[idx_10]
  
  # Settling time (within 2% of final value)
  settled <- abs(y - final_value) <= 0.02 * final_value
  idx_settle <- max(which(!settled))
  settling_time <- ifelse(idx_settle == -Inf, t[1], t[idx_settle])
  
  # Peak overshoot
  peak_value <- max(y)
  overshoot_pct <- (peak_value - final_value) / final_value * 100
  
  # Peak time
  peak_time <- t[which.max(y)]
  
  return(list(
    rise_time = rise_time,
    settling_time = settling_time,
    overshoot_pct = overshoot_pct,
    peak_value = peak_value,
    peak_time = peak_time
  ))
}

###############################################################################
# PART 3: FREQUENCY-DOMAIN ANALYSIS (BODE PLOTS)
###############################################################################

#' Compute frequency response (Bode plot data)
#'
#' @param tf Transfer function (list with num and den)
#' @param freq_range Frequency range (rad/time)
#' @return Data frame with frequency, magnitude (dB), and phase (degrees)
frequency_response <- function(tf, freq_range = 10^seq(-3, 3, length.out = 500)) {
  s <- 1i * freq_range  # s = j*omega
  
  # Evaluate transfer function at each frequency
  H <- sapply(s, function(si) {
    num_val <- polyval(tf$num, si)
    den_val <- polyval(tf$den, si)
    return(num_val / den_val)
  })
  
  # Compute magnitude (dB) and phase (degrees)
  magnitude_dB <- 20 * log10(abs(H))
  phase_deg <- atan2(Im(H), Re(H)) * 180 / pi
  
  # Unwrap phase (handle discontinuities)
  phase_deg <- unwrap_phase(phase_deg)
  
  return(data.frame(
    freq = freq_range,
    magnitude_dB = magnitude_dB,
    phase_deg = phase_deg
  ))
}

#' Unwrap phase to handle -180/+180 discontinuities
unwrap_phase <- function(phase) {
  diff_phase <- diff(phase)
  # Find jumps greater than 180 degrees
  jumps <- which(abs(diff_phase) > 180)
  
  for (idx in jumps) {
    if (diff_phase[idx] > 180) {
      phase[(idx+1):length(phase)] <- phase[(idx+1):length(phase)] - 360
    } else if (diff_phase[idx] < -180) {
      phase[(idx+1):length(phase)] <- phase[(idx+1):length(phase)] + 360
    }
  }
  
  return(phase)
}

#' Plot Bode diagram
#'
#' @param bode_data Data frame from frequency_response()
#' @param title Plot title
#' @return ggplot object
plot_bode <- function(bode_data, title = "Bode Plot") {
  # Magnitude plot
  p1 <- ggplot(bode_data, aes(x = freq, y = magnitude_dB)) +
    geom_line(color = "steelblue", linewidth = 1) +
    scale_x_log10() +
    labs(x = "Frequency (rad/min)", y = "Magnitude (dB)",
         title = paste(title, "- Magnitude")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Phase plot
  p2 <- ggplot(bode_data, aes(x = freq, y = phase_deg)) +
    geom_line(color = "darkorange", linewidth = 1) +
    scale_x_log10() +
    labs(x = "Frequency (rad/min)", y = "Phase (degrees)",
         title = paste(title, "- Phase")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Combine plots
  grid.arrange(p1, p2, ncol = 1)
}

#' Extract key metrics from Bode plot
#'
#' @param bode_data Data frame from frequency_response()
#' @return List of metrics
bode_metrics <- function(bode_data) {
  # DC gain
  dc_gain_dB <- bode_data$magnitude_dB[1]
  
  # Cutoff frequency (-3 dB point)
  target_mag <- dc_gain_dB - 3
  idx_cutoff <- which(bode_data$magnitude_dB <= target_mag)[1]
  freq_cutoff <- bode_data$freq[idx_cutoff]
  
  # Bandwidth
  bandwidth <- freq_cutoff
  
  # Phase margin (at gain crossover frequency where magnitude = 0 dB)
  idx_gc <- which.min(abs(bode_data$magnitude_dB))
  phase_margin <- 180 + bode_data$phase_deg[idx_gc]
  
  return(list(
    dc_gain_dB = dc_gain_dB,
    cutoff_freq = freq_cutoff,
    bandwidth = bandwidth,
    phase_margin = phase_margin
  ))
}

###############################################################################
# PART 4: ODE-BASED NETWORK SIMULATION
###############################################################################

#' Define ODE system for network of second-order elements
#'
#' @param t Time
#' @param x State vector [S1, C1, S2, C2, ...]
#' @param parms Parameter list
#' @return Rate of change vector
network_ode <- function(t, x, parms) {
  N <- parms$N  # Number of elements
  
  dxdt <- numeric(2 * N)
  
  for (i in 1:N) {
    # Extract state variables for element i
    S_i <- x[2*i - 1]  # Substrate/hormone
    C_i <- x[2*i]      # Controller
    
    # Substrate compartment equation: dS/dt = Ra - k3*C
    dxdt[2*i - 1] <- parms$Ra[i] - parms$k3[i] * C_i
    
    # Controller compartment equation: dC/dt = k4*S - k2*C
    dxdt[2*i] <- parms$k4[i] * S_i - parms$k2[i] * C_i
    
    # Add feedback from other elements if specified
    if (!is.null(parms$feedback_matrix)) {
      for (j in 1:N) {
        if (i != j && parms$feedback_matrix[i, j] != 0) {
          # Add feedback term: coupling * S_j
          dxdt[2*i] <- dxdt[2*i] + parms$feedback_matrix[i, j] * x[2*j - 1]
        }
      }
    }
  }
  
  return(list(dxdt))
}

#' Simulate network response
#'
#' @param elements List of elements
#' @param times Time vector
#' @param Ra Production rates (vector, one per element)
#' @param initial_conditions Initial state vector (optional)
#' @param feedback_matrix Feedback coupling matrix (optional)
#' @return Data frame with time and state variables
simulate_network <- function(elements, times, Ra, 
                             initial_conditions = NULL,
                             feedback_matrix = NULL) {
  N <- length(elements)
  
  # Default initial conditions (all zeros)
  if (is.null(initial_conditions)) {
    initial_conditions <- rep(0, 2 * N)
  }
  
  # Extract rate constants
  k2 <- sapply(elements, function(e) e$k2)
  k3 <- sapply(elements, function(e) e$k3)
  k4 <- sapply(elements, function(e) e$k4)
  
  # Parameters
  parms <- list(
    N = N,
    Ra = Ra,
    k2 = k2,
    k3 = k3,
    k4 = k4,
    feedback_matrix = feedback_matrix
  )
  
  # Solve ODE
  out <- ode(y = initial_conditions, times = times, func = network_ode, parms = parms)
  
  # Convert to data frame with meaningful names
  out_df <- as.data.frame(out)
  colnames(out_df) <- c("time", paste0(rep(c("S", "C"), N), rep(1:N, each = 2)))
  
  return(out_df)
}

###############################################################################
# PART 5: EXAMPLE SIMULATIONS
###############################################################################

# --- EXAMPLE 1: Single Element Step Response ---
cat("\n=== EXAMPLE 1: Single Insulin Element ===\n")

# Load parameters (assuming CSV is available, otherwise use hardcoded values)
# params <- read.csv("substrates_hormones.csv")
# insulin_params <- params[params$name == "Insulin", ]
# omega_n_insulin <- insulin_params$omega_n_rad_per_min
# zeta_insulin <- insulin_params$zeta

# Hardcoded for demonstration
omega_n_insulin <- 0.1386  # rad/min
zeta_insulin <- 0.5

# Create element
insulin <- create_element(omega_n = omega_n_insulin, zeta = zeta_insulin, type = "G")

# Time vector (60 minutes)
time <- seq(0, 60, by = 0.1)

# Analytical step response
response <- step_response_analytical(time, omega_n_insulin, zeta_insulin, amplitude = 1)

# Plot
df_insulin <- data.frame(time = time, response = response)
p_insulin <- ggplot(df_insulin, aes(x = time, y = response)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  labs(title = "Insulin Step Response (Analytical Solution)",
       subtitle = bquote(omega[n] == .(round(omega_n_insulin, 4)) ~ "rad/min, " ~ 
                          zeta == .(zeta_insulin)),
       x = "Time (minutes)", y = "Normalized Response") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5)

print(p_insulin)

# Compute metrics
metrics <- compute_metrics(time, response)
cat("\nPerformance Metrics:\n")
print(metrics)

# --- EXAMPLE 2: Cascade of Two Elements (Glucose → Insulin) ---
cat("\n\n=== EXAMPLE 2: Glucose-Insulin Cascade ===\n")

# Parameters
omega_n_glucose <- 0.0667  # rad/min
zeta_glucose <- 0.7
omega_n_insulin <- 0.1386
zeta_insulin <- 0.5

# Create elements
glucose <- create_element(omega_n_glucose, zeta_glucose, type = "G")
insulin <- create_element(omega_n_insulin, zeta_insulin, type = "G")

# Create cascade
cascade <- cascade_tf(list(glucose, insulin))

# Frequency response
bode_data <- frequency_response(cascade, freq_range = 10^seq(-3, 2, length.out = 500))

# Plot Bode diagram
plot_bode(bode_data, title = "Glucose-Insulin Cascade")

# Extract metrics
metrics_bode <- bode_metrics(bode_data)
cat("\nBode Plot Metrics:\n")
print(metrics_bode)

# --- EXAMPLE 3: ODE Simulation of Glucose-Insulin Network ---
cat("\n\n=== EXAMPLE 3: Glucose-Insulin Network (ODE Simulation) ===\n")

# Create elements
elem_glucose <- create_element(omega_n_glucose, zeta_glucose, type = "G")
elem_insulin <- create_element(omega_n_insulin, zeta_insulin, type = "G")

# Time vector
times <- seq(0, 120, by = 0.5)  # 120 minutes

# Production rates (constant)
Ra <- c(1.0, 0.5)  # Glucose and insulin production

# Initial conditions (start at zero)
initial <- rep(0, 4)  # [S_glucose, C_glucose, S_insulin, C_insulin]

# Simulate
out <- simulate_network(
  elements = list(elem_glucose, elem_insulin),
  times = times,
  Ra = Ra,
  initial_conditions = initial
)

# Plot results
out_long <- pivot_longer(out, cols = -time, names_to = "variable", values_to = "value")

p_network <- ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) +
  labs(title = "Glucose-Insulin Network Dynamics",
       x = "Time (minutes)", y = "Concentration (arbitrary units)") +
  theme_minimal() +
  scale_color_manual(values = c("S1" = "steelblue", "C1" = "lightblue",
                                  "S2" = "darkorange", "C2" = "peachpuff"),
                     labels = c("S1" = "Glucose (S)", "C1" = "Glucose Controller (C)",
                               "S2" = "Insulin (S)", "C2" = "Insulin Controller (C)"))

print(p_network)

# --- EXAMPLE 4: Three-Element Network with Feedback ---
cat("\n\n=== EXAMPLE 4: Glucose-Insulin-Glucagon with Feedback ===\n")

# Parameters
omega_n_glucagon <- 0.0866  # rad/min
zeta_glucagon <- 0.5

# Create elements
elem_glucose <- create_element(omega_n_glucose, zeta_glucose, type = "G")
elem_insulin <- create_element(omega_n_insulin, zeta_insulin, type = "G")
elem_glucagon <- create_element(omega_n_glucagon, zeta_glucagon, type = "H")

# Feedback matrix (simplified regulatory network)
# Row i, col j: effect of element j on element i
feedback <- matrix(0, nrow = 3, ncol = 3)
feedback[1, 2] <- -0.5  # Insulin suppresses glucose
feedback[1, 3] <- 0.3   # Glucagon stimulates glucose
feedback[2, 1] <- 0.8   # Glucose stimulates insulin
feedback[3, 1] <- -0.4  # Glucose suppresses glucagon

# Time vector
times <- seq(0, 200, by = 0.5)

# Production rates
Ra <- c(1.0, 0.3, 0.2)

# Initial conditions
initial <- rep(0, 6)

# Simulate
out_3elem <- simulate_network(
  elements = list(elem_glucose, elem_insulin, elem_glucagon),
  times = times,
  Ra = Ra,
  initial_conditions = initial,
  feedback_matrix = feedback
)

# Plot
out_3elem_long <- pivot_longer(out_3elem, cols = -time, names_to = "variable", values_to = "value")

# Filter to show only S variables (substrates/hormones)
out_substrates <- out_3elem_long %>% 
  filter(variable %in% c("S1", "S2", "S3"))

p_3network <- ggplot(out_substrates, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Glucose-Insulin-Glucagon Network with Feedback",
       x = "Time (minutes)", y = "Concentration (arbitrary units)") +
  theme_minimal() +
  scale_color_manual(values = c("S1" = "steelblue", "S2" = "darkorange", "S3" = "darkgreen"),
                     labels = c("S1" = "Glucose", "S2" = "Insulin", "S3" = "Glucagon"))

print(p_3network)

###############################################################################
# PART 6: COMPARISON OF G-TYPE VS H-TYPE
###############################################################################

cat("\n\n=== EXAMPLE 5: G-type vs H-type Comparison ===\n")

# Same omega_n and zeta, different architecture
omega_n <- 0.1
zeta <- 0.7

elem_G <- create_element(omega_n, zeta, type = "G")
elem_H <- create_element(omega_n, zeta, type = "H")

# Step response
time <- seq(0, 80, by = 0.1)
response_G <- step_response_analytical(time, omega_n, zeta)
response_H <- step_response_analytical(time, omega_n, zeta)

# Note: For H-type, the DC gain is different, so we need to scale
# Actually, for proper comparison, let's use frequency response

# Frequency response
bode_G <- frequency_response(list(num = elem_G$num, den = elem_G$den))
bode_H <- frequency_response(list(num = elem_H$num, den = elem_H$den))

bode_G$type <- "G-type"
bode_H$type <- "H-type"
bode_comparison <- rbind(bode_G, bode_H)

# Plot magnitude comparison
p_comp_mag <- ggplot(bode_comparison, aes(x = freq, y = magnitude_dB, color = type)) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() +
  labs(title = "G-type vs H-type: Magnitude Response",
       x = "Frequency (rad/min)", y = "Magnitude (dB)") +
  theme_minimal() +
  scale_color_manual(values = c("G-type" = "steelblue", "H-type" = "darkorange"))

print(p_comp_mag)

# Plot phase comparison
p_comp_phase <- ggplot(bode_comparison, aes(x = freq, y = phase_deg, color = type)) +
  geom_line(linewidth = 1.2) +
  scale_x_log10() +
  labs(title = "G-type vs H-type: Phase Response",
       x = "Frequency (rad/min)", y = "Phase (degrees)") +
  theme_minimal() +
  scale_color_manual(values = c("G-type" = "steelblue", "H-type" = "darkorange"))

print(p_comp_phase)

cat("\n=== Simulation Complete ===\n")
cat("All core functions are now loaded and examples have been run.\n")
cat("You can now use these functions to build your own networks!\n\n")

cat("Key functions available:\n")
cat("  - create_element(omega_n, zeta, type)\n")
cat("  - cascade_tf(elements)\n")
cat("  - step_response_analytical(t, omega_n, zeta)\n")
cat("  - frequency_response(tf, freq_range)\n")
cat("  - plot_bode(bode_data, title)\n")
cat("  - simulate_network(elements, times, Ra, ...)\n")
cat("  - compute_metrics(t, y)\n")
cat("  - bode_metrics(bode_data)\n")
