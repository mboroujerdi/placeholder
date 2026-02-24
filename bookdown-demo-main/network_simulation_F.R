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
# OUTPUT FOLDERS: Create directories for saving plots and metrics
###############################################################################

# Define folder paths
plots_dir   <- "output_plots"
metrics_dir <- "output_metrics"

# Create folders if they do not already exist
if (!dir.exists(plots_dir))   dir.create(plots_dir)
if (!dir.exists(metrics_dir)) dir.create(metrics_dir)

cat("Output folders ready:\n")
cat("  Plots   ->", plots_dir, "\n")
cat("  Metrics ->", metrics_dir, "\n\n")

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
  k2 <- 2 * zeta * omega_n
  k3 <- omega_n
  k4 <- omega_n

  if (type == "G") {
    num <- c(k3, k3 * k2)
    den <- c(1, k2, k3 * k4)
  } else if (type == "H") {
    num <- c(k2 * k4)
    den <- c(1, k2, k3 * k4)
  } else {
    stop("Type must be 'G' or 'H'")
  }

  return(list(
    num = num, den = den,
    omega_n = omega_n, zeta = zeta, type = type,
    k2 = k2, k3 = k3, k4 = k4
  ))
}

#' Polynomial evaluation at a point
polyval <- function(p, x) {
  n <- length(p)
  sum(p * x^((n-1):0))
}

#' Polynomial multiplication (convolution)
polymult <- function(a, b) {
  convolve(a, rev(b), type = "open")
}

#' Create cascade transfer function from multiple elements
cascade_tf <- function(elements) {
  num <- elements[[1]]$num
  den <- elements[[1]]$den

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
step_response_analytical <- function(t, omega_n, zeta, amplitude = 1) {
  y <- numeric(length(t))

  if (zeta < 1) {
    omega_d <- omega_n * sqrt(1 - zeta^2)
    phi     <- atan(sqrt(1 - zeta^2) / zeta)
    y <- amplitude * (1 - exp(-zeta * omega_n * t) / sqrt(1 - zeta^2) *
                        sin(omega_d * t + phi))
  } else if (abs(zeta - 1) < 1e-10) {
    y <- amplitude * (1 - exp(-omega_n * t) * (1 + omega_n * t))
  } else {
    sqrt_term <- sqrt(zeta^2 - 1)
    s1 <- -omega_n * (zeta + sqrt_term)
    s2 <- -omega_n * (zeta - sqrt_term)
    y  <- amplitude * (1 - ((zeta + sqrt_term) * exp(s1 * t) -
                              (zeta - sqrt_term) * exp(s2 * t)) / (2 * sqrt_term))
  }

  return(y)
}

#' Compute performance metrics from step response
compute_metrics <- function(t, y, final_value = 1) {
  idx_10   <- which(y >= 0.1 * final_value)[1]
  idx_90   <- which(y >= 0.9 * final_value)[1]
  rise_time <- t[idx_90] - t[idx_10]

  settled      <- abs(y - final_value) <= 0.02 * final_value
  idx_settle   <- max(which(!settled))
  settling_time <- ifelse(idx_settle == -Inf, t[1], t[idx_settle])

  peak_value    <- max(y)
  overshoot_pct <- (peak_value - final_value) / final_value * 100
  peak_time     <- t[which.max(y)]

  return(list(
    rise_time     = rise_time,
    settling_time = settling_time,
    overshoot_pct = overshoot_pct,
    peak_value    = peak_value,
    peak_time     = peak_time
  ))
}

###############################################################################
# PART 3: FREQUENCY-DOMAIN ANALYSIS (BODE PLOTS)
###############################################################################

#' Compute frequency response (Bode plot data)
frequency_response <- function(tf, freq_range = 10^seq(-3, 3, length.out = 500)) {
  s <- 1i * freq_range

  H <- sapply(s, function(si) {
    polyval(tf$num, si) / polyval(tf$den, si)
  })

  magnitude_dB <- 20 * log10(abs(H))
  phase_deg    <- unwrap_phase(atan2(Im(H), Re(H)) * 180 / pi)

  return(data.frame(freq = freq_range, magnitude_dB = magnitude_dB, phase_deg = phase_deg))
}

#' Unwrap phase to handle -180/+180 discontinuities
unwrap_phase <- function(phase) {
  diff_phase <- diff(phase)
  jumps <- which(abs(diff_phase) > 180)
  for (idx in jumps) {
    if (diff_phase[idx] > 180) {
      phase[(idx+1):length(phase)] <- phase[(idx+1):length(phase)] - 360
    } else {
      phase[(idx+1):length(phase)] <- phase[(idx+1):length(phase)] + 360
    }
  }
  return(phase)
}

#' Plot Bode diagram and optionally save to file
#'
#' @param bode_data Data frame from frequency_response()
#' @param title Plot title
#' @param save_file Optional filename (without path) to save into plots_dir
#' @return Invisible: arrangeGrob object
plot_bode <- function(bode_data, title = "Bode Plot", save_file = NULL) {

  # Magnitude plot
  p1 <- ggplot(bode_data, aes(x = freq, y = magnitude_dB)) +
    geom_line(color = "steelblue",
              linewidth = 1) +   # <-- CHANGE linewidth HERE to adjust magnitude line thickness
    scale_x_log10() +
    labs(x = "Frequency (rad/min)", y = "Magnitude (dB)",
         title = paste(title, "- Magnitude")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # Phase plot
  p2 <- ggplot(bode_data, aes(x = freq, y = phase_deg)) +
    geom_line(color = "darkorange",
              linewidth = 1) +   # <-- CHANGE linewidth HERE to adjust phase line thickness
    scale_x_log10() +
    labs(x = "Frequency (rad/min)", y = "Phase (degrees)",
         title = paste(title, "- Phase")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  combined <- arrangeGrob(p1, p2, ncol = 1)
  grid.arrange(combined)

  # Save if filename provided
  if (!is.null(save_file)) {
    ggsave(filename = file.path(plots_dir, save_file),
           plot = combined, width = 8, height = 6, dpi = 300)
    cat("Bode plot saved to:", file.path(plots_dir, save_file), "\n")
  }

  invisible(combined)
}

#' Extract key metrics from Bode plot
bode_metrics <- function(bode_data) {
  dc_gain_dB  <- bode_data$magnitude_dB[1]
  target_mag  <- dc_gain_dB - 3
  idx_cutoff  <- which(bode_data$magnitude_dB <= target_mag)[1]
  freq_cutoff <- bode_data$freq[idx_cutoff]
  bandwidth   <- freq_cutoff
  idx_gc      <- which.min(abs(bode_data$magnitude_dB))
  phase_margin <- 180 + bode_data$phase_deg[idx_gc]

  return(list(
    dc_gain_dB   = dc_gain_dB,
    cutoff_freq  = freq_cutoff,
    bandwidth    = bandwidth,
    phase_margin = phase_margin
  ))
}

###############################################################################
# PART 4: ODE-BASED NETWORK SIMULATION
###############################################################################

#' ODE system for network of second-order elements
network_ode <- function(t, x, parms) {
  N    <- parms$N
  dxdt <- numeric(2 * N)

  for (i in 1:N) {
    S_i <- x[2*i - 1]
    C_i <- x[2*i]

    dxdt[2*i - 1] <- parms$Ra[i] - parms$k3[i] * C_i
    dxdt[2*i]     <- parms$k4[i] * S_i - parms$k2[i] * C_i

    if (!is.null(parms$feedback_matrix)) {
      for (j in 1:N) {
        if (i != j && parms$feedback_matrix[i, j] != 0) {
          dxdt[2*i] <- dxdt[2*i] + parms$feedback_matrix[i, j] * x[2*j - 1]
        }
      }
    }
  }

  return(list(dxdt))
}

#' Simulate network response
simulate_network <- function(elements, times, Ra,
                             initial_conditions = NULL,
                             feedback_matrix    = NULL) {
  N <- length(elements)
  if (is.null(initial_conditions)) initial_conditions <- rep(0, 2 * N)

  parms <- list(
    N               = N,
    Ra              = Ra,
    k2              = sapply(elements, function(e) e$k2),
    k3              = sapply(elements, function(e) e$k3),
    k4              = sapply(elements, function(e) e$k4),
    feedback_matrix = feedback_matrix
  )

  out    <- ode(y = initial_conditions, times = times, func = network_ode, parms = parms)
  out_df <- as.data.frame(out)
  colnames(out_df) <- c("time", paste0(rep(c("S", "C"), N), rep(1:N, each = 2)))

  return(out_df)
}

###############################################################################
# HELPER: Save metrics to CSV
###############################################################################

#' Save a named list of metrics to a CSV file in metrics_dir
#'
#' @param metrics Named list of scalar metric values
#' @param filename Filename (without path), e.g. "insulin_metrics.csv"
save_metrics <- function(metrics, filename) {
  df   <- as.data.frame(metrics)
  path <- file.path(metrics_dir, filename)
  write.csv(df, path, row.names = FALSE)
  cat("Metrics saved to:", path, "\n")
}

###############################################################################
# PART 5: EXAMPLE SIMULATIONS
###############################################################################

# --- EXAMPLE 1: Single Element Step Response ---
cat("\n=== EXAMPLE 1: Single Insulin Element ===\n")

omega_n_insulin <- 0.1386
zeta_insulin    <- 0.5

insulin  <- create_element(omega_n = omega_n_insulin, zeta = zeta_insulin, type = "G")
time     <- seq(0, 60, by = 0.1)
response <- step_response_analytical(time, omega_n_insulin, zeta_insulin, amplitude = 1)

df_insulin <- data.frame(time = time, response = response)

p_insulin <- ggplot(df_insulin, aes(x = time, y = response)) +
  geom_line(color = "steelblue",
            linewidth = 1.2) +   # <-- CHANGE linewidth HERE to adjust step response line thickness
  geom_hline(yintercept = 1, linetype = "dashed",
             linewidth  = 0.5,   # <-- CHANGE linewidth HERE to adjust reference line thickness
             color = "red", alpha = 0.5) +
  labs(title    = "Insulin Step Response (Analytical Solution)",
       subtitle = bquote(omega[n] == .(round(omega_n_insulin, 4)) ~ "rad/min, " ~
                           zeta == .(zeta_insulin)),
       x = "Time (minutes)", y = "Normalized Response") +
  theme_minimal()

print(p_insulin)

# Save plot
ggsave(filename = file.path(plots_dir, "example1_insulin_step_response.png"),
       plot = p_insulin, width = 8, height = 5, dpi = 300)
cat("Plot saved: example1_insulin_step_response.png\n")

# Compute and save metrics
metrics_ex1 <- compute_metrics(time, response)
cat("\nPerformance Metrics:\n"); print(metrics_ex1)
save_metrics(metrics_ex1, "example1_insulin_metrics.csv")


# --- EXAMPLE 2: Cascade of Two Elements (Glucose → Insulin) ---
cat("\n\n=== EXAMPLE 2: Glucose-Insulin Cascade ===\n")

omega_n_glucose <- 0.0667
zeta_glucose    <- 0.7
omega_n_insulin <- 0.1386
zeta_insulin    <- 0.5

glucose <- create_element(omega_n_glucose, zeta_glucose, type = "G")
insulin <- create_element(omega_n_insulin, zeta_insulin, type = "G")
cascade <- cascade_tf(list(glucose, insulin))

bode_data <- frequency_response(cascade, freq_range = 10^seq(-3, 2, length.out = 500))

# Plot and save Bode diagram (line thickness controlled inside plot_bode())
plot_bode(bode_data, title = "Glucose-Insulin Cascade",
          save_file = "example2_glucose_insulin_bode.png")

# Extract and save Bode metrics
metrics_bode_ex2 <- bode_metrics(bode_data)
cat("\nBode Plot Metrics:\n"); print(metrics_bode_ex2)
save_metrics(metrics_bode_ex2, "example2_glucose_insulin_bode_metrics.csv")


# --- EXAMPLE 3: ODE Simulation of Glucose-Insulin Network ---
cat("\n\n=== EXAMPLE 3: Glucose-Insulin Network (ODE Simulation) ===\n")

elem_glucose <- create_element(omega_n_glucose, zeta_glucose, type = "G")
elem_insulin <- create_element(omega_n_insulin, zeta_insulin, type = "G")

times   <- seq(0, 120, by = 0.5)
Ra      <- c(1.0, 0.5)
initial <- rep(0, 4)

out <- simulate_network(
  elements           = list(elem_glucose, elem_insulin),
  times              = times,
  Ra                 = Ra,
  initial_conditions = initial
)

out_long <- pivot_longer(out, cols = -time, names_to = "variable", values_to = "value")

p_network <- ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) +   # <-- CHANGE linewidth HERE to adjust all network trajectory line thicknesses
  labs(title = "Glucose-Insulin Network Dynamics",
       x = "Time (minutes)", y = "Concentration (arbitrary units)") +
  theme_minimal() +
  scale_color_manual(
    values = c("S1" = "steelblue", "C1" = "lightblue",
               "S2" = "darkorange", "C2" = "peachpuff"),
    labels = c("S1" = "Glucose (S)", "C1" = "Glucose Controller (C)",
               "S2" = "Insulin (S)",  "C2" = "Insulin Controller (C)")
  )

print(p_network)

# Save plot and simulation data
ggsave(filename = file.path(plots_dir, "example3_glucose_insulin_ode.png"),
       plot = p_network, width = 8, height = 5, dpi = 300)
cat("Plot saved: example3_glucose_insulin_ode.png\n")

write.csv(out, file.path(metrics_dir, "example3_glucose_insulin_ode_data.csv"), row.names = FALSE)
cat("ODE simulation data saved: example3_glucose_insulin_ode_data.csv\n")


# --- EXAMPLE 4: Three-Element Network with Feedback ---
cat("\n\n=== EXAMPLE 4: Glucose-Insulin-Glucagon with Feedback ===\n")

omega_n_glucagon <- 0.0866
zeta_glucagon    <- 0.5

elem_glucose  <- create_element(omega_n_glucose,  zeta_glucose,  type = "G")
elem_insulin  <- create_element(omega_n_insulin,  zeta_insulin,  type = "G")
elem_glucagon <- create_element(omega_n_glucagon, zeta_glucagon, type = "H")

feedback        <- matrix(0, nrow = 3, ncol = 3)
feedback[1, 2]  <- -0.5   # Insulin suppresses glucose
feedback[1, 3]  <-  0.3   # Glucagon stimulates glucose
feedback[2, 1]  <-  0.8   # Glucose stimulates insulin
feedback[3, 1]  <- -0.4   # Glucose suppresses glucagon

times   <- seq(0, 200, by = 0.5)
Ra      <- c(1.0, 0.3, 0.2)
initial <- rep(0, 6)

out_3elem <- simulate_network(
  elements           = list(elem_glucose, elem_insulin, elem_glucagon),
  times              = times,
  Ra                 = Ra,
  initial_conditions = initial,
  feedback_matrix    = feedback
)

out_substrates <- pivot_longer(out_3elem, cols = -time,
                               names_to = "variable", values_to = "value") %>%
  filter(variable %in% c("S1", "S2", "S3"))

p_3network <- ggplot(out_substrates, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1.2) +   # <-- CHANGE linewidth HERE to adjust three-element network line thicknesses
  labs(title = "Glucose-Insulin-Glucagon Network with Feedback",
       x = "Time (minutes)", y = "Concentration (arbitrary units)") +
  theme_minimal() +
  scale_color_manual(
    values = c("S1" = "steelblue", "S2" = "darkorange", "S3" = "darkgreen"),
    labels = c("S1" = "Glucose", "S2" = "Insulin", "S3" = "Glucagon")
  )

print(p_3network)

# Save plot and simulation data
ggsave(filename = file.path(plots_dir, "example4_three_element_feedback.png"),
       plot = p_3network, width = 8, height = 5, dpi = 300)
cat("Plot saved: example4_three_element_feedback.png\n")

write.csv(out_3elem, file.path(metrics_dir, "example4_three_element_ode_data.csv"), row.names = FALSE)
cat("ODE simulation data saved: example4_three_element_ode_data.csv\n")


###############################################################################
# PART 6: COMPARISON OF G-TYPE VS H-TYPE
###############################################################################

cat("\n\n=== EXAMPLE 5: G-type vs H-type Comparison ===\n")

omega_n <- 0.1
zeta    <- 0.7

elem_G <- create_element(omega_n, zeta, type = "G")
elem_H <- create_element(omega_n, zeta, type = "H")

bode_G       <- frequency_response(list(num = elem_G$num, den = elem_G$den))
bode_H       <- frequency_response(list(num = elem_H$num, den = elem_H$den))
bode_G$type  <- "G-type"
bode_H$type  <- "H-type"
bode_comparison <- rbind(bode_G, bode_H)

# Magnitude comparison
p_comp_mag <- ggplot(bode_comparison, aes(x = freq, y = magnitude_dB, color = type)) +
  geom_line(linewidth = 1.2) +   # <-- CHANGE linewidth HERE to adjust G vs H magnitude comparison line thicknesses
  scale_x_log10() +
  labs(title = "G-type vs H-type: Magnitude Response",
       x = "Frequency (rad/min)", y = "Magnitude (dB)") +
  theme_minimal() +
  scale_color_manual(values = c("G-type" = "steelblue", "H-type" = "darkorange"))

print(p_comp_mag)

# Phase comparison
p_comp_phase <- ggplot(bode_comparison, aes(x = freq, y = phase_deg, color = type)) +
  geom_line(linewidth = 1.2) +   # <-- CHANGE linewidth HERE to adjust G vs H phase comparison line thicknesses
  scale_x_log10() +
  labs(title = "G-type vs H-type: Phase Response",
       x = "Frequency (rad/min)", y = "Phase (degrees)") +
  theme_minimal() +
  scale_color_manual(values = c("G-type" = "steelblue", "H-type" = "darkorange"))

print(p_comp_phase)

# Save both comparison plots
ggsave(filename = file.path(plots_dir, "example5_GH_magnitude_comparison.png"),
       plot = p_comp_mag, width = 8, height = 5, dpi = 300)
ggsave(filename = file.path(plots_dir, "example5_GH_phase_comparison.png"),
       plot = p_comp_phase, width = 8, height = 5, dpi = 300)
cat("Comparison plots saved: example5_GH_magnitude_comparison.png and example5_GH_phase_comparison.png\n")

# Save Bode metrics for both types
save_metrics(bode_metrics(bode_G), "example5_Gtype_bode_metrics.csv")
save_metrics(bode_metrics(bode_H), "example5_Htype_bode_metrics.csv")


###############################################################################
# SUMMARY
###############################################################################

cat("\n=== Simulation Complete ===\n")
cat("Plots saved to   : '", plots_dir,   "/' folder\n", sep = "")
cat("Metrics saved to : '", metrics_dir, "/' folder\n", sep = "")
cat("\nKey functions available:\n")
cat("  - create_element(omega_n, zeta, type)\n")
cat("  - cascade_tf(elements)\n")
cat("  - step_response_analytical(t, omega_n, zeta)\n")
cat("  - frequency_response(tf, freq_range)\n")
cat("  - plot_bode(bode_data, title, save_file)\n")
cat("  - simulate_network(elements, times, Ra, ...)\n")
cat("  - compute_metrics(t, y)\n")
cat("  - bode_metrics(bode_data)\n")
cat("  - save_metrics(metrics, filename)\n")
