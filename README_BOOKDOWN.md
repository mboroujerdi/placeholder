# Network Extension Project - Bookdown Setup

## Files Provided

### 1. **intro.md** 
Introduction chapter covering:
- Motivation and background
- Control theory in biological systems
- Second-order negative feedback framework
- Network extension to endocrine systems
- Scope and objectives
- Relationship to previous work
- Significance

### 2. **methods.md**
Detailed methodology including:
- Mathematical framework (differential equations, state-space, transfer functions)
- G-type and H-type formulations
- Network composition rules
- Parameter assignment from MCR
- **Time-domain analytical solutions** (inverse Laplace transform for step inputs)
- Frequency-domain analysis (Bode plots)
- Stability analysis (Routh-Hurwitz, eigenvalues, Nyquist)
- **R implementation algorithms** (complete pseudocode)
- Validation approach
- Data sources
- Reporting standards

### 3. **substrates_hormones.csv**
Comprehensive parameter database with 35 hormones and substrates including:
- **name**: Hormone/substrate identifier
- **category**: substrate or hormone
- **MCR_min**: Metabolic clearance rate (min⁻¹)
- **MCR_h**: Metabolic clearance rate (h⁻¹)
- **omega_n_rad_per_min**: Natural frequency (rad/min) for use in simulations
- **zeta**: Suggested damping ratio
- **suggested_type**: G or H architecture recommendation
- **halflife_min**: Biological half-life (minutes)
- **reference_key**: BibTeX key for citation
- **notes**: Additional context

**Major substrates:**
- Glucose, FFA, Amino Acids, Lactate, Ketone Bodies

**Key hormones:**
- Insulin, Glucagon, Growth Hormone, Cortisol, ACTH
- T4, T3, TSH (thyroid axis)
- Testosterone, Estradiol, Progesterone, LH, FSH (reproductive)
- And 20+ more

### 4. **references.bib**
Complete BibTeX bibliography with ~40 references covering:
- Your foundational work (Boroujerdi 2013, 2024)
- Control theory in biology (Del Vecchio, Khammash)
- Pharmacokinetic studies for each hormone/substrate
- R packages (deSolve, ggplot2, bookdown)
- Endocrine physiology

## Setting Up Your Bookdown Project

### Basic Setup

1. **Create project structure:**
```bash
project/
├── index.Rmd          # Title page and setup
├── 01-intro.Rmd       # Copy from intro.md
├── 02-methods.Rmd     # Copy from methods.md
├── 03-results.Rmd     # Your results
├── 04-discussion.Rmd  # Your discussion
├── _bookdown.yml      # Bookdown configuration
├── _output.yml        # Output formats
├── book.bib           # Copy from references.bib
├── data/
│   └── substrates_hormones.csv
└── R/
    └── functions.R    # Your R functions
```

2. **Rename files:**
```bash
# In your RStudio project
cp intro.md 01-intro.Rmd
cp methods.md 02-methods.Rmd
cp references.bib book.bib
cp substrates_hormones.csv data/
```

3. **Edit `_bookdown.yml`:**
```yaml
book_filename: "network_extension"
delete_merged_file: true
language:
  ui:
    chapter_name: "Chapter "
```

4. **Edit `_output.yml`:**
```yaml
bookdown::gitbook:
  css: style.css
  config:
    toc:
      before: |
        <li><a href="./">Network Extension Framework</a></li>
      after: |
        <li><a href="https://github.com/rstudio/bookdown" target="blank">Published with bookdown</a></li>
    download: ["pdf", "epub"]
bookdown::pdf_book:
  includes:
    in_header: preamble.tex
  latex_engine: xelatex
  citation_package: natbib
  keep_tex: yes
bookdown::epub_book: default
```

5. **Create `index.Rmd`:**
```yaml
--- 
title: "Control-Theoretic Framework for Endocrine Network Dynamics"
author: "Massoud Boroujerdi"
date: "`r Sys.Date()`"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
description: "Network extension of second-order negative feedback model with applications to endocrine and metabolic systems"
---

# Preface {-}

This work extends our previous frequency-domain analysis framework to networks of hormones and metabolic substrates. Each element is characterized by its metabolic clearance rate (MCR), and network dynamics emerge from the interplay of individual element properties and network topology.

```

## Using the CSV Data in R

### Load the data:
```r
# Load substrate and hormone parameters
params <- read.csv("data/substrates_hormones.csv", stringsAsFactors = FALSE)

# View structure
str(params)

# Filter by category
substrates <- params[params$category == "substrate", ]
hormones <- params[params$category == "hormone", ]

# Get specific hormone
insulin <- params[params$name == "Insulin", ]
omega_n_insulin <- insulin$omega_n_rad_per_min
zeta_insulin <- insulin$zeta
```

### Create network from data:
```r
# Example: Glucose-Insulin-Glucagon network
network_elements <- c("Glucose", "Insulin", "Glucagon")
network_params <- params[params$name %in% network_elements, ]

# Build parameter vectors
omega_n <- network_params$omega_n_rad_per_min
zeta <- network_params$zeta
architecture <- network_params$suggested_type

# Compute rate constants
k2 <- 2 * zeta * omega_n
k3 <- omega_n
k4 <- omega_n
```

## Time-Domain Simulations

The methods section includes analytical solutions for step inputs. To simulate:

```r
# Analytical solution for underdamped second-order system
step_response <- function(t, omega_n, zeta, amplitude = 1) {
  if (zeta < 1) {
    # Underdamped
    omega_d <- omega_n * sqrt(1 - zeta^2)
    phi <- atan(sqrt(1 - zeta^2) / zeta)
    y <- amplitude * (1 - exp(-zeta * omega_n * t) / sqrt(1 - zeta^2) * 
                        sin(omega_d * t + phi))
  } else if (zeta == 1) {
    # Critically damped
    y <- amplitude * (1 - exp(-omega_n * t) * (1 + omega_n * t))
  } else {
    # Overdamped
    s1 <- -omega_n * (zeta + sqrt(zeta^2 - 1))
    s2 <- -omega_n * (zeta - sqrt(zeta^2 - 1))
    y <- amplitude * (1 - (zeta + sqrt(zeta^2 - 1)) * exp(s1 * t) / (2 * sqrt(zeta^2 - 1)) +
                        (zeta - sqrt(zeta^2 - 1)) * exp(s2 * t) / (2 * sqrt(zeta^2 - 1)))
  }
  return(y)
}

# Example: Insulin response to glucose step
library(ggplot2)

time <- seq(0, 60, by = 0.1)  # 60 minutes
insulin_params <- params[params$name == "Insulin", ]
omega_n <- insulin_params$omega_n_rad_per_min
zeta <- insulin_params$zeta

response <- step_response(time, omega_n, zeta)

ggplot(data.frame(time = time, response = response), 
       aes(x = time, y = response)) +
  geom_line(linewidth = 1, color = "steelblue") +
  labs(title = "Insulin Response to Step Input",
       subtitle = paste("ω_n =", round(omega_n, 4), "rad/min, ζ =", zeta),
       x = "Time (minutes)", y = "Normalized Response") +
  theme_minimal()
```

## Building Bookdown

### In RStudio:

1. Open your project
2. Install bookdown: `install.packages("bookdown")`
3. Build menu → Build Book → Or use:

```r
bookdown::render_book("index.Rmd", "bookdown::gitbook")
bookdown::render_book("index.Rmd", "bookdown::pdf_book")
```

### From command line:

```bash
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
```

## Adding Your Results Chapters

Create new Rmd files following this structure:

**03-glucose-insulin.Rmd:**
```markdown
# Glucose-Insulin Network {#glucose-insulin}

## Network Architecture

The glucose-insulin regulatory system consists of...

## Parameter Values

```{r glucose-insulin-params}
params <- read.csv("data/substrates_hormones.csv")
gi_params <- params[params$name %in% c("Glucose", "Insulin"), ]
knitr::kable(gi_params, caption = "Glucose-Insulin Parameters")
```

## Simulation Results

```{r glucose-insulin-sim}
# Your simulation code
```
```

## Citation Examples

In your Rmd files:

```markdown
The second-order negative feedback framework [@Boroujerdi2013; @Boroujerdi2024]
predicts that glucose clearance rate [@DeFronzo1979] determines...

Insulin has a half-life of approximately 5 minutes [@Polonsky1986], 
corresponding to MCR = 0.1386 min⁻¹.
```

## Tips for Collaborators

1. **Consistent naming:** Use the exact names from `substrates_hormones.csv`
2. **Parameter lookup:** Always reference the CSV rather than hardcoding values
3. **Add new hormones:** Follow the CSV format and add references to `book.bib`
4. **Cite properly:** Use `[@key]` for citations
5. **Cross-reference:** Use `\@ref(methods)` to reference other chapters

## Next Steps

1. ✅ Copy files to your bookdown project
2. ✅ Set up _bookdown.yml and _output.yml
3. ✅ Create index.Rmd with your title/author
4. ✅ Rename intro.md and methods.md to .Rmd
5. ⬜ Create results chapters (03-*.Rmd)
6. ⬜ Add your simulation code in R/ directory
7. ⬜ Build and review
8. ⬜ Share with collaborators

## Contact

For questions about the mathematical framework or R implementation, refer to:
- Methods chapter for detailed algorithms
- CSV file for all parameter values and references
- BibTeX file for full citations

Good luck with your project! The foundation is all here - now you can focus on the exciting science.
