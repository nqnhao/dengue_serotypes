# -----------------------------------------
# Time-Varying Multi-Cutoff FOI Simulation
# -----------------------------------------

library(tidyverse)
library(rstan)

# --- Parameters ---
lambda1_values <- c(0.05, 0.03, 0.01, 0.04, 0.02)  # Multiple FOI values for different phases
lambda2_values <- c(0.04, 0.02, 0.01, 0.03, 0.02)  # Multiple FOI values for different phases
sigma12 <- 15  # Cross-protection between serotypes
sigma21 <- 0.1  # Cross-protection between serotypes

# Define time periods (e.g., every 5 years)
time_periods <- seq(0, 60, by = 5)  # Time periods from 0 to 60 years (13 phases)

# Current year
current_year <- 2025  

# Define multiple cut-off years (e.g., different phases)
cutoff_years <- c(2000, 2010, 2020)  # Example: Each cut-off year represents a new FOI phase

# --- Model functions ---
# Probability of remaining susceptible over time
s <- function(t, lambda1_values, lambda2_values) {
  phase <- findInterval(t, time_periods)  # Determine the phase based on the time
  lambda1 <- lambda1_values[phase]  # Force of infection for current phase (lambda1)
  lambda2 <- lambda2_values[phase]  # Force of infection for current phase (lambda2)
  exp(-t * (lambda1 + lambda2))
}

# Probability of being infected with serotype 1
x1 <- function(t, lambda1_values, lambda2_values, sigma12) {
  phase <- findInterval(t, time_periods)
  lambda1 <- lambda1_values[phase]
  lambda2 <- lambda2_values[phase]
  num <- exp(-t * lambda2 * sigma12) * (-1 + exp(t * (-lambda1 - lambda2 + lambda2 * sigma12))) * lambda1
  denom <- (-lambda1 - lambda2 + lambda2 * sigma12)
  return(num / denom)
}

# Probability of being infected with serotype 2
x2 <- function(t, lambda1_values, lambda2_values, sigma21) {
  phase <- findInterval(t, time_periods)
  lambda1 <- lambda1_values[phase]
  lambda2 <- lambda2_values[phase]
  num <- exp(-t * lambda1 * sigma21) * (-1 + exp(t * (-lambda1 - lambda2 + lambda1 * sigma21))) * lambda2
  denom <- (-lambda1 - lambda2 + lambda1 * sigma21)
  return(num / denom)
}

# Probability of being infected with both serotypes
x12 <- function(t, lambda1_values, lambda2_values, sigma12, sigma21) {
  phase <- findInterval(t, time_periods)
  lambda1 <- lambda1_values[phase]
  lambda2 <- lambda2_values[phase]
  expPart <- exp(-t * lambda2 * sigma12 - t * lambda1 * sigma21)
  term <- exp(t * lambda1 * sigma21) * lambda1^2 -
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1^2 +
    exp(t * lambda2 * sigma12) * lambda1 * lambda2 +
    exp(t * lambda1 * sigma21) * lambda1 * lambda2 -
    2 * exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 +
    exp(t * lambda2 * sigma12) * lambda2^2 -
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda2^2 +
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma12 -
    exp(t * (-lambda1 - lambda2) + t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma12 -
    exp(t * lambda2 * sigma12) * lambda2^2 * sigma12 +
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda2^2 * sigma12 -
    exp(t * lambda1 * sigma21) * lambda1^2 * sigma21 +
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1^2 * sigma21 +
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma21 -
    exp(t * (-lambda1 - lambda2) + t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma21 -
    exp(t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma12 * sigma21 +
    exp(t * (-lambda1 - lambda2) + t * lambda2 * sigma12 + t * lambda1 * sigma21) * lambda1 * lambda2 * sigma12 * sigma21
  
  denom <- (-lambda1 - lambda2 + lambda2 * sigma12) * (lambda1 + lambda2 - lambda1 * sigma21)
  return(expPart * term / denom)
}

# --- Age structure ---
age_bins <- seq(5, 65, by = 5)
age_midpoints <- age_bins[-length(age_bins)] + diff(age_bins) / 2

# --- Data simulation ---
set.seed(234)

sim_data <- map_dfr(age_midpoints, function(age) {
  # Time elapsed for each individual (age)
  t <- age  # Time since birth is just their age
  
  # Determine phase based on cutoff years
  phase <- "piecewise"
  
  # Simulate using the time-varying functions for each phase
  p_s <- 0
  p_x1 <- 0
  p_x2 <- 0
  p_x12 <- 0
  
  # Loop through cutoff years
  for (i in 1:(length(cutoff_years) + 1)) {
    cutoff_year <- ifelse(i <= length(cutoff_years), cutoff_years[i], current_year)
    t_phase <- min(age, cutoff_year)  # Time until the next cutoff
    
    if (i == 1) {
      # First phase, before first cutoff
      p_s <- s(t_phase, lambda1_values, lambda2_values)  # Calculate probabilities for the first phase
      p_x1 <- x1(t_phase, lambda1_values, lambda2_values, sigma12)
      p_x2 <- x2(t_phase, lambda1_values, lambda2_values, sigma21)
      p_x12 <- x12(t_phase, lambda1_values, lambda2_values, sigma12, sigma21)
    } else {
      # Later phases, after cutoffs
      p_s <- p_s * s(t_phase, lambda1_values, lambda2_values)  # Update probabilities after each cutoff
      p_x1 <- p_x1 + p_s * x1(t_phase, lambda1_values, lambda2_values, sigma12)
      p_x2 <- p_x2 + p_s * x2(t_phase, lambda1_values, lambda2_values, sigma21)
      p_x12 <- p_x12 + p_s * x12(t_phase, lambda1_values, lambda2_values, sigma12, sigma21)
    }
  }
  
  # Ensure probabilities are valid (avoid NA or zero values)
  probs <- c(p_s, p_x1, p_x2, p_x12)
  
  # Normalize probabilities to ensure the sum is 1
  if (any(is.na(probs)) | any(probs <= 0)) {
    probs[is.na(probs) | probs <= 0] <- 1e-6  # Set invalid probabilities to a small value (near-zero)
  }
  
  probs <- probs / sum(probs)  # Normalize the probabilities to sum to 1
  
  # Simulate the number of individuals in each state
  n <- sample(1000:2000, 1)  # Random sample size for each age group
  counts <- rmultinom(1, size = n, prob = probs)
  
  tibble(
    age = age,
    phase = phase,
    n_tested = n,
    n_s = counts[1],
    n_denv1 = counts[2],
    n_denv2 = counts[3],
    n_denv12 = counts[4],
    p_s = probs[1],
    p_x1 = probs[2],
    p_x2 = probs[3],
    p_x12 = probs[4],
    sum = sum(p_s, p_x1, p_x2, p_x12)  ## should sum up to 1
  )
})

# Check the data
sim_data %>%
  mutate(diff_from_1 = sum - 1) %>%
  arrange(desc(diff_from_1))

# View first few rows of simulation results
head(sim_data)



# Prepare stan_data with correct y format
stan_data <- list(
  N = nrow(sim_data),  # Number of age groups
  K = length(cutoff_years),  # Number of phases (cutoff years)
  age = sim_data$age,
  n_tested = sim_data$n_tested,
  n_s = sim_data$n_s,
  n_denv1 = sim_data$n_denv1,
  n_denv2 = sim_data$n_denv2,
  n_denv12 = sim_data$n_denv12,
  cutoff_years = cutoff_years,  # Cutoff years for phases
  n_age_fine = length(age_fine),  # Length of the age_fine vector
  age_fine = age_fine,  # Fine age intervals
  y = cbind(sim_data$n_s, sim_data$n_denv1, sim_data$n_denv2, sim_data$n_denv12)  # counts for each compartment
)


# Compile the Stan model
stan_model <- stan_model("stan/time_varying_FOI.stan")

# Fit the model with simplified parameters for debugging
fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4, seed = 234)

# Check the results
print(fit, pars = c("lambda1", "lambda2", "sigma12", "sigma21"))

