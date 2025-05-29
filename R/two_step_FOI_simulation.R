# -----------------------------------------
# Forward Time Two-Step FOI Simulation
# -----------------------------------------

library(tidyverse)

# --- Parameters ---
lambda1_a <- 0.05
lambda1_b <- 0.01
lambda2   <- 0.04
sigma12 <- 15
sigma21 <- 0.1

current_year <- 2025
cutoff_year <- 2000
cutoff_age <- current_year - cutoff_year

# --- Model functions (unchanged) ---
s <- function(t, lambda1, lambda2) {
  exp(-t * (lambda1 + lambda2))
}

x1 <- function(t, lambda1, lambda2, sigma12) {
  num <- exp(-t * lambda2 * sigma12) * (-1 + exp(t * (-lambda1 - lambda2 + lambda2 * sigma12))) * lambda1
  denom <- (-lambda1 - lambda2 + lambda2 * sigma12)
  return(num / denom)
}

x2 <- function(t, lambda1, lambda2, sigma21) {
  num <- exp(-t * lambda1 * sigma21) * (-1 + exp(t * (-lambda1 - lambda2 + lambda1 * sigma21))) * lambda2
  denom <- (-lambda1 - lambda2 + lambda1 * sigma21)
  return(num / denom)
}

x12 <- function(t, lambda1, lambda2, sigma12, sigma21) {
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
  if (age <= cutoff_age) {
    # Fully exposed after the cut-off year
    lambda1 <- lambda1_b
    phase <- "after_2015"
    p_s   <- s(age, lambda1, lambda2)
    p_x1  <- x1(age, lambda1, lambda2, sigma12)
    p_x2  <- x2(age, lambda1, lambda2, sigma21)
    p_x12 <- x12(age, lambda1, lambda2, sigma12, sigma21)
  } else {
    # Two-phase exposure
    phase <- "two_step"
    t1 <- age - cutoff_age  # Years before the cut-off year (lambda1_a)
    t2 <- cutoff_age        # Years after the cut-off year (lambda1_b)
    
    # Step 1: simulate years before the cut-off year under lambda1_a
    S1    <- s(t1, lambda1_a, lambda2)
    X1_1  <- x1(t1, lambda1_a, lambda2, sigma12)
    X2_1  <- x2(t1, lambda1_a, lambda2, sigma21)
    X12_1 <- x12(t1, lambda1_a, lambda2, sigma12, sigma21)
    
    # Step 2: continue (current year - cut-off year) more years under lambda1_b
    p_s <- S1 * s(t2, lambda1_b, lambda2)
    #p_x1 <- X1_1 * s(t2, lambda1_b, lambda2) + S1 * x1(t2, lambda1_b, lambda2, sigma12)
    p_x1 <- X1_1 + S1 * x1(t2, lambda1_b, lambda2, sigma12)
    #p_x2 <- X2_1 * s(t2, lambda1_b, lambda2) + S1 * x2(t2, lambda1_b, lambda2, sigma21)
    p_x2 <- X2_1 + S1 * x2(t2, lambda1_b, lambda2, sigma21)
    p_x12 <- X12_1 +
      X1_1 * x2(t2, lambda1_b, lambda2, sigma21) +
      X2_1 * x1(t2, lambda1_b, lambda2, sigma12) +
      S1 * x12(t2, lambda1_b, lambda2, sigma12, sigma21)
  }
  
  probs <- c(p_s, p_x1, p_x2, p_x12)
  #probs <- pmax(probs, 0)
  #probs <- probs / sum(probs)
  
  n <- sample(1000:2000, 1)
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
    sum = sum(p_s, p_x1, p_x2, p_x12)
  )
})

sim_data

# FITTING
age_fine <- seq(5, max(sim_data$age), by=1)
n_age_fine <- length(age_fine)
stan_data <- list(
  N = nrow(sim_data),
  age = sim_data$age,
  n_tested = sim_data$n_tested,
  n_s = sim_data$n_s,
  n_denv1 = sim_data$n_denv1,
  n_denv2 = sim_data$n_denv2,
  n_denv12 = sim_data$n_denv12,
  y = as.matrix(sim_data[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
  age_fine = age_fine,
  n_age_fine = n_age_fine
)



stan_model <- stan_model("stan/simulation_denv_2_serotypes.stan")
fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4, seed = 234)


print(fit, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), digits =4)

### Posterior predictive check

posterior <- rstan::extract(fit)

y <- stan_data$y                 # observed counts: 13 x 4
n_tested <- stan_data$n_tested  
age <- stan_data$age            

# Predicted category probabilities = count / n_tested
# y_rep: iterations x N x 4
pred_probs <- array(NA, dim = dim(posterior$y_rep))  # 4000 x 13 x 4

for (j in 1:4) {
  pred_probs[,,j] <- sweep(posterior$y_rep[,,j], 2, n_tested, "/")
}

make_pred_df <- function(prob_matrix, age_vec, comp_name, obs_vals, n_vals) {
  ci <- t(apply(prob_matrix, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  tibble(
    age = age_vec,
    median = ci[, 2],
    lower = ci[, 1],
    upper = ci[, 3],
    observed = obs_vals / n_vals,
    lower_obs = observed - 1.96 * sqrt(observed * (1 - observed) / n_vals),
    upper_obs = observed + 1.96 * sqrt(observed * (1 - observed) / n_vals),
    compartment = comp_name
  )
}

prob_s <- 1 - pred_probs[, , 2] - pred_probs[, , 3] - pred_probs[, , 4]
df_s   <- make_pred_df(prob_s, age, "Susceptible", y[, 1], n_tested)
df_d1   <- make_pred_df(pred_probs[,,2], age, "DENV1 only", y[, 2], n_tested)
df_d2   <- make_pred_df(pred_probs[,,3], age, "DENV2 only", y[, 3], n_tested)
df_d12  <- make_pred_df(pred_probs[,,4], age, "DENV1 + DENV2", y[, 4], n_tested)

df_ppc <- bind_rows(df_s, df_d1, df_d2, df_d12)

get_ci_fine <- function(mat) {
  apply(mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
}

ci_s    <- get_ci_fine(posterior$prob_s)
ci_d1   <- get_ci_fine(posterior$prob_x1)
ci_d2   <- get_ci_fine(posterior$prob_x2)
ci_d12  <- get_ci_fine(posterior$prob_x12)

age_fine <- stan_data$age_fine

df_smooth <- tibble(
  age = rep(age_fine, 4),
  median = c(ci_s[2,], ci_d1[2,], ci_d2[2,], ci_d12[2,]),
  lower = c(ci_s[1,], ci_d1[1,], ci_d2[1,], ci_d12[1,]),
  upper = c(ci_s[3,], ci_d1[3,], ci_d2[3,], ci_d12[3,]),
  compartment = rep(c("Susceptible", "DENV1 only", "DENV2 only", "DENV1 + DENV2"), each = length(age_fine))
)

df_smooth <- df_smooth %>%
  mutate(compartment = factor(compartment, levels = c("Susceptible", "DENV1 only", "DENV2 only", "DENV1 + DENV2")))

df_ppc <- df_ppc %>%
  mutate(compartment = factor(compartment, levels = c("Susceptible", "DENV1 only", "DENV2 only", "DENV1 + DENV2")))

ggplot(df_ppc, aes(x = age)) +
  # smooth line from age_fine
  geom_ribbon(data = df_smooth, aes(x = age, ymin = lower, ymax = upper), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = df_smooth, aes(x = age, y = median), color = "red", linewidth = 0.8, inherit.aes = FALSE) +
  
  # observed
  geom_point(aes(y = observed), color = "black", size = 2) +
  geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs), width = 0.4) +
  
  facet_wrap(~ compartment) +
  labs(
    title = "Posterior Predictive Check",
    x = "Age", y = "Probability"
  ) +
  theme_minimal()
