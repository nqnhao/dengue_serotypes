library(tidyverse)
library(rstan)
library(purrr)
library(patchwork)
library(bayesplot)
library(ggplot2)

# Susceptible to both serotypes
s <- function(t, lambda1, lambda2) {
  exp(-t * (lambda1 + lambda2))
}

# Infected with DENV1 only
x1 <- function(t, lambda1, lambda2, sigma12) {
  num <- exp(-t * lambda2 * sigma12) * (-1 + exp(t * (-lambda1 - lambda2 + lambda2 * sigma12))) * lambda1
  denom <- (-lambda1 - lambda2 + lambda2 * sigma12)
  return(num / denom)
}

# Infected with DENV2 only
x2 <- function(t, lambda1, lambda2, sigma21) {
  num <- exp(-t * lambda1 * sigma21) * (-1 + exp(t * (-lambda1 - lambda2 + lambda1 * sigma21))) * lambda2
  denom <- (-lambda1 - lambda2 + lambda1 * sigma21)
  return(num / denom)
}

# Infected with both serotypes
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


# Data simulation

lambda1_fn <- function(age) ifelse(age <= 40, 0.07, 0.03)
lambda2_fn <- function(age) ifelse(age <= 40, 0.08, 0.04)
sigma12 <- 1.5  # Enhancement of DENV2 after DENV1
sigma21 <- 0.7  # Protection of DENV1 after DENV2

age_bins <- seq(5, 65, by = 5)
age_midpoints <- age_bins[-length(age_bins)] + diff(age_bins) / 2


set.seed(123)

sim_data <- map_dfr(age_midpoints, function(age) {
  l1 <- lambda1_fn(age)
  l2 <- lambda2_fn(age)
  
  p_s <- s(age, l1, l2)
  p_x1 <- x1(age, l1, l2, sigma12)
  p_x2 <- x2(age, l1, l2, sigma21)
  p_x12 <- x12(age, l1, l2, sigma12, sigma21)
  
  probs <- c(p_s, p_x1, p_x2, p_x12)
  n_per_bin <- sample(100:200, 1)
  
  counts <- rmultinom(1, size = n_per_bin, prob = probs)
  
  tibble(
    age = age,
    n_tested = n_per_bin,
    n_s = counts[1],
    n_denv1 = counts[2],
    n_denv2 = counts[3],
    n_denv12 = counts[4],
    p_s = p_s,
    p_x1 =  p_x1,
    p_x2 =  p_x2,
    p_x12 =  p_x12,
    lambda1 = l1,
    lambda2 = l2
  )
})


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


print(fit, pars = c("lambda1", "lambda2", "sigma12", "sigma21"))
shinystan::launch_shinystan(fit)

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
  
  # predicted from y_rep
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "skyblue", alpha = 0.3) +
  #geom_line(aes(y = median), color = "blue") +
  
  # observed
  geom_point(aes(y = observed), color = "black", size = 2) +
  geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs), width = 0.4) +
  
  facet_wrap(~ compartment) +
  labs(
    title = "Posterior Predictive Check",
    x = "Age", y = "Probability"
  ) +
  theme_minimal()



########
## Count PPC
y_rep <- rstan::extract(fit)$y_rep
y <- stan_data$y
category_names <- c("Susceptible", "DENV1-only", "DENV2-only", "DENV1+DENV2")

ppc_plots <- map(1:4, function(k) {
  obs_counts <- y[, k]
  ppc_counts <- y_rep[, , k]  # iterations x N
  
  bayesplot::ppc_intervals(
    y = obs_counts,
    yrep = ppc_counts,
    x = stan_data$age,
    prob = 0.8,
    prob_outer = 0.95
  ) +
    labs(
      title = paste("Posterior Predictive Check:", category_names[k]),
      x = "Age",
      y = "Count"
    ) +
    theme_minimal()
})

wrap_plots(ppc_plots, ncol = 2)

## Probability PPC
ppc_prob_plots <- map(1:4, function(k) {
  obs_probs <- y[, k] / stan_data$n_tested
  pred_counts <- as.matrix(y_rep[, , k])
  pred_probs <- sweep(pred_counts, 2, stan_data$n_tested, FUN = "/")
  bayesplot::ppc_intervals(
    y = obs_probs,
    yrep = pred_probs,
    x = stan_data$age,
    prob = 0.8,
    prob_outer = 0.95
  ) +
    labs(
      title = paste("Posterior Predictive Check (Probability):", category_names[k]),
      x = "Age",
      y = "Probability"
    ) +
    theme_minimal()
})

ppc_prob_plots <- map(1:4, function(k) {
  obs_probs <- y[, k] / stan_data$n_tested
  pred_counts <- as.matrix(y_rep[, , k])
  pred_probs <- sweep(pred_counts, 2, stan_data$n_tested, FUN = "/")
  bayesplot::ppc_ribbon(
    y = obs_probs,
    yrep = pred_probs,
    x = stan_data$age,
    prob = 0.8,
    prob_outer = 0.95,
    y_draw = c("points")
  ) +
    labs(
      title = paste("Posterior Predictive Check (Probability):", category_names[k]),
      x = "Age",
      y = "Probability"
    ) +
    theme_minimal()
})

wrap_plots(ppc_prob_plots, ncol = 2)


### Check Fisher test
sim_list_with_p <- map(years, function(y) {
  df <- sim_data %>%
    mutate(
      year = y,
      p_s = n_s / n_tested,
      p_x1 = n_denv1 / n_tested,
      p_x2 = n_denv2 / n_tested,
      p_x12 = n_denv12 / n_tested,
      q_x1 = (n_denv1 + n_denv12) / n_tested,
      q_x2 = (n_denv2 + n_denv12) / n_tested,
      q_x12 = q_x1 * q_x2
    ) %>%
    rowwise() %>%
    mutate(
      expected_infected = round(n_tested * q_x12),
      expected_uninfected = n_tested - expected_infected,
      observed_uninfected = n_tested - n_denv12,
      p_value = fisher.test(
        matrix(c(n_denv12, observed_uninfected,
                 expected_infected, expected_uninfected), nrow = 2)
      )$p.value
    ) %>%
    ungroup()
  return(df)
})

sim_all_pval <- bind_rows(sim_list_with_p)

ggplot(sim_all_pval, aes(x = age)) +
  geom_line(aes(y = p_x12, colour = "Observed p_x12")) +
  geom_line(aes(y = q_x12, colour = "Estimated q_x12")) +
  geom_text(aes(y = p_x12 + 0.02, label = ifelse(p_value < 0.05, "p<0.05", sprintf("p=%.3f", p_value))), size = 2, hjust = 0) +
  theme_minimal() +
  theme(legend.position = "top") + 
  labs(
    x = "Age", y = "Probability",
    colour = "Legend"
  )


ggplot(sim_all_pval, aes(x = factor(age))) +
  geom_point(aes(y = p_x12, colour = "Observed", shape = p_value < 0.05), size = 3) +
  geom_point(aes(y = q_x12, colour = "Expected"), size = 3) +
  geom_text(
    aes(
      y = p_x12 + 0.02,
      label = ifelse(p_value < 0.05, "p<0.05", sprintf("p=%.3f", p_value))
    ),
    size = 3, hjust = 0
  ) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19), name = "p < 0.05") +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x = "Age midpoint", y = "Probability",
    colour = "Seroprevalence of DENV1 + DENV2"
  )