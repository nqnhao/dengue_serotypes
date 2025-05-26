setwd("/Users/nguyenquannhuhao/Documents/RStudio/Dengue_placement/Hao")
library(rstan)
library(tibble)
library(tidyverse)
library(purrr)
#set.seed(123)

## SIMULATING DATA
age_bins = seq(0, 65, by = 5)
n_per_bin = 100
lambda1 = 0.05
lambda2 = 0.04
alpha12 = 1.5
alpha21 = 0.7
bin_midpoints <- age_bins[-length(age_bins)] + diff(age_bins) / 2

sim_data <- map_dfr(bin_midpoints, function(age_mid) {
  n <- n_per_bin
  age <- rep(age_mid, n)
  p1_naive <- 1 - exp(-lambda1 * age)
  p2_naive <- 1 - exp(-lambda2 * age)
  adj_age1 <- age * (alpha21 * p2_naive + (1 - p2_naive))
  adj_age2 <- age * (alpha12 * p1_naive + (1 - p1_naive))
  p1 <- 1 - exp(-lambda1 * adj_age1)
  p2 <- 1 - exp(-lambda2 * adj_age2)
  sero1 <- rbinom(n, 1, p1)
  sero2 <- rbinom(n, 1, p2)
  tibble(age = age, sero1 = sero1, sero2 = sero2)
  
  # Simulate DENV1
  #p1 <- 1 - exp(-lambda1 * age)
  #sero1 <- rbinom(n, 1, p1)
  
  # Simulate DENV2 given DENV1
  #adj_age2 <- ifelse(sero1 == 1, age * alpha12, age)
  #p2 <- 1 - exp(-lambda2 * adj_age2)
  #sero2 <- rbinom(n, 1, p2)
  
  # Re-simulate DENV1 given DENV2
  #adj_age1 <- ifelse(sero2 == 1, age * alpha21, age)
  #p1_adj <- 1 - exp(-lambda1 * adj_age1)
  #sero1_final <- rbinom(n, 1, p1_adj)
  
  #tibble(age = age_mid, sero1 = sero1_final, sero2 = sero2)
  
})

sim_table <- sim_data %>%
  group_by(age) %>%
  summarise(
    n_tested = n(),
    n_denv1 = sum(sero1),
    n_denv2 = sum(sero2),
    .groups = "drop"
  )

stan_data <- list(
  N = nrow(sim_table),
  age = sim_table$age,
  n_tested = sim_table$n_tested,
  n_denv1 = sim_table$n_denv1,
  n_denv2 = sim_table$n_denv2
)

## FITTING
stanmodel <- stan_model("denv_2serotype_simulation.stan")
fit <- sampling(stanmodel, data = stan_data, iter = 2000, chains = 4)
print(fit, pars = c("lambda1", "lambda2", "sigma12", "sigma21"))




posterior <- rstan::extract(fit)

ci_denv1 <- t(apply(posterior$pred_denv1, 2, quantile, probs = c(0.025, 0.5, 0.975)))
ci_denv2 <- t(apply(posterior$pred_denv2, 2, quantile, probs = c(0.025, 0.5, 0.975)))

data_plot <- sim_table %>%
  mutate(
    pred_denv1 = ci_denv1[, 2],
    lower1 = ci_denv1[, 1],
    upper1 = ci_denv1[, 3],
    pred_denv2 = ci_denv2[, 2],
    lower2 = ci_denv2[, 1],
    upper2 = ci_denv2[, 3]
  )

ggplot(data_plot, aes(x = age)) +
  geom_point(aes(y = n_denv1), color = "black", size = 2) +
  geom_line(aes(y = pred_denv1), color = "red") +
  geom_ribbon(aes(ymin = lower1, ymax = upper1), fill = "red", alpha = 0.2) +
  labs(title = "Posterior Predictive Check: DENV1", y = "Seropositives", x = "Age")

ggplot(data_plot, aes(x = age)) +
  geom_point(aes(y = n_denv2), color = "black", size = 2) +
  geom_line(aes(y = pred_denv2), color = "blue") +
  geom_ribbon(aes(ymin = lower2, ymax = upper2), fill = "blue", alpha = 0.2) +
  labs(title = "Posterior Predictive Check: DENV2", y = "Seropositives", x = "Age")
