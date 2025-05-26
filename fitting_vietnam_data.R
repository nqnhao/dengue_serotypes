
library(readxl)
library(rstan)
library(table1)
library(tidyverse)
library(ggplot2)
library(purrr)

### VIETNAM
setwd("/Users/nguyenquannhuhao/Documents/RStudio/Dengue_placement/Hao")


data_vietnam <- read_excel("/Users/nguyenquannhuhao/Documents/RStudio/Dengue_placement/cleaned_data_dengue_vietnam.xlsx")

### Serotype 1 & 2
data_vietnam_total <- data_vietnam %>%
  select("Age", "DV1", "DV2", "DV3", "DV4") %>% 
  rename("age" = "Age") %>% 
  mutate(
    denv1_status = case_when(
      DV1 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv2_status = case_when(
      DV2 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv3_status = case_when(
      DV3 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv4_status = case_when(
      DV4 > 5 ~ 1,
      TRUE ~ 0
    ),
    age = case_when(
      age > 0 & age <= 5 ~ 2.5,
      age > 5  & age <= 10 ~ 7.5,
      age > 10 & age <= 15 ~ 12.5,
      age > 15 & age <= 20 ~ 17.5,
      age > 20 & age <= 25 ~ 22.5,
      age > 25 & age <= 30 ~ 27.5,
      age > 30 & age <= 35 ~ 32.5,
      age > 35 & age <= 40 ~ 37.5,
      age > 40 & age <= 45 ~ 42.5,
      age > 45 & age <= 50 ~ 47.5,
      age > 50 & age <= 55 ~ 52.5,
      age > 55 & age <= 60 ~ 57.5,
      age > 60 ~ 62.5,
      TRUE ~ NA_real_
    ),
    only_denv1 = as.numeric(denv1_status == 1 & denv2_status == 0),
    only_denv2 = as.numeric(denv1_status == 0 & denv2_status == 1),
    both_denv12 = as.numeric(denv1_status == 1 & denv2_status == 1)
  ) %>%
  filter(denv3_status %in% c(0, "0", "NA"),
         denv4_status %in% c(0, "0", "NA")) %>%
  group_by(age) %>%
  summarise(
    n_tested = n(),
    n_denv1 = sum(only_denv1, na.rm = TRUE),
    n_denv2 = sum(only_denv2, na.rm = TRUE),
    n_denv12 = sum(both_denv12, na.rm = TRUE),
    n_s = n_tested - n_denv1 - n_denv2 - n_denv12,
    .groups = "drop"
  )


data_vietnam_HCM <- data_vietnam %>%
  filter(Site == "HC") %>% 
  select("Age", "DV1", "DV2", "DV3", "DV4") %>% 
  rename("age" = "Age") %>% 
  mutate(
    denv1_status = case_when(
      DV1 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv2_status = case_when(
      DV2 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv3_status = case_when(
      DV3 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv4_status = case_when(
      DV4 > 5 ~ 1,
      TRUE ~ 0
    ),
    age = case_when(
      age > 0 & age <= 5 ~ 2.5,
      age > 5  & age <= 10 ~ 7.5,
      age > 10 & age <= 15 ~ 12.5,
      age > 15 & age <= 20 ~ 17.5,
      age > 20 & age <= 25 ~ 22.5,
      age > 25 & age <= 30 ~ 27.5,
      age > 30 & age <= 35 ~ 32.5,
      age > 35 & age <= 40 ~ 37.5,
      age > 40 & age <= 45 ~ 42.5,
      age > 45 & age <= 50 ~ 47.5,
      age > 50 & age <= 55 ~ 52.5,
      age > 55 & age <= 60 ~ 57.5,
      age > 60 ~ 62.5,
      TRUE ~ NA_real_
    ),
    only_denv1 = as.numeric(denv1_status == 1 & denv2_status == 0),
    only_denv2 = as.numeric(denv1_status == 0 & denv2_status == 1),
    both_denv12 = as.numeric(denv1_status == 1 & denv2_status == 1)
  ) %>%
  filter(denv3_status %in% c(0, "0", "NA"),
         denv4_status %in% c(0, "0", "NA")) %>%
  group_by(age) %>%
  summarise(
    n_tested = n(),
    n_denv1 = sum(only_denv1, na.rm = TRUE),
    n_denv2 = sum(only_denv2, na.rm = TRUE),
    n_denv12 = sum(both_denv12, na.rm = TRUE),
    n_s = n_tested - n_denv1 - n_denv2 - n_denv12,
    .groups = "drop"
  )

data_vietnam_KH <- data_vietnam %>%
  filter(Site == "KH") %>% 
  select("Age", "DV1", "DV2", "DV3", "DV4") %>% 
  rename("age" = "Age") %>% 
  mutate(
    denv1_status = case_when(
      DV1 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv2_status = case_when(
      DV2 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv3_status = case_when(
      DV3 > 5 ~ 1,
      TRUE ~ 0
    ),
    denv4_status = case_when(
      DV4 > 5 ~ 1,
      TRUE ~ 0
    ),
    age = case_when(
      age > 0 & age <= 5 ~ 2.5,
      age > 5  & age <= 10 ~ 7.5,
      age > 10 & age <= 15 ~ 12.5,
      age > 15 & age <= 20 ~ 17.5,
      age > 20 & age <= 25 ~ 22.5,
      age > 25 & age <= 30 ~ 27.5,
      age > 30 & age <= 35 ~ 32.5,
      age > 35 & age <= 40 ~ 37.5,
      age > 40 & age <= 45 ~ 42.5,
      age > 45 & age <= 50 ~ 47.5,
      age > 50 & age <= 55 ~ 52.5,
      age > 55 & age <= 60 ~ 57.5,
      age > 60 ~ 62.5,
      TRUE ~ NA_real_
    ),
    only_denv1 = as.numeric(denv1_status == 1 & denv2_status == 0),
    only_denv2 = as.numeric(denv1_status == 0 & denv2_status == 1),
    both_denv12 = as.numeric(denv1_status == 1 & denv2_status == 1)
  ) %>%
  filter(denv3_status %in% c(0, "0", "NA"),
         denv4_status %in% c(0, "0", "NA")) %>%
  group_by(age) %>%
  summarise(
    n_tested = n(),
    n_denv1 = sum(only_denv1, na.rm = TRUE),
    n_denv2 = sum(only_denv2, na.rm = TRUE),
    n_denv12 = sum(both_denv12, na.rm = TRUE),
    n_s = n_tested - n_denv1 - n_denv2 - n_denv12,
    .groups = "drop"
  )

fit <- function(data){
  age_fine <- seq(0, max(data$age), by=1)
  n_age_fine <- length(age_fine)
  stan_data <- list(
    N = nrow(data),
    age = data$age,
    n_tested = data$n_tested,
    n_s = data$n_s,
    n_denv1 = data$n_denv1,
    n_denv2 = data$n_denv2,
    n_denv12 = data$n_denv12,
    y = as.matrix(data[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
    age_fine = age_fine,
    n_age_fine = n_age_fine
  )
  stan_model <- stan_model("simulation_denv_2_serotypes.stan")
  fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4)
}

fit_Vietnam_total <- fit(data_vietnam_total)
print(fit_Vietnam_total, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits = 4)

fit_Vietnam_HCM <- fit(data_vietnam_HCM)
print(fit_Vietnam_HCM, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits = 4)

fit_Vietnam_KH <- fit(data_vietnam_KH)
print(fit_Vietnam_KH, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits = 4)

## Data visualisation

data_vietnam <- data_vietnam %>%
  mutate(
    p_s = n_s/n_tested,
    p_x1 =  n_denv1/n_tested,
    p_x2 =  n_denv2/n_tested,
    p_x12 =  n_denv12/n_tested,
    q_x1 = (n_denv1 + n_denv12) / n_tested,
    q_x2 = (n_denv2 + n_denv12) / n_tested,
    q_x12 = q_x1 * q_x2,
    p_1_2 = p_x12 / q_x2,
    p_2_1 = p_x12 / q_x1,
    lower_q_x1 = q_x1 - 1.96 * sqrt(q_x1 * (1 - q_x1) / n_tested),
    upper_q_x1 = q_x1 + 1.96 * sqrt(q_x1 * (1 - q_x1) / n_tested),
    lower_q_x2 = q_x2 - 1.96 * sqrt(q_x2 * (1 - q_x2) / n_tested),
    upper_q_x2 = q_x2 + 1.96 * sqrt(q_x2 * (1 - q_x2) / n_tested),
  )
ggplot(data=data_vietnam) + geom_line(aes(x = age, y = p_s, colour = "p_s")) + geom_line(aes(x = age, y = p_x1, colour = "p_x1")) + geom_line(aes(x= age, y = p_x2, colour = "p_x2")) + geom_line(aes(x = age, y = p_x12, colour = "p_x12"))+ theme_minimal()


vietnam_long <- data_vietnam %>%
  pivot_longer(cols = c(p_s, p_x1, p_x2, p_x12, q_x1, q_x2, q_x12, p_1_2, p_2_1, lower_q_x1, upper_q_x1, lower_q_x2, upper_q_x2),
               names_to = "prob_type",
               values_to = "probability") %>% 
  #filter(prob_type %in% c("p_x12", "q_x1", "q_x2", "q_x12")) #%>% 
  filter(prob_type %in% c("q_x1", "q_x2","lower_q_x1", "upper_q_x1", "lower_q_x2", "upper_q_x2"))
#filter(prob_type %in% c("p_s", "p_x1", "p_x2", "p_x12"))
#filter(prob_type %in% c("p_1_2", "p_2_1"))

ggplot(data_vietnam, aes(x = factor(age))) +
  geom_point(aes(y = q_x1, colour = "DENV1"), size = 2) +
  geom_errorbar(aes(ymin = lower_q_x1, ymax = upper_q_x1, colour = "DENV1"), width = 0.2) +
  geom_point(aes(y = q_x2, colour = "DENV2"), size = 2) +
  geom_errorbar(aes(ymin = lower_q_x2, ymax = upper_q_x2, colour = "DENV2"), width = 0.2) +
  theme_bw() +
  theme(legend.position = "top") + 
  labs(
    x = "Age", y = "Probability",
    colour = "Seroprevalence"
  )

## PPC - Total data
age_fine <- seq(0, max(data_vietnam_total$age), by=1)
n_age_fine <- length(age_fine)
stan_data <- list(
  N = nrow(data_vietnam_total),
  age = data_vietnam_total$age,
  n_tested = data_vietnam_total$n_tested,
  n_s = data_vietnam_total$n_s,
  n_denv1 = data_vietnam_total$n_denv1,
  n_denv2 = data_vietnam_total$n_denv2,
  n_denv12 = data_vietnam_total$n_denv12,
  y = as.matrix(data_vietnam_total[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
  age_fine = age_fine,
  n_age_fine = n_age_fine
)



posterior <- rstan::extract(fit_Vietnam_total)

y <- stan_data$y                 # observed counts: 13 x 4
n_tested <- stan_data$n_tested  
age <- stan_data$age            
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
  geom_ribbon(data = df_smooth, aes(x = age, ymin = lower, ymax = upper), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = df_smooth, aes(x = age, y = median), color = "red", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(y = observed), color = "black", size = 2) +
  geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs), width = 0.4) +
  
  facet_wrap(~ compartment) +
  labs(
    title = "Posterior Predictive Check",
    x = "Age", y = "Probability"
  ) +
  theme_minimal()

## PPC - HCM data
age_fine <- seq(0, max(data_vietnam_HCM$age), by=1)
n_age_fine <- length(age_fine)
stan_data <- list(
  N = nrow(data_vietnam_HCM),
  age = data_vietnam_HCM$age,
  n_tested = data_vietnam_HCM$n_tested,
  n_s = data_vietnam_HCM$n_s,
  n_denv1 = data_vietnam_HCM$n_denv1,
  n_denv2 = data_vietnam_HCM$n_denv2,
  n_denv12 = data_vietnam_HCM$n_denv12,
  y = as.matrix(data_vietnam_HCM[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
  age_fine = age_fine,
  n_age_fine = n_age_fine
)



posterior <- rstan::extract(fit_Vietnam_HCM)

y <- stan_data$y                 # observed counts: 13 x 4
n_tested <- stan_data$n_tested  
age <- stan_data$age            
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
  geom_ribbon(data = df_smooth, aes(x = age, ymin = lower, ymax = upper), fill = "red", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = df_smooth, aes(x = age, y = median), color = "red", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(y = observed), color = "black", size = 2) +
  geom_errorbar(aes(ymin = lower_obs, ymax = upper_obs), width = 0.4) +
  
  facet_wrap(~ compartment) +
  labs(
    title = "Posterior Predictive Check",
    x = "Age", y = "Probability"
  ) +
  theme_minimal()



## Check Indepence of 2 serotypes

vietnam_list_with_p <- data_vietnam %>%
    mutate(
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

vietnam_all_pval <- bind_rows(vietnam_list_with_p)

ggplot(vietnam_all_pval, aes(x = age)) +
  geom_line(aes(y = p_x12, colour = "Observed p_x12")) +
  geom_line(aes(y = q_x12, colour = "Estimated q_x12")) +
  geom_text(aes(y = p_x12 + 0.02, label = ifelse(p_value < 0.05, "p<0.05", sprintf("p=%.3f", p_value))), size = 2, hjust = 0) +
  theme_minimal() +
  theme(legend.position = "top") + 
  labs(
    x = "Age", y = "Probability",
    colour = "Legend"
  )


ggplot(vietnam_all_pval, aes(x = factor(age))) +
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
