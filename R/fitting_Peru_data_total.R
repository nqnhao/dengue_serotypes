library(readxl)
library(rstan)
library(table1)
library(dplyr)
library(ggplot2)
library(purrr)
#library(writexl)

## DATA CLEANING: SEROTYPES 1 & 2
data_peru <- read_excel("data/cohorts_data.xlsx")



data_cleaning <- function(data, year) {
  data %>%
    filter(sampleyear == year) %>%
    filter(denv3_status %in% c(0, "0", "NA"),
           denv4_status %in% c(0, "0", "NA")) %>%
    select(agegrp, denv1_status, denv2_status, denv3_status, denv4_status) %>%
    filter(!is.na(denv1_status), !is.na(denv2_status)) %>%
    mutate(
      denv1_status = case_when(
        denv1_status %in% c(1, "1") ~ 1,
        denv1_status %in% c(0, "0") ~ 0,
        TRUE ~ NA_real_
      ),
      denv2_status = case_when(
        denv2_status %in% c(1, "1") ~ 1,
        denv2_status %in% c(0, "0") ~ 0,
        TRUE ~ NA_real_
      ),
      age = case_when(
        #agegrp == "age05" ~ 5,
        #agegrp == "age05-10" ~ 7.5,
        agegrp %in% c("age05", "age05-10") ~ 7.5,
        agegrp == "age10-15" ~ 12.5,
        agegrp == "age15-20" ~ 17.5,
        agegrp == "age20-25" ~ 22.5,
        agegrp == "age25-30" ~ 27.5,
        agegrp == "age30-35" ~ 32.5,
        agegrp == "age35-40" ~ 37.5,
        agegrp == "age40-45" ~ 42.5,
        agegrp == "age45-50" ~ 47.5,
        agegrp == "age50-55" ~ 52.5,
        agegrp == "age55-60" ~ 57.5,
        agegrp == "age60+" ~ 62.5,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(age), !is.na(denv1_status), !is.na(denv2_status)) %>%
    mutate(
      only_denv1 = as.numeric(denv1_status == 1 & denv2_status == 0),
      only_denv2 = as.numeric(denv1_status == 0 & denv2_status == 1),
      both_denv12 = as.numeric(denv1_status == 1 & denv2_status == 1)
    ) %>%
    group_by(age) %>%
    summarise(
      n_tested = n(),
      n_denv1 = sum(only_denv1, na.rm = TRUE),
      n_denv2 = sum(only_denv2, na.rm = TRUE),
      n_denv12 = sum(both_denv12, na.rm = TRUE),
      n_s = n_tested - n_denv1 - n_denv2 - n_denv12,
      .groups = "drop"
    )
}


peru_2010 <- data_cleaning (data_peru, 2010)
peru_2008 <- data_cleaning (data_peru, 2008)
peru_2006 <- data_cleaning (data_peru, 2006)
peru_2004 <- data_cleaning (data_peru, 2004)
peru_2002 <- data_cleaning (data_peru, 2002)
peru_2001 <- data_cleaning (data_peru, 2001)
peru_1999 <- data_cleaning (data_peru, 1999)
peru_1996 <- data_cleaning (data_peru, 1996)
peru_1995 <- data_cleaning (data_peru, 1995)
peru_1994 <- data_cleaning (data_peru, 1994)
peru_1993 <- data_cleaning (data_peru, 1993)

peru_2010 <- peru_2010 %>%
  mutate(
    p_s = n_s/n_tested,
    p_x1 =  n_denv1/n_tested,
    p_x2 =  n_denv2/n_tested,
    p_x12 =  n_denv12/n_tested
  )
ggplot(data=peru_2010) + geom_line(aes(x = age, y = p_s, colour = "p_s")) + geom_line(aes(x = age, y = p_x1, colour = "p_x1")) + geom_line(aes(x= age, y = p_x2, colour = "p_x2")) + geom_line(aes(x = age, y = p_x12, colour = "p_x12")) + theme_minimal()


years <- c(2010,2008, 2006, 2004, 2002, 2001, 1999, 1996, 1995, 1994, 1993)
#years <- c(2010,2008, 2006, 2004, 2002, 2001, 1999)
peru_list <- map(years, function(y) {
  df <- get(paste0("peru_", y))
  df %>%
    mutate(
      year = y,
      p_s = n_s / n_tested,
      p_x1 = n_denv1 / n_tested,
      p_x2 = n_denv2 / n_tested,
      p_x12 = n_denv12 / n_tested,
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
})

peru_all <- bind_rows(peru_list) #%>% 
  #filter(year==2010)

peru_long <- peru_all %>%
  pivot_longer(cols = c(p_s, p_x1, p_x2, p_x12, q_x1, q_x2, q_x12, p_1_2, p_2_1, lower_q_x1, upper_q_x1, lower_q_x2, upper_q_x2),
               names_to = "prob_type",
               values_to = "probability") %>% 
  #filter(prob_type %in% c("p_x12", "q_x1", "q_x2", "q_x12")) #%>% 
  filter(prob_type %in% c("q_x1", "q_x2","lower_q_x1", "upper_q_x1", "lower_q_x2", "upper_q_x2"))
  #filter(prob_type %in% c("p_s", "p_x1", "p_x2", "p_x12"))
  #filter(prob_type %in% c("p_1_2", "p_2_1"))

ggplot(peru_all, aes(x = factor(age))) +
  geom_point(aes(y = q_x1, colour = "DENV1"), size = 2) +
  geom_errorbar(aes(ymin = lower_q_x1, ymax = upper_q_x1, colour = "DENV1"), width = 0.2) +
  geom_point(aes(y = q_x2, colour = "DENV2"), size = 2) +
  geom_errorbar(aes(ymin = lower_q_x2, ymax = upper_q_x2, colour = "DENV2"), width = 0.2) +
  facet_wrap(~ year) +
  theme_bw() +
  theme(legend.position = "top") + 
  labs(
    x = "Age", y = "Probability",
    colour = "Seroprevalence"
  )

ggplot(peru_long, aes(x = age, y = probability, colour = prob_type)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Age", y = "Probability", colour = "Probability Type") +
  facet_wrap(~ year)


peru_list_with_p <- map(years, function(y) {
  df <- get(paste0("peru_", y)) %>%
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

peru_all_pval <- bind_rows(peru_list_with_p)


ggplot(peru_all_pval, aes(x = age)) +
  geom_line(aes(y = p_x12, colour = "Observed p_x12")) +
  geom_line(aes(y = q_x12, colour = "Estimated q_x12")) +
  geom_text(aes(y = p_x12 + 0.02, label = ifelse(p_value < 0.05, "p<0.05", sprintf("p=%.3f", p_value))), size = 2, hjust = 0) +
  facet_wrap(~ year) +
  theme_minimal() +
  theme(legend.position = "top") + 
  labs(
    x = "Age", y = "Probability",
    colour = "Legend"
  )


ggplot(peru_all_pval, aes(x = factor(age))) +
  geom_point(aes(y = p_x12, colour = "Observed", shape = p_value < 0.05), size = 2.5) +
  geom_point(aes(y = q_x12, colour = "Expected"), size = 2) +
  geom_text(
    aes(
      y = p_x12 + 0.02,
      label = ifelse(p_value < 0.05, "p<0.05", sprintf("p=%.3f", p_value))
    ),
    size = 2, hjust = 0
  ) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19), name = "p < 0.05") +
  facet_wrap(~ year) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x = "Age midpoint", y = "Probability",
    colour = "Seroprevalence of DENV1 + DENV2"
  )



## FITTING DATA EACH YEAR


fit <- function(data){
  age_fine <- seq(5, max(data$age), by=1)
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
  stan_model <- stan_model("stan/simulation_denv_2_serotypes.stan")
  fit <- sampling(stan_model, data = stan_data, iter = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
}



fit_peru_2010 <- fit(peru_2010)
print(fit_peru_2010, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits = 3)

shinystan::launch_shinystan(fit_peru_2010)

#### Fit all datasets

fit_peru_2008 <- fit(peru_2008)
print(fit_peru_2008, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_2006 <- fit(peru_2006)
print(fit_peru_2006, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_2004 <- fit(peru_2004)
print(fit_peru_2004, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_2002 <- fit(peru_2002)
print(fit_peru_2002, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_2001 <- fit(peru_2001)
print(fit_peru_2001, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_1999 <- fit(peru_1999)
print(fit_peru_1999, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_1996 <- fit(peru_1996)
print(fit_peru_1996, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_1995 <- fit(peru_1995)
print(fit_peru_1995, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_1994 <- fit(peru_1994)
print(fit_peru_1994, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

fit_peru_1993 <- fit(peru_1993)
print(fit_peru_1993, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975), digits=3)

# List of model fits and associated years
fits <- list(
  `2010` = fit_peru_2010,
  `2008` = fit_peru_2008,
  `2006` = fit_peru_2006,
  `2004` = fit_peru_2004,
  `2002` = fit_peru_2002,
  `2001` = fit_peru_2001,
  `1999` = fit_peru_1999,
  `1996` = fit_peru_1996,
  `1995` = fit_peru_1995,
  `1994` = fit_peru_1994,
  `1993` = fit_peru_1993
)

# Function to extract summary for selected parameters
extract_summary <- function(fit, year) {
  summ <- summary(fit, pars = c("lambda1", "lambda2", "sigma12", "sigma21"), probs = c(0.025, 0.5, 0.975))$summary
  tibble(
    year = year,
    parameter = rownames(summ),
    mean = summ[, "mean"],
    q025 = summ[, "2.5%"],
    median = summ[, "50%"],
    q975 = summ[, "97.5%"]
  )
}

# Combine all into one tidy data frame
summary_all_years <- imap_dfr(fits, extract_summary)

summary_all_years

ggplot(summary_all_years, aes(x = factor(year), y = mean, ymin = q025, ymax = q975, color = parameter)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.4) +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(title = "Posterior Mean and 95% CI by Year", y = "Estimate", x = "Year") +
  theme_bw()

##PPC peru_2010
age_fine <- seq(5, max(peru_2010$age), by=1)
n_age_fine <- length(age_fine)
stan_data <- list(
  N = nrow(peru_2010),
  age = peru_2010$age,
  n_tested = peru_2010$n_tested,
  n_s = peru_2010$n_s,
  n_denv1 = peru_2010$n_denv1,
  n_denv2 = peru_2010$n_denv2,
  n_denv12 = peru_2010$n_denv12,
  y = as.matrix(peru_2010[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
  age_fine = age_fine,
  n_age_fine = n_age_fine
)



posterior <- rstan::extract(fit_peru_2010)

y <- stan_data$y
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

##PPC peru_1996
age_fine <- seq(5, max(peru_1996$age), by=1)
n_age_fine <- length(age_fine)
stan_data <- list(
  N = nrow(peru_1996),
  age = peru_1996$age,
  n_tested = peru_1996$n_tested,
  n_s = peru_1996$n_s,
  n_denv1 = peru_1996$n_denv1,
  n_denv2 = peru_1996$n_denv2,
  n_denv12 = peru_1996$n_denv12,
  y = as.matrix(peru_1996[, c("n_s", "n_denv1", "n_denv2", "n_denv12")]),
  age_fine = age_fine,
  n_age_fine = n_age_fine
)



posterior <- rstan::extract(fit_peru_1996)

y <- stan_data$y  
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



