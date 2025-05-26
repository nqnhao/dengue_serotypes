functions {
  real cumulative_hazard(real t, vector lambda, vector breaks, int K) {
    real h = 0;
    for (k in 1:K) {
      real t_start = breaks[k];
      real t_end = breaks[k + 1];
      if (t > t_end) {
        h += lambda[k] * (t_end - t_start);
      } else if (t > t_start) {
        h += lambda[k] * (t - t_start);
        break;
      }
    }
    return h;
  }

  real p_s(real t, vector lambda1, vector lambda2, vector breaks, int K) {
    return exp(-cumulative_hazard(t, lambda1, breaks, K) - cumulative_hazard(t, lambda2, breaks, K));
  }

  real p_x1(real t, vector lambda1, vector lambda2, real sigma12, vector breaks, int K) {
    real L1 = cumulative_hazard(t, lambda1, breaks, K);
    real L2 = cumulative_hazard(t, lambda2, breaks, K);
    real L2_sigma = cumulative_hazard(t, lambda2 * sigma12, breaks, K);

    real num = exp(-L2_sigma) * (-1 + exp(-L1 - L2 + L2_sigma)) * L1;
    real den = (-L1 - L2 + L2_sigma);
    return num / den;
  }

  real p_x2(real t, vector lambda1, vector lambda2, real sigma21, vector breaks, int K) {
    real L1 = cumulative_hazard(t, lambda1, breaks, K);
    real L2 = cumulative_hazard(t, lambda2, breaks, K);
    real L1_sigma = cumulative_hazard(t, lambda1 * sigma21, breaks, K);

    real num = exp(-L1_sigma) * (-1 + exp(-L1 - L2 + L1_sigma)) * L2;
    real den = (-L1 - L2 + L1_sigma);
    return num / den;
  }

  real p_x12(real t, vector lambda1, vector lambda2, real sigma12, real sigma21, vector breaks, int K) {
    return 1 - p_s(t, lambda1, lambda2, breaks, K)
             - p_x1(t, lambda1, lambda2, sigma12, breaks, K)
             - p_x2(t, lambda1, lambda2, sigma21, breaks, K);
  }
}

data {
  int<lower=1> N;
  vector[N] age;
  int<lower=1> n_age_fine;
  vector[n_age_fine] age_fine;
  int<lower=0> n_tested[N];
  int<lower=0> y[N, 4];
  int<lower=1> n_intervals;
  vector[n_intervals + 1] change_points;
}

parameters {
  vector<lower=0>[n_intervals] lambda1;
  vector<lower=0>[n_intervals] lambda2;
  real<lower=0> sigma12;
  real<lower=0> sigma21;
}

transformed parameters {
  vector[N] prob_s;
  vector[N] prob_x1;
  vector[N] prob_x2;
  vector[N] prob_x12;

  for (i in 1:N) {
    real t = age[i];
    prob_s[i]   = p_s(t, lambda1, lambda2, change_points, n_intervals);
    prob_x1[i]  = p_x1(t, lambda1, lambda2, sigma12, change_points, n_intervals);
    prob_x2[i]  = p_x2(t, lambda1, lambda2, sigma21, change_points, n_intervals);
    prob_x12[i] = p_x12(t, lambda1, lambda2, sigma12, sigma21, change_points, n_intervals);
  }
}

model {
  lambda1 ~ exponential(1);
  lambda2 ~ exponential(1);
  sigma12 ~ uniform(0.0, 20.0);
  sigma21 ~ uniform(0.0, 20.0);

  for (i in 1:N) {
    vector[4] probs;
probs[1] = prob_s[i];
probs[2] = prob_x1[i];
probs[3] = prob_x2[i];
probs[4] = prob_x12[i];

    y[i] ~ multinomial(probs);
  }
}

generated quantities {
  int y_rep[N, 4];
  vector[n_age_fine] p_s_grid;
  vector[n_age_fine] p_x1_grid;
  vector[n_age_fine] p_x2_grid;
  vector[n_age_fine] p_x12_grid;
  vector[n_age_fine] seroprev_grid;

  for (i in 1:N) {
   vector[4] probs;
probs[1] = prob_s[i];
probs[2] = prob_x1[i];
probs[3] = prob_x2[i];
probs[4] = prob_x12[i];

    y_rep[i] = multinomial_rng(probs, n_tested[i]);
  }

  for (i in 1:n_age_fine) {
    real t = age_fine[i];
    real s = p_s(t, lambda1, lambda2, change_points, n_intervals);
    real x1 = p_x1(t, lambda1, lambda2, sigma12, change_points, n_intervals);
    real x2 = p_x2(t, lambda1, lambda2, sigma21, change_points, n_intervals);
    real x12 = p_x12(t, lambda1, lambda2, sigma12, sigma21, change_points, n_intervals);

    p_s_grid[i] = s;
    p_x1_grid[i] = x1;
    p_x2_grid[i] = x2;
    p_x12_grid[i] = x12;
    seroprev_grid[i] = 1 - s;
  }
}
