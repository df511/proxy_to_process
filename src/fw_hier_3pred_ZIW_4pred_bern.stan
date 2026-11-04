// Hierarchical Zero-Inflated Weibull Model: three Weibull predictors, four zero-inflation predictors
// Identical to fw_hier_3pred_ZIW_1pred_bern.stan, plus P2-P4 / delta2-delta4 in the zero-inflation part.

data {
  int<lower=1> N;
  int<lower=1> J;                     // Number of sites
  int<lower=1,upper=J> site_id[N];    // Site ID for each grid cell
  vector[N] X1;
  vector[N] X2;
  vector[N] X3;
  vector[N] P1;
  vector[N] P2;
  vector[N] P3;
  vector[N] P4;
  real<lower=0> fw_score_weighted[N];
}

parameters {
  real<lower=0> alpha;           // Weibull shape
  real mu0;                      // Global mean of site intercepts
  real<lower=0> tau_mu;          // Std dev of site-level intercepts
  vector[J] mu_raw;              // Raw intercepts for sites
  real beta1;
  real beta2;
  real beta3;

  real gamma;                    // Intercept for zero-inflation
  real delta1;
  real delta2;
  real delta3;
  real delta4;
}

transformed parameters {
  vector[J] mu = mu0 + tau_mu * mu_raw;  // Site-specific intercepts
  vector[N] log_scale;
  vector<lower=0>[N] sigma;
  vector<lower=0, upper=1>[N] pi;

  for (n in 1:N) {
    log_scale[n] = mu[site_id[n]] + beta1 * X1[n] + beta2 * X2[n] + beta3 * X3[n];
    sigma[n] = exp(log_scale[n]);
    pi[n] = inv_logit(gamma + delta1 * P1[n] + delta2 * P2[n] + delta3 * P3[n] + delta4 * P4[n]);
  }
}

model {
  // Hyperpriors
  mu0 ~ normal(0, 1);
  tau_mu ~ normal(0, 1);
  mu_raw ~ normal(0, 1);
  // Priors
  alpha ~ normal(0.5, 0.5);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);

  gamma ~ normal(1, 1);
  delta1 ~ normal(0, 1);
  delta2 ~ normal(0, 1);
  delta3 ~ normal(0, 1);
  delta4 ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    if (fw_score_weighted[n] == 0) {
      target += bernoulli_lpmf(1 | pi[n]);
    } else {
      target += bernoulli_lpmf(0 | pi[n]) +
                weibull_lpdf(fw_score_weighted[n] | alpha, sigma[n]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] fw_score_pred;

  for (n in 1:N) {

    // Pointwise log-likelihood for LOO/WAIC
    if (fw_score_weighted[n] == 0) {
      log_lik[n] = bernoulli_lpmf(1 | pi[n]);
    } else {
      log_lik[n] =
        bernoulli_lpmf(0 | pi[n]) +
        weibull_lpdf(fw_score_weighted[n] | alpha, sigma[n]);
    }

    // Posterior predictive draw
    if (bernoulli_rng(pi[n]) == 1) {
      fw_score_pred[n] = 0;
    } else {
      fw_score_pred[n] = weibull_rng(alpha, sigma[n]);
    }
  }
}
