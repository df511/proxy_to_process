
// stan code for a hierarchical model (predicts per site)
// with 1 predictor of a Tweedie distribution

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=J> site_id[N];

  vector[N] X1;              // single predictor
  real<lower=0> fw_score_weighted[N];      // outcome (must be >= 0)
}

parameters {
  // hierarchical intercept
  real mu0;
  real<lower=0> tau_mu;
  vector[J] mu_raw;

  // slope
  real beta1;

  // Tweedie parameters
  real<lower=0> phi;         // dispersion
  real<lower=1, upper=2> p;  // Tweedie power

  // latent Poisson rate (for compound representation)
  vector<lower=0>[N] lambda;
}

transformed parameters {
  vector[J] mu_site;
  vector[N] eta;
  vector<lower=0>[N] mu;

  mu_site = mu0 + tau_mu * mu_raw;

  for (n in 1:N) {
    eta[n] = mu_site[site_id[n]] + beta1 * X1[n];
    mu[n] = exp(eta[n]);   // mean on positive scale
  }
}

model {
  // priors
  mu0 ~ normal(0, 1);
  tau_mu ~ normal(0, 1);
  mu_raw ~ normal(0, 1);

  beta1 ~ normal(0, 1);

  phi ~ normal(1, 0.5);
  p ~ uniform(1, 2);

  // IMPORTANT: compound Poisson construction
  for (n in 1:N) {

    // latent Poisson count
    lambda[n] ~ gamma(phi, phi / mu[n]);

    if (fw_score_weighted[n] == 0) {
      target += poisson_lpmf(0 | lambda[n]);
    } else {
      target += poisson_lpmf(1 | lambda[n]) +
                gamma_lpdf(fw_score_weighted[n] | 1 / (p - 1),
                                 (p - 1) / (phi * pow(mu[n], 2 - p)));
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] lambda_rep;

  for (n in 1:N) {

    real lp;

    lambda_rep[n] = poisson_rng(gamma_rng(phi, phi / mu[n]));

    if (fw_score_weighted[n] == 0) {

      lp = poisson_lpmf(0 | lambda[n]);

      y_rep[n] = 0;

    } else {

      lp =
        poisson_lpmf(1 | lambda[n]) +
        gamma_lpdf(fw_score_weighted[n] |
                   1 / (p - 1),
                   (p - 1) / (phi * pow(mu[n], 2 - p)));

      y_rep[n] = gamma_rng(
        1 / (p - 1),
        (p - 1) / (phi * pow(mu[n], 2 - p))
      );
    }

    log_lik[n] = lp;
  }
}