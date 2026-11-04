## fit_presence_models.R
## Let the socio-political factors affect the PRESENCE of food waste as well as its amount.
##
##   F123_D   : F123 + Deprivation in the zero-inflation part (diagnostic check D, refit and saved)
##   F123_all : F123 + all three factors in the zero-inflation part
##
## Weibull (amount) part, both models: beta1 = Housing, beta2 = Deprivation, beta3 = Enfranchisement
## Zero-inflation part (probability of NO waste; negative = waste more likely present):
##   delta1 = service roads, delta2 = Deprivation,
##   F123_all only: delta3 = Housing, delta4 = Enfranchisement
##
## Compares both by LOO with F12, F13, F123 and the proxy. Saves each fit, then frees memory.
## Run from the proxy_to_process project root, in a fresh R session. About 10-15 minutes.

library(here)
library(sf)
library(rstan)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(here("data", "factor_analysis_16var.Rdata"))   # Lambda
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub

pick <- function(v) which.max(abs(Lambda[v, ]))
fac  <- function(k) dat_final_16var[[paste0("Factor", k)]]
Hous <- fac(pick("housing_highdensity"))
Dep  <- fac(pick("pct_rent_subsidized"))
Enf  <- fac(pick("employed_pct"))

base <- list(
  N = nrow(dat_no_ub),
  J = length(unique(dat_no_ub$site)),
  site_id = as.integer(factor(dat_no_ub$site)),
  X1 = Hous, X2 = Dep, X3 = Enf,
  P1 = (dat_no_ub$service_road_weight - mean(dat_no_ub$service_road_weight)) /
       sd(dat_no_ub$service_road_weight),
  fw_score_weighted = dat_no_ub$fw_score_weighted
)

pars_amount <- c("mu0", "tau_mu", "alpha", "beta1", "beta2", "beta3", "gamma", "delta1")

fit_save <- function(name, stanfile, data, pars) {
  cat("\n\n=====", name, "=====\n")
  fit <- stan(here("src", stanfile), data = data, chains = 4, iter = 2000, seed = 2026)
  assign(paste0("fit_", name), fit)
  save(list = paste0("fit_", name), file = here("data", paste0("fit_hier_ZIW_16var_", name, ".Rdata")))
  print(fit, pars = pars, probs = c(0.025, 0.5, 0.975))
  check_hmc_diagnostics(fit)
  ll <- extract_log_lik(fit, merge_chains = FALSE)
  rm(list = c("fit", paste0("fit_", name))); gc()
  lo <- loo(ll, r_eff = relative_eff(exp(ll)))
  print(lo)
  lo
}

new_loos <- list()
new_loos$F123_D <- fit_save("F123_D", "fw_hier_3pred_ZIW_2pred_bern.stan",
                            c(base, list(P2 = Dep)),
                            c(pars_amount, "delta2"))
new_loos$F123_all <- fit_save("F123_all", "fw_hier_3pred_ZIW_4pred_bern.stan",
                              c(base, list(P2 = Dep, P3 = Hous, P4 = Enf)),
                              c(pars_amount, "delta2", "delta3", "delta4"))

load(here("data", "loo_16var_models.Rdata"))   # loos: F12, F123, F13, proxy
loos_all <- c(loos, new_loos)
cat("\n\n===== LOO comparison, all models =====\n")
print(loo_compare(loos_all))
save(loos_all, file = here("data", "loo_16var_models_all.Rdata"))
