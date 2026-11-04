## fit_F13.R
## Fit the original manuscript pairing on the corrected scores:
##   F13: Housing-tenure-density + Enfranchisement
##   beta1 = Housing-tenure-density, beta2 = Enfranchisement (same coding as the original model)
## Same data, priors and settings as refit_factor_models.R. Saves the fit and frees memory.
##
## Run from the proxy_to_process project root, in a fresh R session, before
## resume_after_fits.R (which picks the F13 fit up automatically).

library(here)
library(sf)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(here("data", "factor_analysis_16var.Rdata"))   # Lambda, ...
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var (corrected Factor1-3)
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub

pick <- function(v) which.max(abs(Lambda[v, ]))
k_housing <- pick("housing_highdensity")
k_enfr    <- pick("employed_pct")
cat("Housing = Factor", k_housing, "| Enfranchisement = Factor", k_enfr, "\n")

stan_F13 <- list(
  N = nrow(dat_no_ub),
  J = length(unique(dat_no_ub$site)),
  site_id = as.integer(factor(dat_no_ub$site)),
  X1 = dat_final_16var[[paste0("Factor", k_housing)]],
  X2 = dat_final_16var[[paste0("Factor", k_enfr)]],
  P1 = (dat_no_ub$service_road_weight - mean(dat_no_ub$service_road_weight)) /
       sd(dat_no_ub$service_road_weight),
  fw_score_weighted = dat_no_ub$fw_score_weighted
)

fit_F13 <- stan(here("src", "fw_hier_2pred_ZIW_1pred_bern.stan"), data = stan_F13,
                chains = 4, iter = 2000, seed = 2026)
save(fit_F13, file = here("data", "fit_hier_ZIW_16var_F13.Rdata"))

print(fit_F13, pars = c("mu0", "tau_mu", "alpha", "beta1", "beta2", "gamma", "delta1"),
      probs = c(0.025, 0.5, 0.975))
check_hmc_diagnostics(fit_F13)

rm(fit_F13); gc()
