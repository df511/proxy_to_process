## refit_factor_models.R
## Refit the process-based model with the corrected 3-factor, 16-variable factor scores:
##   F12  : Housing-tenure-density + Deprivation
##   F123 : Housing-tenure-density + Deprivation + Enfranchisement
##   F13  : Housing-tenure-density + Enfranchisement (the manuscript's original pairing;
##          optional, off by default: set fit_F13 <- TRUE to include it)
## Compare against the proxy model with LOO, and remake the city-wide prediction maps.
##
## Run AFTER src/factor_analysis_revision.R, from the proxy_to_process project root.
## Requires: src/fw_hier_2pred_ZIW_1pred_bern.stan and src/fw_hier_3pred_ZIW_1pred_bern.stan
## The proxy (Tweedie) fit does not use the factors, so it is NOT refit; only its city map
## is redone with the transect standardization of distance to food retail.

library(here)
library(dplyr)
library(sf)
library(rstan)
library(loo)
library(ggplot2)
options(mc.cores = parallel::detectCores())

fit_F13 <- FALSE   # set TRUE to also refit the original Factor 1 + Factor 3 pairing
rstan_options(auto_write = TRUE)

## ---------------------------------------------------------------
## 1. Training data and corrected factor scores
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # fa_result, Lambda, W, train_means, train_sds, fa_vars
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var (unscaled + corrected Factor1-3)
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub (same rows/order)
stopifnot(identical(dat_final_16var$fw_score_weighted, dat_no_ub$fw_score_weighted))

## Identify factors by content rather than by number (factanal orders by variance explained)
pick <- function(v) which.max(abs(Lambda[v, ]))
k_housing <- pick("housing_highdensity")
k_enfr    <- pick("employed_pct")
k_depriv  <- pick("pct_rent_subsidized")
stopifnot(length(unique(c(k_housing, k_enfr, k_depriv))) == 3)
cat("Housing-tenure-density = Factor", k_housing,
    "| Enfranchisement = Factor", k_enfr,
    "| Deprivation = Factor", k_depriv, "\n")
print(round(Lambda[, c(k_housing, k_enfr, k_depriv)], 2))

## Factor signs are arbitrary. Check these read the way the manuscript describes,
## e.g. high Housing-tenure-density = affluent / low density (income positive, density negative).
cat("\nCorrelation of each factor with key raw variables (transects):\n")
print(round(cor(st_drop_geometry(dat_no_ub)[, c("household_income", "housing_highdensity",
                                                 "rent_pct", "employed_pct", "pct_rent_subsidized")],
                st_drop_geometry(dat_final_16var)[, paste0("Factor", c(k_housing, k_enfr, k_depriv))]), 2))

sr_mean <- mean(dat_no_ub$service_road_weight)
sr_sd   <- sd(dat_no_ub$service_road_weight)

fac <- function(k) dat_final_16var[[paste0("Factor", k)]]

stan_base <- list(
  N = nrow(dat_no_ub),
  J = length(unique(dat_no_ub$site)),
  site_id = as.integer(factor(dat_no_ub$site)),   # same coding as dat_final$site_id
  P1 = (dat_no_ub$service_road_weight - sr_mean) / sr_sd,
  fw_score_weighted = dat_no_ub$fw_score_weighted
)
## beta1 = Housing-tenure-density, beta2 = Deprivation, beta3 = Enfranchisement
stan_F12  <- c(stan_base, list(X1 = fac(k_housing), X2 = fac(k_depriv)))
stan_F123 <- c(stan_base, list(X1 = fac(k_housing), X2 = fac(k_depriv), X3 = fac(k_enfr)))
## F13 keeps the original manuscript coding: beta1 = Housing, beta2 = Enfranchisement
stan_F13  <- c(stan_base, list(X1 = fac(k_housing), X2 = fac(k_enfr)))

## ---------------------------------------------------------------
## 2. Fit the models
## ---------------------------------------------------------------
## Each fit is saved as soon as it finishes, so a crash later doesn't lose it.
## To skip a fit that already finished, load() its .Rdata file and comment out the stan() call.
fits <- list()

fits$F12 <- stan(here("src", "fw_hier_2pred_ZIW_1pred_bern.stan"), data = stan_F12,
                 chains = 4, iter = 2000, seed = 2026)
fit_F12 <- fits$F12
save(fit_F12, file = here("data", "fit_hier_ZIW_16var_F12.Rdata"))

fits$F123 <- stan(here("src", "fw_hier_3pred_ZIW_1pred_bern.stan"), data = stan_F123,
                  chains = 4, iter = 2000, seed = 2026)
fit_F123 <- fits$F123
save(fit_F123, file = here("data", "fit_hier_ZIW_16var_F123.Rdata"))

if (fit_F13) {
  fits$F13 <- stan(here("src", "fw_hier_2pred_ZIW_1pred_bern.stan"), data = stan_F13,
                   chains = 4, iter = 2000, seed = 2026)
  fit_F13_obj <- fits$F13
  save(fit_F13_obj, file = here("data", "fit_hier_ZIW_16var_F13.Rdata"))
}

pars <- c("mu0", "tau_mu", "alpha", "beta1", "beta2", "gamma", "delta1")
for (m in names(fits)) {
  cat("\n=====", m, "=====\n")
  print(fits[[m]], pars = if (m == "F123") c(pars, "beta3") else pars,
        probs = c(0.025, 0.5, 0.975))
  check_hmc_diagnostics(fits[[m]])   # want: all Rhat < 1.01, no divergences
}

## ---------------------------------------------------------------
## 3. LOO: process models vs proxy (same rows and response)
## ---------------------------------------------------------------
load(here("data", "fit_hier_1pred_tweedie.Rdata"))
loo_of <- function(fit) {
  ll <- extract_log_lik(fit, merge_chains = FALSE)
  loo(ll, r_eff = relative_eff(exp(ll)))
}
loos <- lapply(fits, loo_of)
loos$proxy <- loo_of(fit_hier_1pred_tweedie)
for (m in names(loos)) { cat("\n=====", m, "=====\n"); print(loos[[m]]) }   # check Pareto k
print(loo_compare(loos))
save(loos, file = here("data", "loo_16var_models.Rdata"))

## ---------------------------------------------------------------
## 4. City-wide predictions
##    Posterior median of the expected value per cell; transect intercepts set to mu0.
##    Computed in chunks of cells so memory stays bounded.
## ---------------------------------------------------------------
load(here("data", "van_dat_final_750_16var.Rdata"))  # van_dat_fa

chunked_median <- function(n_cells, f, chunk = 5000) {
  out <- numeric(n_cells)
  for (start in seq(1, n_cells, by = chunk)) {
    idx <- start:min(start + chunk - 1, n_cells)
    out[idx] <- apply(f(idx), 1, median)   # f returns cells x draws
  }
  out
}

city_expected_ziw <- function(fit, X, P1) {
  post <- rstan::extract(fit)
  B <- cbind(post$beta1, post$beta2, if (!is.null(post$beta3)) post$beta3)
  shape_term <- gamma(1 + 1 / post$alpha)               # one value per draw
  chunked_median(nrow(X), function(idx) {
    log_sigma <- sweep(X[idx, , drop = FALSE] %*% t(B), 2, post$mu0, "+")
    pi0 <- plogis(outer(P1[idx], post$delta1) + rep(post$gamma, each = length(idx)))
    (1 - pi0) * sweep(exp(log_sigma), 2, shape_term, "*")
  })
}

city <- st_drop_geometry(van_dat_fa)
## Column order must match X1, X2, X3 in the Stan data above
x_cols <- list(F12  = c(k_housing, k_depriv),
               F123 = c(k_housing, k_depriv, k_enfr),
               F13  = c(k_housing, k_enfr))
for (m in names(fits)) {
  X <- as.matrix(city[, paste0("Factor", x_cols[[m]])])
  van_dat_fa[[paste0("pred_", m)]] <- city_expected_ziw(fits[[m]], X, city$service_road_weight_z)
}

## Proxy model, mapped with the transect standardization of distance to food retail
post_tw <- rstan::extract(fit_hier_1pred_tweedie)
van_dat_fa$pred_proxy <- chunked_median(nrow(city), function(idx)
  exp(sweep(outer(city$nearest_food_retail_z[idx], post_tw$beta1), 2, post_tw$mu0, "+")))

save(van_dat_fa, k_housing, k_enfr, k_depriv,
     file = here("data", "van_dat_predictions_16var.Rdata"))
print(summary(st_drop_geometry(van_dat_fa)[, c(paste0("pred_", names(fits)), "pred_proxy")]))

## ---------------------------------------------------------------
## 5. Maps for review: each factor and each model's predictions
## ---------------------------------------------------------------
map_var <- function(var, file, title, log10 = FALSE) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    scale_fill_viridis_c(option = "inferno", trans = if (log10) "log10" else "identity",
                         name = NULL) +
    labs(title = title) +
    theme_void()
  ggsave(here("figs", file), p, width = 12, height = 8)
}
map_var(paste0("Factor", k_housing), "city_factor_housing_16var.pdf",        "Housing-tenure-density factor")
map_var(paste0("Factor", k_enfr),    "city_factor_enfranchisement_16var.pdf", "Enfranchisement factor")
map_var(paste0("Factor", k_depriv),  "city_factor_deprivation_16var.pdf",     "Deprivation factor")
model_titles <- c(F12 = "Housing + Deprivation", F123 = "Housing + Deprivation + Enfranchisement",
                  F13 = "Housing + Enfranchisement")
for (m in names(fits)) {
  map_var(paste0("pred_", m), paste0("city_pred_", m, "_log10.pdf"),
          paste("Expected FW score:", model_titles[[m]]), log10 = TRUE)
}
map_var("pred_proxy", "city_pred_proxy_log10.pdf", "Expected FW score: proxy model", log10 = TRUE)
