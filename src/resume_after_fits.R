## resume_after_fits.R
## Finishes refit_factor_models.R from the saved fits, WITHOUT refitting:
## LOO comparison, city-wide predictions, and the review maps.
## Memory-light: loads one fit at a time, keeps only what it needs, then frees it.
##
## Needs: data/fit_hier_ZIW_16var_F12.Rdata, data/fit_hier_ZIW_16var_F123.Rdata,
##        data/fit_hier_1pred_tweedie.Rdata, data/factor_analysis_16var.Rdata,
##        data/van_dat_final_750_16var.Rdata (from factor_analysis_revision.R)
## Run from the proxy_to_process project root, in a fresh R session.

library(here)
library(sf)
library(rstan)
library(loo)
library(ggplot2)

## ---------------------------------------------------------------
## 0. Which factor is which (same rule as refit_factor_models.R)
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # Lambda, ...
pick <- function(v) which.max(abs(Lambda[v, ]))
k_housing <- pick("housing_highdensity")
k_enfr    <- pick("employed_pct")
k_depriv  <- pick("pct_rent_subsidized")
cat("Housing = Factor", k_housing, "| Deprivation = Factor", k_depriv,
    "| Enfranchisement = Factor", k_enfr, "\n")

## Column order must match X1, X2, X3 used in the fits
models <- list(
  F12  = list(file = "fit_hier_ZIW_16var_F12.Rdata",  cols = c(k_housing, k_depriv)),
  F123 = list(file = "fit_hier_ZIW_16var_F123.Rdata", cols = c(k_housing, k_depriv, k_enfr))
)
## Original pairing (Housing + Enfranchisement), included if fit_F13.R has been run
if (file.exists(here("data", "fit_hier_ZIW_16var_F13.Rdata"))) {
  models$F13 <- list(file = "fit_hier_ZIW_16var_F13.Rdata", cols = c(k_housing, k_enfr))
}
cat("Models:", paste(names(models), collapse = ", "), "\n")

## Load one saved fit and return just the object (frees the rest when it goes out of scope)
load_fit <- function(file) {
  env <- new.env()
  load(here("data", file), envir = env)
  get(ls(env)[1], envir = env)
}

## ---------------------------------------------------------------
## 1. LOO, one model at a time
## ---------------------------------------------------------------
loo_from_file <- function(file) {
  fit <- load_fit(file)
  ll <- extract_log_lik(fit, merge_chains = FALSE)
  rm(fit); gc()
  out <- loo(ll, r_eff = relative_eff(exp(ll)))
  rm(ll); gc()
  out
}

loos <- list()
for (m in names(models)) {
  cat("\nLOO:", m, "\n")
  loos[[m]] <- loo_from_file(models[[m]]$file)
  print(loos[[m]])
}
cat("\nLOO: proxy\n")
loos$proxy <- loo_from_file("fit_hier_1pred_tweedie.Rdata")
print(loos$proxy)

cat("\nLOO comparison:\n")
print(loo_compare(loos))
save(loos, file = here("data", "loo_16var_models.Rdata"))

## ---------------------------------------------------------------
## 2. Posterior draws needed for prediction (small), one model at a time
## ---------------------------------------------------------------
draws <- list()
for (m in names(models)) {
  fit <- load_fit(models[[m]]$file)
  pars <- c("mu0", "alpha", "gamma", "delta1", paste0("beta", seq_along(models[[m]]$cols)))
  draws[[m]] <- rstan::extract(fit, pars = pars)
  rm(fit); gc()
}
fit_tw <- load_fit("fit_hier_1pred_tweedie.Rdata")
draws$proxy <- rstan::extract(fit_tw, pars = c("mu0", "beta1"))
rm(fit_tw); gc()

## ---------------------------------------------------------------
## 3. City-wide predictions: posterior median of the expected value per cell,
##    transect intercepts set to mu0; computed in chunks of cells
## ---------------------------------------------------------------
load(here("data", "van_dat_final_750_16var.Rdata"))  # van_dat_fa
city <- st_drop_geometry(van_dat_fa)

chunked_median <- function(n_cells, f, chunk = 2000) {
  out <- numeric(n_cells)
  for (start in seq(1, n_cells, by = chunk)) {
    idx <- start:min(start + chunk - 1, n_cells)
    out[idx] <- apply(f(idx), 1, median)   # f returns cells x draws
  }
  out
}

for (m in names(models)) {
  post <- draws[[m]]
  B <- do.call(cbind, post[paste0("beta", seq_along(models[[m]]$cols))])
  X <- as.matrix(city[, paste0("Factor", models[[m]]$cols)])
  P1 <- city$service_road_weight_z
  shape_term <- gamma(1 + 1 / post$alpha)
  cat("Predicting", m, "...\n")
  van_dat_fa[[paste0("pred_", m)]] <- chunked_median(nrow(X), function(idx) {
    log_sigma <- sweep(X[idx, , drop = FALSE] %*% t(B), 2, post$mu0, "+")
    pi0 <- plogis(outer(P1[idx], post$delta1) + rep(post$gamma, each = length(idx)))
    (1 - pi0) * sweep(exp(log_sigma), 2, shape_term, "*")
  })
  gc()
}

cat("Predicting proxy ...\n")
van_dat_fa$pred_proxy <- chunked_median(nrow(city), function(idx)
  exp(sweep(outer(city$nearest_food_retail_z[idx], draws$proxy$beta1), 2, draws$proxy$mu0, "+")))

save(van_dat_fa, k_housing, k_enfr, k_depriv,
     file = here("data", "van_dat_predictions_16var.Rdata"))

cat("\nQuantiles of predicted expected FW score:\n")
print(signif(sapply(c(paste0("pred_", names(models)), "pred_proxy"), function(v)
  quantile(van_dat_fa[[v]], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))), 3))

## ---------------------------------------------------------------
## 4. Review maps: factors and predictions
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
map_var(paste0("Factor", k_depriv),  "city_factor_deprivation_16var.pdf",     "Deprivation factor")
map_var(paste0("Factor", k_enfr),    "city_factor_enfranchisement_16var.pdf", "Enfranchisement factor")
model_titles <- c(F12 = "Housing + Deprivation",
                  F123 = "Housing + Deprivation + Enfranchisement",
                  F13 = "Housing + Enfranchisement")
for (m in names(models)) {
  map_var(paste0("pred_", m), paste0("city_pred_", m, "_log10.pdf"),
          paste("Expected FW score:", model_titles[[m]]), log10 = TRUE)
}
map_var("pred_proxy", "city_pred_proxy_log10.pdf", "Expected FW score: proxy model", log10 = TRUE)
cat("\nDone.\n")
