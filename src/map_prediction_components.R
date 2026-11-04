## map_prediction_components.R
## City-wide maps that show the socio-political gradient, which the full expected-value maps
## hide because service roads (the zero-inflation part) dominate them.
##
##   1. Amount where waste is present: the Weibull mean only, sigma * Gamma(1 + 1/alpha),
##      which depends only on the factor scores (plug-in posterior medians).
##   2. Full expected FW score averaged to 750 m cells (UTM 10N), which smooths out
##      individual lanes and shows the neighbourhood-scale pattern.
##   3. Full expected FW score with the colour ramp capped at the 95th percentile.
##
## Needs data/van_dat_predictions_16var.Rdata (from resume_after_fits.R or
## refit_factor_models.R) and the saved F12 / F123 fits.
## Run from the proxy_to_process project root.

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)

load(here("data", "van_dat_predictions_16var.Rdata"))  # van_dat_fa, k_housing, k_enfr, k_depriv
## Fits are loaded one at a time below, to keep memory use down
load_fit <- function(file) {
  env <- new.env()
  load(here("data", file), envir = env)
  get(ls(env)[1], envir = env)
}

city <- st_drop_geometry(van_dat_fa)

## Column order must match X1, X2, X3 in refit_factor_models.R
models <- list(
  F12  = list(file = "fit_hier_ZIW_16var_F12.Rdata",  cols = c(k_housing, k_depriv),
              title = "Housing + Deprivation"),
  F123 = list(file = "fit_hier_ZIW_16var_F123.Rdata", cols = c(k_housing, k_depriv, k_enfr),
              title = "Housing + Deprivation + Enfranchisement")
)
if (file.exists(here("data", "fit_hier_ZIW_16var_F13.Rdata")) && "pred_F13" %in% names(city)) {
  models$F13 <- list(file = "fit_hier_ZIW_16var_F13.Rdata", cols = c(k_housing, k_enfr),
                     title = "Housing + Enfranchisement")
}

map_theme <- theme_void() +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(size = 10))

save_map <- function(p, file) ggsave(here("figs", file), p, width = 12, height = 8)

## ---------------------------------------------------------------
## 1. Amount where waste is present (Weibull part only)
## ---------------------------------------------------------------
for (m in names(models)) {
  pars <- c("mu0", "alpha", paste0("beta", seq_along(models[[m]]$cols)))
  fit <- load_fit(models[[m]]$file)
  post <- rstan::extract(fit, pars = pars)
  rm(fit); gc()
  b <- c(median(post$beta1), median(post$beta2), if (!is.null(post$beta3)) median(post$beta3))
  X <- as.matrix(city[, paste0("Factor", models[[m]]$cols)])
  van_dat_fa[[paste0("amount_", m)]] <-
    as.numeric(exp(median(post$mu0) + X %*% b) * gamma(1 + 1 / median(post$alpha)))
}

cat("Quantiles of amount where present:\n")
print(signif(sapply(paste0("amount_", names(models)), function(v)
  quantile(van_dat_fa[[v]], c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1))), 3))

lims_amount <- quantile(c(van_dat_fa$amount_F12, van_dat_fa$amount_F123), c(0.01, 0.99))
for (m in names(models)) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[paste0("amount_", m)]]), color = NA) +
    scale_fill_viridis_c(option = "inferno", limits = lims_amount, oob = squish,
                         name = "Expected amount\nwhere present") +
    labs(title = paste("Amount of food waste where present:", models[[m]]$title),
         subtitle = "Weibull part only (service roads excluded); limits at 1st-99th percentile") +
    map_theme
  save_map(p, paste0("city_amount_", m, ".pdf"))
}

## ---------------------------------------------------------------
## 2. Full expected value averaged to 750 m cells
## ---------------------------------------------------------------
city_utm <- st_transform(van_dat_fa, 26910)
grid750 <- st_sf(geometry = st_make_grid(city_utm, cellsize = 750))
grid750$cell <- seq_len(nrow(grid750))

pts <- suppressWarnings(st_centroid(city_utm[, paste0("pred_", names(models))]))
pts <- st_join(pts, grid750, join = st_within)

pred_cols <- paste0("pred_", names(models))
cell_means <- st_drop_geometry(pts) |>
  group_by(cell) |>
  summarise(across(all_of(pred_cols), mean), n_sub = n(), .groups = "drop")
grid750 <- inner_join(grid750, cell_means, by = "cell")

lims_750 <- quantile(unlist(st_drop_geometry(grid750)[, pred_cols]), c(0.01, 0.99))
for (m in names(models)) {
  p <- ggplot(grid750) +
    geom_sf(aes(fill = .data[[paste0("pred_", m)]]), color = NA) +
    scale_fill_viridis_c(option = "inferno", limits = lims_750, oob = squish,
                         name = "Mean expected\nFW score") +
    labs(title = paste("Expected FW score, 750 m average:", models[[m]]$title),
         subtitle = "Full model (incl. service roads), averaged within 750 m cells; limits at 1st-99th percentile") +
    map_theme
  save_map(p, paste0("city_pred_", m, "_750m.pdf"))
}

## ---------------------------------------------------------------
## 3. Full expected value, ramp capped at the 95th percentile
## ---------------------------------------------------------------
lims_cap <- quantile(unlist(city[, pred_cols]), c(0.01, 0.95))
for (m in names(models)) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[paste0("pred_", m)]]), color = NA) +
    scale_fill_viridis_c(option = "inferno", limits = lims_cap, oob = squish,
                         name = "Expected\nFW score") +
    labs(title = paste("Expected FW score:", models[[m]]$title),
         subtitle = sprintf("Ramp capped at the 95th percentile (%.3g); higher values drawn in the top colour",
                            lims_cap[2])) +
    map_theme
  save_map(p, paste0("city_pred_", m, "_cap95.pdf"))
}
