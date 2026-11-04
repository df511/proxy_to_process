## remap_predictions.R
## Redraw the city-wide prediction maps from data/van_dat_predictions_16var.Rdata
## (written by refit_factor_models.R) with a more readable colour scale. No refitting.
##
## Colour scale: log10, with limits at the 1st and 99th percentiles of the process-model
## predictions; cells beyond the limits are drawn in the end colours ("squished").
## F12 and F123 share the same limits so they can be compared directly.
##
## Also prints a unit check on the distance variables (the proxy map came out flat).
## Run from the proxy_to_process project root.

library(here)
library(sf)
library(ggplot2)
library(scales)

load(here("data", "van_dat_predictions_16var.Rdata"))  # van_dat_fa, k_housing, k_enfr, k_depriv
load(here("data", "dat_no_ub_750.Rdata"))              # dat_no_ub (transect cells)

## ---------------------------------------------------------------
## 1. Unit check: are the distances on the same scale in the transect and city data?
## ---------------------------------------------------------------
cat("nearest_food_retail, transects:\n");  print(summary(dat_no_ub$nearest_food_retail))
cat("nearest_food_retail, city:\n");       print(summary(van_dat_fa$nearest_food_retail))
cat("nearest_shelter_dist, transects:\n"); print(summary(dat_no_ub$nearest_shelter_dist))
cat("nearest_shelter_dist, city:\n");      print(summary(van_dat_fa$nearest_shelter_dist))
cat("nearest_food_retail_z, city (transect standardization):\n")
print(summary(van_dat_fa$nearest_food_retail_z))

## ---------------------------------------------------------------
## 2. Distribution of the predictions
## ---------------------------------------------------------------
preds <- intersect(c("pred_F12", "pred_F123", "pred_F13", "pred_proxy"), names(van_dat_fa))
probs <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
cat("\nQuantiles of predicted expected FW score:\n")
print(signif(sapply(preds, function(v) quantile(van_dat_fa[[v]], probs, na.rm = TRUE)), 3))

## ---------------------------------------------------------------
## 3. Maps
## ---------------------------------------------------------------
map_pred <- function(var, title, file, lims) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    scale_fill_viridis_c(
      option = "inferno",
      trans = "log10",
      limits = lims,
      oob = squish,
      labels = label_number(drop0trailing = TRUE),
      name = "Expected\nFW score"
    ) +
    labs(title = title,
         subtitle = sprintf("Colour limits: 1st-99th percentile (%.3g to %.3g); values beyond are squished",
                            lims[1], lims[2])) +
    theme_void()
  ggsave(here("figs", file), p, width = 12, height = 8)
}

process_preds <- setdiff(preds, "pred_proxy")
lims_process <- quantile(unlist(lapply(process_preds, function(v) van_dat_fa[[v]])),
                         c(0.01, 0.99), na.rm = TRUE)

titles <- c(pred_F12 = "Housing + Deprivation",
            pred_F123 = "Housing + Deprivation + Enfranchisement",
            pred_F13 = "Housing + Enfranchisement")
for (v in process_preds) {
  map_pred(v, paste("Expected FW score:", titles[[v]]),
           paste0("city_", v, "_log10_p01p99.pdf"), lims_process)
}

lims_proxy <- quantile(van_dat_fa$pred_proxy, c(0.01, 0.99), na.rm = TRUE)
if (diff(lims_proxy) > 0) {
  map_pred("pred_proxy", "Expected FW score: proxy model",
           "city_pred_proxy_log10_p01p99.pdf", lims_proxy)
} else {
  cat("\nProxy predictions are constant across the 1st-99th percentile; map skipped.\n")
}
