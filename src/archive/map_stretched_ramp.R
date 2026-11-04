## map_stretched_ramp.R
## City-wide expected FW score maps over the FULL range of values (nothing capped or hidden),
## with the colour ramp stretched so variation between 0.25 and 1 is easy to see:
##   minimum to 0.25 : 20% of the colour ramp
##   0.25 to 1       : 70%
##   1 to maximum    : 10%
## All models, including the proxy, share the same scale so the maps compare directly.
## The proxy is also mapped on its own scale, since all its values fall below 0.25.
##
## Needs data/van_dat_predictions_16var.Rdata (from resume_after_fits.R).
## Run from the proxy_to_process project root.

library(here)
library(sf)
library(ggplot2)
library(scales)
library(viridisLite)

load(here("data", "van_dat_predictions_16var.Rdata"))  # van_dat_fa

## ---- Colour stretch: data breakpoints and the share of the ramp they get ----
mid_lo <- 0.25
mid_hi <- 1
ramp_at <- c(0, 0.20, 0.90, 1)   # ramp position at min, mid_lo, mid_hi, max

titles <- c(pred_F13  = "Housing + Enfranchisement",
            pred_F12  = "Housing + Deprivation",
            pred_F123 = "Housing + Deprivation + Enfranchisement",
            pred_proxy = "Proxy model (distance to food retail)")
preds <- intersect(names(titles), names(van_dat_fa))

rng <- range(unlist(lapply(preds, function(v) van_dat_fa[[v]])), na.rm = TRUE)
cat("Shared range across models:", signif(rng, 3), "\n")
stopifnot(rng[1] < mid_lo, rng[2] > mid_hi)

## Place 101 inferno colours along the ramp, then convert each ramp position back to a data
## value with the piecewise-linear stretch; gradientn needs those as positions in [0, 1].
n_col <- 101
ramp_pos <- seq(0, 1, length.out = n_col)
data_at <- approx(x = ramp_at, y = c(rng[1], mid_lo, mid_hi, rng[2]), xout = ramp_pos)$y
col_values <- rescale(data_at, from = rng)

brks <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
brks <- c(brks[brks > rng[1]], rng[2])
labs_brks <- c(format(brks[-length(brks)], drop0trailing = TRUE),
               paste0(signif(rng[2], 3), " (max)"))

for (v in preds) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[v]]), color = NA) +
    scale_fill_gradientn(
      colours = inferno(n_col),
      values = col_values,
      limits = rng,
      breaks = brks,
      labels = labs_brks,
      name = "Expected\nFW score"
    ) +
    guides(fill = guide_colourbar(barheight = unit(12, "lines"))) +
    labs(title = paste("Expected FW score:", titles[[v]]),
         subtitle = "Full range shown; colour ramp stretched so 0.25–1 takes 70% of the colours") +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))
  ggsave(here("figs", paste0("city_", v, "_stretched.pdf")), p, width = 12, height = 8)
}

cat("\nShare of city cells in each band:\n")
print(round(sapply(preds, function(v) {
  x <- van_dat_fa[[v]]
  c(below_0.25 = mean(x < mid_lo), `0.25_to_1` = mean(x >= mid_lo & x <= mid_hi), above_1 = mean(x > mid_hi))
}), 3))

## ---- Proxy on its own scale (all its values are below 0.25) ----
if ("pred_proxy" %in% names(van_dat_fa)) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = pred_proxy), color = NA) +
    scale_fill_viridis_c(option = "inferno", name = "Expected\nFW score") +
    guides(fill = guide_colourbar(barheight = unit(12, "lines"))) +
    labs(title = "Expected FW score: Proxy model (distance to food retail)",
         subtitle = sprintf("Own colour scale: %.3g to %.3g", min(van_dat_fa$pred_proxy),
                            max(van_dat_fa$pred_proxy))) +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))
  ggsave(here("figs", "city_pred_proxy_ownscale.pdf"), p, width = 12, height = 8)
}
