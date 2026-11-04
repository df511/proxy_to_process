## map_laneways.R
## City-wide maps that keep the full range of values visible but give most of the colour
## ramp to the high-value cells (mostly lanes/service roads), so variation between
## 0.25 and 1 is easier to see.
##
## Colour allocation: the darkest 30% of the inferno palette covers 0 to 0.25;
## the remaining 70% covers 0.25 to 1; values above 1 take the top colour.
##
## Needs data/van_dat_predictions_16var.Rdata (from resume_after_fits.R).
## Run from the proxy_to_process project root.

library(here)
library(sf)
library(ggplot2)
library(scales)
library(viridisLite)

load(here("data", "van_dat_predictions_16var.Rdata"))  # van_dat_fa

split_at  <- 0.25   # value where the detailed part of the ramp starts
top       <- 1      # values above this take the top colour
low_share <- 0.30   # share of the palette used for 0 to split_at

titles <- c(pred_F13   = "Housing + Enfranchisement",
            pred_F12   = "Housing + Deprivation",
            pred_F123  = "Housing + Deprivation + Enfranchisement",
            pred_proxy = "Proxy model")
preds <- intersect(names(titles), names(van_dat_fa))

cat("Share of city cells with expected FW score >= ", split_at, ":\n", sep = "")
print(round(sapply(preds, function(v) mean(van_dat_fa[[v]] >= split_at)), 3))

## Palette positions: n_low colours spread over [0, split_at], the rest over [split_at, top]
pal   <- inferno(256)
n_low <- round(low_share * length(pal))
split_pos <- split_at / top
values <- c(seq(0, split_pos, length.out = n_low),
            seq(split_pos, 1, length.out = length(pal) - n_low + 1)[-1])

for (v in preds) {
  p <- ggplot(van_dat_fa) +
    geom_sf(aes(fill = .data[[v]]), color = NA) +
    scale_fill_gradientn(
      colours = pal,
      values = values,
      limits = c(0, top),
      oob = squish,
      breaks = c(0, 0.1, 0.25, 0.4, 0.6, 0.8, 1),
      labels = c("0", "0.1", "0.25", "0.4", "0.6", "0.8", "≥ 1"),
      name = "Expected\nFW score",
      guide = guide_colourbar(barheight = unit(12, "lines"))
    ) +
    labs(title = paste("Expected FW score:", titles[[v]]),
         subtitle = sprintf("Darkest %d%% of the colour ramp covers 0 to %.2g; the rest covers %.2g to %.2g; values above %.2g in the top colour",
                            round(100 * low_share), split_at, split_at, top, top)) +
    theme_void() +
    theme(plot.title = element_text(face = "bold"))
  ggsave(here("figs", paste0("city_", v, "_split_ramp.pdf")), p, width = 12, height = 8)
}
