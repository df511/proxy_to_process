## final_figures_F123_all.R
## Figures 3 and 4 for the main model, F123_all (all three factors in both model parts),
## with the same layout, sizes and styling as manuscript_main_figs.Rmd. Writes NEW files;
## the existing Figure3.pdf / Figure4.pdf are not touched.
##
##   figs/Figure3_F123_all.pdf  (a) posterior coefficients, (b) expected FW score map
##   figs/Figure4_F123_all.pdf  (a) process model, (b) proxy on the same scale,
##                              (c) proxy on its own scale, (d) log10 ratio process/proxy
##   figs/FigureS_F123_all_750m.pdf  expected FW score averaged to 750 m cells, dropping
##                              edge cells less than half covered by the city grid
##
## Maps use posterior-median predictions from data/van_dat_predictions_16var.Rdata
## (pred_F123_all from predict_F123_all.R; pred_proxy from resume_after_fits.R).
## Expected-FW-score maps share one "stretched" colour ramp: minimum-0.25 gets 20% of the
## colours, 0.25-1 gets 70%, 1-maximum gets 10%; nothing is capped.
## Run from the proxy_to_process project root, in a fresh R session. A few minutes.

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(scales)
library(viridisLite)
library(cowplot)
library(scico)
library(rstan)

## ---------------------------------------------------------------
## Inputs
## ---------------------------------------------------------------
env <- new.env()
load(here("data", "fit_hier_ZIW_16var_F123_all.Rdata"), envir = env)
post <- rstan::extract(get(ls(env)[1], envir = env),
                       pars = c("beta1", "beta2", "beta3", "delta1", "delta2", "delta3", "delta4"))
rm(env); gc()

load(here("data", "van_dat_predictions_16var.Rdata"))   # van_dat_fa
stopifnot(all(c("pred_F123_all", "pred_proxy") %in% names(van_dat_fa)))

## Cells in masked DAs (domain_mask.R: population density below the threshold) are removed
mask_col <- "grey85"
if ("in_domain" %in% names(van_dat_fa)) {
  ## Masked DAs are removed from the maps entirely (drawn as empty space)
  cat("Removing", sum(!van_dat_fa$in_domain), "masked cells from the maps\n")
  van_dat_fa <- van_dat_fa[van_dat_fa$in_domain, ]
} else {
  cat("No in_domain column found: run src/domain_mask.R to mask parkland/unpopulated cells.\n")
}

## ---------------------------------------------------------------
## Shared stretched colour scale for expected FW score
## ---------------------------------------------------------------
rng <- range(c(van_dat_fa$pred_F123_all, van_dat_fa$pred_proxy), na.rm = TRUE)
data_at <- approx(c(0, 0.20, 0.90, 1), c(rng[1], 0.25, 1, rng[2]),
                  xout = seq(0, 1, length.out = 101))$y
brks <- c(c(0.1, 0.25, 0.5, 0.75, 1)[c(0.1, 0.25, 0.5, 0.75, 1) > rng[1]], rng[2])
stretched_fill <- function(name = "Expected FW Score") {
  scale_fill_gradientn(
    colours = inferno(101), values = rescale(data_at, from = rng), limits = rng,
    na.value = mask_col,
    breaks = brks,
    labels = c(format(brks[-length(brks)], drop0trailing = TRUE), signif(rng[2], 3)),
    name = name
  )
}
cat("Shared colour range:", signif(rng, 3), "\n")

## =========================================================
## Figure 3
## =========================================================

## (a) Ridge plot: quantity effects (purple) and presence effects (orange).
##
## Signs are re-expressed for readability; no refit is needed, because the model is linear in
## each factor and every coefficient prior is symmetric around 0, so reversing a scale gives the
## same model with the coefficient's sign reversed.
##   * Presence: the model's zero-inflation part is P(NO waste); the log-odds of presence are
##     minus those of absence, so every delta is multiplied by -1.
##   * Housing-tenure-density is reoriented so HIGH = denser, more renter-occupied
##     (the factor as fitted has high = affluent, low-density): multiply its effects by -1.
##   * Enfranchisement is reoriented as DISenfranchisement, HIGH = less employment:
##     multiply its effects by -1.
##   * Deprivation keeps its orientation (HIGH = more deprived).
## Net multipliers on the fitted parameters:
##   beta1 (Housing, quantity)   -1     delta1 (roads, presence)              -1
##   beta2 (Deprivation, quant.) +1     delta2 (Deprivation, presence)        -1
##   beta3 (Enfranch., quantity) -1     delta3 (Housing, presence)  -1 x -1 = +1
##                                      delta4 (Enfranch., presence) -1 x -1 = +1
## Reorient Housing and Enfranchisement the same way in any factor maps and loading plots.
sign_flip <- c(beta1 = -1, beta2 = 1, beta3 = -1,
               delta1 = -1, delta2 = -1, delta3 = 1, delta4 = 1)
## One row per scale; quantity (purple) and presence (orange) overlaid, explained by a legend.
row_of <- c(beta1 = "Housing-tenure\u2014density", beta2 = "Deprivation", beta3 = "Disenfranchisement",
            delta1 = "Service roads", delta3 = "Housing-tenure\u2014density",
            delta2 = "Deprivation", delta4 = "Disenfranchisement")
part_of <- c(beta1 = "Effect on quantity", beta2 = "Effect on quantity", beta3 = "Effect on quantity",
             delta1 = "Effect on presence", delta2 = "Effect on presence",
             delta3 = "Effect on presence", delta4 = "Effect on presence")
row_order <- c("Service roads", "Housing-tenure\u2014density", "Deprivation", "Disenfranchisement")
part_colors <- c("Effect on quantity" = "purple", "Effect on presence" = "orange")

pars_a <- names(sign_flip)
posterior_long <- as.data.frame(mapply(function(p, f) p * f, post[pars_a], sign_flip[pars_a],
                                       SIMPLIFY = FALSE)) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  mutate(scale_row = factor(row_of[parameter], levels = rev(row_order)),   # first row on top
         part = factor(part_of[parameter], levels = names(part_colors)))

cat("\nFigure 3a, re-expressed effects: median [95% interval]\n")
print(posterior_long |>
        group_by(scale_row, part) |>
        summarise(median = round(median(value), 2),
                  lo = round(quantile(value, 0.025), 2),
                  hi = round(quantile(value, 0.975), 2), .groups = "drop") |>
        arrange(desc(scale_row), part))

plot_3a <- ggplot(posterior_long,
                  aes(x = value, y = scale_row, fill = part, group = interaction(scale_row, part))) +
  geom_density_ridges(scale = 0.9, alpha = 0.6, color = "white", rel_min_height = 0.01) +
  theme_minimal(base_size = 10) +
  scale_fill_manual(values = part_colors, name = NULL) +
  labs(x = "Coefficient Value", y = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  theme(
    axis.text.y = element_text(face = "bold", size = 8),
    legend.position = "top",
    legend.justification = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

## (b) Expected FW score map
map_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

plot_3b <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_F123_all), color = NA) +
  stretched_fill() +
  map_theme +
  theme(plot.margin = margin(2, 2, 2, 2),
        legend.justification = "center",
        legend.box.margin = margin(0, 0, 0, -10)) +
  coord_sf(expand = FALSE)

## Layout: share of the figure height given to the map (b); panel (a) gets the rest.
## Increase map_share to squish (a) further.
fig3_height <- 6.0    # inches
map_share   <- 0.64
fig3 <- ggdraw() +
  draw_plot(plot_3a, x = 0, y = map_share, width = 1, height = 1 - map_share - 0.02) +
  draw_plot(plot_3b, x = 0, y = 0,         width = 1, height = map_share) +
  draw_plot_label(label = c("(a)", "(b)"), x = c(0.015, 0.015), y = c(0.995, map_share + 0.005),
                  size = 14, fontface = "bold")

ggsave(here("figs", "Figure3_F123_all.pdf"), fig3, width = 4.875, height = fig3_height, dpi = 300)

## =========================================================
## Figure 4
## =========================================================
pred_diff <- log10(van_dat_fa$pred_F123_all) - log10(van_dat_fa$pred_proxy)
cat("Log10 ratio (process / proxy): ", paste(signif(range(pred_diff), 3), collapse = " to "), "\n")

plot_4a <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_F123_all), color = NA) +
  stretched_fill() + map_theme

plot_4b <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_proxy), color = NA) +
  stretched_fill() + map_theme

plot_4c <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_proxy), color = NA) +
  scale_fill_viridis_c(option = "inferno", name = "Expected FW Score",
                       labels = function(x) sprintf("%.2f", x), na.value = mask_col) +
  map_theme

van_dat_fa$pred_diff <- pred_diff

## Even (symmetric) diverging scale with a neutral band around 0. The limits are +/- the larger
## of |1st percentile| and |99th percentile|; values beyond take the end colours. Legend in ratios.
q <- quantile(pred_diff, c(0.01, 0.99), na.rm = TRUE)
L  <- max(abs(q), 0.01)            # symmetric limits: same range on both sides of 0
lo <- -L; hi <- L
cat("Log10 ratio quantiles:\n"); print(signif(quantile(pred_diff, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 3))
## Neutral band: cells where the models are within this factor of each other get the
## neutral colour (1.25 means ratios between 0.8x and 1.25x). Widen or narrow as needed.
neutral_ratio <- 1.25
band <- log10(neutral_ratio)
L <- max(L, 2 * band); lo <- -L; hi <- L

## vik (Crameri Scientific Colour Maps): perceptually uniform and colour-vision-deficiency
## friendly. begin/end drop the darkest 15% at each end so the extremes are less dark.
n_side <- 50
vik_cols <- scico(2 * n_side + 1, palette = "vik", begin = 0.15, end = 0.85)
mid_col  <- vik_cols[n_side + 1]
cols_d <- c(vik_cols[1:n_side], mid_col, mid_col, vik_cols[(n_side + 2):(2 * n_side + 1)])
vals_d <- c(seq(lo, -band, length.out = n_side + 1)[-(n_side + 1)],   # blue ramp up to -band
            -band, band,                                              # flat neutral band
            seq(band, hi, length.out = n_side + 1)[-1])               # red ramp from +band
## Legend ticks: evenly spaced in log, symmetric around 1x (the neutral band is not labelled)
ratio_brks <- c(1/10, 1/5, 1/2, 1, 2, 5, 10)
ratio_brks <- ratio_brks[log10(ratio_brks) >= lo & log10(ratio_brks) <= hi]
ratio_lab  <- ifelse(ratio_brks < 1, paste0("1/", round(1 / ratio_brks), "\u00d7"),
                     paste0(ratio_brks, "\u00d7"))

plot_4d <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_diff), color = NA) +
  scale_fill_gradientn(colours = cols_d, values = rescale(vals_d, from = c(lo, hi)),
                       limits = c(lo, hi), oob = squish, na.value = mask_col,
                       breaks = log10(ratio_brks), labels = ratio_lab,
                       name = "Ratio") +
  guides(fill = guide_colourbar(barheight = unit(10, "lines"))) +
  map_theme

## Build each panel on its own first, so any error stops the script with its real message
## (cowplot otherwise draws a blank space in place of a panel that fails).
for (nm in c("plot_4a", "plot_4b", "plot_4c", "plot_4d")) {
  cat("Building", nm, "...\n")
  invisible(ggplotGrob(get(nm)))
}

## Panel titles as grid labels, so they sit at the same position in every panel
fig4 <- plot_grid(
  plot_4a, plot_4b, plot_4c, plot_4d,
  ncol = 2, align = "hv", axis = "tblr",
  labels = c("(a) Process model", "(b) Proxy model",
             "(c) Proxy model (own scale)", "(d) Process \u00f7 proxy"),
  label_size = 13, hjust = 0, label_x = 0.02, label_y = 0.99
)
ggsave(here("figs", "Figure4_F123_all.pdf"), fig4, device = "pdf",
       width = 10, height = 8.3333334, units = "in")

## =========================================================
## Supplementary: 750 m average, edge cells with < 50% coverage dropped
## =========================================================
city_utm <- st_transform(van_dat_fa["pred_F123_all"], 26910)
grid750 <- st_sf(geometry = st_make_grid(city_utm, cellsize = 750))
grid750$cell <- seq_len(nrow(grid750))
pts <- st_join(suppressWarnings(st_centroid(city_utm)), grid750, join = st_within)
cell_means <- st_drop_geometry(pts) |>
  group_by(cell) |>
  summarise(pred = mean(pred_F123_all, na.rm = TRUE), n_sub = n(),
            n_valid = sum(!is.na(pred_F123_all)), .groups = "drop") |>
  mutate(coverage = n_sub / max(n_sub), valid_share = n_valid / n_sub)
keep750 <- cell_means$coverage >= 0.5 & cell_means$valid_share >= 0.5
cat("750 m cells:", nrow(cell_means), "| dropped (< 50% inside the city, or < 50% within",
    "sampled conditions):", sum(!keep750), "\n")
grid750 <- inner_join(grid750, cell_means[keep750, ], by = "cell")

lims750 <- quantile(grid750$pred, c(0.01, 0.99))
p750 <- ggplot(grid750) +
  geom_sf(aes(fill = pred), color = NA) +
  scale_fill_viridis_c(option = "inferno", limits = lims750, oob = squish,
                       name = "Mean expected\nFW Score") +
  map_theme
ggsave(here("figs", "FigureS_F123_all_750m.pdf"), p750, width = 8, height = 6)

cat("\nDone: Figure3_F123_all.pdf, Figure4_F123_all.pdf, FigureS_F123_all_750m.pdf\n")
