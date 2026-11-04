## supplement_factor_maps.R
## City-wide maps of the three socio-political scales, oriented so that a HIGHER score is the
## less privileged condition (same orientation as Figure 3):
##   Housing-tenure-density : -Factor(housing)   high = more renters, lower income, denser
##   Deprivation            : +Factor(deprivation) high = more social and material deprivation
##   Disenfranchisement     : -Factor(enfranchisement) high = less employment, lower female share
##
## Shared, symmetric diverging scale centred on 0 (the transect mean), vik palette
## (colour-vision-deficiency friendly), limits at the largest |1st/99th percentile| across the
## three scales; values beyond are drawn in the end colours. Study sites outlined and labelled.
##
## Outputs:
##   figs/FigureS_factor_maps.pdf                 three panels, one column
##   figs/FigureS_factor_map_housing.pdf, _deprivation.pdf, _disenfranchisement.pdf
##   data/factor_loadings_16var_oriented.csv      loadings with the same orientation
## Run from the proxy_to_process project root. A few minutes (three 410k-cell maps).

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)
library(scico)
library(cowplot)

## ---------------------------------------------------------------
## 1. Oriented factor scores and loadings
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))     # Lambda
load(here("data", "van_dat_final_750_16var.Rdata"))   # van_dat_fa (city cells, corrected scores)

pick <- function(v) unname(which.max(abs(Lambda[v, ])))   # unname: keep "housing", not "housing.Factor1"
k <- c(housing = pick("housing_highdensity"),
       deprivation = pick("pct_rent_subsidized"),
       disenfranchisement = pick("employed_pct"))
orient <- c(housing = -1, deprivation = 1, disenfranchisement = -1)
titles <- c(housing = "Housing-tenure—density",
            deprivation = "Deprivation",
            disenfranchisement = "Disenfranchisement")
subtitles <- c(housing = "Higher = more renters, lower income, denser housing",
               deprivation = "Higher = more social and material deprivation",
               disenfranchisement = "Higher = less employment, lower female share")

for (s in names(k)) van_dat_fa[[s]] <- orient[[s]] * van_dat_fa[[paste0("Factor", k[[s]])]]

## Cells in masked DAs (domain_mask.R: population density below the threshold) are removed
mask_col <- "grey85"
if ("in_domain" %in% names(van_dat_fa)) {
  ## Masked DAs are removed from the maps entirely (drawn as empty space)
  cat("Removing", sum(!van_dat_fa$in_domain), "masked cells from the maps\n")
  van_dat_fa <- van_dat_fa[van_dat_fa$in_domain, ]
} else {
  cat("No in_domain column found: run src/domain_mask.R to mask parkland/unpopulated cells.\n")
}

## Sanity check: oriented scores vs. key raw variables (city cells)
cat("Correlations of oriented scores with raw variables (city):\n")
raw <- st_drop_geometry(van_dat_fa)[, c("rent_pct", "household_income", "housing_highdensity",
                                        "pct_rent_subsidized", "one_parent_pct",
                                        "employed_pct", "female_pct")]
print(round(cor(raw, st_drop_geometry(van_dat_fa)[, names(k)], use = "complete.obs"), 2))

L_oriented <- sweep(Lambda[, k], 2, orient, "*")
colnames(L_oriented) <- titles[names(k)]
write.csv(round(L_oriented, 3), here("data", "factor_loadings_16var_oriented.csv"))
cat("\nOriented loadings:\n"); print(round(L_oriented, 2))

## ---------------------------------------------------------------
## 2. Study sites and labels (same placement as the supplementary covariate maps)
## ---------------------------------------------------------------
load(here("data", "sites_3857.Rdata"))
sites_3857 <- sites_3857[-c(8), ]
site_labels <- c("Collingwood", "McKenzie", "Shaugnessy", "Queen Elizabeth",
                 "False Creek", "Downtown", "West End", "43 & Churchill",
                 "14 & Spruce", "19 & Yukon", "Strathcona", "2nd Beach")
sites_3857 <- st_transform(sites_3857, st_crs(van_dat_fa))

label_coords <- st_coordinates(suppressWarnings(st_centroid(sites_3857)))
label_center <- c(mean(label_coords[, 1]), min(label_coords[, 2]) - 500)
dx <- label_coords[, 1] - label_center[1]
dy <- label_coords[, 2] - label_center[2]
dist <- sqrt(dx^2 + dy^2)
push_dist <- 1000
label_df <- data.frame(X = label_coords[, 1] + dx / dist * push_dist,
                       Y = label_coords[, 2] + dy / dist * push_dist,
                       label = site_labels)
i_col <- which(site_labels == "Collingwood"); i_mck <- which(site_labels == "McKenzie")
label_df$Y[i_col] <- label_df$Y[i_col] + 120; label_df$X[i_col] <- label_df$X[i_col] - 120
label_df$Y[i_mck] <- label_df$Y[i_mck] + 300; label_df$X[i_mck] <- label_df$X[i_mck] + 300

offsets <- expand.grid(dx = c(-0.5, 0, 0.5), dy = c(-0.5, 0, 0.5))
offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), ]
outline_df <- do.call(rbind, lapply(seq_len(nrow(offsets)), function(i)
  transform(label_df, X = X + offsets$dx[i], Y = Y + offsets$dy[i])))

## ---------------------------------------------------------------
## 3. Shared symmetric diverging scale
## ---------------------------------------------------------------
L <- max(abs(unlist(lapply(names(k), function(s)
  quantile(van_dat_fa[[s]], c(0.01, 0.99), na.rm = TRUE)))))
cat("\nShared colour limits: ±", signif(L, 3), "\n")
fill_scale <- scale_fill_gradientn(
  colours = scico(101, palette = "vik", begin = 0.15, end = 0.85),
  limits = c(-L, L), oob = squish, na.value = mask_col, name = "Score"
)

map_one <- function(s, label_size = 3) {
  ggplot() +
    geom_sf(data = van_dat_fa, aes(fill = .data[[s]]), color = NA) +
    geom_sf(data = sites_3857, fill = NA, color = "black", linewidth = 1.2) +
    geom_sf(data = sites_3857, fill = NA, color = "white", linewidth = 0.5) +
    geom_text(data = outline_df, aes(X, Y, label = label),
              size = label_size + 0.05, fontface = "bold", color = "black") +
    geom_text(data = label_df, aes(X, Y, label = label),
              size = label_size, fontface = "bold", color = "white") +
    fill_scale +
    guides(fill = guide_colourbar(barheight = unit(7, "lines"))) +
    labs(title = titles[[s]], subtitle = subtitles[[s]]) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8.5),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
}

maps <- lapply(setNames(names(k), names(k)), map_one)
for (s in names(maps)) {
  ggsave(here("figs", paste0("FigureS_factor_map_", s, ".pdf")), maps[[s]], width = 8, height = 6)
}

fig <- plot_grid(plotlist = maps, ncol = 1, labels = c("(a)", "(b)", "(c)"),
                 label_size = 12, hjust = 0, label_x = 0.01)
ggsave(here("figs", "FigureS_factor_maps.pdf"), fig, width = 7, height = 15)

cat("\nDone: FigureS_factor_maps.pdf and one PDF per scale.\n")
