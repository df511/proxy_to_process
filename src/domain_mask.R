## domain_mask.R
## Flag city cells whose census-derived scales do not describe local conditions, decided at the
## DISSEMINATION AREA (DA) level, since every cell takes its DA's census values.
##
## DAs are reconstructed from DA-level attributes that every cell of a DA shares
## (DA area, population density, income, % renters, deprivation indices). For each DA:
##   veg_share = mean NaturalVegSum over the DA's cells  (share of the DA that is natural vegetation)
## A whole DA is MASKED if its population density is low:
##   pop_km < pop_thresh                       (main rule)
## and, optionally (use_veg = TRUE), also if it is mostly natural vegetation:
##   veg_share >= veg_thresh
## except that cells inside a study-site outline are always kept (sampled areas, e.g. the
## 2nd Beach concession edge).
## Adds logical in_domain to van_dat_fa in data/van_dat_predictions_16var.Rdata and
## data/van_dat_final_750_16var.Rdata; draws figs/domain_mask_check.pdf;
## writes data/domain_mask_DAs.csv (one row per DA).
## Adjust the thresholds below and rerun until the check map looks right.
## Run from the proxy_to_process project root.

library(here)
library(sf)
library(dplyr)
library(ggplot2)

load(here("data", "dat_no_ub_750.Rdata"))   # dat_no_ub (transect cells), for reference values

pop_thresh <- 650     # people per km2; mask DAs below this density (main rule); 650 chosen by inspecting figs/domain_mask_check.pdf
use_veg    <- FALSE   # also mask DAs that are mostly natural vegetation?
veg_thresh <- 0.5     # DA share of natural vegetation (only used if use_veg = TRUE)

## ---------------------------------------------------------------
## 1. Reconstruct DAs
## ---------------------------------------------------------------
load(here("data", "van_dat_predictions_16var.Rdata"))   # van_dat_fa, k_housing, k_enfr, k_depriv
city <- st_drop_geometry(van_dat_fa)

da_cols <- intersect(c("Area (sq km)", "pop_km", "household_income", "rent_pct",
                       "SCOREMAT", "SCORESOC"), names(city))
cat("Reconstructing DAs from:", paste(da_cols, collapse = ", "), "\n")
da_key <- do.call(paste, c(lapply(city[da_cols], function(x) signif(x, 8)), sep = "|"))
da_id  <- as.integer(factor(da_key))
cat("Reconstructed DAs:", max(da_id),
    "(the City of Vancouver has roughly 1,000 DAs; a number in that range means the grouping worked)\n")

da <- data.frame(da_id = da_id, veg = city$NaturalVegSum, pop_km = city$pop_km,
                 area_km2 = city[["Area (sq km)"]]) |>
  group_by(da_id) |>
  summarise(n_cells = n(),
            veg_share = mean(veg, na.rm = TRUE),
            pop_km = first(pop_km),
            area_km2 = first(area_km2),
            .groups = "drop") |>
  mutate(population = round(pop_km * area_km2),
         mask_veg = use_veg & veg_share >= veg_thresh,
         mask_pop = pop_km < pop_thresh,
         masked_DA = mask_veg | mask_pop)

## ---------------------------------------------------------------
## 2. Cell-level flag
## ---------------------------------------------------------------
load(here("data", "sites_3857.Rdata"))
sites <- st_transform(sites_3857[-c(8), ], st_crs(van_dat_fa))
in_site <- lengths(st_intersects(suppressWarnings(st_centroid(st_geometry(van_dat_fa))), sites)) > 0

cell_masked_DA <- da$masked_DA[match(da_id, da$da_id)]
masked <- cell_masked_DA & !in_site
van_dat_fa$in_domain <- !masked

cat("\nDA population density (per km2), quantiles across DAs:\n")
print(round(quantile(da$pop_km, c(0, .01, .02, .05, .1, .25, .5))))
cat("Transect cells, population density quantiles:\n")
print(round(quantile(dat_no_ub$pop_km, c(0, .01, .05, .1, .5))))
cat("\nRule: mask DAs with population density <", pop_thresh, "per km2",
    if (use_veg) paste("or natural vegetation share >=", veg_thresh) else "(vegetation rule off)", "\n")
cat("DAs masked:", sum(da$masked_DA), "of", nrow(da),
    "(vegetation:", sum(da$mask_veg), "| low population:", sum(da$mask_pop),
    "| both:", sum(da$mask_veg & da$mask_pop), ")\n")
cat("Cells masked:", sum(masked), sprintf("(%.1f%% of the city)", 100 * mean(masked)),
    "| kept because inside a study site:", sum(cell_masked_DA & in_site), "\n")

cat("\nLowest-density DAs (area in km2; veg_share = share natural vegetation):\n")
print(da |> arrange(pop_km) |>
        transmute(da_id, n_cells, area_km2 = round(area_km2, 3), population,
                  pop_km = round(pop_km), veg_share = round(veg_share, 2), masked_DA) |>
        head(20))
cat("\nDAs with the highest natural vegetation share:\n")
print(da |> arrange(desc(veg_share)) |>
        transmute(da_id, n_cells, area_km2 = round(area_km2, 3), population,
                  pop_km = round(pop_km), veg_share = round(veg_share, 2), masked_DA) |>
        head(15))
write.csv(da, here("data", "domain_mask_DAs.csv"), row.names = FALSE)

save(van_dat_fa, k_housing, k_enfr, k_depriv, file = here("data", "van_dat_predictions_16var.Rdata"))

## ---------------------------------------------------------------
## 3. Check map
## ---------------------------------------------------------------
reason <- ifelse(!cell_masked_DA, "kept",
          ifelse(in_site, "kept (inside study site)",
          ifelse(da$mask_veg[match(da_id, da$da_id)] & da$mask_pop[match(da_id, da$da_id)],
                 "DA: vegetation and low population",
          ifelse(da$mask_veg[match(da_id, da$da_id)], "DA: mostly natural vegetation",
                 "DA: low population"))))
van_dat_fa$mask_reason <- factor(reason, levels = c("kept", "kept (inside study site)",
  "DA: mostly natural vegetation", "DA: low population", "DA: vegetation and low population"))

p <- ggplot() +
  geom_sf(data = van_dat_fa, aes(fill = mask_reason), color = NA) +
  geom_sf(data = sites, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("kept" = "grey85", "kept (inside study site)" = "#4393c3",
                               "DA: mostly natural vegetation" = "#1b7837",
                               "DA: low population" = "#e08214",
                               "DA: vegetation and low population" = "#762a83"),
                    name = NULL, drop = FALSE) +
  labs(title = "Masked dissemination areas (coloured) and kept cells (grey / blue)",
       subtitle = paste0("Mask a DA if population density < ", pop_thresh, " per km²",
                         if (use_veg) paste0(" or ≥ ", 100 * veg_thresh, "% natural vegetation") else "",
                         "; study sites always kept")) +
  theme_void()
ggsave(here("figs", "domain_mask_check.pdf"), p, width = 10, height = 7)

## Same flag on the factor-score file (same cells and order as the predictions file)
mask <- van_dat_fa$in_domain
load(here("data", "van_dat_final_750_16var.Rdata"))
stopifnot(nrow(van_dat_fa) == length(mask))
van_dat_fa$in_domain <- mask
save(van_dat_fa, file = here("data", "van_dat_final_750_16var.Rdata"))

cat("\nSaved in_domain to both files. Check figs/domain_mask_check.pdf.\n")
