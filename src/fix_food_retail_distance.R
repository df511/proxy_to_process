## fix_food_retail_distance.R
## The city grid's nearest_food_retail is in DEGREES: compile_data_YVRreproject.R
## (lines ~675-688) transforms to EPSG:4326 and runs FNN::get.knnx on lon/lat coordinates.
## The transect grid's nearest_food_retail (compile_data_adjustfw.R, line ~1058) is
## st_distance() in EPSG:3857, i.e. metres.
##
## This script rebuilds the food-retail locations from the same raw files, computes the
## distance the transect way (cell polygon to nearest food location, EPSG:3857) for BOTH
## grids, checks that it reproduces the transect values the models were fit on, and saves
## a corrected city file.
##
## Run from the proxy_to_process project root. Output:
##   data/van_dat_demo_lc_retail_shelters_roads_mods_750_fooddist.Rdata  (van_dat, corrected)

library(here)
library(sf)
library(dplyr)
library(jsonlite)

to3857 <- function(x) st_geometry(st_transform(x, 3857))
pts_4326 <- function(lon, lat) {
  ok <- !is.na(lon) & !is.na(lat)
  st_as_sf(data.frame(lon = lon[ok], lat = lat[ok]), coords = c("lon", "lat"), crs = 4326)
}

## ---------------------------------------------------------------
## 1. Food-retail locations (same sources and filters as the compile scripts)
## ---------------------------------------------------------------
stores <- fromJSON(here("data", "raw", "storefronts-inventory.json"))
stores <- stores[stores$retail_category %in%
                   c("Food & Beverage", "Convenience Goods", "Entertainment and Leisure"), ]
g_stores <- to3857(pts_4326(stores$geo_point_2d$lon, stores$geo_point_2d$lat))

vendors <- fromJSON(here("data", "raw", "food-vendors.json"))
g_vendors <- to3857(pts_4326(vendors$geo_point_2d$lon, vendors$geo_point_2d$lat))

lowcost <- fromJSON(here("data", "raw", "free-and-low-cost-food-programs.json"))
g_lowcost <- to3857(pts_4326(as.numeric(lowcost$longitude), as.numeric(lowcost$latitude)))

g_rest <- to3857(st_read(here("data", "raw", "restaurants_fastfood_CoV_quickOSM.geojson"),
                         quiet = TRUE))

cat("Food locations: storefronts", length(g_stores), "| vendors", length(g_vendors),
    "| low-cost programs", length(g_lowcost), "| OSM restaurants", length(g_rest), "\n")

sources <- list(
  all            = c(g_stores, g_vendors, g_lowcost, g_rest),
  no_restaurants = c(g_stores, g_vendors, g_lowcost)
)

nearest_dist <- function(cells, food) {
  idx <- st_nearest_feature(cells, food)
  as.numeric(st_distance(cells, food[idx], by_element = TRUE))
}

## ---------------------------------------------------------------
## 2. Check: does this reproduce the transect values used in the models?
## ---------------------------------------------------------------
load(here("data", "dat_no_ub_750.Rdata"))   # dat_no_ub
transect_cells <- st_transform(st_geometry(dat_no_ub), 3857)
old_t <- dat_no_ub$nearest_food_retail

check <- sapply(sources, function(food) {
  new_t <- nearest_dist(transect_cells, food)
  c(correlation = cor(new_t, old_t),
    median_abs_diff_m = median(abs(new_t - old_t)),
    max_abs_diff_m = max(abs(new_t - old_t)))
})
cat("\nRecomputed vs. original transect distances:\n")
print(round(check, 4))

best <- names(sources)[which.min(check["median_abs_diff_m", ])]
cat("\nSource set that best reproduces the transect values:", best, "\n")
if (check["max_abs_diff_m", best] > 1) {
  cat("WARNING: no source set reproduces the transect values exactly (max difference",
      round(check["max_abs_diff_m", best], 1), "m). Stop and check before continuing.\n")
}

## ---------------------------------------------------------------
## 3. Corrected city distances
## ---------------------------------------------------------------
load(here("data", "van_dat_demo_lc_retail_shelters_roads_mods_750.Rdata"))  # van_dat
van_dat$nearest_food_retail_deg <- van_dat$nearest_food_retail   # keep the old (degree) value
van_dat$nearest_food_retail <- nearest_dist(st_transform(st_geometry(van_dat), 3857),
                                            sources[[best]])

cat("\nCity nearest_food_retail, corrected (m):\n");  print(summary(van_dat$nearest_food_retail))
cat("Transect nearest_food_retail (m):\n");           print(summary(old_t))
cat("Correlation of corrected city distance with old degree values:",
    round(cor(van_dat$nearest_food_retail, van_dat$nearest_food_retail_deg), 3), "\n")

save(van_dat, file = here("data", "van_dat_demo_lc_retail_shelters_roads_mods_750_fooddist.Rdata"))
