## food_retail_locations.R
## Helper: all food-retail locations used for distance to food retail, as one sfc in EPSG:3857.
## Same sources and filters as compile_data_adjustfw.R (transect grid):
##   City of Vancouver storefronts (Food & Beverage, Convenience Goods, Entertainment and Leisure),
##   food vendors, free and low-cost food programs, and OSM restaurants / fast food.
## fix_food_retail_distance.R confirmed this set reproduces the transect distances exactly.
##
## Usage: source(here("src", "food_retail_locations.R")); food <- food_retail_locations()
##        nearest_food_retail_m(cells, food)   # metres, cell polygon to nearest location

food_retail_locations <- function() {
  to3857 <- function(x) sf::st_geometry(sf::st_transform(x, 3857))
  pts_4326 <- function(lon, lat) {
    ok <- !is.na(lon) & !is.na(lat)
    sf::st_as_sf(data.frame(lon = lon[ok], lat = lat[ok]), coords = c("lon", "lat"), crs = 4326)
  }
  stores <- jsonlite::fromJSON(here::here("data", "raw", "storefronts-inventory.json"))
  stores <- stores[stores$retail_category %in%
                     c("Food & Beverage", "Convenience Goods", "Entertainment and Leisure"), ]
  vendors <- jsonlite::fromJSON(here::here("data", "raw", "food-vendors.json"))
  lowcost <- jsonlite::fromJSON(here::here("data", "raw", "free-and-low-cost-food-programs.json"))
  rest <- sf::st_read(here::here("data", "raw", "restaurants_fastfood_CoV_quickOSM.geojson"),
                      quiet = TRUE)
  c(to3857(pts_4326(stores$geo_point_2d$lon, stores$geo_point_2d$lat)),
    to3857(pts_4326(vendors$geo_point_2d$lon, vendors$geo_point_2d$lat)),
    to3857(pts_4326(as.numeric(lowcost$longitude), as.numeric(lowcost$latitude))),
    to3857(rest))
}

nearest_food_retail_m <- function(cells, food) {
  cells <- sf::st_transform(sf::st_geometry(cells), 3857)
  idx <- sf::st_nearest_feature(cells, food)
  as.numeric(sf::st_distance(cells, food[idx], by_element = TRUE))
}
