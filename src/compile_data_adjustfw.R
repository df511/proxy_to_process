## load, compile, and clean data for food waste modelling

## Written by Daniel Forrest
## Started on February 27, 2025

####ENSURE THAT YOU ENTER YOUR PERSONAL CENSUSMAPPER API KEY IN LINE 614
####OR THIS SCRIPT WILL NOT RUN

# Load necessary libraries
library(sf)
library(terra)
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(cancensus)
library(viridis)
library(MASS)
library(here)
library(sf)
library(geojsonsf)
library(matrixStats)
library(ggmap)
library(tidyterra)
library(jsonlite)
library(exactextractr)
library(progressr)
library(reshape2)
library(ggrepel)
library(gridExtra)



#### The following script is divided into modules, depending on the data type and intended output, 
#### so as to reduce run time, repeat operations, and organize by source data
#### If you are running for the first time, set all of the following parameters to TRUE
#### If you are running to alter a portion of the data compilation code,
#### set all stages prior to the focal stage to FALSE

### Module 1: gen grids
gen_grids <- TRUE  # Change to FALSE to load data instead
### Module 2: gen food waste
gen_foodwaste <- TRUE  # Change to FALSE to load data instead
### Module 3: gen Canadian census dat (demographics)
gen_demog <- TRUE  # Change to FALSE to load data instead
### Module 4: read in land cover data, calculate proportion of cell containing each of 13 land cover classes
gen_lc <- TRUE
### Module 5: read in food retail location data
gen_foodretail <- TRUE
### Module 6: read in food retail location data
gen_shelters <- TRUE
### Module 7: read in food retail location data
gen_roads <- TRUE
### Module 8: read in property values from tax data
gen_finalmods <- TRUE
### Module 9: read in property values from tax data
gen_center_std <- TRUE
### Module 10: conduct Factor Analysis on final variable set
gen_factanal <- TRUE

#########################################################################################################
############ Module 1: generate grids ###########################

if (gen_grids) {
  ### Set parameters
  
  # Define grid cell size (e.g., 2x2)
  ### size of a single side in meters
  #cell_size <- 2.2360679774998 ### for 5 sq. meters.
  #cell_size <- 5 ## for 25 sq. meters.
  cell_size <- 27.3861278 ### for ~750 sq. meters (I measured several parcels across the region and the range was ~600 to 900 sq. m)
  
  #### build functions
  
  # # Function to create an internally fitting grid for each polygon with its site name
  # create_internal_grid <- function(polygon, cell_size, site_name) {
  #   bbox <- st_bbox(polygon)
  #   
  #   # Generate grid in metric units
  #   grid <- st_make_grid(polygon, cellsize = cell_size, square = TRUE, 
  #                        offset = c(bbox["xmin"], bbox["ymin"]))  
  #   
  #   # Convert to sf object
  #   grid_sf <- st_sf(geometry = grid)
  #   
  #   # Keep only full grid cells within the polygon
  #   grid_inside <- grid_sf[st_within(grid_sf, polygon, sparse = FALSE), ]
  #   
  #   # Add site name as an attribute
  #   grid_inside$name <- site_name
  #   
  #   return(grid_inside)
  # } 
  # 
  
  
  create_internal_grid <- function(polygon, cell_size, site_name, overlap_buffer = 5) {
    # Buffer outward by the overlap distance
    polygon_buffered <- st_buffer(polygon, overlap_buffer)
    
    bbox <- st_bbox(polygon_buffered)
    
    # Generate grid with the buffered polygon
    grid <- st_make_grid(polygon_buffered, cellsize = cell_size, square = TRUE, 
                         offset = c(bbox["xmin"], bbox["ymin"]))  
    
    grid_sf <- st_sf(geometry = grid)
    
    # Option 1: Keep grid cells that intersect (i.e., touch or overlap) original polygon
    grid_overlap <- grid_sf[st_intersects(grid_sf, polygon, sparse = FALSE), ]
    
    # Option 2: OR keep all cells that intersect the buffered polygon (wider coverage)
    # grid_overlap <- grid_sf[st_intersects(grid_sf, polygon_buffered, sparse = FALSE), ]
    
    # Add site name
    grid_overlap$name <- site_name
    
    return(grid_overlap)
  }
  
  
  ###### read files
  
  lancaster <- st_read(here("data","sites","lancaster_sites.kml"))
  weber <- st_read(here("data","sites","weber_site_edits.kml"))
  weber_rep <- weber[1:3,] ## Weber sites I visited and surveyed (drop one outside of CoV borders)
  plot(weber_rep) # check 
  weber_lancaster <- do.call(rbind, list(lancaster, weber_rep))
  #weber_lancaster_no43churchill <- do.call(rbind, list(lancaster, weber_rep2))
  
  plot(weber_lancaster) # check
  
  strathcona <- st_read(here("data","sites","strathcona_final.kml"))
  plot(strathcona)
  second_beach <- st_read(here("data","sites", "2nd_beach_edit.kml"))
  plot(second_beach)
  
  
  sites <- do.call(rbind, list(lancaster, weber_rep, strathcona, second_beach))
  sites_3857 <- st_transform(sites, crs = 3857)
  sites_3857$Name <- c("lancaster_cw","lancaster_mk", "lancaster_sh" ,
                       "lancaster_qe", "lancaster_fc",  "lancaster_dt",
                       "lancaster_we", "lancaster_ub",
                       "43_churchill","14_spruce", "19_yukon",  "strathcona", "2nd_beach")
  # ### plot to check correct naming
  # # Ensure 'name' column exists
  # if (!"Name" %in% colnames(sites_3857)) {
  #   stop("Column 'name' not found in sites_3857!")
  # }
  # 
  # # Create the plot
  # ggplot() +
  #   geom_sf(data = sites_3857, fill = "lightblue", color = "black") +  # Plot polygons
  #   geom_sf_text(data = sites_3857, aes(label = Name), size = 3, color = "red") +  # Add names as text
  #   theme_minimal() +
  #   labs(title = "Sites with Names", x = "Longitude", y = "Latitude")


  save(sites_3857, file = here("data","sites_3857.Rdata"))
  
  
  # Apply function to each polygon, retaining "name" column
  grid_list <- mapply(create_internal_grid, 
                      st_geometry(sites_3857), 
                      cell_size, 
                      sites_3857$Name, 
                      SIMPLIFY = FALSE)
  
  # Combine grids
  final_grid <- do.call(rbind, grid_list)
  
  # ✅ Calculate cell area in square meters
  cell_areas <- st_area(final_grid)
  
  # Print mean cell area (if needed)
  mean_cell_area <- mean(cell_areas)
  print(paste("Mean cell size:", mean_cell_area, "square meters"))
  
 
  # Define output file name and resolution
  pdf(here("figs","grid_plot_750m_buffered.pdf"), width = 80, height = 64)  # Adjust size as needed
  
  # Plot the grid inside polygons
  plot(st_geometry(sites_3857), col = "lightblue", border = "black")
  plot(st_geometry(final_grid), add = TRUE, border = "red")
  
  # Close the PDF device to save the file
  dev.off()
  
  save(final_grid, file = here("data","grid_clipped_750m.Rdata"))
  #######
} else {
  ## file called "final_grid"
  load(file = here("data", "grid_clipped_750m.Rdata"))
}




grid_clipped_750m <- final_grid
st_crs(grid_clipped_750m) <- 3857
# Add a unique ID column
grid_clipped_750m <- grid_clipped_750m %>%
  mutate(grid_id = row_number())


##############################################################################################
############ Module 2: import, clean, and grid-ify food waste data ###########################

if (gen_foodwaste) {

###### import food waste data
foodwaste <- read_csv(here("data", "food_waste_data_09_19_23.csv"))
foodwaste <-  foodwaste %>% fill(c(fid, site, date), .direction = "down")
foodwaste_poly <- st_read(here("data", "waste_receptacles.gpkg"))
#plot(foodwaste_poly)

foodwaste_poly$fid <- seq_len(nrow(foodwaste_poly)) #### replicate FID
merged_data <- merge(foodwaste_poly, foodwaste, by = "fid") ### merge data frames
# filtered_data <- merged_data[merged_data$fid == "155", ] ### Test view
#######remove spaces in column "site" 
merged_data$site <- gsub(" ", "", merged_data$site)

###### import survey data
surv_dat <- read_csv(here("data","survey_info_09_19_23.csv"))
#### filter out trial surveys
surv_dat$pickup_type[surv_dat$pickup_type == "no pickup"] <- "none"
surv_dat <- surv_dat[c(5:48),]
surv_dist <- read_csv(here("data","survey_area_distance.csv"))
surv_dist <- surv_dist[c(1:10),c(1:4)] #### NEED TO REPLACE 19_Yukon WITH SHORTENED SURVEY

### calculate number of site visits
total_visits <- surv_dat %>%
  group_by(site) %>%
  summarize(total_visits = n_distinct(date))


###### merge survey data with food waste data
dat <- left_join(merged_data, surv_dat, by = c("site" = "site", "date" ="date"), relationship = "many-to-many")
dat <- left_join(dat, surv_dist, by = c("site" = "site"))
##Join the unique visits data with the main data
dat <- dat %>%
  left_join(total_visits, by = "site")

##### filter to only req'd columns
dat<- dat[,c(1,9:18, 23:36, 39)]
## rename for interpretability
dat <- dat %>% rename(waste_type = type.y,
                      bin_form = form.y,
                      bin_state = state.y,
                      bin_contents = contents.y,
                      bin_count = count.y,
                      bin_notes = notes.y)

####### format and fill
dat <- as.data.table(dat)

dat$bin_count[is.na(dat$bin_count)] <- 1 #replace NAs with 1 (if missing, always a single bin/receptacle described)

# Filter out rows in "bin_count" with "_spot" in the string
dat <- dat[!grepl("_spot", site)]
#### Fix data entry error
#Identify the row number where the shift should occur
#row_num <- 5375  

# Perform the shift #not performing as expected yet. 
#dat[row_num, `:=`(fl_density = bin_count,
#                  bin_contents = shift(dat[row_num, bin_contents], type = "lead"))]



# Function to convert columns from character to numeric
convert_to_numeric <- function(col) {
  # Check if the column is numeric (including integers and decimals)
  if (all(grepl("^\\d+\\.?\\d*$", col))) {
    as.numeric(col)
  } else {
    col
  }
}

# Identify numeric columns and convert them to numeric
dat[, (names(dat)) := lapply(.SD, convert_to_numeric), .SDcols = names(dat)]


print(dat)


##### convert to numeric
dat$bin_count <- as.numeric(dat$bin_count)
dat$bin_count[is.na(dat$bin_count)] <- 1 #replace NAs with 1 (if missing, always a single bin/receptacle described)
#remove spaces in column "site" 
dat[, site := gsub(" ", "", site)]

# Replace NAs and "no pickup" in the "pickup_type" column with "none"
dat[, pickup_type := ifelse(is.na(pickup_type) | pickup_type == "no pickup", "none", pickup_type)]
# Replace "garbage" in pickup_type with "misc." where site is "lancaster_fc"
dat[site == "lancaster_fc" & pickup_type == "garbage", pickup_type := "misc."]
dat[site == "lancaster_we", pickup_type := "misc."]
# Replace bin_form with "food_litter" where waste_type is "food_litter"
dat[waste_type == "food_litter", bin_form := "food_litter"]
# Replace bin_form with "food_litter" where waste_type is "food_litter"
dat[waste_type == "food_litter", bin_state := "food_litter"]

### correcting duplicate FIDs
dat$fid[c(269,271,272)] <- 452
dat$fid[c(1658)] <- 362
dat$fid[c(2949, 2950)] <- 659
dat$fid[c(5088)] <- 856

# ### create spatial data frame
dat_sf <- st_as_sf(dat, sf_column_name = "geometry")

# Vector of FIDs to remove
fids_to_remove <- c(608, 609, 612,988, 1029, 1064, 1065, 1074, 1150) ### all carcasses, 1150 is apple tree



# Filter out rows where fid is in the removal list
dat_sf <- dat_sf %>%
  filter(!fid %in% fids_to_remove)



dat_sf <- dat_sf %>%
  mutate(
    fl_density = if_else(fid %in% c(1005), "high", fl_density)
  )
# dat_sf <- dat_sf %>%
#   mutate(
#     waste_type = if_else(fid %in% c(1064, 1065, 1029), "food_litter", waste_type),
#     bin_form   = if_else(fid %in% c(1064, 1065, 1029), "food_litter", bin_form),
#     bin_state  = if_else(fid %in% c(1064, 1065, 1029), "food_litter", bin_state),
#     fl_density = if_else(fid %in% c(1029), "low", fl_density)
#   )

dat_sf <- dat_sf %>%
  mutate(
    is_food_litter = waste_type == "other" & !is.na(fl_density),
    waste_type = if_else(is_food_litter, "food_litter", waste_type),
    bin_form   = if_else(is_food_litter, "food_litter", bin_form),
    bin_state  = if_else(is_food_litter, "food_litter", bin_state)
  ) %>%
  select(-is_food_litter)  # clean up



######## Define Weights (can adjust later)
weights_state <- c(
  "food_litter" = 1,
  "open_w_food" = 0.5,
  "open_inbag_food" = 0.25,
  "food_residue" = 0.375,
  "restricted_food" = 0.2,
  "open_unknown" = 0.05,
  "open_inbag_unknown" = 0.05,
  "restricted_unknown" = 0.05,
  "restricted_empty" = 0.01,
  "restricted_nofood" = 0.02,
  "open_no_food" = 0.02,
  "open_empty" = 0.02,
  "open_inbag_nofood" = 0.03,
  "closed" = 0.01,
  "locked" = 0.01,
  "other" = 0.05,
  "NA" = 0
)

weights_form <- c(
  "dumpster" = 1.0,
  "tote" = 0.05,
  "home_bin" = 0.033,
  "plastic_bag" = 0.033,
  #"food_litter" = 0.5,  # will be combined with fld_weight and area
  "other" = 0.066,
  "stationary" = 0.2,
  "NA" = 0.0
)

# # Define weights for each bin_form
# weights_form <- c("dumpster" = 15,
#                   "tote" = 3,
#                   "home_bin" = 1,
#                   "plastic_bag" = 1,
#                   "food_litter" = 1,
#                   "other" = 2,
#                   "stationary" = 0.266, #about the same as a tote
#                   "NA" = 1)

weights_type <- c(
  "food_litter" = 1.0,
  "compost" = 0.33,
  "garbage" = 0.2,
  "recycling-containers" = 0.1,
  "recycling-glass" = 0.05,
  "recycling-papercardboard" = 0.005,
  "other" = 0.1,
  "NA" = 0.0
)

weights_fld <- c(
  "<visible" = 0.0078125,#1/4 of below
  "low" = 0.03125,#1/4 of below
  "medium" = 0.125, #1/4 of below
  "high" = 0.5
)


# Step 1: Compute area for all geometries
dat_sf <- dat_sf %>%
  mutate(area_m2 = as.numeric(st_area(.)))

# Step 2: Vectorized lookup and additive scoring
dat_sf <- dat_sf %>%
  mutate(
    state_weight = weights_state[match(bin_state, names(weights_state))],
    form_weight = weights_form[match(bin_form, names(weights_form))],
    type_weight = weights_type[match(waste_type, names(weights_type))],
    fld_weight = weights_fld[match(fl_density, names(weights_fld))],
    
    # Override form_weight if food_litter
    form_weight = if_else(bin_state == "food_litter", 0, form_weight),
    
    # Area-weighted fld_term ONLY for food litter
    fld_term = if_else(
      bin_state == "food_litter",
      fld_weight * (area_m2)^(1/8),  # log adjust area, add offset so never negative log(area_m2 + 0.0257)
      0
    ),
    
    # # Additive score (weights ~0–1)
    # additive_score = (1 * state_weight)
    # + (1.5 * form_weight)
    # + (1.2 * type_weight)
    # + (2 * fld_term),
    
    # Additive score (weights ~0–1)
    additive_score = (0.5 * state_weight)
    + (1 * form_weight)
    + (1 * type_weight)
    + (1 * fld_term),
    
    # Final weighted value per observation
    weighted_value = bin_count * additive_score
  )





# Ensure UTM for distance calculation
#dat_sf_centroids <- st_transform(dat_sf_centroids, 26910)
grid_clipped_utm <- st_transform(grid_clipped_750m, 26910)




dat_sf_utm <- st_transform(dat_sf, 26910)

dat_sf_centroids <- st_centroid(dat_sf_utm)

# Get site names
sites <- unique(dat_sf$site)

# Initialize the corrected site column
dat_sf_centroids$site_corrected <- NA_character_

# Transform grid for accurate distances
grid_clipped_utm <- st_transform(grid_clipped_750m, 26910)



# Identify which centroids are not inside any grid cell of their site
dat_sf_centroids$in_grid <- FALSE
dat_sf_centroids$geometry_corrected <- st_geometry(dat_sf_centroids)



for (s in unique(dat_sf_centroids$site)) {
  site_centroids <- dat_sf_centroids[dat_sf_centroids$site == s, ]
  site_grid <- grid_clipped_utm[grid_clipped_utm$name == s, ]
  
  if (nrow(site_centroids) == 0 || nrow(site_grid) == 0) next
  
  inside <- st_intersects(site_centroids, site_grid, sparse = FALSE)
  
  not_inside_idx <- which(!rowSums(inside))
  cat("Site:", s, "\n")
  cat("Number of centroids:", nrow(site_centroids), "\n")
  cat("Number inside:", sum(rowSums(inside)), "\n")
  cat("Number outside:", length(not_inside_idx), "\n")
  
  if (length(not_inside_idx) == 0) next
  
  nearest_idx <- st_nearest_feature(site_centroids[not_inside_idx, ], site_grid)
  
  # Fix: Use only geometry
  corrected_geoms <- st_centroid(st_geometry(site_grid)[nearest_idx])
  
  rows_to_correct <- which(dat_sf_centroids$site == s)[not_inside_idx]
  dat_sf_centroids$geometry_corrected[rows_to_correct] <- corrected_geoms
  dat_sf_centroids$in_grid[rows_to_correct] <- TRUE
}


ggplot() +
  geom_sf(data = site_grid, fill = NA, color = "black") +
  geom_sf(data = site_centroids, aes(color = rowSums(inside) > 0)) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
  ggtitle(paste("Inside/outside check for site:", s))


cat("Number of centroids moved:", sum(dat_sf_centroids$in_grid), "\n")


# Update geometry
st_geometry(dat_sf_centroids) <- dat_sf_centroids$geometry_corrected

# Reproject back to EPSG:3857 for compatibility with grid
dat_sf_centroids <- st_transform(dat_sf_centroids, 3857)




sites <- intersect(unique(dat_sf_centroids$site), unique(grid_clipped_750m$name))

dat_sf_joined_list <- lapply(sites, function(s) {
  st_join(
    dat_sf_centroids[dat_sf_centroids$site == s, ],
    grid_clipped_750m[grid_clipped_750m$name == s, ],
    left = FALSE
  )
})

dat_sf_joined <- do.call(rbind, dat_sf_joined_list)

# Step 5: Adjust by pickup type
setDT(dat_sf_joined)



dat_sf_centroids$was_moved <- dat_sf_centroids$in_grid

ggplot() +
  geom_sf(data = grid_clipped_750m, fill = "grey95", color = "black") +
  geom_sf(data = dat_sf_centroids, aes(color = was_moved), size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c(`TRUE` = "red", `FALSE` = "blue"),
    labels = c(`TRUE` = "Moved", `FALSE` = "Original"),
    name = "Centroid Status"
  ) +
  theme_minimal() +
  labs(title = "Centroids Moved to Nearest Grid Cell Within Site")



dat_sf_joined[, adjusted_weighted_value := weighted_value * fifelse(
  pickup_type == "misc.", 1 / total_visits,
  fifelse(pickup_type %in% c("recycling", "garbage"), 3 / 14,
          fifelse(pickup_type == "none", 8 / 14, 1)))
]

daily_sums <- dat_sf_joined[, .(
  daily_sum = sum(adjusted_weighted_value, na.rm = TRUE)
), by = .(grid_id, site, date)]


fw_scores <- daily_sums[, .(
  fw_score_weighted = sum(daily_sum, na.rm = TRUE),  # total weighted score
  fw_score_max = max(daily_sum, na.rm = TRUE)        # max per day per grid
), by = grid_id]


dat_grid <- merge(grid_clipped_750m, fw_scores, by = "grid_id", all.x = TRUE)

dat_grid <- dat_grid %>%
  mutate(fw_score_weighted = replace_na(fw_score_weighted, 0),
         fw_score_max= replace_na(fw_score_max, 0))

# Access the attribute data of the sf object and replace -Inf with NA
dat_grid[] <- lapply(dat_grid[], function(x) {
  if (is.numeric(x)) {
    x[x == -Inf] <- NA
  }
  return(x)
})

dat_grid <- dat_grid %>%
  dplyr::rename(site = name)


#### Save and name file!!
save(dat_grid, file = here("data","dat_grid_750.Rdata"))

#######
} else {
## object called "dat_grid"
load(file = here("data", "dat_grid_750.Rdata"))
}

##############################################################################################
############ Module 3: import, clean, and grid-ify Canadian census data ###########################


if (gen_demog) {

### insert your own Cancensus API key here (sign up here: https://censusmapper.ca/ ; and read the package info here: https://cran.r-project.org/web/packages/cancensus/vignettes/cancensus.html)
options(cancensus.api_key = 'INSERT YOUR API HERE') ### insert personalized API Key, sign-up here: https://censusmapper.ca/users/sign_up


### For how to use the cancensus package: https://cran.r-project.org/web/packages/cancensus/vignettes/cancensus.html


census_vec <- list_census_vectors(dataset ="CA21")

#read in census data via query of census mapper
census_data <- get_census(dataset = 'CA21', regions = list(CSD = "5915022"), 
                          vectors = c("v_CA21_1","v_CA21_906", "v_CA21_6", #population, median income, pop. density
                                      "v_CA21_5808","v_CA21_5865", "v_CA21_4204", "v_CA21_4201",
                                      "v_CA21_4238", "v_CA21_4239", "v_CA21_4237",
                                      "v_CA21_4875","v_CA21_4878", "v_CA21_4881",
                                      "v_CA21_4884", "v_CA21_4887", "v_CA21_4890",
                                      "v_CA21_4893", "v_CA21_4896", "v_CA21_4899",
                                      "v_CA21_4902", "v_CA21_4905", "v_CA21_4908",
                                      "v_CA21_4911","v_CA21_4914",
                                      "v_CA21_9", "v_CA21_10", "v_CA21_8", ## gender
                                      "v_CA21_4297", "v_CA21_4298", "v_CA21_4299", "v_CA21_4300","v_CA21_4314","v_CA21_4315", "v_CA21_4316","v_CA21_4260", "v_CA21_4261","v_CA21_4272","v_CA21_4273","v_CA21_4274", ##housing suitability, repairs...
                                      "v_CA21_434","v_CA21_435","v_CA21_436","v_CA21_437","v_CA21_438", "v_CA21_439", "v_CA21_440", ##housing types; https://www23.statcan.gc.ca/imdb/p3VD.pl?Function=getVD&TVD=144257
                                      "v_CA21_4392", "v_CA21_4401", "v_CA21_4407", "v_CA21_4410", ##immigration... citizens/non-citizens, non-immigrants, immigrants
                                      "v_CA21_4923", "v_CA21_4938", ### ethnicity, English, Chinese, re: Ley 1995
                                      "v_CA21_507", "v_CA21_499", ### total one parent families, total families
                                      "v_CA21_534", "v_CA21_510",## living alone, total persons in private households
                                      "v_CA21_483", "v_CA21_486", "v_CA21_489","v_CA21_453", # Total Separated, Divorced, Widowed, out of all Marital status respondents
                                      "v_CA21_6498", "v_CA21_6495", ## total employed out of all workforce
                                      "v_CA21_5802", "v_CA21_5799", ## total no high school diploma or equivalent out of respondents
                                      "v_CA21_452", # average household size
                                      "v_CA21_4917", #parent vector (ethnic or cultural origin), then 
                                      "v_CA21_5109", "v_CA21_5010", "v_CA21_5205", "v_CA21_4977", "v_CA21_5094", "v_CA21_5202", # Latin, Central or South American, n.o.s., Total African, n.o.s., Total West or Central Asian or Middle Eastern, n.o.s., Total European, n.o.s, Total South Asian, n.o.s,  Total East or Southeast Asian, n.o.s
                                      "v_CA21_4389", "v_CA21_4404","v_CA21_4872", # total_citizen_Q, total_immigrant, total vis minority
                                      "v_CA21_4968",#caucasian
                                      "v_CA21_251", ## 65+, divide by v_CA21_8
                                      "v_CA21_386", "v_CA21_389" #avg. age, median age 
                          ),
                          labels = "detailed", geo_format = "sf", level = "DA")

# 
# var <- colnames(census_data[,43])
# var <- var[1]
# clean_var1 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,44])
# var <- var[1]
# clean_var2 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,45])
# var <- var[1]
# clean_var3 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,46])
# var <- var[1]
# clean_var4 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,47])
# var <- var[1]
# clean_var5 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,48])
# var <- var[1]
# clean_var6 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,49])
# var <- var[1]
# clean_var7 <- gsub("\\\\", "", var)
# # var <- colnames(census_data[,50])
# # var <- var[1]
# # clean_var8 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,50])
# var <- var[1]
# clean_var8 <- gsub("\\\\", "", var)
# var <- colnames(census_data[,51])
# var <- var[1]
# clean_var9 <- gsub("\\\\", "", var)


var <- colnames(census_data[,84])
var <- var[1]
clean_var10 <- gsub("\\\\", "", var)
var <- colnames(census_data[,85])
var <- var[1]
clean_var11 <- gsub("\\\\", "", var)
var <- colnames(census_data[,86])
var <- var[1]
clean_var12 <- gsub("\\\\", "", var)
var <- colnames(census_data[,87])
var <- var[1]
clean_var13 <- gsub("\\\\", "", var)
var <- colnames(census_data[,88])
var <- var[1]
clean_var14 <- gsub("\\\\", "", var)
var <- colnames(census_data[,89])
var <- var[1]
clean_var15 <- gsub("\\\\", "", var)
var <- colnames(census_data[,90])
var <- var[1]
clean_var16 <- gsub("\\\\", "", var)

#rename columns
census_data <- census_data %>% 
  rename(household_income = `v_CA21_906: Median total income of household in 2020 ($)`,
         population = `v_CA21_1: Population, 2021`,
         minority_pop = `v_CA21_4875: Total visible minority population`,
         black_pop = `v_CA21_4884: Black`,
         indig_pop = `v_CA21_4204: Indigenous identity (39)`,
         indig_respondents = `v_CA21_4201: Total - Indigenous identity for the population in private households`,
         pop_km = `v_CA21_6: Population density per square kilometre`,
         owners = `v_CA21_4238: Owner`,
         renters = `v_CA21_4239: Renter`,
         female = `v_CA21_10: Total - Age`,
         male = `v_CA21_9: Total - Age`,
         single_detached = `v_CA21_435: Single-detached house`,
         semi_detached = `v_CA21_436: Semi-detached house`,
         row_houses = `v_CA21_437: Row house`,
         apt_flat = `v_CA21_438: Apartment or flat in a duplex`,
         low_rise = `v_CA21_439: Apartment in a building that has fewer than five storeys`,
         high_rise = `v_CA21_440: Apartment in a building that has five or more storeys`,
         # # thirtyonshelter_notsuitable = !!sym(clean_var1),
         # # thirtyonshelter_majorrepairs = !!sym(clean_var2),
         # # notsuitable_majorrepairs = !!sym(clean_var3),
         # # thirtyonshelter_notsuitable_majorepairs = !!sym(clean_var4),
         # # pct_rent_subsidized = !!sym(clean_var5),
         # # pct_rent_thirty = !!sym(clean_var6),
         # # total_occupied_housing = !!sym(clean_var7),
         # # #pct_housing_notsuitable = !!sym(clean_var7),
         # # #total_occupied_housing = !!sym(clean_var8),
         # # minor_repairs = !!sym(clean_var8),
         # # major_repairs = !!sym(clean_var9),
         citizens = `v_CA21_4392: Canadian citizens`,
         non_citizens = `v_CA21_4401: Not Canadian citizens`,
         non_immigrants = `v_CA21_4407: Non-immigrants`,
         immigrants = `v_CA21_4410: Immigrants`,
         English = `v_CA21_4923: English`,
         Chinese = `v_CA21_4938: Chinese`,
         LCS_American = !!sym(clean_var10),
         African = !!sym(clean_var11),
         MidEast_WC_Asian = !!sym(clean_var12),
         E_SE_Asian = !!sym(clean_var13),
         European = !!sym(clean_var14),
         Caucasian = !!sym(clean_var15),
         S_Asian = !!sym(clean_var16)
         )

#calculate percemtages of ethnicities per population and create new cols in df
census_data$black_pct <- census_data$black_pop / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$chinese_pct <- census_data$`v_CA21_4881: Chinese` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$southasian_pct <- census_data$`v_CA21_4878: South Asian` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$filipino_pct <- census_data$`v_CA21_4887: Filipino` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$arab_pct <- census_data$`v_CA21_4890: Arab` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$latinamerican_pct <- census_data$`v_CA21_4893: Latin American` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$southeastasian_pct <- census_data$`v_CA21_4896: Southeast Asian` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$westasian_pct <- census_data$`v_CA21_4899: West Asian` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$korean_pct <- census_data$`v_CA21_4902: Korean` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$japanese_pct <- census_data$`v_CA21_4905: Japanese` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$eastasian_pct <- census_data$japanese_pct + census_data$korean_pct + census_data$chinese_pct
census_data$eastsoutheastasian_pct <- census_data$japanese_pct + census_data$korean_pct + census_data$chinese_pct +census_data$filipino_pct + census_data$southeastasian_pct
census_data$visminority_pct <- census_data$`v_CA21_4908: Visible minority, n.i.e.` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$multivisminority_pct <- census_data$`v_CA21_4911: Multiple visible minorities` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$notvisminority_pct <- census_data$`v_CA21_4914: Not a visible minority` / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$indig_pct <- census_data$indig_pop / census_data$indig_respondents
census_data$minority_pct <- census_data$minority_pop / census_data$`v_CA21_4872: Total - Visible minority for the population in private households`
census_data$own_pct <- census_data$owners / census_data$`v_CA21_4237: Total - Private households by tenure`
census_data$rent_pct <- census_data$renters / census_data$`v_CA21_4237: Total - Private households by tenure`
census_data$male_pct <- census_data$male / census_data$`v_CA21_8: Total - Age`
census_data$female_pct <- census_data$female / census_data$`v_CA21_8: Total - Age`
census_data$english_pct <- census_data$English / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$chinese_origin_pct <- census_data$Chinese / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$LCS_American_origin_pct <- census_data$LCS_American/ census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$African_origin_pct <- census_data$African / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$MidEast_WC_Asian_origin_pct <- census_data$MidEast_WC_Asian / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$E_SE_Asian_origin_pct <- census_data$E_SE_Asian / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$European_origin_pct <- census_data$European / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$S_Asian_origin_pct <- census_data$S_Asian / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$Caucasian_origin_pct <- census_data$Caucasian / census_data$`v_CA21_4917: Total - Ethnic or cultural origin for the population in private households`
census_data$citizen_pct <- census_data$citizens/ census_data$`v_CA21_4389: Total - Citizenship for the population in private households`
census_data$noncitizen_pct <- census_data$non_citizens / census_data$`v_CA21_4389: Total - Citizenship for the population in private households`
census_data$nonimmigrant_pct <- census_data$non_immigrants/ census_data$`v_CA21_4404: Total - Immigrant status and period of immigration for the population in private households`
census_data$immigrant_pct <- census_data$immigrants / census_data$`v_CA21_4404: Total - Immigrant status and period of immigration for the population in private households`
census_data$separated_divorced_widowed_pct <- (census_data$`v_CA21_483: Separated` + census_data$`v_CA21_489: Widowed`+census_data$`v_CA21_486: Divorced`)/census_data$`v_CA21_453: Marital status for the total population aged 15 years and over`
census_data$employed_pct <- census_data$`v_CA21_6498: Employed`/census_data$`v_CA21_6495: In the labour force`
census_data$no_diploma_pct <- census_data$`v_CA21_5802: No high school diploma or equivalency certificate`/census_data$`v_CA21_5799: Total - Secondary (high) school diploma or equivalency certificate for the population aged 15 years and over in private households`
census_data$live_alone_pct <- census_data$`v_CA21_534: Living alone`/census_data$`v_CA21_510: Persons in private households`
census_data$elderly_pct = census_data$`v_CA21_251: 65 years and over`/census_data$`v_CA21_8: Total - Age`
census_data <- census_data[, -c(2, 4,6,7,9, 11:13)]  # Drops columns 2, 5, and 7

MSDI <- read_excel(here("data","raw", "INDQ_MSDI_Canada", "A-MSDIData_Can2021_en", "1. EquivalenceTableCanada2021_en.xlsx")) ### read about the indices here: https://www.inspq.qc.ca/sites/default/files/2024-04/3476-material-social-deprivation-index-guide-2021.pdf

# Filter the dataframe
MSDI_yvr <- MSDI %>% 
  filter(MUNIC == "5915022")

MSDI_yvr <- MSDI_yvr %>% mutate(DA = as.character(DA))

census_data <- census_data %>%
  left_join(MSDI_yvr %>% select(DA, SCOREMAT, SCORESOC), by = c("name" = "DA"))


census_data <- st_transform(census_data, st_crs(dat_grid)) ### ensure the same CRS



# Perform a spatial join
dat_grid <- st_join(dat_grid, census_data)



### save in file, indicating contents (food waste + demographic vars)
save(dat_grid, file = here("data","dat_fw_demo_750.Rdata"))
} else {
  ## object called "dat_grid"
  load(file = here("data", "dat_fw_demo_750.Rdata"))
}


##############################################################################################
############ Module 4: import, clean, and grid-ify Metro Van land cover data ###########################

### expect this module to take ~3 minutes to run

if (gen_lc) {
file_path <- here("data","raw","MetroVan_LCC_2020_final.tif")

# Force reading as a raster with assumed dimensions
lcc2020 <- rast(file_path)

lcc2020 <- project(lcc2020, "EPSG:3857")
van_ext<- c(-13715667.2184339, -13702478.8456052, 6313777.37423673, 6325346.85547016)
lcc2020_van <- crop(lcc2020, van_ext)
lcc2020_van <- as.int(lcc2020_van)

lc_prop <- exact_extract(lcc2020_van, dat_grid, fun = "frac")
#### values from Report (No #14, Snow/Ice) https://metrovancouver.org/services/regional-planning/Documents/mv-land-cover-classification-sei-update-2022.pdf
colnames(lc_prop) <- c("Buildings", "Paved", "OtherBuilt", "Barren", "Soil","Conifer","Deciduous","Shrub","ModGrassHerb","NatGrassHerb","NonphotoVeg","Water","Shadow")

lc_prop$grid_id <- dat_grid$grid_id

### check for identical duplicates (yes, all 1s returned)
lc_prop %>%
  group_by(grid_id) %>%
  filter(n() > 1) %>%
  summarise(across(where(is.numeric), ~ length(unique(.)))) %>%
  summarise(across(everything(), max))

# Check if rows sum to 1
row_sums <- rowSums(lc_prop[, 1:13], na.rm = TRUE)
# Quick summary of sums
summary(row_sums)
# Optional: flag any rows that don’t sum to (almost) 1
which(abs(row_sums - 1) > 1e-6)

lc_prop_clean <- lc_prop %>%
  distinct(grid_id, .keep_all = TRUE)


#lc_prop$grid_id <- 1:nrow(lc_prop)
dat_grid <- left_join(dat_grid, lc_prop_clean, by = "grid_id")
rm(lcc2020,lcc2020_van, lc_prop)
### save file
save(dat_grid, file = here("data","dat_fw_demo_lc_750.Rdata"))

} else {
  ## object called "dat_grid"
  load(file = here("data", "dat_fw_demo_lc_750.Rdata"))
}




#########################################################################################################
############## Module 5: read in food retail location data #################


if (gen_foodretail) {

### Storefronts! ####
### Read storefronts file, metadata available here: https://opendata.vancouver.ca/explore/dataset/storefronts-inventory/information/
stores <- fromJSON(here("data","raw" ,"storefronts-inventory.json"))

### Extract lat/lon from the nested column
stores$longitude <- stores$geo_point_2d$lon
stores$latitude <- stores$geo_point_2d$lat
### Convert to sf object
stores_sf <- st_as_sf(stores, coords = c("longitude", "latitude"), crs = 4326)
### filter to only food related businesses
unique(stores_sf$retail_category)
categories <- c("Food & Beverage", "Convenience Goods", "Entertainment and Leisure")

# Filter the sf dataframe
food_stores <- stores_sf %>%
  filter(retail_category %in% categories)
rm(stores,stores_sf)

# Check results
print(food_stores)

ggplot(data = food_stores) +
  geom_sf(aes(color = retail_category, fill = retail_category), 
          shape = 21, size = 3, alpha = 0.6) +  # Shape 21 allows fill + border
  scale_fill_viridis_d() +  # Use a nice color palette
  scale_color_viridis_d() + 
  theme_minimal() +
  labs(title = "Retail Store Categories")


######## food vendors

### Read vendors file, metadata available here: https://opendata.vancouver.ca/explore/dataset/food-vendors/information/
vendors <- fromJSON(here("data","raw" ,"food-vendors.json"))

### Extract lat/lon from the nested column
vendors$longitude <- vendors$geo_point_2d$lon
vendors$latitude <- vendors$geo_point_2d$lat
### Convert to sf object
vendors_sf <- st_as_sf(vendors, coords = c("longitude", "latitude"), crs = 4326)


## object called "dat_grid"
load(file = here("data", "dat_fw_demo_lc_750.Rdata"))

food_stores <- food_stores %>%
  mutate(
    lon = geo_point_2d$lon,
    lat = geo_point_2d$lat
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 3857)

food_stores <- food_stores[, -c(1:4,8,9:11)]

food_stores <- st_transform(food_stores, st_crs(dat_grid))



vendors_sf <- vendors_sf %>%
  mutate(
    lon = geo_point_2d$lon,
    lat = geo_point_2d$lat
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 3857)
vendors_sf <- st_transform(vendors_sf, st_crs(dat_grid))

vendors_sf <- vendors_sf %>%
  rename(retail_category = vendor_type)

vendors_sf <- vendors_sf[, -c(1,3,5:9,12)]
vendors_sf$year_recorded <- "2020"

food_stores <- bind_rows(food_stores, vendors_sf)


## metadata available here: https://opendata.vancouver.ca/explore/dataset/free-and-low-cost-food-programs/table/?sort=program_name

lowcost_food <- fromJSON(here("data","raw" ,"free-and-low-cost-food-programs.json"))
lowcost_food <- lowcost_food[!is.na(lowcost_food$latitude),]
### Convert to sf object
lowcost_food_sf <- st_as_sf(lowcost_food, coords = c("longitude", "latitude"), crs = 4326)
lowcost_food_sf <- lowcost_food_sf %>%
  mutate(lon = st_coordinates(.)[,1],  # X = longitude
         lat = st_coordinates(.)[,2])  # Y = latitude
lowcost_food_sf <- st_transform(lowcost_food_sf, st_crs(dat_grid))

lowcost_food_sf$retail_category <- "free_lowcost_food"
lowcost_food_sf <- lowcost_food_sf %>%
  rename(business_name = program_name)

lowcost_food_sf  <- lowcost_food_sf  %>%
  mutate(year = substr(as.character(last_update_date), 1, 4)) %>%  # Extract year
  select(-last_update_date)  # Optionally drop original column

lowcost_food_sf  <- lowcost_food_sf[, c(1,24:28)]
lowcost_food_sf$retail_category <- "free_lowcost_food"
lowcost_food_sf <- lowcost_food_sf %>%
  rename(year_recorded = year)

food_stores <- bind_rows(food_stores, lowcost_food_sf)

### read in restaurants and fast food options obtained via quickOSM QGIS plugin, amenities: restaurants, fast food
restaurants_CoV <- st_read(here("data","raw","restaurants_fastfood_CoV_quickOSM.geojson"))

restaurants_CoV <- restaurants_CoV %>% 
  select(amenity, shop, name)

restaurants_CoV <- restaurants_CoV %>% 
  rename(retail_category =amenity,
         business_name = name)
restaurants_CoV <- restaurants_CoV %>%
  mutate(retail_category = ifelse(is.na(retail_category), shop, retail_category)) %>%
  select(-shop)


restaurants_CoV <- restaurants_CoV %>%
  mutate(lon = st_coordinates(.)[,1],  # X = longitude
         lat = st_coordinates(.)[,2]) 

restaurants_CoV$year_recorded <- "2025"


food_stores <- bind_rows(food_stores, restaurants_CoV)



# ggplot() +
#   geom_sf(data = food_stores, aes(color = retail_category, fill = retail_category), 
#           shape = 21, size = 3, alpha = 0.6) +  # Shape 21 allows fill + border
#   scale_fill_viridis_d() +  # Use a nice color palette
#   scale_color_viridis_d() + 
#   geom_sf(data = sites_3857, color = "orange", fill = NA) +
#   theme_minimal() +
#   labs(title = "Retail Store Categories")


# Perform spatial join: Assigns food store points to grid cells
joined <- st_join(dat_grid, food_stores, left = TRUE) 




# Convert to data.table
dt_joined <- as.data.table(joined)

# Convert `retail_category` to a factor for faster comparison
dt_joined[, retail_category := as.factor(retail_category)]

retail_categories <- unique(food_stores$retail_category)
retail_categories <- retail_categories[!retail_categories %in% c("NA", NA)]


# # Enable progress bar
handlers(global = TRUE)
dat_grid_final <- with_progress({
  p <- progressor(along = unique(dt_joined$grid_id))  # Progress for unique grid IDs
  
  # Identify grid cells containing at least one relevant food retail location
  food_grids <- dt_joined[, .(
    Food_Retail = as.integer(any(retail_category %in% retail_categories, na.rm = TRUE))
  ), by = .(grid_id)]
  
  p()  # Update progress
  
  # Merge back with original data to retain all columns
  result <- merge(dt_joined, food_grids, by = "grid_id", all.x = TRUE)
  
  st_as_sf(result)  # Convert back to sf
  
})

dat_grid_final <- dat_grid_final[,-c(143:147)]

#### gen nearest food retailer column
dat_grid_final <- dat_grid_final %>%
  mutate(nearest_food_retail = apply(st_distance(geometry, food_stores), 1, min))

### test plot
ggplot(dat_grid_final) +
  geom_sf(aes(fill = factor(Food_Retail)), color = NA) +  # Remove borders
  scale_fill_manual(values = c("0" = "gray90", "1" = "red")) +
  theme_minimal() +
  labs(title = "Food & Beverage Locations", fill = "Presence")





save(dat_grid_final, file = here("data","dat_fw_demo_lc_retail_750.Rdata"))



} else {
  ## object called "dat_grid_final"
  load(file = here("data", "dat_fw_demo_lc_retail_750.Rdata"))
}


#########################################################################################################
############## Module 6: read in nearest homeless shelter #################

if (gen_shelters) {
### data sourced from: https://opendata.vancouver.ca/explore/dataset/homeless-shelter-locations/information/
shelters <- st_read(here("data","raw","homeless-shelter-locations"))

shelters <- shelters %>%
  mutate(lon = st_coordinates(.)[,1],  # X = longitude
         lat = st_coordinates(.)[,2])  # Y = latitude
shelters <- st_transform(shelters, st_crs(dat_grid_final))


shelters$year <- "2025"

shelters <- shelters[, c(1,9:11)]

# Compute the minimum distance from each grid cell to the nearest shelter
dat_grid_final <- dat_grid_final %>%
  mutate(nearest_shelter_dist = apply(st_distance(geometry, shelters), 1, min))
save(dat_grid_final, file = here("data","dat_fw_demo_lc_retail_shelters_750.Rdata"))



} else {
  ## object called "dat_grid_final"
  load(file = here("data", "dat_fw_demo_lc_retail_shelters_750.Rdata"))
}


if (gen_roads) {
  ### Set parameters
  roads <- st_read(here("data","raw","british-columbia-latest-free.shp","gis_osm_roads_free_1.shp")) #### download here: https://download.geofabrik.de/north-america/canada/british-columbia.html
  roads <- st_transform(roads, st_crs(dat_grid_final)) ### metadata here: https://download.geofabrik.de/osm-data-in-gis-formats-free.pdf
  roads_cropped <- st_intersection(roads, dat_grid_final)
  #plot(roads_cropped$geometry)
  #roads_service <- roads_cropped[roads_cropped$fclass == c("service","unclassified", "cycleway") ]
  roads_service <- roads_cropped[roads_cropped$fclass %in% c("service", "unclassified", "cycleway"), ]
  
  # Assign weights based on road type
  roads_service$road_weight <- dplyr::case_when(
    roads_service$fclass %in% c("service", "unclassified") ~ 1,
    roads_service$fclass == "cycleway" ~ 0.25
  )
  
  
  
  #test <- roads_cropped[roads_cropped$fclass == "cycleway", ]
  plot(roads_service$geometry)
  #plot(test$geometry)
  
  roads_in_grid <- st_intersection(roads_service, dat_grid_final)
  # Then calculate lengths in meters
  roads_in_grid$length_m <- st_length(roads_in_grid)
  lengths_by_grid <- roads_in_grid %>%
    st_drop_geometry() %>%
    group_by(grid_id) %>%
    summarise(service_road_length_m = sum(as.numeric(length_m)))
  
  dat_grid_final <- dat_grid_final %>%
    left_join(lengths_by_grid, by = "grid_id")
  
  # If no roads intersected a grid cell, it will have NA — replace with 0:
  dat_grid_final$service_road_length_m[is.na(dat_grid_final$service_road_length_m)] <- 0
  
  
  # Get centroids of each grid cell to measure from center
  grid_centroids <- st_centroid(dat_grid_final)
  
  # Compute distance from centroids to all service roads
  dist_to_road <- st_distance(grid_centroids, roads_service)
  
  # Find index of the nearest road for each centroid
  min_dist_index <- apply(dist_to_road, 1, which.min)
  
  # Get the minimum distance values
  min_dist_m <- apply(dist_to_road, 1, min)
  
  # Store minimum distances in your data
  dat_grid_final$dist_to_service_road_m <- as.numeric(min_dist_m)
  
  # Identify the road type (fclass) of the nearest road for each centroid
  nearest_fclass <- roads_service$fclass[min_dist_index]
  
  # Compute distance-based weight
  dist_weight <- case_when(
    dat_grid_final$dist_to_service_road_m <= 3 ~ 1,
    dat_grid_final$dist_to_service_road_m <= 20 ~ {
      d <- (dat_grid_final$dist_to_service_road_m - 3) / (20 - 3)
      (1 - d^2)
    },
    TRUE ~ 0
  )
  
  # Modify by fclass: cycleway = 0.33 multiplier
  road_type_multiplier <- ifelse(nearest_fclass == "cycleway", 0.33, 1)
  
  # Final weighted score
  dat_grid_final$service_road_weight <- dist_weight * road_type_multiplier
  
  dat_grid_final$service_road_bool <- case_when(
    dat_grid_final$dist_to_service_road_m <= 25 ~ 1,
    TRUE ~ 0
  )
  
  save(dat_grid_final, file = here("data","dat_fw_demo_lc_retail_shelters_roads_750.Rdata"))
  
} else {
  ## object called "dat_grid_final"
  load(file = here("data", "dat_fw_demo_lc_retail_shelters_roads_750.Rdata"))
}
  
  


#########################################################################################################
############## Module 6: Make final cleaning modifications for units (e.g., percentages), missing data, sums, etc.  #################

if (gen_finalmods) {

dat_grid_final$NaturalVegSum2 <- dat_grid_final$Conifer +dat_grid_final$Deciduous+dat_grid_final$Shrub +dat_grid_final$NatGrassHerb +dat_grid_final$Soil +dat_grid_final$Barren +dat_grid_final$NonphotoVeg+dat_grid_final$OtherBuilt
dat_grid_final$NaturalVegSum <- dat_grid_final$Conifer +dat_grid_final$Deciduous+dat_grid_final$Shrub +dat_grid_final$NatGrassHerb
dat_grid_final$BuiltSum <- dat_grid_final$Paved+dat_grid_final$Buildings+dat_grid_final$OtherBuilt
dat_grid_final$GrassSoilSum <- dat_grid_final$ModGrassHerb + dat_grid_final$Soil
dat_grid_final$housing_lowdensity <- (dat_grid_final$single_detached + dat_grid_final$semi_detached + dat_grid_final$row_houses)/dat_grid_final$`v_CA21_434: Occupied private dwellings by structural type of dwelling data` ### rowhouses may need to be removed, because Vancouver doesn't uniformly service ALL of these
dat_grid_final$housing_highdensity <- (dat_grid_final$apt_flat + dat_grid_final$low_rise + dat_grid_final$high_rise)/dat_grid_final$`v_CA21_434: Occupied private dwellings by structural type of dwelling data`
dat_grid_final$housing_suitability <- dat_grid_final$`v_CA21_4261: Suitable`/dat_grid_final$`v_CA21_4260: Total - Private households by housing suitability`
dat_grid_final$pct_unsuitable_housing <- 1 - dat_grid_final$housing_suitability
dat_grid_final$pct_unsuitable_housing[dat_grid_final$pct_unsuitable_housing == 1] <- 0 ### replaces Pacific Spirit and Ocean near 2nd beach with 0s (i.e., no housing at all)
dat_grid_final$pct_rent_thirty <- dat_grid_final$`v_CA21_4315: % of tenant households spending 30% or more of its income on shelter costs (55)`/100
dat_grid_final$pct_rent_subsidized <- dat_grid_final$`v_CA21_4314: % of tenant households in subsidized housing (61)`/100
#dat_grid_final$poorservice_housing <- dat_grid_final$rent_pct + dat_grid_final$pct_housing_notsuitable + dat_grid_final$pct_rent_thirty
dat_grid_final$rent_subsidized_interaction <- dat_grid_final$pct_rent_subsidized*dat_grid_final$rent_pct
dat_grid_final$pct_minor_repairs <- dat_grid_final$`v_CA21_4273: Only regular maintenance and minor repairs needed`/dat_grid_final$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
dat_grid_final$pct_major_repairs <- dat_grid_final$`v_CA21_4274: Major repairs needed`/dat_grid_final$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
dat_grid_final$total_repairs <- (dat_grid_final$`v_CA21_4273: Only regular maintenance and minor repairs needed` + dat_grid_final$`v_CA21_4274: Major repairs needed`)/dat_grid_final$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
dat_grid_final$one_parent_pct <- dat_grid_final$`v_CA21_507: Total one-parent families`/dat_grid_final$`v_CA21_499: Total number of census families in private households`
dat_grid_final$avg_household_size <- dat_grid_final$`v_CA21_452: Average household size`
# logit <- function(p) {
#   log(p / (1 - p))
# }

#dat_grid_final$service_road_weight_logit <- logit(dat_grid_final$service_road_weight)

dat_grid_final <- dat_grid_final %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

dat_grid_final <- dat_grid_final %>%
  select(-name) ### remove "name" column. No use for it and is NA for all grid cells in pacific spirit ("lancaster_ub") transect. Causes problems.



save(dat_grid_final, file = here("data","dat_fw_demo_lc_retail_shelters_roads_mods_750.Rdata"))


} else {
  ## object called "dat_grid_final"
  load(file = here("data", "dat_fw_demo_lc_retail_shelters_roads_mods_750.Rdata"))
}


if (gen_center_std) {

  
dat_clean <- na.omit(dat_grid_final)
# Extract geometry and data separately
geom <- st_geometry(dat_clean)
data_only <- st_drop_geometry(dat_clean)
# Identify columns to scale (exclude ID columns etc.)
cols_to_scale <- setdiff(names(data_only), c("site", "grid_id"))  # adjust as needed
data_only <- as.data.frame(st_drop_geometry(dat_clean))
# Scale those columns safely
data_only[cols_to_scale] <- lapply(data_only[cols_to_scale], function(x) as.numeric(scale(x)))
# Recombine scaled data with geometry
dat_scaled <- st_sf(data_only, geometry = st_geometry(dat_clean))
dat_scaled$site_id <- as.integer(factor(dat_scaled$site))
save(dat_scaled, file = here("data", "dat_scaled_750.Rdata"))

dat_no_ub <- dat_clean %>%
  filter(site != "lancaster_ub")
geom <- st_geometry(dat_no_ub)
data_only <- st_drop_geometry(dat_no_ub)
# Identify columns to scale (exclude ID columns etc.)
cols_to_scale <- setdiff(names(data_only), c("site", "grid_id"))  # adjust as needed
data_only <- as.data.frame(st_drop_geometry(dat_no_ub))
# Scale those columns safely
data_only[cols_to_scale] <- lapply(data_only[cols_to_scale], function(x) as.numeric(scale(x)))
# Recombine scaled data with geometry
dat_no_ub_scaled <- st_sf(data_only, geometry = st_geometry(dat_no_ub))
dat_no_ub_scaled$site_id <- as.integer(factor(dat_no_ub_scaled$site))
save(dat_no_ub_scaled, file = here("data", "dat_no_ub_scaled_750.Rdata"))
save(dat_no_ub, file = here("data", "dat_no_ub_750.Rdata"))


} else {
  ## object called "dat_scaled"
  load(file = here("data", "dat_no_ub_scaled_750.Rdata"))
}


if (gen_factanal) {
## Factor Analysis of Housing Vars

my_data <- dat_no_ub_scaled %>% select(where(is.numeric))
my_data <- as.data.frame(my_data)
my_data <- my_data[, colSums(is.na(my_data)) == 0]  # Remove NA columns



# Define the variables to keep in the heatmap
selected_vars <- c("rent_pct","nearest_shelter_dist","nearest_food_retail", "pop_km","household_income","employed_pct" ,"no_diploma_pct","separated_divorced_widowed_pct",
                   "female_pct","one_parent_pct", "avg_household_size","housing_highdensity",
                   "pct_major_repairs", "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing")


# # Define the variables to keep in the heatmap
# selected_vars <- c("rent_pct","nearest_shelter_dist","nearest_food_retail", "pop_km","household_income","employed_pct" ,"no_diploma_pct","separated_divorced_widowed_pct",
#                    "female_pct","one_parent_pct", "avg_household_size","housing_highdensity",
#                    "pct_major_repairs", "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing", "v_CA21_386: Average age")

# # Define the variables to keep in the heatmap
# selected_vars <- c("rent_pct","nearest_shelter_dist","nearest_food_retail", "pop_km","household_income","employed_pct" ,"no_diploma_pct","separated_divorced_widowed_pct",
#                    "female_pct","one_parent_pct", "avg_household_size","housing_highdensity",
#                    "pct_major_repairs", "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing", "elderly_pct")




# Ensure the selected variables exist in the dataset
selected_vars <- selected_vars[selected_vars %in% colnames(my_data)]

# Subset data to only include selected variables
my_data_subset <- my_data[, selected_vars]

# Ensure no missing values
my_data_subset <- na.omit(my_data_subset)

# Compute eigenvalues of the correlation matrix
eigen_vals <- eigen(cor(my_data_subset))$values

# Scree plot
plot(eigen_vals, type = "b", pch = 19,
     main = "Scree Plot",
     xlab = "Factor Number", ylab = "Eigenvalue")
abline(h = 1, col = "red", lty = 2)  # Kaiser criterion line


# Set number of factors to extract
## due to Kaiser Criterion: https://psycnet.apa.org/fulltext/2016-15750-001.html
num_factors <- 3
#num_factors <- 4 # for avg age/eldery

# Perform factor analysis
fa_result <- factanal(my_data_subset, factors = num_factors, rotation = "varimax", scores = "regression")


loadings_matrix <- fa_result$loadings  # or fa_model$weights for regression method


fa_scores <-as.data.frame(fa_result$scores)

# View summary of the result
print(fa_result)

# Extract factor loadings
loadings_df <- data.frame(fa_result$loadings[, 1:num_factors])
loadings_df$Variable <- rownames(loadings_df)
loadings_melted <- melt(loadings_df, id.vars = "Variable", variable.name = "Factor", value.name = "Loading")

# Plot factor loadings
fa_loadings_plot <- ggplot(loadings_melted, aes(x = Variable, y = Loading, fill = Factor)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Factor Loadings", x = "Variable", y = "Loading") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(fa_loadings_plot)



# Convert to data frame for plotting
loadings_df <- as.data.frame(unclass(fa_result$loadings))
loadings_df$Variable <- rownames(loadings_df)

# Reshape to long format manually
loadings_long <- reshape(loadings_df, 
                         varying = names(loadings_df)[1:num_factors], 
                         v.names = "Loading",
                         timevar = "Factor",
                         times = paste0("Factor", 1:num_factors),
                         direction = "long")

# Plot (using base R barplot)
par(mar = c(10, 4, 4, 2))  # increase bottom margin
with(loadings_long, {
  barplot(Loading, names.arg = Variable, las = 2,
          main = "Factor Loadings", col = as.factor(Factor))
})

# Extract and scale loadings
loadings <- fa_result$loadings[, 1:3]


# Convert to data frame
loadings_df <- as.data.frame(loadings)
loadings_df$Variable <- rownames(loadings_df)

# Optional: scale loadings for visibility
scale_factor <- 1.5
loadings_df$Factor1 <- loadings_df[,1] * scale_factor
loadings_df$Factor2 <- loadings_df[,2] * scale_factor
loadings_df$Factor3 <- loadings_df[,3] * scale_factor
#loadings_df$Factor4 <- loadings_df[,4] * scale_factor

# Save plot to PDF
pdf(here("figs","FA12.pdf"), width = 10, height = 8)


# Plot
ggplot(loadings_df, aes(x = 0, y = 0, xend = Factor1, yend = Factor2)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text_repel(aes(x = Factor1, y = Factor2, label = Variable), size = 4) +
  coord_equal() +
  xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  labs(title = "Factor Analysis Loadings Biplot",
       x = "Factor 1", y = "Factor 2") +
  theme_minimal()

dev.off()

# Save plot to PDF
pdf(here("figs","FA13.pdf"), width = 10, height = 8)

# Plot
ggplot(loadings_df, aes(x = 0, y = 0, xend = Factor1, yend = Factor3)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text_repel(aes(x = Factor1, y = Factor3, label = Variable), size = 4) +
  coord_equal() +
  xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  labs(title = "Factor Analysis Loadings Biplot",
       x = "Factor 1", y = "Factor 3") +
  theme_minimal()

dev.off()


# Save plot to PDF
pdf(here("figs","FA23.pdf"), width = 10, height = 8)

# Plot
ggplot(loadings_df, aes(x = 0, y = 0, xend = Factor2, yend = Factor3)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text_repel(aes(x = Factor2, y = Factor3, label = Variable), size = 4) +
  coord_equal() +
  xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  labs(title = "Factor Analysis Loadings Biplot",
       x = "Factor 2", y = "Factor 3") +
  theme_minimal()

dev.off()

# 
# # Save plot to PDF
# pdf(here("figs","FA14.pdf"), width = 10, height = 8)
# 
# # Plot
# ggplot(loadings_df, aes(x = 0, y = 0, xend = Factor1, yend = Factor4)) +
#   geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = "red") +
#   geom_text_repel(aes(x = Factor1, y = Factor4, label = Variable), size = 4) +
#   coord_equal() +
#   xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
#   labs(title = "Factor Analysis Loadings Biplot",
#        x = "Factor 1", y = "Factor 3") +
#   theme_minimal()
# 
# dev.off()

dat_final <- cbind(dat_no_ub_scaled, fa_scores)
save(dat_final, file = here("data", "dat_final_750.Rdata"))
save(loadings_matrix, file = here("data", "loadings_matrix.Rdata"))

} else {
  ## object called "dat_final"
  load(file = here("data", "dat_final_750.Rdata"))
}


