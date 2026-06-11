## load, compile, and clean data for food waste modelling

## Written by Daniel Forrest
## Started on February 27, 2025

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
library(FNN)
library(reshape2)
library(ggrepel)
library(gridExtra)
library(spatstat)



#### The following script is divided into modules, depending on the data type and intended output, 
#### so as to reduce run time, repeat operations, and organize by source data
#### If you are running for the first time, set all of the following parameters to TRUE
#### If you are running to alter a portion of the data compilation code,
#### set all stages prior to the focal stage to FALSE

### Module 1: gen grids
gen_grids <- FALSE  # Change to FALSE to load data instead
### Module 2: gen food waste
gen_foodwaste <- FALSE  # Change to FALSE to load data instead
### Module 3: gen Canadian census dat (demographics)
gen_demog <- FALSE # Change to FALSE to load data instead
### Module 4: read in land cover data, calculate proportion of cell containing each of 13 land cover classes
gen_lc <- FALSE
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
### Module 10: conduct a factor analysis
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
  
  # 1. Make a bbox polygon in WGS84
  bbox_wgs84 <- st_bbox(c(xmin = -123.2247, ymin = 49.1982, 
                          xmax = -123.0235, ymax = 49.3167), 
                        crs = st_crs(4326))
  
  bbox_poly_wgs84 <- st_as_sfc(bbox_wgs84)
  
  # 2. Transform to EPSG:3857
  bbox_poly_3857 <- st_transform(bbox_poly_wgs84, crs = 3857)
  
  # 3. Get the new bounding box
  bbox_3857 <- st_bbox(bbox_poly_3857)
  
  # Convert bbox to an sf polygon
  bbox_vancouver_poly <- st_as_sfc(bbox_3857)
  
  
  
  # Now use the bbox polygon with your grid function
  grid_vancouver <- create_internal_grid(
    polygon = bbox_vancouver_poly, 
    cell_size = cell_size,        # adjust this for your desired grid size (~degrees)
    site_name = "Vancouver", 
    overlap_buffer = 0        # optional: here we set buffer to 0 since bbox is exact
  )
  
  # 
  # ggplot() +
  #   geom_sf(data = grid_vancouver, fill = NA, color = "blue") +
  #   geom_sf(data = bbox_vancouver_poly, fill = NA, color = "red") +
  #   theme_minimal() +
  #   labs(title = "Grid Over Vancouver Bounding Box")
  # 
  # 
  # # ✅ Calculate cell area in square meters
  # cell_areas <- st_area(grid_vancouver)
  # 
  # # Print mean cell area (if needed)
  # mean_cell_area <- mean(cell_areas)
  # print(paste("Mean cell size:", mean_cell_area, "square meters"))
  # 
  # 
  # # Define output file name and resolution
  # pdf(here("figs","grid_plot_750m_vancouver.pdf"), width = 80, height = 64)  # Adjust size as needed
  # 
  # # Plot the grid inside polygons
  # plot(st_geometry(bbox_vancouver_poly), col = "lightblue", border = "black")
  # plot(st_geometry(grid_vancouver), add = TRUE, border = "red")
  # 
  # # Close the PDF device to save the file
  # dev.off()
  
  st_crs(grid_vancouver) <- 3857
  # Add a unique ID column
  grid_vancouver <- grid_vancouver %>%
    mutate(grid_id = row_number())
  
  save(grid_vancouver, file = here("data","grid_750m_vancouver.Rdata"))
  #######
} else {
  ## file called "final_grid"
  load(file = here("data","grid_750m_vancouver.Rdata"))
}




if (gen_demog) {
  
  ### insert your own Google Maps key here
  options(cancensus.api_key = 'CensusMapper_e5809a95267db864f07fd91906dfdc24') ### insert personalized API Key, sign-up here: https://censusmapper.ca/users/sign_up
  
  
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
                                        "v_CA21_4968" #caucasian
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
  census_data <- census_data[, -c(2, 4,6,7,9, 11:13)]  # Drops columns 2, 5, and 7
  
  MSDI <- read_excel(here("data","raw", "INDQ_MSDI_Canada", "A-MSDIData_Can2021_en", "1. EquivalenceTableCanada2021_en.xlsx")) ### read about the indices here: https://www.inspq.qc.ca/sites/default/files/2024-04/3476-material-social-deprivation-index-guide-2021.pdf
  
  # Filter the dataframe
  MSDI_yvr <- MSDI %>% 
    filter(MUNIC == "5915022")
  
  MSDI_yvr <- MSDI_yvr %>% mutate(DA = as.character(DA))
  
  census_data <- census_data %>%
    left_join(MSDI_yvr %>% select(DA, SCOREMAT, SCORESOC), by = c("name" = "DA"))
  
  
  census_data <- st_transform(census_data, st_crs(grid_vancouver)) ### ensure the same CRS
  
  # Loop through all columns with missing values in the census_data dataset
  for (col in names(census_data)) {
    # Check if the column contains any NAs
    if (any(is.na(census_data[[col]]))) {
      # Calculate the median (ignoring NAs) and replace NAs with the median
      median_value <- median(census_data[[col]], na.rm = TRUE)
      census_data[[col]][is.na(census_data[[col]])] <- median_value
    }
  }
  
  
  # #### test for spatial autocorrelation (yes!)
  # 
  # # Create neighbors list using Queen contiguity
  # nb <- poly2nb(census_data[!is.na(census_data$employed_pct),], queen = TRUE)
  # 
  # # Create spatial weights (row-standardized)
  # lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  # 
  # # Moran’s I test
  # moran_test <- moran.test(census_data$employed_pct[!is.na(census_data$employed_pct)], lw, zero.policy = TRUE)
  # 
  # # Print results
  # print(moran_test)
  # 
  # 
  
  # Perform a spatial join
  van_dat <- st_join(grid_vancouver, census_data)
  
  
  
  # van_dat <- van_dat %>%
  #   filter(!is.na(van_dat$household_income),
  #          !is.na(van_dat$employed_pct),
  #          !is.na(van_dat$female_pct),
  #          !is.na(van_dat$`v_CA21_507: Total one-parent families`),
  #          !is.na(van_dat$no_diploma_pct),
  #          !is.na(van_dat$`v_CA21_4314: % of tenant households in subsidized housing (61)`),
  #          !is.na(van_dat$pop_km),
  #          !is.na(van_dat$`v_CA21_4274: Major repairs needed`),
  #          !is.na(van_dat$`v_CA21_4260: Total - Private households by housing suitability`),
  #          !is.na(van_dat$rent_pct),
  #          !is.na(van_dat$`v_CA21_4315: % of tenant households spending 30% or more of its income on shelter costs (55)`),
  #          !is.na(van_dat$`v_CA21_452: Average household size`),
  #          !is.na(van_dat$separated_divorced_widowed_pct))  # Removes rows empty of census data
  # 
  # 
  ### save in file, indicating contents (food waste + demographic vars)
  save(van_dat, file = here("data","van_dat_750.Rdata"))
} else {
  ## object called "dat_grid"
  load(file = here("data","van_dat_750.Rdata"))
}


##############################################################################################
############ Module 3: import, clean, and grid-ify Metro Van land cover data ###########################

### expect this module to take ~3 minutes to run

if (gen_lc) {
  file_path <- here("data","raw","MetroVan_LCC_2020_final.tif")
  
  # Force reading as a raster with assumed dimensions
  lcc2020 <- rast(file_path)
  
  lcc2020 <- project(lcc2020, "EPSG:3857")
  van_ext <- ext(project(vect(van_dat), "EPSG:3857"))
  #van_ext<- c(-13715667.2184339, -13702478.8456052, 6313777.37423673, 6325346.85547016)
  lcc2020_van <- crop(lcc2020, van_ext)
  lcc2020_van <- as.int(lcc2020_van)
  
  lc_prop <- exact_extract(lcc2020_van, van_dat, fun = "frac")
  #### values from Report (No #14, Snow/Ice) https://metrovancouver.org/services/regional-planning/Documents/mv-land-cover-classification-sei-update-2022.pdf
  colnames(lc_prop) <- c("Buildings", "Paved", "OtherBuilt", "Barren", "Soil","Conifer","Deciduous","Shrub","ModGrassHerb","NatGrassHerb","NonphotoVeg","Water","Shadow")
  
  lc_prop$grid_id <- van_dat$grid_id
  
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
  van_dat <- left_join(van_dat, lc_prop_clean, by = "grid_id")
  rm(lcc2020,lcc2020_van, lc_prop)
  
  van_dat <- van_dat %>%
    filter(!is.na(van_dat$rent_pct)) ### filter that removes cells without any census data
  
  
  port_shape <- st_read(here("data", "port_removal.kml"))
  port_shape <- st_transform(port_shape, crs = st_crs(van_dat))
  
  # Check intersection: returns logical vector
  intersects <- st_intersects(van_dat, port_shape, sparse = FALSE)[, 1]
  sum(intersects)
  
  plot(st_geometry(van_dat), col = "grey")
  plot(st_geometry(port_shape), border = "red", add = TRUE)
  
  
  # Remove features that intersect the drawn shape
  van_dat <- van_dat[!intersects, ]
  
  ### save file
  save(van_dat, file = here("data","van_dat_demo_lc_750.Rdata"))
  
} else {
  ## object called "dat_grid"
  load(file = here("data", "van_dat_demo_lc_750.Rdata"))
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
  
  
  food_stores <- food_stores %>%
    mutate(
      lon = geo_point_2d$lon,
      lat = geo_point_2d$lat
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 3857)
  
  food_stores <- food_stores[, -c(1:4,8,9:11)]
  
  food_stores <- st_transform(food_stores, st_crs(van_dat))
  
  
  
  vendors_sf <- vendors_sf %>%
    mutate(
      lon = geo_point_2d$lon,
      lat = geo_point_2d$lat
    ) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 3857)
  vendors_sf <- st_transform(vendors_sf, st_crs(van_dat))
  
  vendors_sf <- vendors_sf %>%
    rename(retail_category = vendor_type)
  
  vendors_sf <- vendors_sf[, -c(1,3,5:9,12)]
  vendors_sf$year_recorded <- "2020"
  
  food_stores <- bind_rows(food_stores, vendors_sf)
  
  
  
  lowcost_food <- fromJSON(here("data","raw" ,"free-and-low-cost-food-programs.json"))
  lowcost_food <- lowcost_food[!is.na(lowcost_food$latitude),]
  ### Convert to sf object
  lowcost_food_sf <- st_as_sf(lowcost_food, coords = c("longitude", "latitude"), crs = 4326)
  lowcost_food_sf <- lowcost_food_sf %>%
    mutate(lon = st_coordinates(.)[,1],  # X = longitude
           lat = st_coordinates(.)[,2])  # Y = latitude
  lowcost_food_sf <- st_transform(lowcost_food_sf, st_crs(van_dat))
  
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
  
  
  restaurants_CoV <- st_read(here("data","restaurants_fastfood_CoV_quickOSM.geojson"))
  
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
  joined <- st_join(van_dat, food_stores, left = TRUE) 
  
  
  
  
  # Convert to data.table
  dt_joined <- as.data.table(joined)
  
  # Convert `retail_category` to a factor for faster comparison
  dt_joined[, retail_category := as.factor(retail_category)]
  
  retail_categories <- unique(food_stores$retail_category)
  retail_categories <- retail_categories[!retail_categories %in% c("NA", NA)]
  
  
  # # Enable progress bar
  handlers(global = TRUE)
  van_dat <- with_progress({
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
  
  van_dat <- van_dat[,-c(140:144)]

  
  # Ensure both datasets are in WGS84
  van_dat <- st_transform(van_dat, 4326)
  food_stores <- st_transform(food_stores, 4326)
  
  van_centroids <- st_centroid(van_dat)
  
  # Then extract coordinates
  van_coords <- st_coordinates(van_centroids)[, 1:2]
  food_coords <- st_coordinates(food_stores)[, 1:2]
  
  # Use FNN for fast nearest neighbor
  nn_result <- get.knnx(data = food_coords, query = van_coords, k = 1)
  
  # Add nearest distance to your data (in meters)
  van_dat$nearest_food_retail <- nn_result$nn.dist[,1]
  
  
  
  ### test plot
  ggplot(van_dat) +
    geom_sf(aes(fill = factor(Food_Retail)), color = NA) +  # Remove borders
    scale_fill_manual(values = c("0" = "gray90", "1" = "red")) +
    theme_minimal() +
    labs(title = "Food & Beverage Locations", fill = "Presence")
  
  
  van_dat <- st_transform(van_dat, 3857)
  
  
  save(van_dat, file = here("data","van_dat_demo_lc_retail_750.Rdata"))
  
  
  
} else {
  ## object called "van_dat"
  load(file = here("data", "van_dat_demo_lc_retail_750.Rdata"))
}


#########################################################################################################
############## Module 6: read in nearest homeless shelter #################

if (gen_shelters) {
  ### data sourced from: https://opendata.vancouver.ca/explore/dataset/homeless-shelter-locations/information/
  shelters <- st_read(here("data","raw","homeless-shelter-locations"))
  
  shelters <- shelters %>%
    mutate(lon = st_coordinates(.)[,1],  # X = longitude
           lat = st_coordinates(.)[,2])  # Y = latitude
  shelters <- st_transform(shelters, st_crs(van_dat))
  
  
  shelters$year <- "2025"
  
  shelters <- shelters[, c(1,9:11)]
  
  
  # Compute the minimum distance from each grid cell to the nearest shelter
  van_dat <- van_dat %>%
    mutate(nearest_shelter_dist = apply(st_distance(geometry, shelters), 1, min))
  
  
  save(van_dat, file = here("data","van_dat_demo_lc_retail_shelters_750.Rdata"))
  
  
  
} else {
  ## object called "van_dat"
  load(file = here("data", "van_dat_demo_lc_retail_shelters_750.Rdata"))
}


if (gen_roads) {
  ### Set parameters
  roads <- st_read(here("data","raw","british-columbia-latest-free.shp","gis_osm_roads_free_1.shp")) #### download here: https://download.geofabrik.de/north-america/canada/british-columbia.html
  roads <- st_transform(roads, st_crs(van_dat)) ### metadata here: https://download.geofabrik.de/osm-data-in-gis-formats-free.pdf
  roads_cropped <- st_intersection(roads, van_dat)
  #plot(roads_cropped$geometry)
  #roads_service <- roads_cropped[roads_cropped$fclass == c("service","unclassified", "cycleway") ]
  roads_service <- roads_cropped[roads_cropped$fclass %in% c("service", "unclassified", "cycleway"), ]
  
  # Assign weights based on road type
  roads_service$road_weight <- dplyr::case_when(
    roads_service$fclass %in% c("service", "unclassified") ~ 1,
    roads_service$fclass == "cycleway" ~ 0.25
  )
  
  
  
  #test <- roads_cropped[roads_cropped$fclass == "cycleway", ]
  #plot(roads_service$geometry)
  #plot(test$geometry)
  
  roads_in_grid <- st_intersection(roads_service, van_dat)
  # Then calculate lengths in meters
  roads_in_grid$length_m <- st_length(roads_in_grid)
  lengths_by_grid <- roads_in_grid %>%
    st_drop_geometry() %>%
    group_by(grid_id) %>%
    summarise(service_road_length_m = sum(as.numeric(length_m)))
  
  van_dat <- van_dat %>%
    left_join(lengths_by_grid, by = "grid_id")
  
  # If no roads intersected a grid cell, it will have NA — replace with 0:
  van_dat$service_road_length_m[is.na(van_dat$service_road_length_m)] <- 0
  # 1. Sample midpoints along roads, cast to LINESTRING if needed
  roads_service_lines <- roads_service[
    st_geometry_type(roads_service) %in% c("LINESTRING", "MULTILINESTRING"), ]
  roads_service_lines <- st_cast(roads_service_lines, "LINESTRING")
  road_points <- st_line_sample(roads_service_lines, sample = 0.5)
  road_points_sf <- st_sf(geometry = road_points)
  road_points_sf$fclass <- roads_service_lines$fclass  # carry over fclass
  
  # 2. Ensure everything is projected (EPSG:3857 for meter-based distance)
  van_dat <- st_transform(van_dat, 3857)
  road_points_sf <- st_transform(road_points_sf, 3857)
  
  # 3. Get coordinates
  van_coords <- st_coordinates(st_centroid(van_dat))[, 1:2]
  road_coords <- st_coordinates(road_points_sf)[, 1:2]
  
  # 4. Use FNN to find nearest road point and distance
  nn_result <- get.knnx(data = road_coords, query = van_coords, k = 1)
  
  # 5. Add nearest distance to data
  van_dat$dist_to_service_road_m <- nn_result$nn.dist[, 1]
  
  # 6. Get index of nearest road and extract road type (fclass)
  nearest_idx <- nn_result$nn.index[, 1]
  nearest_fclass <- road_points_sf$fclass[nearest_idx]
  
  # 7. Compute distance-based weight
  dist_weight <- case_when(
    van_dat$dist_to_service_road_m <= 3 ~ 1,
    van_dat$dist_to_service_road_m <= 20 ~ {
      d <- (van_dat$dist_to_service_road_m - 3) / (20 - 3)
      (1 - d^2)
    },
    TRUE ~ 0
  )
  
  # 8. Apply road type multiplier
  road_type_multiplier <- ifelse(nearest_fclass == "cycleway", 0.33, 1)
  
  # 9. Final weighted score
  van_dat$service_road_weight <- dist_weight * road_type_multiplier
  
  # 10. Binary flag for within 25m
  van_dat$service_road_bool <- ifelse(van_dat$dist_to_service_road_m <= 25, 1, 0)
  
  save(van_dat, file = here("data","van_dat_demo_lc_retail_shelters_roads_750.Rdata"))
  
} else {
  ## object called "van_dat"
  load(file = here("data", "van_dat_demo_lc_retail_shelters_roads_750.Rdata"))
}




#########################################################################################################
############## Module 6: read in nearest homeless shelter #################

if (gen_finalmods) {
  
  van_dat$NaturalVegSum2 <- van_dat$Conifer +van_dat$Deciduous+van_dat$Shrub +van_dat$NatGrassHerb +van_dat$Soil +van_dat$Barren +van_dat$NonphotoVeg+van_dat$OtherBuilt
  van_dat$NaturalVegSum <- van_dat$Conifer +van_dat$Deciduous+van_dat$Shrub +van_dat$NatGrassHerb
  van_dat$BuiltSum <- van_dat$Paved+van_dat$Buildings+van_dat$OtherBuilt
  van_dat$GrassSoilSum <- van_dat$ModGrassHerb + van_dat$Soil
  van_dat$housing_lowdensity <- (van_dat$single_detached + van_dat$semi_detached + van_dat$row_houses)/van_dat$`v_CA21_434: Occupied private dwellings by structural type of dwelling data` ### rowhouses may need to be removed, because vancouver doesn't uniformly service ALL of these
  van_dat$housing_highdensity <- (van_dat$apt_flat + van_dat$low_rise + van_dat$high_rise)/van_dat$`v_CA21_434: Occupied private dwellings by structural type of dwelling data`
  van_dat$housing_suitability <- van_dat$`v_CA21_4261: Suitable`/van_dat$`v_CA21_4260: Total - Private households by housing suitability`
  van_dat$pct_unsuitable_housing <- 1 - van_dat$housing_suitability
  van_dat$pct_unsuitable_housing[van_dat$pct_unsuitable_housing == 1] <- 0 ### replaces Pacific Spirit and Ocean near 2nd beach with 0s (i.e., no housing at all)
  van_dat$pct_rent_thirty <- van_dat$`v_CA21_4315: % of tenant households spending 30% or more of its income on shelter costs (55)`/100
  van_dat$pct_rent_subsidized <- van_dat$`v_CA21_4314: % of tenant households in subsidized housing (61)`/100
  #van_dat$poorservice_housing <- van_dat$rent_pct + van_dat$pct_housing_notsuitable + van_dat$pct_rent_thirty
  van_dat$rent_subsidized_interaction <- van_dat$pct_rent_subsidized*van_dat$rent_pct
  van_dat$pct_minor_repairs <- van_dat$`v_CA21_4273: Only regular maintenance and minor repairs needed`/van_dat$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
  van_dat$pct_major_repairs <- van_dat$`v_CA21_4274: Major repairs needed`/van_dat$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
  van_dat$total_repairs <- (van_dat$`v_CA21_4273: Only regular maintenance and minor repairs needed` + van_dat$`v_CA21_4274: Major repairs needed`)/van_dat$`v_CA21_4272: Total - Occupied private dwellings by dwelling condition`
  van_dat$one_parent_pct <- van_dat$`v_CA21_507: Total one-parent families`/van_dat$`v_CA21_499: Total number of census families in private households`
  van_dat$avg_household_size <- van_dat$`v_CA21_452: Average household size`
  # logit <- function(p) {
  #   log(p / (1 - p))
  # }
  
  #van_dat$service_road_weight_logit <- logit(van_dat$service_road_weight)
  
  van_dat <- van_dat %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  van_dat <- van_dat %>% select(-name.y, -business_name, -retail_category)
  
  save(van_dat, file = here("data","van_dat_demo_lc_retail_shelters_roads_mods_750.Rdata"))
  
  
} else {
  ## object called "van_dat"
  load(file = here("data", "van_dat_demo_lc_retail_shelters_roads_mods_750.Rdata"))
}


if (gen_center_std) {
  
  
  van_dat_clean <- na.omit(van_dat)
  # Extract geometry and data separately
  geom <- st_geometry(van_dat_clean)
  data_only <- st_drop_geometry(van_dat_clean)
  # Identify columns to scale (exclude ID columns etc.)
  cols_to_scale <- sapply(data_only, is.numeric)
  data_only <- as.data.frame(st_drop_geometry(van_dat_clean))
  # Scale those columns safely
  data_only[cols_to_scale] <- lapply(data_only[cols_to_scale], function(x) as.numeric(scale(x)))
  # Recombine scaled data with geometry
  van_dat_scaled <- st_sf(data_only, geometry = st_geometry(van_dat_clean))
  save(van_dat_scaled, file = here("data", "van_dat_scaled_750.Rdata"))
  
  
  
} else {
  ## object called "dat_scaled"
  load(file = here("data", "van_dat_scaled_750.Rdata"))
}



if (gen_factanal) {
  ## Factor Analysis of Housing Vars
  
  load(file = here("data", "loadings_matrix.Rdata"))
  loadings_matrix <- unclass(loadings_matrix)
  
  # Define the variables to keep in the heatmap
  selected_vars <- c("rent_pct","nearest_shelter_dist","nearest_food_retail", "pop_km","household_income","employed_pct" ,"no_diploma_pct","separated_divorced_widowed_pct",
                     "female_pct","one_parent_pct", "avg_household_size","housing_highdensity",
                     "pct_major_repairs", "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing")
  
  
  
  
  # Ensure the selected variables exist in the dataset
  van_dat_factors <- van_dat_scaled %>% select(any_of(selected_vars))
  van_dat_matrix <- as.matrix(van_dat_factors)
  van_dat_matrix <- van_dat_matrix[,1:16]
  storage.mode(van_dat_matrix) <- "double"  # Just to be sure
  
  # Then project:
  factor_scores <- van_dat_matrix %*% loadings_matrix
  
  van_dat <- cbind(van_dat_scaled, factor_scores)
  
  winsorize_base <- function(x, lower = 0.01, upper = 0.99) {
    q <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
    x_wins <- pmin(pmax(x, q[1]), q[2])
    return(x_wins)
  }
  
  # Winsorize Factor3 in van_dat at 1% and 99% levels
  van_dat$Factor3_wins <- winsorize_base(van_dat$Factor3, lower = 0.001, upper = 0.999)
  # Winsorize Factor3 in van_dat at 1% and 99% levels
  van_dat$Factor1_wins <- winsorize_base(van_dat$Factor1, lower = 0.001, upper = 0.999)
  
  
  save(van_dat, file = here("data", "van_dat_final_750.Rdata"))
  
  
} else {
  ## object called "dat_final"
  load(file = here("data", "van_dat_final_750.Rdata"))
}



if (gen_factanal2) {
  ## Factor Analysis of Housing Vars
  
  load(file = here("data", "loadings_matrix_bycollection.Rdata"))
  loadings_matrix <- unclass(loadings_matrix)
  
  # Define the variables to keep in the heatmap
  selected_vars <- c("rent_pct","nearest_shelter_dist","nearest_food_retail", "pop_km","household_income","employed_pct" ,"no_diploma_pct","separated_divorced_widowed_pct",
                     "female_pct","one_parent_pct", "avg_household_size","housing_highdensity",
                     "pct_major_repairs", "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing")
  
  
  
  
  # Ensure the selected variables exist in the dataset
  van_dat_factors <- van_dat_scaled %>% select(any_of(selected_vars))
  van_dat_matrix <- as.matrix(van_dat_factors)
  van_dat_matrix <- van_dat_matrix[,1:16]
  storage.mode(van_dat_matrix) <- "double"  # Just to be sure
  
  # Then project:
  factor_scores <- van_dat_matrix %*% loadings_matrix
  
  van_dat <- cbind(van_dat_scaled, factor_scores)
  
  winsorize_base <- function(x, lower = 0.01, upper = 0.99) {
    q <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
    x_wins <- pmin(pmax(x, q[1]), q[2])
    return(x_wins)
  }
  
  # Winsorize Factor3 in van_dat at 1% and 99% levels
  van_dat$Factor3_wins <- winsorize_base(van_dat$Factor3, lower = 0.001, upper = 0.999)
  # Winsorize Factor3 in van_dat at 1% and 99% levels
  van_dat$Factor2_wins <- winsorize_base(van_dat$Factor2, lower = 0.001, upper = 0.999)
  van_dat$Factor1_wins <- winsorize_base(van_dat$Factor1, lower = 0.001, upper = 0.999)
  
  save(van_dat, file = here("data", "van_dat_final_750_2.Rdata"))
  
  
} else {
  ## object called "dat_final"
  load(file = here("data", "van_dat_final_750_2.Rdata"))
}

