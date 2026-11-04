## factor_analysis_revision.R
## Recompute the manuscript's factor analysis and score the City of Vancouver grid
## consistently with the transect (training) data.
##
## Fixes two problems in compile_data_adjustfw.R / compile_data_YVRreproject.R:
##   1. City-wide factor scores were computed as Z_city %*% loadings. The transect scores
##      came from factanal(scores = "regression"), whose weights are solve(cor) %*% loadings.
##      Here both are scored with the same regression weights W.
##   2. City covariates were standardized with the city's own means and SDs, but the
##      models were fit on covariates standardized with the transect means and SDs.
##      Here everything is standardized with the transect means and SDs.
## It also fixes the factor analysis itself: the saved dat_final_750.Rdata and
## loadings_matrix.Rdata came from a 4-factor run that included average age and dropped
## female_pct. This script reruns the 3-factor, 16-variable analysis described in the
## manuscript.
##
## Run AFTER src/fix_food_retail_distance.R, from the proxy_to_process project root.
## Outputs:
##   data/factor_analysis_16var.Rdata     fa_result, Lambda, W, train_means, train_sds, fa_vars
##   data/factor_loadings_16var.csv       loadings table (for the manuscript)
##   data/dat_final_750_16var.Rdata       transect cells with corrected Factor1-3 scores
##   data/van_dat_final_750_16var.Rdata   city cells (van_dat_fa) with corrected scores and
##                                        transect-standardized service_road_weight_z and
##                                        nearest_food_retail_z

library(here)
library(sf)
library(dplyr)

## ---------------------------------------------------------------
## 1. Data (both UNSCALED)
## ---------------------------------------------------------------
load(here("data", "dat_no_ub_750.Rdata"))                                    # dat_no_ub: transect cells
## City cells (van_dat). nearest_food_retail must be in metres. compile_data_YVRreproject.R now
## writes it that way; files built before the Nov 2026 fix have it in degrees, in which case the
## corrected copy written by fix_food_retail_distance.R is used instead.
load(here("data", "van_dat_demo_lc_retail_shelters_roads_mods_750.Rdata"))  # van_dat
if (max(van_dat$nearest_food_retail, na.rm = TRUE) < 1) {
  cat("City nearest_food_retail is in degrees (file built before the fix); using the corrected copy.\n")
  load(here("data", "van_dat_demo_lc_retail_shelters_roads_mods_750_fooddist.Rdata"))
}
stopifnot(max(van_dat$nearest_food_retail, na.rm = TRUE) > 1)

fa_vars <- c("rent_pct", "nearest_shelter_dist", "nearest_food_retail", "pop_km",
             "household_income", "employed_pct", "no_diploma_pct",
             "separated_divorced_widowed_pct", "female_pct", "one_parent_pct",
             "avg_household_size", "housing_highdensity", "pct_major_repairs",
             "pct_rent_thirty", "pct_rent_subsidized", "pct_unsuitable_housing")

stopifnot(all(fa_vars %in% names(dat_no_ub)), all(fa_vars %in% names(van_dat)))

## ---------------------------------------------------------------
## 2. Factor analysis on the transect cells
## ---------------------------------------------------------------
X_train <- as.matrix(st_drop_geometry(dat_no_ub)[, fa_vars])
stopifnot(!anyNA(X_train))

train_means <- colMeans(X_train)
train_sds   <- apply(X_train, 2, sd)
Z_train     <- scale(X_train, center = train_means, scale = train_sds)

eigen_vals <- eigen(cor(Z_train))$values
cat("Eigenvalues:", round(eigen_vals, 2), "\n")
cat("Factors with eigenvalue > 1 (Kaiser):", sum(eigen_vals > 1), "\n")

fa_result <- factanal(Z_train, factors = 3, rotation = "varimax", scores = "regression")
print(fa_result, cutoff = 0.3)

Lambda <- unclass(fa_result$loadings)
W <- solve(cor(Z_train), Lambda)   # regression (Thomson) scoring weights, as factanal uses

scores_train <- Z_train %*% W
colnames(scores_train) <- paste0("Factor", 1:3)
stopifnot(max(abs(scores_train - fa_result$scores)) < 1e-6)  # identical to factanal's scores

write.csv(round(Lambda, 3), here("data", "factor_loadings_16var.csv"))

## ---------------------------------------------------------------
## 3. Score the city cells with the SAME standardization and weights
## ---------------------------------------------------------------
city_df <- st_drop_geometry(van_dat)
keep <- complete.cases(city_df[, c(fa_vars, "service_road_weight")])
cat("City cells:", nrow(van_dat), "| complete for scoring:", sum(keep), "\n")

van_dat_fa <- van_dat[keep, ]
X_city <- as.matrix(st_drop_geometry(van_dat_fa)[, fa_vars])
Z_city <- scale(X_city, center = train_means, scale = train_sds)
scores_city <- Z_city %*% W

for (k in 1:3) van_dat_fa[[paste0("Factor", k)]] <- as.numeric(scores_city[, k])

## Non-factor predictors, standardized with the transect means and SDs
## (matches dat_final$service_road_weight and dat_final$nearest_food_retail used in the fits)
van_dat_fa$service_road_weight_z <-
  (van_dat_fa$service_road_weight - mean(dat_no_ub$service_road_weight)) / sd(dat_no_ub$service_road_weight)
van_dat_fa$nearest_food_retail_z <-
  (van_dat_fa$nearest_food_retail - mean(dat_no_ub$nearest_food_retail)) / sd(dat_no_ub$nearest_food_retail)

## ---------------------------------------------------------------
## 4. Diagnostics
## ---------------------------------------------------------------
cat("\nFactor score correlations, transects (should be ~0):\n")
print(round(cor(scores_train), 2))
cat("\nFactor score correlations, city (corrected scoring):\n")
print(round(cor(scores_city), 2))

## Compare with the old city scores, if that file is present
if (file.exists(here("data", "van_dat_final_750.Rdata"))) {
  old_env <- new.env()
  load(here("data", "van_dat_final_750.Rdata"), envir = old_env)   # van_dat (old scoring)
  old <- st_drop_geometry(old_env$van_dat)[, paste0("Factor", 1:3)]
  cat("\nFactor score correlations, city (old scoring, as mapped):\n")
  print(round(cor(old, use = "complete.obs"), 2))
  rm(old_env)
}

## Extrapolation: how much of the city lies outside the range sampled by the transects?
outside <- sweep(Z_city, 2, apply(Z_train, 2, min), "<") | sweep(Z_city, 2, apply(Z_train, 2, max), ">")
cat("\nShare of city cells outside the transect range, by variable:\n")
print(round(colMeans(outside), 3))
cat("Share of city cells outside the transect range on ANY variable:",
    round(mean(rowSums(outside) > 0), 3), "\n")

rng_train <- apply(scores_train, 2, range)
cat("\nShare of city cells with factor scores outside the transect range:\n")
print(round(sapply(1:3, function(k)
  mean(scores_city[, k] < rng_train[1, k] | scores_city[, k] > rng_train[2, k])), 3))

## ---------------------------------------------------------------
## 5. Save
## ---------------------------------------------------------------
dat_final_16var <- dat_no_ub
for (k in 1:3) dat_final_16var[[paste0("Factor", k)]] <- as.numeric(scores_train[, k])

save(fa_result, Lambda, W, train_means, train_sds, fa_vars,
     file = here("data", "factor_analysis_16var.Rdata"))
save(dat_final_16var, file = here("data", "dat_final_750_16var.Rdata"))
save(van_dat_fa, file = here("data", "van_dat_final_750_16var.Rdata"))
