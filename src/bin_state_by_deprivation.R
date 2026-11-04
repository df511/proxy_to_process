## bin_state_by_deprivation.R
## Are receptacles in more deprived (and denser) areas more often contained (closed, locked,
## restricted) and less often open with food? Descriptive tables only, no model fitting.
##
## Re-runs Module 2 of src/compile_data_adjustfw.R (unchanged) to get one row per observation,
## links each observation to its grid cell's factor scores, and tabulates bin state and form
## by quintile of Deprivation and of Housing-tenure-density.
##
## Output: data/fw_observations_750.Rdata (observation table, reusable), printed tables.
## Run from the proxy_to_process project root. A couple of minutes.

library(here)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)

## ---------------------------------------------------------------
## 1. Observations (Module 2), or reuse the saved table
## ---------------------------------------------------------------
obs_file <- here("data", "fw_observations_750.Rdata")
if (file.exists(obs_file)) {
  load(obs_file)   # obs
} else {
  load(here("data", "grid_clipped_750m.Rdata"))
  grid_clipped_750m <- final_grid
  st_crs(grid_clipped_750m) <- 3857
  grid_clipped_750m <- grid_clipped_750m %>% mutate(grid_id = row_number())

  src   <- readLines(here("src", "compile_data_adjustfw.R"))
  start <- grep("^if \\(gen_foodwaste\\) \\{", src)
  stop  <- grep('save\\(dat_grid, file = here\\("data", ?"dat_grid_750.Rdata"\\)\\)', src)
  stopifnot(length(start) == 1, length(stop) == 1, stop > start)
  invisible(capture.output(eval(parse(text = src[(start + 1):(stop - 1)]))))

  obs <- as.data.table(dat_sf_joined)
  obs <- obs[, .(grid_id, site, date, fid, waste_type, bin_form, bin_state, bin_count,
                 weighted_value, adjusted_weighted_value)]
  save(obs, file = obs_file)
}

## ---------------------------------------------------------------
## 2. Link to cell factor scores (modelled cells only)
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # Lambda
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub
pick <- function(v) which.max(abs(Lambda[v, ]))
cell_fac <- data.table(
  grid_id = dat_no_ub$grid_id,
  Dep  = dat_final_16var[[paste0("Factor", pick("pct_rent_subsidized"))]],
  Hous = dat_final_16var[[paste0("Factor", pick("housing_highdensity"))]]
)
## Quintiles over all modelled cells (ntile handles tied values; many cells share census values)
q_lab <- function(x, lo, hi) factor(paste0("Q", x, ifelse(x == 1, lo, ifelse(x == 5, hi, ""))),
                                    levels = paste0("Q", 1:5, c(lo, "", "", "", hi)))
cell_fac[, Dep_q  := q_lab(ntile(Dep, 5),  " (least deprived)", " (most deprived)")]
cell_fac[, Hous_q := q_lab(ntile(Hous, 5), " (dense/rental)", " (affluent/low-density)")]

d <- merge(obs, cell_fac, by = "grid_id")   # drops observations outside modelled cells
d[, n := fifelse(is.na(bin_count), 1, bin_count)]

## Group bin states into accessibility classes
d[, access := fcase(
  bin_state == "food_litter",                                   "food litter",
  bin_state %in% c("open_w_food", "open_inbag_food", "food_residue", "restricted_food"),
                                                                "open/restricted WITH food",
  bin_state %in% c("open_no_food", "open_empty", "open_inbag_nofood"), "open, no food",
  bin_state %in% c("open_unknown", "open_inbag_unknown", "restricted_unknown", "other"),
                                                                "open/other, contents unknown",
  bin_state %in% c("closed", "locked", "restricted_empty", "restricted_nofood"),
                                                                "contained (closed/locked/restricted)",
  default = "missing"
)]
cat("Observations linked to modelled cells:", nrow(d), "(receptacles/patches:", sum(d$n), ")\n")

## ---------------------------------------------------------------
## 3. Tables
## ---------------------------------------------------------------
pct_table <- function(by, what) {
  t <- d[, .(n = sum(n)), by = c(by, what)]
  t[, pct := round(100 * n / sum(n), 1), by = by]
  dcast(t, as.formula(paste(by, "~", what)), value.var = "pct", fill = 0)
}

cat("\n===== Accessibility class (% of receptacles/patches) by DEPRIVATION quintile =====\n")
print(pct_table("Dep_q", "access"))
cat("\n===== Accessibility class by HOUSING-TENURE-DENSITY quintile =====\n")
print(pct_table("Hous_q", "access"))
cat("\n===== Receptacle form by DEPRIVATION quintile =====\n")
print(pct_table("Dep_q", "bin_form"))
cat("\n===== Waste stream by DEPRIVATION quintile =====\n")
print(pct_table("Dep_q", "waste_type"))

## Receptacles per cell and mean score per receptacle
per_cell <- d[, .(items = sum(n), score = sum(adjusted_weighted_value, na.rm = TRUE)),
              by = .(grid_id, Dep_q)]
cat("\n===== Per cell with any record, by DEPRIVATION quintile =====\n")
print(per_cell[, .(cells_with_records = .N,
                   mean_items_per_cell = round(mean(items), 1),
                   mean_score_per_item = round(sum(score) / sum(items), 4)),
               by = Dep_q][order(Dep_q)])
cat("Share of ALL modelled cells with any record, by deprivation quintile:\n")
print(cell_fac[, .(share_with_records = round(mean(grid_id %in% per_cell$grid_id), 3)),
               by = Dep_q][order(Dep_q)])
