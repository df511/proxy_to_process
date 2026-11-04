## run_all_fixes.R
## Runs the factor-analysis and city-distance fixes in order:
##   1. fix_food_retail_distance.R   city distance to food retail in metres (was degrees)
##   2. factor_analysis_revision.R   3-factor, 16-variable FA; consistent city scoring
##   3. refit_factor_models.R        F12 and F123 fits, LOO vs proxy, city predictions (~20 min)
##   4. map_prediction_components.R  amount-where-present, 750 m and capped maps
## Stops after step 1 if the recomputed transect distances don't match the originals.
##
## Run from the proxy_to_process project, in a fresh R session (Session > Restart R):
##   source("src/run_all_fixes.R")

library(here)

run_step <- function(n, file) {
  cat("\n\n#################################################\n")
  cat("## Step", n, ":", file, "\n")
  cat("## Started", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("#################################################\n\n")
  source(here("src", file), echo = FALSE)
  cat("\n## Step", n, "finished", format(Sys.time(), "%H:%M:%S"), "\n")
}

## Keep a copy of everything printed, to send back afterwards
log_file <- here("run_all_fixes_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

run_step(1, "fix_food_retail_distance.R")
if (check["max_abs_diff_m", best] > 1) {
  stop("Recomputed transect distances don't match the originals (max difference ",
       round(check["max_abs_diff_m", best], 1),
       " m). Stopping before the factor analysis; send the log to check.")
}

run_step(2, "factor_analysis_revision.R")
run_step(3, "refit_factor_models.R")
run_step(4, "map_prediction_components.R")

cat("\nAll steps finished. Output log:", log_file, "\n")
