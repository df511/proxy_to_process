## run_revised_analysis.R
## Run order for the revised (November 2026) analysis: corrected factor analysis and city
## scoring, the model fits, city predictions, the population-density mask, diagnostics, and the
## manuscript figures and tables.
##
## Each step runs in its own environment. The model-fitting and prediction steps use a lot of
## memory; on a laptop it is safest to run this script one section at a time and restart R
## between sections (Session > Restart R). Switch sections on/off with the flags below.
## Run from the proxy_to_process project root.

library(here)
step <- function(file) {
  cat("\n\n#################### ", file, " (", format(Sys.time(), "%H:%M"), ")\n", sep = "")
  source(here("src", file), local = new.env())
  invisible(gc())
}

run_data        <- FALSE  # 1. Rebuild data from raw sources (slow; needs a CensusMapper API key)
run_scoring     <- TRUE   # 2. Factor analysis and city scoring (fast)
run_fits        <- FALSE  # 3. Fit all models (each fit is saved to data/fit_hier_ZIW_16var_*.Rdata)
run_predictions <- TRUE   # 4. LOO comparison, city predictions, mask
run_diagnostics <- FALSE  # 5. Deprivation checks, bins vs litter, bin states, weight sensitivity
run_figures     <- TRUE   # 6. Figures and tables

## 1. Data ------------------------------------------------------------------------------------
if (run_data) {
  step("compile_data_adjustfw.R")       # transect grid: food waste scores, covariates
  step("compile_data_YVRreproject.R")   # city grid; its factor-analysis step runs
                                        # factor_analysis_revision.R (so step 2 is then redundant)
}

## 2. Factor analysis and city scoring ---------------------------------------------------------
if (run_scoring && !run_data) {
  ## City files built before the November 2026 fix have distance to food retail in degrees;
  ## fix_food_retail_distance.R writes a corrected copy (in metres), which
  ## factor_analysis_revision.R uses automatically when needed.
  step("fix_food_retail_distance.R")
  step("factor_analysis_revision.R")    # 3 factors, 16 variables; transect weights and scaling
}

## 3. Model fits -------------------------------------------------------------------------------
## The proxy (Tweedie) model, fit_hier_1pred_tweedie.Rdata, comes from final_model_configs_stan.Rmd
## and does not depend on the factors.
if (run_fits) {
  step("refit_factor_models.R")         # F12 and F123 (then LOO and predictions; heavy)
  step("fit_F13.R")                     # F13
  step("fit_presence_models.R")         # F123_D and F123_all (main model); LOO of all models
}

## 4. LOO, city predictions, mask --------------------------------------------------------------
if (run_predictions) {
  step("resume_after_fits.R")           # LOO (F12, F123, F13, proxy) and city predictions
  step("predict_F123_all.R")            # city predictions for the main model
  step("domain_mask.R")                 # mask DAs with population density < 650 per km2
}

## 5. Diagnostics (supplementary analyses) -----------------------------------------------------
if (run_diagnostics) {
  step("deprivation_checks.R")          # site structure, within-site effect, influential cells
  step("split_bins_litter.R")           # F123_all refit to bins and to food litter separately
  step("bin_state_by_deprivation.R")    # receptacle state and form by deprivation (descriptive)
  step("weight_sensitivity.R")          # alternative receptacle-form weights
}

## 6. Figures and tables -----------------------------------------------------------------------
## All main figures (2-4) come from manuscript_main_figs.Rmd and all supplementary figures
## (S1-S16) and tables from manuscript_supplement_figs.Rmd. In RStudio you can instead open each
## Rmd and use Run All (Ctrl/Cmd + Alt + R).
if (run_figures) {
  for (rmd in c("manuscript_main_figs.Rmd", "manuscript_supplement_figs.Rmd")) {
    cat("\n\n#################### ", rmd, " (", format(Sys.time(), "%H:%M"), ")\n", sep = "")
    rmarkdown::render(here("src", rmd), envir = new.env(), quiet = TRUE)
    invisible(gc())
  }
}

cat("\nDone.\n")
