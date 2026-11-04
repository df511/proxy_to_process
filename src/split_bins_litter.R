## split_bins_litter.R
## Is the "more frequent, smaller amounts" deprivation pattern driven by waste in bins
## or by loose food litter?
##
## 1. Re-runs ONLY Module 2 (food waste scoring) of src/compile_data_adjustfw.R, unchanged,
##    to get observation-level weighted values (dat_sf_joined).
## 2. Sums them per grid cell separately for food litter (bin_state == "food_litter") and
##    for everything else (bins/receptacles). The original score is a plain sum, so
##    bins + litter must equal fw_score_weighted exactly; the script checks this.
## 3. Fits the all-factors model (fw_hier_3pred_ZIW_4pred_bern.stan) to each part.
##
## Amount part: beta1 = Housing, beta2 = Deprivation, beta3 = Enfranchisement
## Presence part (probability of NO waste; negative = more likely present):
##   delta1 = service roads, delta2 = Deprivation, delta3 = Housing, delta4 = Enfranchisement
##
## Run from the proxy_to_process project root, in a fresh R session. About 10-20 minutes.

library(here)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ---------------------------------------------------------------
## 1. Re-run Module 2 of compile_data_adjustfw.R
## ---------------------------------------------------------------
load(here("data", "grid_clipped_750m.Rdata"))   # final_grid (Module 1 output)
grid_clipped_750m <- final_grid
st_crs(grid_clipped_750m) <- 3857
grid_clipped_750m <- grid_clipped_750m %>% mutate(grid_id = row_number())

src   <- readLines(here("src", "compile_data_adjustfw.R"))
start <- grep("^if \\(gen_foodwaste\\) \\{", src)
stop  <- grep('save\\(dat_grid, file = here\\("data", ?"dat_grid_750.Rdata"\\)\\)', src)
stopifnot(length(start) == 1, length(stop) == 1, stop > start)
cat("Running Module 2: lines", start + 1, "to", stop - 1, "of compile_data_adjustfw.R\n")
eval(parse(text = src[(start + 1):(stop - 1)]))
## Creates dat_sf_joined (one row per observation, with grid_id and adjusted_weighted_value)

## ---------------------------------------------------------------
## 2. Split the score
## ---------------------------------------------------------------
obs <- as.data.table(dat_sf_joined)
obs[, part := fifelse(!is.na(bin_state) & bin_state == "food_litter", "litter", "bins")]
cat("\nObservations by part:\n"); print(obs[, .N, by = part])

split <- obs[, .(score = sum(adjusted_weighted_value, na.rm = TRUE)), by = .(grid_id, part)]
split <- dcast(split, grid_id ~ part, value.var = "score", fill = 0)

load(here("data", "dat_no_ub_750.Rdata"))   # dat_no_ub (the modelled cells)
cells <- data.frame(grid_id = dat_no_ub$grid_id) %>%
  left_join(as.data.frame(split), by = "grid_id") %>%
  mutate(across(c(bins, litter), ~ replace_na(.x, 0)))

check <- max(abs(cells$bins + cells$litter - dat_no_ub$fw_score_weighted))
cat("\nMax |bins + litter - fw_score_weighted|:", signif(check, 3), "\n")
if (check > 1e-6) stop("Split scores don't add up to the original score; stopping.")

cat("\nShare of total score from litter:",
    round(sum(cells$litter) / sum(dat_no_ub$fw_score_weighted), 3), "\n")
cat("Cells with any bin score:", sum(cells$bins > 0),
    "| any litter score:", sum(cells$litter > 0), "of", nrow(cells), "\n")
save(cells, file = here("data", "fw_score_bins_litter_750.Rdata"))

## ---------------------------------------------------------------
## 3. Fit the all-factors model to each part
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # Lambda
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var
pick <- function(v) which.max(abs(Lambda[v, ]))
fac  <- function(k) dat_final_16var[[paste0("Factor", k)]]
Hous <- fac(pick("housing_highdensity"))
Dep  <- fac(pick("pct_rent_subsidized"))
Enf  <- fac(pick("employed_pct"))

base <- list(
  N = nrow(dat_no_ub),
  J = length(unique(dat_no_ub$site)),
  site_id = as.integer(factor(dat_no_ub$site)),
  X1 = Hous, X2 = Dep, X3 = Enf,
  P1 = (dat_no_ub$service_road_weight - mean(dat_no_ub$service_road_weight)) /
       sd(dat_no_ub$service_road_weight),
  P2 = Dep, P3 = Hous, P4 = Enf
)
pars <- c("mu0", "tau_mu", "alpha", "beta1", "beta2", "beta3",
          "gamma", "delta1", "delta2", "delta3", "delta4")

summaries <- list()
for (part in c("bins", "litter")) {
  cat("\n\n=====", toupper(part), "=====\n")
  fit <- stan(here("src", "fw_hier_3pred_ZIW_4pred_bern.stan"),
              data = c(base, list(fw_score_weighted = cells[[part]])),
              chains = 4, iter = 2000, seed = 2026, refresh = 0)
  s <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  print(round(s[, c("50%", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
  cat("Divergences:", sum(get_divergent_iterations(fit)), "\n")
  summaries[[part]] <- s
  rm(fit); gc()
}

cat("\n\n===== Deprivation, bins vs litter =====\n")
print(round(rbind(
  `amount (beta2), bins`     = summaries$bins["beta2",   c("50%", "2.5%", "97.5%")],
  `amount (beta2), litter`   = summaries$litter["beta2", c("50%", "2.5%", "97.5%")],
  `presence (delta2), bins`  = summaries$bins["delta2",   c("50%", "2.5%", "97.5%")],
  `presence (delta2), litter`= summaries$litter["delta2", c("50%", "2.5%", "97.5%")]
), 3))
save(summaries, file = here("data", "bins_litter_fits_summary.Rdata"))
