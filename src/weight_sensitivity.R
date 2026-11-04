## weight_sensitivity.R
## Do the F123_all results depend on the receptacle-form weights in the food waste score?
## (Current weights: dumpster 1.0, stationary 0.2, other 0.066, tote 0.05,
##  home_bin 0.033, plastic_bag 0.033, i.e. a dumpster counts ~30x a home bin.)
##
## Rebuilds the score from the saved observations (data/fw_observations_750.Rdata, written by
## bin_state_by_deprivation.R), changing ONLY the form weight:
##   additive_score = 0.5 * state + form + type + litter term   (compile_data_adjustfw.R)
##   so new additive = old additive - old form weight + new form weight.
## Collection-schedule adjustment and aggregation are unchanged.
##
##   original : current weights (reproduces fw_score_weighted exactly; checked)
##   equal    : every receptacle form gets the same weight (0.05); no dumpster premium
##   older    : the commented-out weights in compile_data_adjustfw.R (dumpster 15, tote 3,
##              home bin 1, plastic bag 1, other 2, stationary 0.266), divided by 15
## Food litter has no form weight, so its score is identical in all three.
##
## Refits F123_all (fw_hier_3pred_ZIW_4pred_bern.stan) for "equal" and "older" and prints the
## coefficients next to the original F123_all fit.
## Run from the proxy_to_process project root, in a fresh R session. About 5-10 minutes.

library(here)
library(sf)
library(dplyr)
library(data.table)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(here("data", "fw_observations_750.Rdata"))   # obs
load(here("data", "dat_no_ub_750.Rdata"))          # dat_no_ub

## ---------------------------------------------------------------
## 1. Form-weight schemes
## ---------------------------------------------------------------
forms <- c("dumpster", "tote", "home_bin", "plastic_bag", "other", "stationary")
schemes <- list(
  original = c(dumpster = 1.0, tote = 0.05, home_bin = 0.033, plastic_bag = 0.033,
               other = 0.066, stationary = 0.2),
  equal    = setNames(rep(0.05, length(forms)), forms),
  older    = c(dumpster = 15, tote = 3, home_bin = 1, plastic_bag = 1,
               other = 2, stationary = 0.266) / 15
)
print(round(do.call(rbind, schemes), 4))

form_w <- function(w) {
  out <- unname(w[match(obs$bin_form, names(w))])          # NA for unmatched forms, as in the original
  out[!is.na(obs$bin_state) & obs$bin_state == "food_litter"] <- 0
  out
}

## ---------------------------------------------------------------
## 2. Rebuild cell scores under each scheme
## ---------------------------------------------------------------
cnt    <- fifelse(is.na(obs$bin_count), 1, obs$bin_count)
add0   <- obs$weighted_value / cnt                                 # original additive score
adjfac <- fifelse(obs$weighted_value != 0, obs$adjusted_weighted_value / obs$weighted_value, 0)
f0     <- form_w(schemes$original)

cell_score <- function(w) {
  adj <- cnt * (add0 - f0 + form_w(w)) * adjfac
  s <- data.table(grid_id = obs$grid_id, v = adj)[, .(score = sum(v, na.rm = TRUE)), by = grid_id]
  out <- s$score[match(dat_no_ub$grid_id, s$grid_id)]
  out[is.na(out)] <- 0
  out
}
scores <- lapply(schemes, cell_score)

chk <- max(abs(scores$original - dat_no_ub$fw_score_weighted))
cat("\nCheck, original weights vs fw_score_weighted: max abs difference =", signif(chk, 3), "\n")
if (chk > 1e-6) stop("Rebuilt original scores don't match fw_score_weighted; stopping.")

cat("\nCorrelation of cell scores across schemes:\n")
print(round(cor(do.call(cbind, scores)), 3))
cat("\nShare of total score from dumpsters:\n")
dump <- !is.na(obs$bin_form) & obs$bin_form == "dumpster"
print(round(sapply(schemes, function(w) {
  adj <- cnt * (add0 - f0 + form_w(w)) * adjfac
  sum(adj[dump], na.rm = TRUE) / sum(adj, na.rm = TRUE)
}), 3))

## ---------------------------------------------------------------
## 3. Refit F123_all under the alternative schemes
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
pars <- c("beta1", "beta2", "beta3", "delta1", "delta2", "delta3", "delta4", "alpha", "tau_mu")
lab  <- c(beta1 = "amount: Housing", beta2 = "amount: Deprivation", beta3 = "amount: Enfranchisement",
          delta1 = "presence: service roads", delta2 = "presence: Deprivation",
          delta3 = "presence: Housing", delta4 = "presence: Enfranchisement",
          alpha = "Weibull shape", tau_mu = "site SD")

ci <- function(s) sprintf("%6.2f [%5.2f, %5.2f]", s[, "50%"], s[, "2.5%"], s[, "97.5%"])

## Original F123_all, from the saved fit
env <- new.env()
load(here("data", "fit_hier_ZIW_16var_F123_all.Rdata"), envir = env)
s_orig <- summary(get(ls(env)[1], envir = env), pars = pars, probs = c(0.025, 0.5, 0.975))$summary
rm(env); gc()

res <- list(original = s_orig)
for (sch in c("equal", "older")) {
  cat("\nFitting F123_all with", sch, "form weights...\n")
  fit <- stan(here("src", "fw_hier_3pred_ZIW_4pred_bern.stan"),
              data = c(base, list(fw_score_weighted = scores[[sch]])),
              chains = 4, iter = 2000, seed = 2026, refresh = 0)
  res[[sch]] <- summary(fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  cat("Divergences:", sum(get_divergent_iterations(fit)),
      "| max Rhat:", round(max(res[[sch]][, "Rhat"]), 3), "\n")
  rm(fit); gc()
}

cat("\n===== F123_all coefficients by form-weight scheme: median [95% interval] =====\n")
tab <- data.frame(parameter = lab[pars],
                  original = ci(res$original[pars, ]),
                  equal    = ci(res$equal[pars, ]),
                  older    = ci(res$older[pars, ]),
                  row.names = NULL)
print(tab, right = FALSE)

save(scores, res, file = here("data", "weight_sensitivity_F123_all.Rdata"))
