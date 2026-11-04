## deprivation_checks.R
## Why does Deprivation lower the amount of food waste (F123 beta2 = -0.18)?
## Quick diagnostic fits, each summarised and then freed from memory:
##
##   A. Site table: site intercepts vs. site-mean deprivation, housing and FW score (no fit)
##   B. Pooled F123: no site intercepts (does the effect depend on the hierarchy?)
##   C. Within-site F123: deprivation centred on each site's mean (within-site gradient only)
##   D. F123 + Deprivation in the zero-inflation part (less waste at all, or smaller amounts?)
##   E. F123 without the most deprived site
##   F. F123 without the most deprived 2.5% of cells
##
## In F123: beta1 = Housing-tenure-density, beta2 = Deprivation, beta3 = Enfranchisement,
##          delta1 = service roads (zero-inflation), delta2 = Deprivation (zero-inflation, D only).
## Needs: data/factor_analysis_16var.Rdata, data/dat_final_750_16var.Rdata,
##        data/dat_no_ub_750.Rdata, data/fit_hier_ZIW_16var_F123.Rdata,
##        data/loo_16var_models.Rdata, src/fw_hier_3pred_ZIW_1pred_bern.stan,
##        src/fw_hier_3pred_ZIW_2pred_bern.stan
## Run from the proxy_to_process project root, in a fresh R session. Roughly 10-30 minutes.

library(here)
library(sf)
library(dplyr)
library(rstan)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## ---------------------------------------------------------------
## Data
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # Lambda
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub

pick <- function(v) which.max(abs(Lambda[v, ]))
fac  <- function(k) dat_final_16var[[paste0("Factor", k)]]
Hous <- fac(pick("housing_highdensity"))
Dep  <- fac(pick("pct_rent_subsidized"))
Enf  <- fac(pick("employed_pct"))

site_f <- factor(dat_no_ub$site)
y  <- dat_no_ub$fw_score_weighted
P1 <- (dat_no_ub$service_road_weight - mean(dat_no_ub$service_road_weight)) /
      sd(dat_no_ub$service_road_weight)

make_data <- function(keep = rep(TRUE, length(y)), X2 = Dep, pooled = FALSE, P2 = NULL) {
  s <- droplevels(site_f[keep])
  d <- list(
    N = sum(keep),
    J = if (pooled) 1L else nlevels(s),
    site_id = if (pooled) rep(1L, sum(keep)) else as.integer(s),
    X1 = Hous[keep], X2 = X2[keep], X3 = Enf[keep],
    P1 = P1[keep],
    fw_score_weighted = y[keep]
  )
  if (!is.null(P2)) d$P2 <- P2[keep]
  d
}

pars_main <- c("mu0", "tau_mu", "alpha", "beta1", "beta2", "beta3", "gamma", "delta1")

## Fit, print a compact summary, optionally compute LOO, then free the fit
fit_check <- function(label, stanfile, data, pars = pars_main, want_loo = FALSE) {
  cat("\n\n=====", label, "=====  (N =", data$N, ", sites =", data$J, ")\n")
  f <- stan(here("src", stanfile), data = data, chains = 4, iter = 2000,
            seed = 2026, refresh = 0)
  s <- summary(f, pars = pars, probs = c(0.025, 0.5, 0.975))$summary
  print(round(s[, c("50%", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
  cat("Divergences:", sum(get_divergent_iterations(f)), "\n")
  out <- list(summary = s)
  if (want_loo) {
    ll <- extract_log_lik(f, merge_chains = FALSE)
    out$loo <- loo(ll, r_eff = relative_eff(exp(ll)))
  }
  rm(f); gc()
  out
}

results <- list()

## ---------------------------------------------------------------
## A. Site table
## ---------------------------------------------------------------
load_fit <- function(file) {
  env <- new.env(); load(here("data", file), envir = env); get(ls(env)[1], envir = env)
}
fit <- load_fit("fit_hier_ZIW_16var_F123.Rdata")
mu <- rstan::extract(fit, pars = "mu")$mu       # draws x sites, same site order as levels(site_f)
rm(fit); gc()

site_tab <- data.frame(
  site          = levels(site_f),
  n_cells       = as.vector(table(site_f)),
  mean_dep      = as.vector(tapply(Dep,  site_f, mean)),
  mean_housing  = as.vector(tapply(Hous, site_f, mean)),
  mean_fw       = as.vector(tapply(y,    site_f, mean)),
  prop_zero     = as.vector(tapply(y == 0, site_f, mean)),
  intercept_med = apply(mu, 2, median),
  intercept_lo  = apply(mu, 2, quantile, 0.025),
  intercept_hi  = apply(mu, 2, quantile, 0.975)
)
cat("\n===== A. Site intercepts vs. site means (sorted by mean deprivation) =====\n")
print(site_tab |> arrange(desc(mean_dep)) |> mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("Correlation, site-mean deprivation vs. site intercept:",
    round(cor(site_tab$mean_dep, site_tab$intercept_med), 2), "\n")
cat("Share of deprivation variance that is between sites:",
    round(var(ave(Dep, site_f)) / var(Dep), 2), "\n")
results$A <- site_tab

## ---------------------------------------------------------------
## B. Pooled: no site intercepts
## ---------------------------------------------------------------
results$B <- fit_check("B. Pooled F123 (no site intercepts)",
                       "fw_hier_3pred_ZIW_1pred_bern.stan", make_data(pooled = TRUE),
                       want_loo = TRUE)

## ---------------------------------------------------------------
## C. Within-site deprivation only
## ---------------------------------------------------------------
Dep_within <- Dep - ave(Dep, site_f)
results$C <- fit_check("C. F123 with deprivation centred within site",
                       "fw_hier_3pred_ZIW_1pred_bern.stan", make_data(X2 = Dep_within))

## ---------------------------------------------------------------
## D. Deprivation also in the zero-inflation (presence) part
## ---------------------------------------------------------------
results$D <- fit_check("D. F123 + Deprivation in zero-inflation (delta2; positive = more zeros)",
                       "fw_hier_3pred_ZIW_2pred_bern.stan", make_data(P2 = Dep),
                       pars = c(pars_main, "delta2"), want_loo = TRUE)

## ---------------------------------------------------------------
## E. Without the most deprived site
## ---------------------------------------------------------------
top_site <- site_tab$site[which.max(site_tab$mean_dep)]
results$E <- fit_check(paste("E. F123 without the most deprived site:", top_site),
                       "fw_hier_3pred_ZIW_1pred_bern.stan",
                       make_data(keep = dat_no_ub$site != top_site))

## ---------------------------------------------------------------
## F. Without the most deprived 2.5% of cells
## ---------------------------------------------------------------
cut <- quantile(Dep, 0.975)
results$F <- fit_check("F. F123 without the most deprived 2.5% of cells",
                       "fw_hier_3pred_ZIW_1pred_bern.stan", make_data(keep = Dep < cut))

## ---------------------------------------------------------------
## Summary
## ---------------------------------------------------------------
load(here("data", "loo_16var_models.Rdata"))   # loos (F12, F123, F13, proxy)
cat("\n\n===== LOO: original F123 vs. pooled (B) and deprivation-in-presence (D) =====\n")
print(loo_compare(list(F123 = loos$F123, B_pooled = results$B$loo, D_dep_presence = results$D$loo)))

beta2_row <- function(r) r$summary["beta2", c("50%", "2.5%", "97.5%")]
cat("\n===== Deprivation effect on amount (beta2), across checks =====\n")
orig <- load_fit("fit_hier_ZIW_16var_F123.Rdata")
orig_b2 <- summary(orig, pars = "beta2", probs = c(0.025, 0.5, 0.975))$summary[1, c("50%", "2.5%", "97.5%")]
rm(orig); gc()
print(round(rbind(
  `F123 (original)`         = orig_b2,
  `B. pooled`               = beta2_row(results$B),
  `C. within-site`          = beta2_row(results$C),
  `D. + dep. in presence`   = beta2_row(results$D),
  `E. without top site`     = beta2_row(results$E),
  `F. without top 2.5%`     = beta2_row(results$F)
), 3))
cat("\nD: Deprivation in zero-inflation (delta2; positive = more likely NO waste):\n")
print(round(results$D$summary["delta2", c("50%", "2.5%", "97.5%")], 3))

results_summaries <- lapply(results[c("B", "C", "D", "E", "F")], function(r) r$summary)
save(site_tab, results_summaries, file = here("data", "deprivation_checks.Rdata"))
