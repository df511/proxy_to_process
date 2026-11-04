## supplement_figures_F123_all.R
## Supplementary figures and tables for the main model, F123_all (all three scales in both the
## presence and quantity parts). Uses the same orientation as Figure 3:
##   Housing-tenure-density : high = more renters, lower income, denser   (fitted factor x -1)
##   Deprivation            : high = more social and material deprivation (as fitted)
##   Disenfranchisement     : high = less employment, lower female share  (fitted factor x -1)
##   Presence effects are on the log-odds of presence (the model's absence coefficients x -1).
##
## Figures (replace the superseded versions in manuscript_supplement_figs.Rmd):
##   figs/FigureS2.pdf   factor-analysis loading biplots (corrected 16-variable analysis)
##   figs/FigureS5.pdf   prior predictive simulations for F123_all (observed covariates)
##   figs/FigureS6.png   posterior predictive checks for F123_all
## Tables (tables/):
##   Table_LOO_all_models.csv          LOO comparison of all models
##   Table_F123_all_coefficients.csv   F123_all posterior summaries (oriented)
##   Table_bins_litter.csv             F123_all refit to bins and to food litter separately
##   Table_weight_sensitivity.csv      F123_all under alternative receptacle-form weights
##   Table_factor_loadings.csv         oriented loadings
## Run from the proxy_to_process project root, in a fresh R session. A few minutes.

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(bayesplot)
library(rstan)

dir.create(here("tables"), showWarnings = FALSE)
dir.create(here("figs"), showWarnings = FALSE)

## ---------------------------------------------------------------
## Common inputs and orientation
## ---------------------------------------------------------------
load(here("data", "factor_analysis_16var.Rdata"))   # fa_result, Lambda, W, ...
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var
load(here("data", "dat_no_ub_750.Rdata"))           # dat_no_ub

pick <- function(v) unname(which.max(abs(Lambda[v, ])))
k <- c(housing = pick("housing_highdensity"),
       deprivation = pick("pct_rent_subsidized"),
       disenfranchisement = pick("employed_pct"))
orient <- c(housing = -1, deprivation = 1, disenfranchisement = -1)
scale_name <- c(housing = "Housing-tenure—density", deprivation = "Deprivation",
                disenfranchisement = "Disenfranchisement")

## Parameters of F123_all, what they are, and the multiplier that gives the reported orientation
par_info <- tribble(
  ~par,     ~component, ~term,                         ~mult,
  "beta1",  "Quantity", "Housing-tenure—density",  -1,
  "beta2",  "Quantity", "Deprivation",                   1,
  "beta3",  "Quantity", "Disenfranchisement",           -1,
  "delta1", "Presence", "Service roads",                -1,
  "delta3", "Presence", "Housing-tenure—density",   1,
  "delta2", "Presence", "Deprivation",                  -1,
  "delta4", "Presence", "Disenfranchisement",            1,
  "mu0",    "Quantity", "Intercept (log Weibull scale)",  1,
  "tau_mu", "Quantity", "Site SD of intercepts",          1,
  "alpha",  "Quantity", "Weibull shape",                  1,
  "gamma",  "Presence", "Intercept (log-odds of presence)", -1
)

## Re-express a rstan/loo-style summary matrix: multiply by mult (swapping the interval when < 0)
orient_summary <- function(s, info = par_info) {
  info <- info[info$par %in% rownames(s), ]
  out <- lapply(seq_len(nrow(info)), function(i) {
    r <- s[info$par[i], ]; m <- info$mult[i]
    lo <- r[["2.5%"]] * m; hi <- r[["97.5%"]] * m
    data.frame(component = info$component[i], term = info$term[i], parameter = info$par[i],
               median = round(r[["50%"]] * m, 3),
               lower_95 = round(min(lo, hi), 3), upper_95 = round(max(lo, hi), 3),
               n_eff = if ("n_eff" %in% names(r)) round(r[["n_eff"]]) else NA,
               Rhat = if ("Rhat" %in% names(r)) round(r[["Rhat"]], 3) else NA)
  })
  do.call(rbind, out)
}

load_fit <- function(file) {
  env <- new.env(); load(here("data", file), envir = env); get(ls(env)[1], envir = env)
}

## =========================================================
## Figure S2: loading biplots (oriented)
## =========================================================
L <- sweep(Lambda[, k], 2, orient, "*")
colnames(L) <- names(k)
write.csv(cbind(variable = rownames(L), as.data.frame(round(L, 3)) |>
                  setNames(scale_name[names(k)])),
          here("tables", "Table_factor_loadings.csv"), row.names = FALSE)

pretty_var <- c(
  rent_pct = "% renters", nearest_shelter_dist = "Distance to shelter",
  nearest_food_retail = "Distance to food retail", pop_km = "Population density",
  household_income = "Household income", employed_pct = "% employed",
  no_diploma_pct = "% no diploma", separated_divorced_widowed_pct = "% separated/divorced/widowed",
  female_pct = "% female", one_parent_pct = "% one-parent families",
  avg_household_size = "Household size", housing_highdensity = "High-density housing",
  pct_major_repairs = "% needing major repairs", pct_rent_thirty = "% rent > 30% of income",
  pct_rent_subsidized = "% subsidized rent", pct_unsuitable_housing = "% unsuitable housing"
)
ld <- data.frame(L, variable = unname(pretty_var[rownames(L)]))

biplot <- function(xs, ys) {
  ggplot(ld, aes(x = 0, y = 0, xend = .data[[xs]], yend = .data[[ys]])) +
    geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
    geom_segment(arrow = arrow(length = unit(0.15, "cm")), color = "firebrick") +
    geom_text_repel(aes(x = .data[[xs]], y = .data[[ys]], label = variable),
                    size = 2.6, max.overlaps = Inf, min.segment.length = 0) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(x = scale_name[[xs]], y = scale_name[[ys]]) +
    theme_minimal(base_size = 10)
}
figS2 <- (biplot("housing", "deprivation") | biplot("housing", "disenfranchisement") |
            biplot("deprivation", "disenfranchisement")) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")
ggsave(here("figs", "FigureS2.pdf"), figS2, width = 14, height = 5.2)

## =========================================================
## Figure S5: prior predictive simulations for F123_all
## Priors as in fw_hier_3pred_ZIW_4pred_bern.stan; covariates are the observed transect values.
## =========================================================
fac <- function(kk) dat_final_16var[[paste0("Factor", kk)]]
Xh <- fac(k[["housing"]]); Xd <- fac(k[["deprivation"]]); Xe <- fac(k[["disenfranchisement"]])
P1 <- (dat_no_ub$service_road_weight - mean(dat_no_ub$service_road_weight)) /
      sd(dat_no_ub$service_road_weight)
site_id <- as.integer(factor(dat_no_ub$site))
N <- length(site_id); J <- max(site_id)

set.seed(1234)
n_sim <- 20
sims <- lapply(seq_len(n_sim), function(i) {
  alpha  <- abs(rnorm(1, 0.5, 0.5))
  mu     <- rnorm(1, 0, 1) + abs(rnorm(1, 0, 1)) * rnorm(J, 0, 1)
  b      <- rnorm(3, 0, 1)
  gam    <- rnorm(1, 1, 1)
  d      <- rnorm(4, 0, 1)
  sigma  <- exp(mu[site_id] + b[1] * Xh + b[2] * Xd + b[3] * Xe)
  pi0    <- plogis(gam + d[1] * P1 + d[2] * Xd + d[3] * Xh + d[4] * Xe)
  y <- ifelse(runif(N) < pi0, 0, rweibull(N, alpha, sigma))
  tibble(sim_id = i, fw_score_weighted = y)
})
sim_df <- bind_rows(sims)
obs_zero_rate <- mean(dat_no_ub$fw_score_weighted == 0)

pp1 <- ggplot(filter(sim_df, sim_id <= 4), aes(x = fw_score_weighted)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ sim_id, scales = "free_x") +
  theme_minimal(base_size = 11) + labs(x = "Food waste score", y = "Count")
pp2 <- sim_df |>
  group_by(sim_id) |> summarise(p_zero = mean(fw_score_weighted == 0)) |>
  ggplot(aes(x = p_zero)) +
  geom_histogram(bins = 40, fill = "firebrick", alpha = 0.8) +
  geom_vline(xintercept = obs_zero_rate, linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 11) + labs(x = "Proportion of zeros (simulation)", y = "Count")
pp3 <- ggplot(filter(sim_df, sim_id <= 4, fw_score_weighted > 0), aes(x = fw_score_weighted)) +
  geom_histogram(bins = 100, fill = "darkgreen", alpha = 0.7) +
  facet_wrap(~ sim_id, scales = "free_x") +
  theme_minimal(base_size = 11) + labs(x = "Positive food waste score", y = "Count")
ggsave(here("figs", "FigureS5.pdf"),
       (pp1 / pp2 / pp3) + plot_annotation(tag_levels = "a", tag_suffix = ")"),
       width = 6, height = 8, device = pdf)

## =========================================================
## Figure S6: posterior predictive checks for F123_all
## =========================================================
fit <- load_fit("fit_hier_ZIW_16var_F123_all.Rdata")
y_rep <- rstan::extract(fit, pars = "fw_score_pred")$fw_score_pred
s_all <- summary(fit, pars = par_info$par, probs = c(0.025, 0.5, 0.975))$summary
rm(fit); gc()

y_obs <- dat_no_ub$fw_score_weighted
set.seed(123)
rows <- sample(nrow(y_rep), 50)
base_theme <- theme_minimal() +
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "none")
shared_legend <- get_legend(ppc_dens_overlay(y_obs, y_rep[rows, ]) +
                              theme_minimal() + theme(legend.position = "right"))
q1 <- ppc_dens_overlay(y_obs, y_rep[rows, ]) +
  labs(title = "Full Range", x = "", y = "Density") + base_theme
q2 <- ppc_dens_overlay(y_obs, y_rep[rows, ]) + xlim(0.000001, 10) +
  labs(title = "Zoom: >0–10", x = "", y = "") + base_theme
q3 <- ppc_dens_overlay(y_obs, y_rep[rows, ]) + coord_cartesian(xlim = c(-0.01, 1)) +
  labs(title = "Zoom: 0–1", x = "", y = "") + base_theme
panel_row <- plot_grid(q1, q2, q3, labels = c("a", "b", "c"), ncol = 3,
                       align = "hv", axis = "tblr", label_size = 16)
main_plot <- plot_grid(panel_row, ggdraw() + draw_label("Food Waste Score", size = 14),
                       ncol = 1, rel_heights = c(1, 0.08))
ggsave(here("figs", "FigureS6.png"),
       plot_grid(main_plot, shared_legend, ncol = 2, rel_widths = c(1, 0.12)),
       width = 12, height = 4, dpi = 300)

cat("\nObserved vs. predicted share of zeros:",
    round(mean(y_obs == 0), 3), "vs.", round(mean(y_rep[rows, ] == 0), 3), "\n")
rm(y_rep); gc()

## =========================================================
## Tables
## =========================================================
## F123_all coefficients
write.csv(orient_summary(s_all), here("tables", "Table_F123_all_coefficients.csv"), row.names = FALSE)

## LOO comparison
load(here("data", "loo_16var_models_all.Rdata"))   # loos_all
model_desc <- c(
  F123_all = "All three scales in presence and quantity (main model)",
  F123_D   = "All three scales in quantity; Deprivation also in presence",
  F123     = "All three scales in quantity only",
  F12      = "Housing + Deprivation in quantity only",
  F13      = "Housing + Disenfranchisement in quantity only",
  proxy    = "Proxy: distance to food retail (Tweedie)"
)
cmp <- loo_compare(loos_all)
loo_tab <- data.frame(
  model = rownames(cmp),
  description = unname(model_desc[rownames(cmp)]),
  elpd_loo = round(cmp[, "elpd_loo"], 1), se_elpd_loo = round(cmp[, "se_elpd_loo"], 1),
  p_loo = round(cmp[, "p_loo"], 1), looic = round(cmp[, "looic"], 1),
  elpd_diff = round(cmp[, "elpd_diff"], 1), se_diff = round(cmp[, "se_diff"], 1),
  row.names = NULL
)
write.csv(loo_tab, here("tables", "Table_LOO_all_models.csv"), row.names = FALSE)
print(loo_tab)

## Bins vs food litter
load(here("data", "bins_litter_fits_summary.Rdata"))   # summaries (bins, litter)
bl <- rbind(cbind(outcome = "Bins and receptacles", orient_summary(summaries$bins)),
            cbind(outcome = "Food litter",          orient_summary(summaries$litter)))
write.csv(bl, here("tables", "Table_bins_litter.csv"), row.names = FALSE)

## Receptacle-form weight sensitivity
load(here("data", "weight_sensitivity_F123_all.Rdata"))   # res (original, equal, older)
ws_labels <- c(original = "Original (dumpster 30x home bin)",
               older = "Older weights (dumpster 15x home bin)",
               equal = "Equal form weights")
ws <- do.call(rbind, lapply(names(ws_labels), function(n)
  cbind(weights = ws_labels[[n]], orient_summary(res[[n]]))))
write.csv(ws, here("tables", "Table_weight_sensitivity.csv"), row.names = FALSE)

cat("\nDone: FigureS2.pdf, FigureS5.pdf, FigureS6.png and tables/ written.\n")
