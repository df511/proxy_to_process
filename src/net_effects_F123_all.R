## net_effects_F123_all.R
## Net effect of each scale on EXPECTED wildlife-accessible food waste in F123_all, combining
## the presence and quantity parts (which can pull in opposite directions).
##
## Expected score: E = P(present) x E(amount | present)
##   P(present)  = 1 - inv_logit(gamma + delta1*Roads + delta2*Dep + delta3*Hous + delta4*Enf)
##   E(amount)   = exp(mu0 + beta1*Hous + beta2*Dep + beta3*Enf) * Gamma(1 + 1/alpha)
## Baseline: every predictor at its transect mean (0), site intercept at mu0.
## For each scale, the comparison point is +1 transect SD in the reported orientation
## (towards the less privileged condition; for service roads, +1 SD more service road).
## Per posterior draw:
##   net ratio      = E(+1 SD) / E(baseline)
##   presence ratio = P(present, +1 SD) / P(present, baseline)
##   quantity ratio = E(amount, +1 SD) / E(amount, baseline)
## (net = presence x quantity). Reported as posterior median and 95% interval.
##
## Outputs: tables/Table_net_effects.csv, figs/FigureS_net_effects.pdf
## Run from the proxy_to_process project root. About a minute.

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(rstan)

load(here("data", "factor_analysis_16var.Rdata"))   # Lambda
load(here("data", "dat_final_750_16var.Rdata"))     # dat_final_16var

pick <- function(v) unname(which.max(abs(Lambda[v, ])))
k <- c(housing = pick("housing_highdensity"),
       deprivation = pick("pct_rent_subsidized"),
       disenfranchisement = pick("employed_pct"))
orient <- c(housing = -1, deprivation = 1, disenfranchisement = -1)
sd_fac <- sapply(k, function(kk) sd(dat_final_16var[[paste0("Factor", kk)]]))
cat("Transect SD of each factor score:", round(sd_fac, 3), "\n")

env <- new.env()
load(here("data", "fit_hier_ZIW_16var_F123_all.Rdata"), envir = env)
post <- rstan::extract(get(ls(env)[1], envir = env),
                       pars = c("mu0", "alpha", "gamma", "beta1", "beta2", "beta3",
                                "delta1", "delta2", "delta3", "delta4"))
rm(env); gc()

## Shifts in the FITTED predictors for "+1 SD towards the less privileged condition"
## (X1 = Housing, X2 = Deprivation, X3 = Enfranchisement; P1 = roads, P2 = Dep, P3 = Hous, P4 = Enf)
shift <- function(scale) {
  z <- c(hous = 0, dep = 0, enf = 0, road = 0)
  if (scale == "roads") z["road"] <- 1
  if (scale == "housing") z["hous"] <- orient[["housing"]] * sd_fac[["housing"]]
  if (scale == "deprivation") z["dep"] <- orient[["deprivation"]] * sd_fac[["deprivation"]]
  if (scale == "disenfranchisement") z["enf"] <- orient[["disenfranchisement"]] * sd_fac[["disenfranchisement"]]
  z
}

p_present <- function(z) 1 - plogis(post$gamma + post$delta1 * z["road"] + post$delta2 * z["dep"] +
                                      post$delta3 * z["hous"] + post$delta4 * z["enf"])
amount <- function(z) exp(post$mu0 + post$beta1 * z["hous"] + post$beta2 * z["dep"] +
                            post$beta3 * z["enf"]) * gamma(1 + 1 / post$alpha)

z0 <- shift("none")
P0 <- p_present(z0); A0 <- amount(z0)
cat("Baseline expected FW score (posterior median):", signif(median(P0 * A0), 3),
    "| P(present):", signif(median(P0), 3), "| amount if present:", signif(median(A0), 3), "\n")

summ <- function(x) c(median = median(x), lower_95 = unname(quantile(x, 0.025)),
                      upper_95 = unname(quantile(x, 0.975)), p_gt_1 = mean(x > 1))

labels <- c(roads = "Service roads",
            housing = "Housing-tenure—density",
            deprivation = "Deprivation",
            disenfranchisement = "Disenfranchisement")
tab <- do.call(rbind, lapply(names(labels), function(s) {
  z <- shift(s)
  P1 <- p_present(z); A1 <- amount(z)
  rbind(
    data.frame(scale = labels[[s]], measure = "Net (expected FW score)", t(summ((P1 * A1) / (P0 * A0)))),
    data.frame(scale = labels[[s]], measure = "Presence",                t(summ(P1 / P0))),
    data.frame(scale = labels[[s]], measure = "Quantity",                t(summ(A1 / A0)))
  )
}))
tab[, c("median", "lower_95", "upper_95", "p_gt_1")] <- round(tab[, c("median", "lower_95", "upper_95", "p_gt_1")], 3)
rownames(tab) <- NULL
dir.create(here("tables"), showWarnings = FALSE)
write.csv(tab, here("tables", "Table_net_effects.csv"), row.names = FALSE)
cat("\nRatio of expected value at +1 SD (towards the less privileged condition) vs. the transect mean:\n")
print(tab)

## Figure: ratios on a log scale; 1 = no change
tab$scale   <- factor(tab$scale, levels = rev(unname(labels)))
tab$measure <- factor(tab$measure, levels = c("Net (expected FW score)", "Presence", "Quantity"))
p <- ggplot(tab, aes(x = median, y = scale, colour = measure, shape = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_95),
                  position = position_dodge(width = 0.6), size = 0.4) +
  scale_x_log10() +
  scale_colour_manual(values = c("Net (expected FW score)" = "black",
                                 "Presence" = "orange", "Quantity" = "purple"), name = NULL) +
  scale_shape_manual(values = c("Net (expected FW score)" = 16, "Presence" = 17, "Quantity" = 15),
                     name = NULL) +
  labs(x = "Ratio at +1 SD vs. transect mean (log scale)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
ggsave(here("figs", "FigureS_net_effects.pdf"), p, width = 6.5, height = 4)

cat("\nDone: tables/Table_net_effects.csv and figs/FigureS_net_effects.pdf\n")
