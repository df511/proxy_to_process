## predict_F123_all.R
## City-wide predictions and maps for F123_all (all three factors in BOTH model parts).
##
## Amount (Weibull): beta1 = Housing, beta2 = Deprivation, beta3 = Enfranchisement
## Presence (zero-inflation, P(no record)): delta1 = service roads, delta2 = Deprivation,
##                                          delta3 = Housing, delta4 = Enfranchisement
##
## Adds pred_F123_all (expected FW score), presence_F123_all (P(any record)) and
## amount_F123_all (expected amount where present) to data/van_dat_predictions_16var.Rdata,
## then draws maps:
##   figs/city_pred_F123_all_stretched.pdf   expected FW score, stretched ramp (shared with
##                                           the other models' *_stretched maps)
##   figs/city_presence_F123_all.pdf         probability anything is recorded
##   figs/city_amount_F123_all.pdf           expected amount where present
##   figs/city_pred_F123_all_750m.pdf        expected FW score averaged to 750 m cells
##
## Posterior median of each quantity per cell; transect intercepts set to mu0.
## Run from the proxy_to_process project root, in a fresh R session. 5-15 minutes.

library(here)
library(sf)
library(dplyr)
library(rstan)
library(ggplot2)
library(scales)
library(viridisLite)

## ---------------------------------------------------------------
## 1. Posterior draws (only the parameters needed)
## ---------------------------------------------------------------
env <- new.env()
load(here("data", "fit_hier_ZIW_16var_F123_all.Rdata"), envir = env)
post <- rstan::extract(get(ls(env)[1], envir = env),
                       pars = c("mu0", "alpha", "beta1", "beta2", "beta3",
                                "gamma", "delta1", "delta2", "delta3", "delta4"))
rm(env); gc()

## ---------------------------------------------------------------
## 2. City cells
## ---------------------------------------------------------------
load(here("data", "van_dat_predictions_16var.Rdata"))  # van_dat_fa, k_housing, k_enfr, k_depriv
city <- st_drop_geometry(van_dat_fa)
Hous <- city[[paste0("Factor", k_housing)]]
Dep  <- city[[paste0("Factor", k_depriv)]]
Enf  <- city[[paste0("Factor", k_enfr)]]
Road <- city$service_road_weight_z

B <- cbind(post$beta1, post$beta2, post$beta3)          # draws x 3
D <- cbind(post$delta1, post$delta2, post$delta3, post$delta4)
X <- cbind(Hous, Dep, Enf)
P <- cbind(Road, Dep, Hous, Enf)
shape_term <- gamma(1 + 1 / post$alpha)

## Returns per-cell posterior medians of: expected value, P(present), amount where present
predict_chunk <- function(idx) {
  amount  <- sweep(exp(sweep(X[idx, , drop = FALSE] %*% t(B), 2, post$mu0, "+")),
                   2, shape_term, "*")                                   # cells x draws
  present <- 1 - plogis(sweep(P[idx, , drop = FALSE] %*% t(D), 2, post$gamma, "+"))
  cbind(expected = apply(present * amount, 1, median),
        present  = apply(present, 1, median),
        amount   = apply(amount, 1, median))
}

n <- nrow(city)
out <- matrix(NA_real_, n, 3)
for (start in seq(1, n, by = 2000)) {
  idx <- start:min(start + 1999, n)
  out[idx, ] <- predict_chunk(idx)
}
van_dat_fa$pred_F123_all     <- out[, 1]
van_dat_fa$presence_F123_all <- out[, 2]
van_dat_fa$amount_F123_all   <- out[, 3]
rm(out); gc()

save(van_dat_fa, k_housing, k_enfr, k_depriv,
     file = here("data", "van_dat_predictions_16var.Rdata"))

cat("Quantiles:\n")
print(signif(sapply(c("pred_F123_all", "presence_F123_all", "amount_F123_all"), function(v)
  quantile(van_dat_fa[[v]], c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1))), 3))

## ---------------------------------------------------------------
## 3. Maps
## ---------------------------------------------------------------
map_theme <- theme_void() + theme(plot.title = element_text(face = "bold"))
bar <- guides(fill = guide_colourbar(barheight = unit(12, "lines")))
save_map <- function(p, file) ggsave(here("figs", file), p, width = 12, height = 8)

## (a) Expected FW score, stretched ramp. Same rule as map_stretched_ramp.R
##     (min-0.25: 20% of colours, 0.25-1: 70%, 1-max: 10%), with the range shared across
##     all models' predictions so this map compares directly with the other *_stretched maps.
pred_cols <- intersect(c("pred_F12", "pred_F123", "pred_F13", "pred_proxy", "pred_F123_all"),
                       names(van_dat_fa))
rng <- range(unlist(lapply(pred_cols, function(v) van_dat_fa[[v]])), na.rm = TRUE)
ramp_at <- c(0, 0.20, 0.90, 1)
data_at <- approx(ramp_at, c(rng[1], 0.25, 1, rng[2]), xout = seq(0, 1, length.out = 101))$y
brks <- c(c(0.05, 0.1, 0.25, 0.5, 0.75, 1)[c(0.05, 0.1, 0.25, 0.5, 0.75, 1) > rng[1]], rng[2])
cat("\nShared range for the stretched ramp:", signif(rng, 3),
    "\n(If this differs from the range printed by map_stretched_ramp.R, rerun that script so",
    "all *_stretched maps use the same scale.)\n")

p <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = pred_F123_all), color = NA) +
  scale_fill_gradientn(colours = inferno(101), values = rescale(data_at, from = rng),
                       limits = rng, breaks = brks,
                       labels = c(format(brks[-length(brks)], drop0trailing = TRUE),
                                  paste0(signif(rng[2], 3), " (max)")),
                       name = "Expected\nFW score") + bar +
  labs(title = "Expected FW score: all factors in presence and amount (F123_all)",
       subtitle = "Full range shown; colour ramp stretched so 0.25–1 takes 70% of the colours") +
  map_theme
save_map(p, "city_pred_F123_all_stretched.pdf")

## (b) Probability that anything is recorded (receptacle or litter)
p <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = presence_F123_all), color = NA) +
  scale_fill_viridis_c(option = "inferno", name = "P(any\nrecord)") + bar +
  labs(title = "Probability of any receptacle or food litter (F123_all, presence part)",
       subtitle = "Service roads plus Housing, Deprivation and Enfranchisement") +
  map_theme
save_map(p, "city_presence_F123_all.pdf")

## (c) Expected amount where present
lims <- quantile(van_dat_fa$amount_F123_all, c(0.01, 0.99))
p <- ggplot(van_dat_fa) +
  geom_sf(aes(fill = amount_F123_all), color = NA) +
  scale_fill_viridis_c(option = "inferno", limits = lims, oob = squish,
                       name = "Expected amount\nwhere present") + bar +
  labs(title = "Amount of food waste where present (F123_all, amount part)",
       subtitle = "Housing, Deprivation and Enfranchisement; limits at 1st-99th percentile") +
  map_theme
save_map(p, "city_amount_F123_all.pdf")

## (d) Expected FW score averaged to 750 m cells (UTM 10N)
city_utm <- st_transform(van_dat_fa["pred_F123_all"], 26910)
grid750 <- st_sf(geometry = st_make_grid(city_utm, cellsize = 750))
grid750$cell <- seq_len(nrow(grid750))
pts <- st_join(suppressWarnings(st_centroid(city_utm)), grid750, join = st_within)
cell_means <- st_drop_geometry(pts) |>
  group_by(cell) |>
  summarise(pred = mean(pred_F123_all), .groups = "drop")
grid750 <- inner_join(grid750, cell_means, by = "cell")
lims750 <- quantile(grid750$pred, c(0.01, 0.99))
p <- ggplot(grid750) +
  geom_sf(aes(fill = pred), color = NA) +
  scale_fill_viridis_c(option = "inferno", limits = lims750, oob = squish,
                       name = "Mean expected\nFW score") + bar +
  labs(title = "Expected FW score, 750 m average (F123_all)",
       subtitle = "Averaged within 750 m cells; limits at 1st-99th percentile") +
  map_theme
save_map(p, "city_pred_F123_all_750m.pdf")

cat("\nDone. Maps written to figs/.\n")
