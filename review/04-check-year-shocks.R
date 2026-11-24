
# Code-review check: are residuals shared across time series within a year
# (i.e., is a year-level random effect missing)? And how much of the midge
# covariate varies among years vs. within years? Run from project root:
#   Rscript review/04-check-year-shocks.R > review/04-check-year-shocks-output.txt 2>&1

source("scripts/00-preamble.R")
options(width = 120, pillar.width = 120)

d <- read_rds(data_rds)
cat("Share of midges_z variance among years (R^2 of midges_z ~ factor(year)):",
    round(summary(lm(midges_z ~ factor(year), data = d))$r.squared, 3), "\n")
cat("... among years after removing plot means (R^2 of year in midges_z ~ plot + year minus plot-only):",
    round(summary(lm(midges_z ~ plot + factor(year), data = d))$r.squared -
              summary(lm(midges_z ~ plot, data = d))$r.squared, 3), "\n")

m <- read_rds(model_rds)
od <- m$orig_data                        # data in the order passed to Stan
obs_per <- m$stan_data$obs_per
ts_id <- rep(seq_along(obs_per), obs_per)
first <- !duplicated(ts_id)

# Posterior-mean process-error deviates (standard normal scale), one per non-first obs:
ze <- colMeans(rstan::extract(m$stan, "ze")[[1]])
yr <- od$time[!first]
tx <- od$taxon[!first]
cat("\nProcess deviates ze (posterior means), n =", length(ze), "\n")
f_all <- lm(ze ~ factor(yr))
cat("R^2 of year (all taxa pooled):", round(summary(f_all)$r.squared, 3),
    "  F-test p =", signif(anova(f_all)[1, "Pr(>F)"], 3), "\n")
f_tx <- lm(ze ~ factor(yr) * tx)
cat("R^2 of year x taxon:", round(summary(f_tx)$r.squared, 3),
    "  F-test p (year x taxon vs none) =",
    signif(anova(lm(ze ~ 1), f_tx)[2, "Pr(>F)"], 3), "\n")
# Expected R^2 under independence ~ (df_model)/(n-1):
cat("Expected R^2 under no shared year effect: year",
    round(10 / (length(ze) - 1), 3), "; year x taxon",
    round((length(coef(f_tx)) - 1) / (length(ze) - 1), 3), "\n")

# Observation-level residuals on log scale: log(y_pred) - (ly_pred + offset)
lyp <- colMeans(rstan::extract(m$stan, "ly_pred")[[1]])
ypr <- colMeans(log(rstan::extract(m$stan, "y_pred")[[1]]))
eo <- ypr - (lyp + m$stan_data$offset)
f_o <- lm(eo ~ factor(od$time) * od$taxon)
cat("\nObservation residuals: R^2 of year x taxon =", round(summary(f_o)$r.squared, 3),
    " p =", signif(anova(lm(eo ~ 1), f_o)[2, "Pr(>F)"], 3), "\n")
cat("Expected under independence:", round((length(coef(f_o)) - 1) / (length(eo) - 1), 3), "\n")

# Year x taxon mean of combined residual (ze * sig_proc + obs residual), to show magnitudes
sp <- median(rstan::extract(m$stan, "sig_proc")[[1]])
comb <- tibble(year = od$time + 2008, taxon = od$taxon, r = eo)
comb$r[!first] <- comb$r[!first] + sp * ze
cat("\nMean combined log-scale residual by year and taxon (16 plots each):\n")
comb |> group_by(taxon, year) |> summarize(r = round(mean(r), 2), .groups = "drop") |>
    pivot_wider(names_from = taxon, values_from = r) |> print(n = 20)
