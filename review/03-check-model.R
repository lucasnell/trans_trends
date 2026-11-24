
# Code-review checks on the fitted model (rds_files/model-fit.rds) and on the
# effect-size calculations in scripts/05-model-summaries.R. Run from project root:
#   Rscript review/03-check-model.R > review/03-check-model-output.txt 2>&1

source("scripts/00-preamble.R")
options(width = 120, pillar.width = 120)

d <- read_rds(data_rds)
m <- read_rds(model_rds)
fit <- m$stan

# ---- 1. Sampler health ----
cat("== Sampler diagnostics ==\n")
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
cat("Divergences:", sum(sapply(sp, \(x) sum(x[, "divergent__"]))), "\n")
cat("Max treedepth hit:", sum(sapply(sp, \(x) sum(x[, "treedepth__"] >= 10))), "\n")
cat("E-BFMI by chain:", round(sapply(sp, \(x) {
    e <- x[, "energy__"]; sum(diff(e)^2) / sum((e - mean(e))^2) }), 3), "\n")
pars <- c("alpha", "sig_beta", "phi", "sig_proc", "sig_obs", "z")
s <- rstan::summary(fit, pars = pars)$summary
cat("Max Rhat (main pars):", round(max(s[, "Rhat"], na.rm = TRUE), 4),
    "  min n_eff:", round(min(s[, "n_eff"], na.rm = TRUE)), "\n")
print(round(s[!grepl("^z\\[", rownames(s)), c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], 3))
sall <- rstan::summary(fit)$summary
cat("Max Rhat (all incl. latent):", round(max(sall[, "Rhat"], na.rm = TRUE), 4),
    "  min n_eff (all):", round(min(sall[, "n_eff"], na.rm = TRUE)), "\n")

# ---- 2. Prior vs posterior for parameters plotted in figures ----
cat("\n== Prior-to-posterior contraction (1 - post var / prior var) ==\n")
prior_var <- c(gamma = 1.5 / 3^2, beta22 = 1 / 20)  # Gamma(1.5, rate 3); Beta(2,2)
pv <- function(p, pvar) {
    x <- rstan::extract(fit, p)[[1]]
    if (is.null(dim(x)) || length(dim(x)) == 1) x <- matrix(x)
    round(1 - apply(x, 2, var) / pvar, 3)
}
cat("sig_beta (int_tax, int_plot, int_trans, midges, time, dist):\n"); print(pv("sig_beta", prior_var[["gamma"]]))
cat("phi (", paste(m$ar_names, collapse = ", "), "):\n"); print(pv("phi", prior_var[["beta22"]]))
cat("sig_proc:", pv("sig_proc", prior_var[["gamma"]]), " sig_obs:", pv("sig_obs", prior_var[["gamma"]]), "\n")

# ---- 3. Is extract(pars = c(...)) returned in the requested order? ----
ex <- rstan::extract(fit, pars = c("sig_proc", "sig_obs", "sig_beta"))
cat("\nextract() order for error SDs:", paste(names(ex), collapse = ", "), "\n")

# ---- 4. Verify re_draws indexing in scripts/05 against Stan's own beta ----
# Rebuild taxon-level coefficient draws as scripts/05 does:
fef <- rstan::extract(fit, "alpha")[[1]]; colnames(fef) <- colnames(m$stan_data$x)
sig_names <- names(fit)[grepl("^sig_beta\\[", names(fit))]
S <- do.call(cbind, lapply(seq_along(sig_names), \(i) {
    Si <- rstan::extract(fit, sig_names[i])[[1]]
    matrix(as.numeric(Si), length(Si), m$stan_data$lev_per_g[i]) }))
Z <- rstan::extract(fit, "z")[[1]]
E <- S * Z
for (i in 1:nrow(m$rnd_lvl_names)) E[, i] <- E[, i] + fef[, m$rnd_lvl_names$Name[i]]
# Stan's beta[ts, coef] for a time series in each taxon (no other groups on slopes):
B <- rstan::extract(fit, "beta")[[1]]   # draws x n_ts x n_coef
ts_taxon <- m$orig_data$taxon[cumsum(m$stan_data$obs_per)]
maxdiff <- 0
for (cf in c("midges_z", "time_z", "dist_z")) {
    for (tx in levels(d$taxon)) {
        j <- which(m$rnd_lvl_names$Name == cf & m$rnd_lvl_names$Groups == "taxon" &
                       m$rnd_lvl_names$Level == tx)
        ts <- which(ts_taxon == tx)[1]
        maxdiff <- max(maxdiff, max(abs(E[, j] - B[, ts, which(colnames(fef) == cf)])))
    }
}
cat("Max |scripts/05 taxon slope - Stan beta|:", signif(maxdiff, 3), "\n")

# ---- 5. What the time and midge coefficients mean under change = TRUE ----
# Model (change = TRUE): mu[t] = x[t] . beta + phi * mu[t-1] + process error,
# mu[1] = x[1] . beta. So the long-run effect of a sustained unit change in a
# predictor on log catch rate is beta / (1 - phi), not beta.
phi <- rstan::extract(fit, "phi")[[1]]
colnames(phi) <- gsub("^taxon", "", m$ar_names)
time_sd <- sd(d$time); midges_sd <- sd(log(d$midges / d$midge_days))
max_time_z <- max(d$time_z)
tt <- sort(unique(d$time_z))

slope <- \(cf, tx) E[, which(m$rnd_lvl_names$Name == cf & m$rnd_lvl_names$Groups == "taxon" &
                                m$rnd_lvl_names$Level == tx)]
int_tx <- \(tx) E[, which(m$rnd_lvl_names$Name == "(Intercept)" &
                             m$rnd_lvl_names$Groups == "taxon" & m$rnd_lvl_names$Level == tx)]
q <- \(x) sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.025), quantile(x, 0.975))

res <- map_dfr(taxa_lvls, \(tx) {
    b_t <- slope("time_z", tx); b_m <- slope("midges_z", tx); p <- phi[, tx]
    b0 <- int_tx(tx)
    # Deterministic implied mean path (no process error), midges_z = dist_z = 0,
    # taxon-level intercept, run through the model's own recursion:
    path <- function(bt) {
        mu <- matrix(NA, length(p), length(tt))
        mu[, 1] <- b0 + bt * tt[1]
        for (k in 2:length(tt)) mu[, k] <- b0 + bt * tt[k] + p^(1) * mu[, k - 1]
        mu
    }
    mu <- path(b_t); mu0 <- path(0 * b_t)
    tibble(taxon = tx,
           phi_med = round(median(p), 3),
           script05_pct_total = q(100 * (exp(b_t * max_time_z) - 1)),
           longrun_pct_total = q(100 * (exp(b_t / (1 - p) * max_time_z) - 1)),
           implied_path_pct = q(100 * (exp(mu[, length(tt)] - mu[, 1]) - 1)),
           path_pct_if_no_trend = q(100 * (exp(mu0[, length(tt)] - mu0[, 1]) - 1)),
           script05_midge10x = sprintf("%.3f", median(10^(b_m / midges_sd))),
           longrun_midge10x = sprintf("%.3f", median(10^(b_m / (1 - p) / midges_sd))),
           P_time_gt0 = round(mean(b_t > 0), 3),
           P_midge_gt0 = round(mean(b_m > 0), 3))
})
cat("\n== Effect sizes: scripts/05 (beta) vs long-run (beta/(1-phi)) vs implied path ==\n")
cat("Cells are posterior median [95% interval]; first-to-last-year % change.\n")
print(res |> select(taxon, phi_med, script05_pct_total, longrun_pct_total,
                    implied_path_pct, path_pct_if_no_trend), width = 200)
print(res |> select(taxon, script05_midge10x, longrun_midge10x, P_time_gt0, P_midge_gt0))
