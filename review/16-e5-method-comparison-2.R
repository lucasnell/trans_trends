
# E5, pairs 1-2: extends review/14-e5-method-comparison.R with
#   (a) gam_ls : same mean as `gam`, but the residual SD varies by plot and
#                season: log(squared GAM residual) ~ plot + s(day-of-year),
#                rescaled so mean predicted variance = GAM residual variance
#   (b) joint scoring: the three responses (big midges, small midges,
#                simuliids) scored with a multivariate normal using the
#                correlation of their cross-validation residuals, instead of
#                as independent (applied to rf_gauss, gam, gam_ls)
#   (c) stacked classifier: logistic regression with no intercept on the
#                per-response log-likelihood-ratio components of rf_gauss,
#                gam, and gam_ls, trained on "true vs. swapped" test pairs.
#                No intercept makes it antisymmetric (swapping the pair flips
#                the sign of the prediction). Evaluated by leave-one-year-out
#                over pairs (all pairs are within one year), so a year's pairs
#                never inform their own classifier.
# Everything else (data, features, folds, test-pair design) follows script 14,
# rerun here so all methods are scored on the same folds and pairs.
#
# Note: the residual correlations for joint scoring are estimated from all
# cross-validation residuals (3 numbers from ~1,000 samples), so they are not
# re-estimated without each test pair; the leakage is negligible.
#
# Run from project root:
#   Rscript review/16-e5-method-comparison-2.R > review/16-e5-method-comparison-2-output.txt 2>&1

.libPaths(c(.libPaths(), "/Users/lucasnell/Library/R/arm64/4.6/library"))  # ranger
suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(ranger)
    library(mgcv)
})
options(width = 160)
set.seed(20261125)
NTREE <- 500
K <- as.integer(Sys.getenv("E5_K", "10"))            # override only for smoke tests
N_PAIRS <- as.integer(Sys.getenv("E5_NPAIRS", "20000"))

# ---- Data ----
key <- c("lakeid", "dist", "trans", "setdate", "coldate")
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
cols <- c(key, "fract_bgch", "bgch", "fract_smch", "smch", "fract_simu", "simu")
i1 <- rd("data/raw/infalls-08-17.csv") |> unique() |> select(all_of(cols))
i2 <- rd("data/raw/infalls-13-19.csv") |> unique() |> select(all_of(cols)) |>
    filter(!(trans == "btl" & dist == 5 & setdate == as.Date("2013-06-05") &
                 coldate == as.Date("2013-08-14")))
inf <- bind_rows(i1, anti_join(i2, i1, by = key)) |>
    filter(lakeid == "myv", trans %in% c("btl", "hag", "kal", "non", "vin"),
           !is.na(setdate), !is.na(coldate), !is.na(dist)) |>
    mutate(row_id = row_number())

pair1 <- inf |> filter(trans == "hag", dist == 5, setdate == as.Date("2012-06-05"),
                       coldate == as.Date("2012-06-19")) |> arrange(smch)
pair2 <- inf |> filter(trans == "non", dist == 50, setdate == as.Date("2012-07-17"),
                       coldate == as.Date("2012-07-31")) |> arrange(fract_smch < 1)
stopifnot(nrow(pair1) == 2, nrow(pair2) == 2,
          pair1$smch[1] == 291, pair2$smch[1] == 360)   # row 1 = line 863 / 942
na_copies <- inf |>
    filter((trans == "non" & dist == 150 & setdate == as.Date("2008-06-04")) |
               (trans == "vin" & dist == 5 & setdate == as.Date("2008-05-30"))) |>
    filter(is.na(smch))
stopifnot(nrow(na_copies) == 2)
inf <- inf |> filter(!row_id %in% c(pair1$row_id, pair2$row_id, na_copies$row_id))

resp_def <- list(big = c("bgch", "fract_bgch"), small = c("smch", "fract_smch"),
                 simu = c("simu", "fract_simu"))
total <- \(row, r) row[[resp_def[[r]][1]]] / row[[resp_def[[r]][2]]]
inf <- inf |> mutate(days = as.numeric(coldate - setdate))
for (r in names(resp_def)) inf[[paste0("y_", r)]] <- log1p(total(inf, r) / inf$days)
cat("Samples:", nrow(inf), "| with small-midge value:", sum(!is.na(inf$y_small)), "\n")

# ---- Features (as in script 13) ----
features <- function(cells, known, r) {
    yk <- paste0("y_", r)
    kn <- known |> filter(!is.na(.data[[yk]])) |>
        select(trans, dist, setdate, coldate, yv = all_of(yk))
    x <- cells |>
        mutate(cell = row_number(), days = as.numeric(coldate - setdate),
               mid_doy = yday(setdate) + days / 2, year = year(coldate),
               ldist = log(dist),
               site = factor(trans, levels = c("btl", "hag", "kal", "non", "vin")))
    oth <- x |> select(cell, trans, dist, setdate, coldate) |>
        inner_join(kn, by = c("trans", "setdate", "coldate"), suffix = c("", ".o"),
                   relationship = "many-to-many") |>
        filter(dist.o != dist) |>
        group_by(cell) |> summarize(other_mean = mean(yv), other_n = n(), .groups = "drop")
    nb <- x |> select(cell, trans, dist, setdate, coldate) |>
        inner_join(kn, by = c("trans", "dist"), suffix = c("", ".k"),
                   relationship = "many-to-many") |>
        filter(year(coldate.k) == year(coldate),
               !(setdate.k == setdate & coldate.k == coldate)) |>
        group_by(cell) |>
        summarize(prev = { i <- which(coldate.k <= setdate + 3)
                           if (length(i)) yv[i][which.max(coldate.k[i])] else NA_real_ },
                  nxt = { i <- which(setdate.k >= coldate - 3)
                          if (length(i)) yv[i][which.min(setdate.k[i])] else NA_real_ },
                  .groups = "drop")
    x |> left_join(oth, by = "cell") |> left_join(nb, by = "cell") |>
        mutate(other_n = replace_na(other_n, 0L))
}
fvars <- c("site", "ldist", "year", "mid_doy", "days", "other_mean", "other_n", "prev", "nxt")

# ---- Model fitting and prediction ----
# Returns, for each row of `newcells`, the predictive mean and SD under each
# method (rf_gauss, gam, gam_ls).
fit_all <- function(train, newcells, r) {
    yk <- paste0("y_", r)
    tr <- features(train, train, r) |> mutate(y = train[[yk]]) |> filter(!is.na(y))
    nd <- features(newcells, train, r)
    # random forest
    m <- ranger(x = as.data.frame(tr[fvars]), y = tr$y, num.trees = NTREE,
                na.action = "na.learn", seed = 1)
    rf_mean <- predict(m, data = as.data.frame(nd[fvars]))$predictions
    # GAMs
    gtr <- tr |> mutate(plot = factor(paste(trans, dist)),
                        sitedate = factor(paste(trans, setdate, coldate)),
                        plotyear = factor(paste(trans, dist, year)))
    f_mu <- y ~ plot + s(mid_doy, by = plot, k = 5) + s(sitedate, bs = "re") +
        s(plotyear, bs = "re")
    g <- bam(f_mu, data = gtr, method = "fREML", discrete = TRUE)
    # Residual-SD model (two-stage; a one-stage gaulss fit fails here because
    # the mean model has about as many coefficients as rows): regress
    # log(squared residual) on plot + smooth of day-of-year, then rescale so
    # the mean predicted variance equals the GAM's residual variance.
    gtr$lr2 <- log(pmax(residuals(g)^2, 1e-6))
    gv <- gam(lr2 ~ plot + s(mid_doy, k = 5), data = gtr)
    vscale <- g$sig2 / mean(exp(fitted(gv)))
    gnd <- nd |> mutate(plot = paste(trans, dist),
                        sitedate = paste(trans, setdate, coldate),
                        plotyear = paste(trans, dist, year))
    ok_plot <- gnd$plot %in% levels(gtr$plot)
    ok_sd <- gnd$sitedate %in% levels(gtr$sitedate)
    ok_py <- gnd$plotyear %in% levels(gtr$plotyear)
    gnd <- gnd |> mutate(plot = factor(ifelse(ok_plot, plot, levels(gtr$plot)[1]), levels(gtr$plot)),
                         sitedate = factor(ifelse(ok_sd, sitedate, levels(gtr$sitedate)[1]), levels(gtr$sitedate)),
                         plotyear = factor(ifelse(ok_py, plotyear, levels(gtr$plotyear)[1]), levels(gtr$plotyear)))
    n <- nrow(gnd)
    gm <- gse <- lsd <- rep(NA_real_, n)
    combos <- paste(ok_sd, ok_py)
    for (cb in unique(combos)) {
        idx <- which(combos == cb)
        ex <- c(if (!ok_sd[idx[1]]) "s(sitedate)", if (!ok_py[idx[1]]) "s(plotyear)")
        ex <- if (length(ex)) ex else NULL
        p <- predict(g, newdata = gnd[idx, ], se.fit = TRUE, exclude = ex)
        gm[idx] <- p$fit; gse[idx] <- p$se.fit
        lsd[idx] <- sqrt(vscale * exp(predict(gv, newdata = gnd[idx, ])))
    }
    gm[!ok_plot] <- NA
    list(rf_gauss = list(mean = rf_mean, sd = rep(sqrt(m$prediction.error), n)),
         gam = list(mean = gm, sd = sqrt(gse^2 + g$sig2)),
         gam_ls = list(mean = gm, sd = sqrt(gse^2 + lsd^2)))
}
methods <- c("rf_gauss", "gam", "gam_ls")

# ---- Folds and test pairs (drawn before any model fitting) ----
folds <- sample(rep(seq_len(K), length.out = nrow(inf)))
cand <- inf |> filter(!is.na(y_small)) |> mutate(year = year(coldate)) |>
    select(row_id, trans, year)
pairs <- cand |> inner_join(cand, by = "year", suffix = c(".a", ".b"),
                            relationship = "many-to-many") |>
    filter(trans.a < trans.b) |> slice_sample(n = N_PAIRS)
cat("Test pairs:", nrow(pairs), "\n")

# ---- Cross-validated predictions: long table row_id x response x method ----
cvp <- map_dfr(seq_len(K), \(k) {
    cat("fold", k, "\n")
    tr <- inf[folds != k, ]; te <- inf[folds == k, ]
    map_dfr(names(resp_def), \(r) {
        P <- fit_all(tr, te, r)
        map_dfr(methods, \(mm) tibble(row_id = te$row_id, response = r, method = mm,
                                      mean = P[[mm]]$mean, sd = P[[mm]]$sd))
    })
})

# Residual correlation among responses, per method (standardized CV residuals)
obs_long <- inf |> select(row_id, starts_with("y_")) |>
    pivot_longer(-row_id, names_to = "response", values_to = "y") |>
    mutate(response = str_remove(response, "^y_"))
R_by_method <- map(set_names(methods), \(mm) {
    w <- cvp |> filter(method == mm) |> left_join(obs_long, by = c("row_id", "response")) |>
        mutate(z = (y - mean) / sd) |> select(row_id, response, z) |>
        pivot_wider(names_from = response, values_from = z)
    cor(w[, names(resp_def)], use = "pairwise.complete.obs")
})
cat("\n===== Correlation of standardized CV residuals among responses =====\n")
walk(methods, \(mm) { cat(mm, ":\n"); print(round(R_by_method[[mm]], 3)) })

# ---- Scoring ----
cvp_list <- split(cvp, paste(cvp$row_id, cvp$method))
get_pred <- function(id, mm) {
    x <- cvp_list[[paste(id, mm)]]
    x[match(names(resp_def), x$response), ]
}
yvec <- function(row, days) map_dbl(names(resp_def), \(r) log1p(total(row, r) / days))
mvn_ld <- function(y, mu, s, R) {
    ok <- !is.na(y) & !is.na(mu)
    if (!any(ok)) return(0)
    y <- y[ok]; mu <- mu[ok]; s <- s[ok]; R <- R[ok, ok, drop = FALSE]
    S <- diag(s, length(s)) %*% R %*% diag(s, length(s))
    d <- y - mu
    as.numeric(-0.5 * (length(d) * log(2 * pi) + determinant(S)$modulus + t(d) %*% solve(S, d)))
}
# Per-response LR components (independent scoring) and joint LR, for one pair
pair_scores <- function(ra, rb, pa, pb, mm) {
    da <- ra$days; db <- rb$days
    ya_a <- yvec(ra, da); yb_b <- yvec(rb, db); yb_a <- yvec(rb, da); ya_b <- yvec(ra, db)
    comp <- map_dbl(seq_along(resp_def), \(k) {
        l <- \(y, p) if (is.na(y) || is.na(p$mean[k])) 0 else dnorm(y, p$mean[k], p$sd[k], log = TRUE)
        (l(ya_a[k], pa) + l(yb_b[k], pb)) - (l(yb_a[k], pa) + l(ya_b[k], pb))
    })
    R <- R_by_method[[mm]]
    joint <- (mvn_ld(ya_a, pa$mean, pa$sd, R) + mvn_ld(yb_b, pb$mean, pb$sd, R)) -
        (mvn_ld(yb_a, pa$mean, pa$sd, R) + mvn_ld(ya_b, pb$mean, pb$sd, R))
    c(set_names(comp, paste0(mm, "_", names(resp_def))), set_names(joint, paste0(mm, "_joint")))
}
rowof <- \(id) inf[match(id, inf$row_id), ]
cat("\nScoring test pairs...\n")
sc <- map_dfr(seq_len(nrow(pairs)), \(q) {
    ra <- rowof(pairs$row_id.a[q]); rb <- rowof(pairs$row_id.b[q])
    v <- unlist(map(methods, \(mm) pair_scores(ra, rb, get_pred(ra$row_id, mm),
                                               get_pred(rb$row_id, mm), mm)))
    as_tibble(as.list(c(q = q, year = pairs$year[q], v)))
})

# LR columns per scoring approach
for (mm in methods) sc[[paste0(mm, "_indep")]] <- rowSums(sc[paste0(mm, "_", names(resp_def))])

# ---- Stacked classifier, leave-one-year-out ----
comp_cols <- unlist(map(methods, \(mm) paste0(mm, "_", names(resp_def))))
stack_fit <- function(dd) {
    # random orientation so the label isn't always 1; features flip sign with it
    flip <- sample(c(-1, 1), nrow(dd), replace = TRUE)
    X <- as.matrix(dd[comp_cols]) * flip
    yv <- as.integer(flip == 1)
    glm(yv ~ X - 1, family = binomial())
}
sc$stacked <- NA_real_
for (yr in unique(sc$year)) {
    fit <- stack_fit(filter(sc, year != yr))
    X <- as.matrix(sc[sc$year == yr, comp_cols])
    sc$stacked[sc$year == yr] <- as.numeric(X %*% coef(fit))
}
stack_all <- stack_fit(sc)
cat("\n===== Stacked classifier weights (fit on all pairs) =====\n")
print(round(setNames(coef(stack_all), comp_cols), 3))

approaches <- c(paste0(methods, "_indep"), paste0(methods, "_joint"), "stacked")
lr_long <- sc |> select(q, all_of(approaches)) |>
    pivot_longer(-q, names_to = "approach", values_to = "lr")

cat("\n===== Overall, per approach =====\n")
lr_long |> group_by(approach) |>
    summarize(accuracy = mean(lr > 0) + 0.5 * mean(lr == 0),
              confident_share = mean(abs(lr) > qlogis(0.95)),
              wrong_with_p_gt_0.95 = mean(lr < -qlogis(0.95)),
              log_score = mean(-log1p(exp(-lr))), .groups = "drop") |>
    arrange(desc(log_score)) |>
    mutate(across(where(is.double), \(x) round(x, 4))) |> as.data.frame() |> print()

cat("\n===== Calibration: stated P(chosen) vs. observed accuracy =====\n")
lr_long |> mutate(p = plogis(abs(lr)),
                  bin = cut(p, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001), right = FALSE)) |>
    group_by(approach, bin) |>
    summarize(n = n(), stated = round(mean(p), 3), observed = round(mean(lr > 0), 3),
              .groups = "drop") |>
    as.data.frame() |> print()

# ---- Real pairs ----
slots1 <- tibble(trans = c("hag", "btl"), dist = c(5, 50),
                 setdate = as.Date(c("2012-06-05", "2012-06-19")),
                 coldate = as.Date(c("2012-06-19", "2012-07-03")))
slots2 <- tibble(trans = c("non", "kal"), dist = c(50, 5),
                 setdate = as.Date(c("2012-07-17", "2012-07-17")),
                 coldate = as.Date(c("2012-07-31", "2012-07-31")))
real <- function(rows, slots) {
    P <- map(set_names(names(resp_def)), \(r) fit_all(inf, slots, r))
    pp <- \(i, mm) tibble(response = names(resp_def),
                          mean = map_dbl(names(resp_def), \(r) P[[r]][[mm]]$mean[i]),
                          sd = map_dbl(names(resp_def), \(r) P[[r]][[mm]]$sd[i]))
    ra <- rows[1, ] |> mutate(days = as.numeric(slots$coldate[1] - slots$setdate[1]))
    rb <- rows[2, ] |> mutate(days = as.numeric(slots$coldate[2] - slots$setdate[2]))
    v <- unlist(map(methods, \(mm) pair_scores(ra, rb, pp(1, mm), pp(2, mm), mm)))
    out <- as_tibble(as.list(v))
    for (mm in methods) out[[paste0(mm, "_indep")]] <- rowSums(out[paste0(mm, "_", names(resp_def))])
    out$stacked <- as.numeric(as.matrix(out[comp_cols]) %*% coef(stack_all))
    out
}
acc_at <- function(a, lr) {
    brks <- c(0, 1, 2, 4, 8, 16, 32, Inf)
    b <- cut(abs(lr), brks, right = FALSE)
    x <- lr_long |> filter(approach == a, cut(abs(lr), brks, right = FALSE) == b)
    sprintf("%.3f (n = %d, |LR| in %s)", mean(x$lr > 0), nrow(x), as.character(b))
}
cat("\n===== Real pairs: LR > 0 favours assignment 1 =====\n")
cat("Pair 1 assignment 1: line 863 -> hag 5 m 06/05-06/19, line 864 -> btl 50 m 06/19-07/03\n")
cat("Pair 2 assignment 1: line 942 -> non 50 m 07/17-07/31, line 944 -> kal 5 m 07/17-07/31\n")
for (pp in list(list("Pair 1", pair1, slots1), list("Pair 2", pair2, slots2))) {
    rr <- real(pp[[2]], pp[[3]])
    for (a in approaches) cat(sprintf("%s %-15s LR = %6.2f | accuracy at this |LR| in test: %s\n",
                                      pp[[1]], a, rr[[a]], acc_at(a, rr[[a]])))
}
