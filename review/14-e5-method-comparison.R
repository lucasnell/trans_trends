
# E5, pairs 1-2: compare three prediction methods on a much larger test set.
#
# Methods (one model per response: big midges, small midges, simuliids; the
# response is log1p(total / days), total = count / subsample fraction):
#   rf_gauss : random forest mean, Gaussian predictive density with constant
#              SD = out-of-bag residual SD (the method in script 13).
#   rf_quant : quantile regression forest; predictive density from the
#              forest's per-tree draws (point mass at zero + kernel density
#              for positive values).
#   gam      : y ~ plot + log-distance + smooth of day-of-year by plot
#              + random intercept per site-date + random intercept per
#              plot-year; Gaussian predictive density with
#              SD = sqrt(se.fit^2 + residual variance).
#
# Data: both infall files combined as in scripts/00-clean-data-round1.R (E6
# row removed), Mývatn, the five analysed sites (btl, hag, kal, non, vin),
# 2008-2019. E5 pairs 1-2 rows removed; all-NA copies of pairs 3-4 removed.
#
# Test: 10-fold cross-validation over samples. For each fold, features and
# models are built from the other nine folds only. Test pairs are random
# pairs of held-out samples from different sites in the same year; each
# pair is scored true vs. swapped by summed log predictive density over the
# three responses, using the slot's days to convert a candidate's total to a
# rate. Accuracy is reported by evidence strength |log-likelihood ratio|.
#
# Run from project root:
#   Rscript review/14-e5-method-comparison.R > review/14-e5-method-comparison-output.txt 2>&1

.libPaths(c(.libPaths(), "/Users/lucasnell/Library/R/arm64/4.6/library"))  # ranger
suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(ranger)
    library(mgcv)
})
options(width = 160)
set.seed(20261124)
NTREE <- 500
K <- 10
N_PAIRS <- 20000

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
# Returns, for each row of `newcells`: rf mean, rf OOB sd, rf quantile draws
# (matrix), gam mean, gam sd.
fit_all <- function(train, newcells, r) {
    yk <- paste0("y_", r)
    tr <- features(train, train, r) |> mutate(y = train[[yk]]) |> filter(!is.na(y))
    nd <- features(newcells, train, r)
    # random forest (quantile forest gives both the mean and per-tree draws)
    m <- ranger(x = as.data.frame(tr[fvars]), y = tr$y, num.trees = NTREE,
                na.action = "na.learn", quantreg = TRUE, keep.inbag = TRUE, seed = 1)
    ndx <- as.data.frame(nd[fvars])
    rf_mean <- predict(m, data = ndx)$predictions
    draws <- predict(m, data = ndx, type = "quantiles", what = \(x) x)$predictions
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = nrow(nd))
    # GAM
    gtr <- tr |> mutate(plot = factor(paste(trans, dist)),
                        sitedate = factor(paste(trans, setdate, coldate)),
                        plotyear = factor(paste(trans, dist, year)))
    g <- bam(y ~ plot + s(mid_doy, by = plot, k = 5) + s(sitedate, bs = "re") +
                 s(plotyear, bs = "re"), data = gtr, method = "fREML", discrete = TRUE)
    gnd <- nd |> mutate(plot = paste(trans, dist),
                        sitedate = paste(trans, setdate, coldate),
                        plotyear = paste(trans, dist, year))
    # unseen levels: substitute a seen level and drop that term from the prediction
    ok_plot <- gnd$plot %in% levels(gtr$plot)
    ok_sd <- gnd$sitedate %in% levels(gtr$sitedate)
    ok_py <- gnd$plotyear %in% levels(gtr$plotyear)
    gnd <- gnd |> mutate(plot = factor(ifelse(ok_plot, plot, levels(gtr$plot)[1]), levels(gtr$plot)),
                         sitedate = factor(ifelse(ok_sd, sitedate, levels(gtr$sitedate)[1]), levels(gtr$sitedate)),
                         plotyear = factor(ifelse(ok_py, plotyear, levels(gtr$plotyear)[1]), levels(gtr$plotyear)))
    gm <- rep(NA_real_, nrow(gnd)); gse <- rep(NA_real_, nrow(gnd))
    combos <- paste(ok_sd, ok_py)
    for (cb in unique(combos)) {
        idx <- which(combos == cb)
        ex <- c(if (!ok_sd[idx[1]]) "s(sitedate)", if (!ok_py[idx[1]]) "s(plotyear)")
        p <- predict(g, newdata = gnd[idx, ], se.fit = TRUE,
                     exclude = if (length(ex)) ex else NULL)
        gm[idx] <- p$fit; gse[idx] <- p$se.fit
    }
    gm[!ok_plot] <- NA  # no plot effect available: no GAM prediction
    list(rf_mean = rf_mean, rf_sd = sqrt(m$prediction.error), draws = draws,
         gam_mean = gm, gam_sd = sqrt(gse^2 + g$sig2))
}

# log predictive density of value y under each method for prediction i
lpd_quant <- function(y, d) {
    d <- d[!is.na(d)]; n <- length(d); n0 <- sum(d == 0)
    p0 <- (n0 + 0.5) / (n + 1)
    if (y == 0) return(log(p0))
    pos <- d[d > 0]
    if (length(pos) < 2) return(log(1 - p0) + dnorm(y, mean(c(pos, y)), 0.5, log = TRUE))
    bw <- max(bw.nrd0(pos), 0.05)
    log(1 - p0) + log(mean(dnorm(y, pos, bw)) + 1e-300)
}
lpd <- function(method, y, P, i) {
    if (is.na(y)) return(0)
    switch(method,
           rf_gauss = dnorm(y, P$rf_mean[i], P$rf_sd, log = TRUE),
           rf_quant = lpd_quant(y, P$draws[i, ]),
           gam = if (is.na(P$gam_mean[i])) 0 else dnorm(y, P$gam_mean[i], P$gam_sd[i], log = TRUE))
}
methods <- c("rf_gauss", "rf_quant", "gam")

# ---- Cross-validated predictions ----
folds <- sample(rep(seq_len(K), length.out = nrow(inf)))
cv <- map(seq_len(K), \(k) {
    cat("fold", k, "\n")
    tr <- inf[folds != k, ]; te <- inf[folds == k, ]
    list(rows = te$row_id,
         P = set_names(map(names(resp_def), \(r) fit_all(tr, te, r)), names(resp_def)))
})
# index: row_id -> (fold, position)
idx <- map_dfr(seq_len(K), \(k) tibble(row_id = cv[[k]]$rows, fold = k,
                                       pos = seq_along(cv[[k]]$rows)))

# ---- Test pairs: different sites, same year, both held out ----
cand <- inf |> filter(!is.na(y_small)) |> mutate(year = year(coldate)) |>
    select(row_id, trans, year)
pairs <- cand |> inner_join(cand, by = "year", suffix = c(".a", ".b"),
                            relationship = "many-to-many") |>
    filter(trans.a < trans.b) |> slice_sample(n = N_PAIRS)
cat("Test pairs:", nrow(pairs), "\n")

rowof <- \(id) inf[match(id, inf$row_id), ]
score_pair <- function(ra, rb, Pa, ia, Pb, ib, method) {
    s <- 0
    for (r in names(resp_def)) {
        ya_at_a <- log1p(total(ra, r) / ra$days); yb_at_b <- log1p(total(rb, r) / rb$days)
        yb_at_a <- log1p(total(rb, r) / ra$days); ya_at_b <- log1p(total(ra, r) / rb$days)
        true <- lpd(method, ya_at_a, Pa[[r]], ia) + lpd(method, yb_at_b, Pb[[r]], ib)
        swap <- lpd(method, yb_at_a, Pa[[r]], ia) + lpd(method, ya_at_b, Pb[[r]], ib)
        s <- s + true - swap
    }
    s
}
lr_tab <- map_dfr(seq_len(nrow(pairs)), \(q) {
    ra <- rowof(pairs$row_id.a[q]); rb <- rowof(pairs$row_id.b[q])
    xa <- idx[idx$row_id == ra$row_id, ]; xb <- idx[idx$row_id == rb$row_id, ]
    Pa <- cv[[xa$fold]]$P; Pb <- cv[[xb$fold]]$P
    tibble(q = q, method = methods,
           lr = map_dbl(methods, \(m) score_pair(ra, rb, Pa, xa$pos, Pb, xb$pos, m)))
})

cat("\n===== Overall, per method (cross-site same-year pairs, all years) =====\n")
lr_tab |> group_by(method) |>
    summarize(n = n(), accuracy = mean(lr > 0) + 0.5 * mean(lr == 0),
              confident_share = mean(abs(lr) > qlogis(0.95)),
              wrong_with_p_gt_0.95 = mean(lr < -qlogis(0.95)),
              log_score = mean(-log1p(exp(-lr))),   # mean log P(true assignment)
              .groups = "drop") |>
    mutate(across(where(is.double), \(x) round(x, 4))) |> as.data.frame() |> print()

cat("\n===== Accuracy by |log-likelihood ratio| =====\n")
lr_tab |> mutate(abs_lr = cut(abs(lr), c(0, 1, 2, 4, 8, 16, 32, Inf), right = FALSE)) |>
    group_by(method, abs_lr) |>
    summarize(n = n(), accuracy = round(mean(lr > 0), 3), .groups = "drop") |>
    pivot_wider(names_from = method, values_from = c(n, accuracy)) |>
    as.data.frame() |> print()

cat("\n===== Calibration: stated P(chosen assignment) vs. observed accuracy =====\n")
lr_tab |> mutate(p_chosen = plogis(abs(lr)),
                 bin = cut(p_chosen, c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001), right = FALSE)) |>
    group_by(method, bin) |>
    summarize(n = n(), stated = round(mean(p_chosen), 3), observed = round(mean(lr > 0), 3),
              .groups = "drop") |>
    as.data.frame() |> print()

# ---- The real pairs ----
slots1 <- tibble(trans = c("hag", "btl"), dist = c(5, 50),
                 setdate = as.Date(c("2012-06-05", "2012-06-19")),
                 coldate = as.Date(c("2012-06-19", "2012-07-03")))
slots2 <- tibble(trans = c("non", "kal"), dist = c(50, 5),
                 setdate = as.Date(c("2012-07-17", "2012-07-17")),
                 coldate = as.Date(c("2012-07-31", "2012-07-31")))
real <- function(rows, slots) {
    P <- set_names(map(names(resp_def), \(r) fit_all(inf, slots, r)), names(resp_def))
    ra <- rows[1, ] |> mutate(days = as.numeric(slots$coldate[1] - slots$setdate[1]))
    rb <- rows[2, ] |> mutate(days = as.numeric(slots$coldate[2] - slots$setdate[2]))
    map_dbl(set_names(methods), \(m) score_pair(ra, rb, P, 1, P, 2, m))
}
# accuracy of each method at the evidence strength it shows for a real pair
acc_at <- function(m, lr) {
    b <- cut(abs(lr), c(0, 1, 2, 4, 8, 16, 32, Inf), right = FALSE)
    x <- lr_tab |> filter(method == m) |>
        filter(cut(abs(lr), c(0, 1, 2, 4, 8, 16, 32, Inf), right = FALSE) == b)
    sprintf("%.3f (n = %d, bin %s)", mean(x$lr > 0), nrow(x), as.character(b))
}
cat("\n===== Real pairs: LR > 0 favours assignment 1 =====\n")
cat("Pair 1 assignment 1: line 863 -> hag 5 m 06/05-06/19, line 864 -> btl 50 m 06/19-07/03\n")
cat("Pair 2 assignment 1: line 942 -> non 50 m 07/17-07/31, line 944 -> kal 5 m 07/17-07/31\n")
for (pp in list(list("Pair 1", pair1, slots1), list("Pair 2", pair2, slots2))) {
    lr <- real(pp[[2]], pp[[3]])
    for (m in methods) cat(sprintf("%s %-9s LR = %6.2f | accuracy at this |LR| in test: %s\n",
                                   pp[[1]], m, lr[[m]], acc_at(m, lr[[m]])))
}
