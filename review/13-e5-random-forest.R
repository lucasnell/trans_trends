
# E5, pairs 1-2: does a random forest assign the duplicated 2012 infall rows
# more reliably than the per-site multiplicative model in
# review/12-e5-assign-pairs.R?
#
# Data: both infall files combined as in scripts/00-clean-data-round1.R (with
# the E6 row removed), Mývatn only, all years 2008-2019. The four rows of E5
# pairs 1-2 are removed; for E5 pairs 3-4 only the all-NA copies are removed.
#
# Response (one forest per response): log1p(total / days) for big midges,
# small midges, and simuliids, where total = count / subsample fraction.
#
# Predictors for a sample (plot = site x distance):
#   site, log(distance), year, midpoint day-of-year, days out;
#   mean response at the other distances of the same site and dates, and how
#   many there are; response of the previous and next sample at the same plot
#   (NA allowed: ranger's na.action = "na.learn").
# Features are rebuilt after removing held-out rows, so a held-out sample's
# own value never appears in any training row's features.
#
# Scoring an assignment: Gaussian log-likelihood of each candidate row's
# response at its slot, with SD = out-of-bag residual SD of that forest;
# summed over responses. Calibration: as in script 12, every pair of known
# 2012 samples (one from each site) is held out and scored true vs. swapped.
#
# `ranger` isn't in this project's renv library; it's loaded from the global
# R library (appended to the library path, so renv packages take precedence).
# Review-only dependency; nothing is installed.
#
# Run from project root:
#   Rscript review/13-e5-random-forest.R > review/13-e5-random-forest-output.txt 2>&1

.libPaths(c(.libPaths(), "/Users/lucasnell/Library/R/arm64/4.6/library"))
suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
    library(ranger)
})
options(width = 160)
set.seed(20261123)
NTREE <- 500

key <- c("lakeid", "dist", "trans", "setdate", "coldate")
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
cols <- c(key, "fract_bgch", "bgch", "fract_smch", "smch", "fract_simu", "simu")
i1 <- rd("data/raw/infalls-08-17.csv") |> unique() |> select(all_of(cols))
i2 <- rd("data/raw/infalls-13-19.csv") |> unique() |> select(all_of(cols)) |>
    filter(!(trans == "btl" & dist == 5 & setdate == as.Date("2013-06-05") &
                 coldate == as.Date("2013-08-14")))
inf <- bind_rows(i1, anti_join(i2, i1, by = key)) |>
    filter(lakeid == "myv", !is.na(setdate), !is.na(coldate), !is.na(dist)) |>
    mutate(row_id = row_number())

# E5 rows
pair1 <- inf |> filter(trans == "hag", dist == 5, setdate == as.Date("2012-06-05"),
                       coldate == as.Date("2012-06-19"))
pair2 <- inf |> filter(trans == "non", dist == 50, setdate == as.Date("2012-07-17"),
                       coldate == as.Date("2012-07-31"))
stopifnot(nrow(pair1) == 2, nrow(pair2) == 2)
pair1 <- arrange(pair1, smch); pair2 <- arrange(pair2, fract_smch < 1)  # row 1 = line 863 / 942
na_copies <- inf |>
    filter((trans == "non" & dist == 150 & setdate == as.Date("2008-06-04")) |
               (trans == "vin" & dist == 5 & setdate == as.Date("2008-05-30"))) |>
    filter(is.na(smch))
stopifnot(nrow(na_copies) == 2)
inf <- inf |> filter(!row_id %in% c(pair1$row_id, pair2$row_id, na_copies$row_id))

resp_def <- list(big = c("bgch", "fract_bgch"), small = c("smch", "fract_smch"),
                 simu = c("simu", "fract_simu"))
add_y <- function(dd) {
    dd <- dd |> mutate(days = as.numeric(coldate - setdate))
    for (r in names(resp_def)) {
        y <- resp_def[[r]][1]; f <- resp_def[[r]][2]
        dd[[paste0("y_", r)]] <- log1p((dd[[y]] / dd[[f]]) / dd$days)
    }
    dd
}
inf <- add_y(inf)

# Build features for `cells` (rows to describe) from `known` (rows whose
# responses may be used). Cells need: trans, dist, setdate, coldate.
features <- function(cells, known, r) {
    yk <- paste0("y_", r)
    kn <- known |> filter(!is.na(.data[[yk]])) |>
        select(trans, dist, setdate, coldate, yv = all_of(yk))
    cells |>
        mutate(cell = row_number(),
               days = as.numeric(coldate - setdate),
               mid_doy = yday(setdate) + days / 2,
               year = year(coldate),
               ldist = log(dist),
               site = factor(trans, levels = sort(unique(known$trans)))) |>
        (\(x) {
            # other distances, same site and dates
            oth <- x |> select(cell, trans, dist, setdate, coldate) |>
                inner_join(kn, by = c("trans", "setdate", "coldate"),
                           suffix = c("", ".o"), relationship = "many-to-many") |>
                filter(dist.o != dist) |>
                group_by(cell) |> summarize(other_mean = mean(yv), other_n = n(),
                                            .groups = "drop")
            # previous and next sample at the same plot (by dates, excluding itself)
            nb <- x |> select(cell, trans, dist, setdate, coldate) |>
                inner_join(kn, by = c("trans", "dist"), suffix = c("", ".k"),
                           relationship = "many-to-many") |>
                filter(year(coldate.k) == year(coldate),
                       !(setdate.k == setdate & coldate.k == coldate)) |>
                group_by(cell) |>
                summarize(prev = {
                    i <- which(coldate.k <= setdate + 3)
                    if (length(i)) yv[i][which.max(coldate.k[i])] else NA_real_
                }, nxt = {
                    i <- which(setdate.k >= coldate - 3)
                    if (length(i)) yv[i][which.min(setdate.k[i])] else NA_real_
                }, .groups = "drop")
            x |> left_join(oth, by = "cell") |> left_join(nb, by = "cell") |>
                mutate(other_n = replace_na(other_n, 0L))
        })()
}
fvars <- c("site", "ldist", "year", "mid_doy", "days", "other_mean", "other_n",
           "prev", "nxt")

# Train one forest on `train` (features built from train itself) and predict
# `newcells`. Returns predictions and OOB residual SD.
fit_predict <- function(train, newcells, r, seed = 1) {
    yk <- paste0("y_", r)
    tr <- features(train, train, r) |>
        mutate(y = train[[yk]]) |> filter(!is.na(y))
    m <- ranger(x = as.data.frame(tr[fvars]), y = tr$y, num.trees = NTREE,
                na.action = "na.learn", seed = seed)
    nd <- features(newcells, train, r)
    nd$site <- factor(as.character(nd$site), levels = levels(tr$site))
    list(pred = predict(m, data = as.data.frame(nd[fvars]))$predictions,
         sd = sqrt(m$prediction.error))
}

# Gaussian log density of a candidate row's response at a slot
cand_y <- function(row, r, days) {
    y <- row[[resp_def[[r]][1]]]; f <- row[[resp_def[[r]][2]]]
    log1p((y / f) / days)
}
ldn <- function(y, mu, s) if (is.na(y) || is.na(mu)) 0 else dnorm(y, mu, s, log = TRUE)

score <- function(rows, slots, train, seed = 1) {
    map_dfr(names(resp_def), \(r) {
        fp <- fit_predict(train, slots, r, seed)
        d <- as.numeric(slots$coldate - slots$setdate)
        ya <- c(cand_y(rows[1, ], r, d[1]), cand_y(rows[2, ], r, d[2]))  # assignment 1
        yb <- c(cand_y(rows[2, ], r, d[1]), cand_y(rows[1, ], r, d[2]))  # assignment 2
        tibble(response = r, pred_slot1 = fp$pred[1], pred_slot2 = fp$pred[2],
               oob_sd = fp$sd,
               row1_y = cand_y(rows[1, ], r, d[1]), row2_y = cand_y(rows[2, ], r, d[2]),
               loglik_assign1 = ldn(ya[1], fp$pred[1], fp$sd) + ldn(ya[2], fp$pred[2], fp$sd),
               loglik_assign2 = ldn(yb[1], fp$pred[1], fp$sd) + ldn(yb[2], fp$pred[2], fp$sd))
    })
}

slot_rows <- function(...) tribble(~trans, ~dist, ~setdate, ~coldate, ...)
s1 <- tibble(trans = c("hag", "btl"), dist = c(5, 50),
             setdate = as.Date(c("2012-06-05", "2012-06-19")),
             coldate = as.Date(c("2012-06-19", "2012-07-03")))
s2 <- tibble(trans = c("non", "kal"), dist = c(50, 5),
             setdate = as.Date(c("2012-07-17", "2012-07-17")),
             coldate = as.Date(c("2012-07-31", "2012-07-31")))

for (pp in list(list("Pair 1: slot 1 = hag 5 m 06/05-06/19, slot 2 = btl 50 m 06/19-07/03",
                     pair1, s1),
                list("Pair 2: slot 1 = non 50 m 07/17-07/31, slot 2 = kal 5 m 07/17-07/31",
                     pair2, s2))) {
    cat("\n=====", pp[[1]], "=====\n")
    cat("Row 1: big", pp[[2]]$bgch[1], "small", pp[[2]]$smch[1], "fract", pp[[2]]$fract_smch[1],
        "simu", pp[[2]]$simu[1], "\n")
    cat("Row 2: big", pp[[2]]$bgch[2], "small", pp[[2]]$smch[2], "fract", pp[[2]]$fract_smch[2],
        "simu", pp[[2]]$simu[2], "\n")
    cat("Assignment 1 = row 1 -> slot 1, row 2 -> slot 2. Values are log1p(per-day rate).\n")
    sc <- score(pp[[2]], pp[[3]], inf)
    print(sc |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> as.data.frame())
    lr <- sum(sc$loglik_assign1) - sum(sc$loglik_assign2)
    cat(sprintf("Combined log-likelihood ratio (assignment 1 vs 2): %.2f; plug-in P(assignment 1) = %.3f\n",
                lr, plogis(lr)))
}

# ---- Calibration on known 2012 samples ----
# Leave-one-out prediction for each known 2012 sample at the four sites, then
# score every cross-site pair true vs. swapped.
cat("\n===== Calibration on known 2012 samples (leave-one-out forests) =====\n")
test <- inf |> filter(year(coldate) == 2012, trans %in% c("hag", "btl", "non", "kal"),
                      !is.na(y_small))
loo <- map_dfr(seq_len(nrow(test)), \(i) {
    tr <- filter(inf, row_id != test$row_id[i])
    cell <- test[i, c("trans", "dist", "setdate", "coldate")]
    map_dfr(names(resp_def), \(r) {
        fp <- fit_predict(tr, cell, r)
        tibble(row_id = test$row_id[i], response = r, pred = fp$pred, sd = fp$sd)
    })
})
loo_w <- loo |> left_join(select(test, row_id, trans, dist, setdate, coldate,
                                 all_of(unname(unlist(resp_def)))), by = "row_id")

cal <- function(site1, site2) {
    a <- filter(test, trans == site1); b <- filter(test, trans == site2)
    res <- expand_grid(i = seq_len(nrow(a)), j = seq_len(nrow(b))) |>
        mutate(lr = map2_dbl(i, j, \(i, j) {
            ra <- a[i, ]; rb <- b[j, ]
            da <- as.numeric(ra$coldate - ra$setdate); db <- as.numeric(rb$coldate - rb$setdate)
            s <- 0
            for (r in names(resp_def)) {
                pa <- filter(loo, row_id == ra$row_id, response == r)
                pb <- filter(loo, row_id == rb$row_id, response == r)
                true <- ldn(cand_y(ra, r, da), pa$pred, pa$sd) + ldn(cand_y(rb, r, db), pb$pred, pb$sd)
                swap <- ldn(cand_y(rb, r, da), pa$pred, pa$sd) + ldn(cand_y(ra, r, db), pb$pred, pb$sd)
                s <- s + true - swap
            }
            s
        }))
    binned <- res |>
        mutate(abs_lr = cut(abs(lr), c(0, 1, 2, 4, 8, 16, 32, Inf), right = FALSE)) |>
        group_by(abs_lr) |> summarize(n = n(), accuracy = mean(lr > 0), .groups = "drop")
    cat("\n", site1, " + ", site2, ": accuracy by |log-likelihood ratio|\n", sep = "")
    print(as.data.frame(binned))
    tibble(sites = paste(site1, site2, sep = " + "), n_pairs = nrow(res),
           correct = mean(res$lr > 0),
           confident_share = mean(plogis(res$lr) > 0.95 | plogis(res$lr) < 0.05),
           wrong_with_p_gt_0.95 = mean(plogis(res$lr) < 0.05))
}
summ <- bind_rows(cal("hag", "btl"), cal("non", "kal"))
cat("\n")
print(summ |> mutate(across(where(is.double), \(x) round(x, 3))) |> as.data.frame())

cat("\nLeave-one-out prediction error (RMSE, log1p scale) on these 2012 samples:\n")
loo_w |> group_by(response) |>
    summarize(rmse = sqrt(mean((pred - case_when(
        response == "big" ~ log1p(bgch / fract_bgch / as.numeric(coldate - setdate)),
        response == "small" ~ log1p(smch / fract_smch / as.numeric(coldate - setdate)),
        TRUE ~ log1p(simu / fract_simu / as.numeric(coldate - setdate))))^2, na.rm = TRUE)),
        oob_sd = mean(sd), .groups = "drop") |>
    as.data.frame() |> print()

# ---- Stability of the two pair results across forest seeds ----
cat("\n===== Pair log-likelihood ratios across 10 forest seeds =====\n")
stab <- map_dfr(1:10, \(sd_) {
    tibble(seed = sd_,
           pair1 = with(score(pair1, s1, inf, sd_), sum(loglik_assign1) - sum(loglik_assign2)),
           pair2 = with(score(pair2, s2, inf, sd_), sum(loglik_assign1) - sum(loglik_assign2)))
})
print(stab |> mutate(across(c(pair1, pair2), \(x) round(x, 2))) |> as.data.frame())

# ---- How typical are the candidates at the hag 5 m slot? ----
cat("\n===== Pair 1 candidates at the hag 5 m slot: (observed - predicted) / OOB SD =====\n")
sc1 <- score(pair1, s1, inf)
print(sc1 |> transmute(response,
                       row1_at_hag5 = round((row1_y - pred_slot1) / oob_sd, 2),
                       row2_at_hag5 = round((row2_y - pred_slot1) / oob_sd, 2),
                       row1_at_btl50 = round((row1_y - pred_slot2) / oob_sd, 2),
                       row2_at_btl50 = round((row2_y - pred_slot2) / oob_sd, 2)) |>
          as.data.frame())

# ---- Is hag 5 m generally over-predicted? Leave-one-out z for its known 2012 samples ----
cat("\n===== Known 2012 hag 5 m samples: (observed - LOO predicted) / OOB SD =====\n")
loo_w |> filter(trans == "hag", dist == 5) |>
    mutate(days = as.numeric(coldate - setdate),
           obs = case_when(response == "big" ~ log1p(bgch / fract_bgch / days),
                           response == "small" ~ log1p(smch / fract_smch / days),
                           TRUE ~ log1p(simu / fract_simu / days)),
           z = round((obs - pred) / sd, 2)) |>
    select(setdate, response, z) |>
    pivot_wider(names_from = response, values_from = z) |>
    as.data.frame() |> print()
