
# E5, pairs 1-2: use other distances at the same site and dates to judge which
# copy of each duplicated 2012 infall sample belongs where.
#
# Pair 1: two rows labelled hag 5 m, 2012/06/05-06/19. The curation note says
#   one may be btl 50 m, 2012/06/19-07/03 (that sample is missing).
# Pair 2: two rows labelled non 50 m, 2012/07/17-07/31. The curation note says
#   one may be kal 5 m, 2012/07/17-07/31 (that sample is missing).
#
# Model (fit separately per site, 2012 only, per response: big midges `bgch`,
# small midges `smch`, simuliids `simu`):
#   raw count ~ NegBin(mean = fract * exp(a_distance + b_period), theta)
# i.e. within a site, each distance is a constant multiple of the others across
# sampling periods. `fract` is the subsample fraction for that response.
# Both rows of the ambiguous pair are excluded from fitting. Each candidate
# assignment is scored by the plug-in NegBin log-likelihood of the candidate
# rows at the two empty cells; posterior probability assumes equal prior odds.
#
# Calibration: the same procedure applied to pairs of *known* cells (one from
# each of the two sites), each predicted leave-one-out, shows how often the
# method picks the correct assignment and how often it is confidently wrong.
#
# Run from project root:
#   Rscript review/12-e5-assign-pairs.R > review/12-e5-assign-pairs-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(MASS, exclude = "select")
})
options(width = 160)
set.seed(20261123)

d <- suppressMessages(read_csv("data/raw/infalls-08-17.csv", show_col_types = FALSE,
                               guess_max = 1e5)) |>
    unique() |>
    mutate(line = row_number()) |>   # identifies rows within this script only
    filter(lakeid == "myv", format(coldate, "%Y") == "2012",
           trans %in% c("btl", "hag", "kal", "non")) |>
    mutate(period = format(setdate, "%m/%d"))

resp <- list(big = c("bgch", "fract_bgch"), small = c("smch", "fract_smch"),
             simu = c("simu", "fract_simu"))

# Fit the per-site model for one response; return function giving the
# predicted total (before subsampling) for a distance/period cell.
fit_site <- function(dd, y, f) {
    dd <- dd |> filter(!is.na(.data[[y]]), !is.na(.data[[f]])) |>
        mutate(yy = round(.data[[y]]), off = log(.data[[f]]),
               distf = factor(dist), periodf = factor(period))
    m <- tryCatch(suppressWarnings(glm.nb(yy ~ distf + periodf + offset(off), data = dd)),
                  error = \(e) NULL)
    if (is.null(m)) return(NULL)
    list(model = m, theta = m$theta, data = dd)
}
pred_total <- function(fit, dist, period) {
    nd <- tibble(distf = factor(dist, levels = levels(fit$data$distf)),
                 periodf = factor(period, levels = levels(fit$data$periodf)),
                 off = 0)
    if (anyNA(nd$distf) || anyNA(nd$periodf)) return(NA_real_)
    as.numeric(predict(fit$model, newdata = nd, type = "response"))
}
# log density of a candidate row's raw count at a cell with predicted total mu
ld <- function(count, fract, mu, theta) {
    if (is.na(count) || is.na(fract) || is.na(mu)) return(0)  # no information
    dnbinom(round(count), size = theta, mu = max(mu * fract, 0.01), log = TRUE)
}

# Score assignments for an ambiguous pair.
# rows: the two ambiguous rows; slots: list of two (site, dist, period).
score_pair <- function(rows, slots, data_excl) {
    out <- map_dfr(names(resp), \(r) {
        y <- resp[[r]][1]; f <- resp[[r]][2]
        fits <- map(slots, \(s) fit_site(filter(data_excl, trans == s$site), y, f))
        if (any(map_lgl(fits, is.null))) return(NULL)
        mu <- map2_dbl(fits, slots, \(ft, s) pred_total(ft, s$dist, s$period))
        th <- map_dbl(fits, "theta")
        # assignment 1: row 1 -> slot 1, row 2 -> slot 2; assignment 2: swapped
        l1 <- ld(rows[[y]][1], rows[[f]][1], mu[1], th[1]) +
            ld(rows[[y]][2], rows[[f]][2], mu[2], th[2])
        l2 <- ld(rows[[y]][2], rows[[f]][2], mu[1], th[1]) +
            ld(rows[[y]][1], rows[[f]][1], mu[2], th[2])
        tibble(response = r,
               pred_slot1 = mu[1], pred_slot2 = mu[2],
               row1_total = rows[[y]][1] / rows[[f]][1],
               row2_total = rows[[y]][2] / rows[[f]][2],
               loglik_assign1 = l1, loglik_assign2 = l2)
    })
    out
}

report <- function(name, rows, slots, slot_labels) {
    data_excl <- anti_join(d, rows, by = "line")
    sc <- score_pair(rows, slots, data_excl)
    cat("\n=====", name, "=====\n")
    cat("Row 1:", sprintf("big %s, small %s (fract %s), simu %s",
                          rows$bgch[1], rows$smch[1], rows$fract_smch[1], rows$simu[1]), "\n")
    cat("Row 2:", sprintf("big %s, small %s (fract %s), simu %s",
                          rows$bgch[2], rows$smch[2], rows$fract_smch[2], rows$simu[2]), "\n")
    cat("Slot 1:", slot_labels[1], "| Slot 2:", slot_labels[2], "\n")
    cat("Assignment 1 = row 1 -> slot 1, row 2 -> slot 2; assignment 2 = swapped.\n")
    cat("pred_slot = model's predicted total (after undoing subsampling) at each empty cell.\n\n")
    print(sc |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> as.data.frame())
    lr <- sum(sc$loglik_assign1) - sum(sc$loglik_assign2)
    cat(sprintf("\nCombined log-likelihood ratio (assignment 1 vs 2): %.2f\n", lr))
    cat(sprintf("P(assignment 1) with equal prior odds (plug-in, ignores parameter uncertainty): %.3f\n",
                plogis(lr)))
    for (r in sc$response) {
        x <- filter(sc, response == r)
        cat(sprintf("  %-5s alone: P(assignment 1) = %.3f\n", r,
                    plogis(x$loglik_assign1 - x$loglik_assign2)))
    }
    invisible(lr)
}

# ---- Pair 1 ----
p1 <- d |> filter(trans == "hag", dist == 5, period == "06/05") |> arrange(line)
stopifnot(nrow(p1) == 2)
report("Pair 1: hag 5 m 06/05-06/19 vs btl 50 m 06/19-07/03", p1,
       list(list(site = "hag", dist = 5, period = "06/05"),
            list(site = "btl", dist = 50, period = "06/19")),
       c("hag 5 m 06/05-06/19", "btl 50 m 06/19-07/03"))

# ---- Pair 2 ----
p2 <- d |> filter(trans == "non", dist == 50, period == "07/17") |> arrange(line)
stopifnot(nrow(p2) == 2)
report("Pair 2: non 50 m 07/17-07/31 vs kal 5 m 07/17-07/31", p2,
       list(list(site = "non", dist = 50, period = "07/17"),
            list(site = "kal", dist = 5, period = "07/17")),
       c("non 50 m 07/17-07/31", "kal 5 m 07/17-07/31"))

# ---- Calibration on known cells ----
# For every pair of known cells (one at each site of the pair), hide both,
# refit, and score the true vs. swapped assignment exactly as above.
calibrate <- function(site1, site2, excl_rows) {
    base <- anti_join(d, excl_rows, by = "line")
    c1 <- filter(base, trans == site1); c2 <- filter(base, trans == site2)
    # leave-one-out log densities per cell, per response, for its own row and
    # for each candidate row from the other site
    loo <- function(cells, site) {
        map(seq_len(nrow(cells)), \(i) {
            fits <- map(names(resp), \(r) {
                y <- resp[[r]][1]; f <- resp[[r]][2]
                ft <- fit_site(filter(base, trans == site, line != cells$line[i]), y, f)
                if (is.null(ft)) return(NULL)
                list(mu = pred_total(ft, cells$dist[i], cells$period[i]), theta = ft$theta)
            })
            set_names(fits, names(resp))
        })
    }
    L1 <- loo(c1, site1); L2 <- loo(c2, site2)
    res <- expand_grid(i = seq_len(nrow(c1)), j = seq_len(nrow(c2))) |>
        mutate(lr = map2_dbl(i, j, \(i, j) {
            s <- 0
            for (r in names(resp)) {
                y <- resp[[r]][1]; f <- resp[[r]][2]
                a <- L1[[i]][[r]]; b <- L2[[j]][[r]]
                if (is.null(a) || is.null(b)) next
                true <- ld(c1[[y]][i], c1[[f]][i], a$mu, a$theta) +
                    ld(c2[[y]][j], c2[[f]][j], b$mu, b$theta)
                swap <- ld(c2[[y]][j], c2[[f]][j], a$mu, a$theta) +
                    ld(c1[[y]][i], c1[[f]][i], b$mu, b$theta)
                s <- s + true - swap
            }
            s
        })) |>
        mutate(p_true = plogis(lr))
    binned <- res |>
        mutate(abs_lr = cut(abs(lr), c(0, 1, 2, 4, 8, 16, 32, Inf), right = FALSE)) |>
        group_by(abs_lr) |>
        summarize(n = n(), accuracy = mean(lr > 0), .groups = "drop")
    cat("\n", site1, " + ", site2, ": accuracy by |log-likelihood ratio|\n", sep = "")
    print(as.data.frame(binned))
    tibble(sites = paste(site1, site2, sep = " + "), n_pairs = nrow(res),
           correct = mean(res$p_true > 0.5),
           confident_share = mean(res$p_true > 0.95 | res$p_true < 0.05),
           correct_when_confident = mean(res$p_true[res$p_true > 0.95 | res$p_true < 0.05] > 0.5),
           wrong_with_p_gt_0.95 = mean(res$p_true < 0.05))
}
cat("\n===== Calibration on known 2012 cells (leave-one-out) =====\n")
cat("correct = share of known pairs assigned correctly;\n",
    "confident_share = share with P > 0.95 for one assignment;\n",
    "correct_when_confident = accuracy among those;\n",
    "wrong_with_p_gt_0.95 = share of all pairs confidently wrong.\n", sep = "")
cal <- bind_rows(calibrate("hag", "btl", p1), calibrate("non", "kal", p2))
print(cal |> mutate(across(where(is.double), \(x) round(x, 3))) |> as.data.frame())

# ---- Does the choice matter? Annual midge catch rate per plot-year ----
# Rate = sum(big/fract_big + small/fract_small) / sum(days) over samples with
# non-NA counts (i.e. with the E3 fix), for the four affected 2012 plots.
cat("\n===== Annual 2012 midge rate (midges/day) under each choice =====\n")
rate <- function(dd) {
    dd <- dd |> mutate(m = bgch / fract_bgch + smch / fract_smch,
                       days = as.numeric(coldate - setdate)) |>
        filter(!is.na(m))
    sum(dd$m) / sum(dd$days)
}
base <- anti_join(d, bind_rows(p1, p2), by = "line")
relabel <- function(r, site, new_dist, set, col) {
    mutate(r, trans = site, dist = new_dist, setdate = as.Date(set), coldate = as.Date(col))
}
choices <- list(
    "assignment 1" = bind_rows(
        relabel(p1[1, ], "hag", 5, "2012-06-05", "2012-06-19"),
        relabel(p1[2, ], "btl", 50, "2012-06-19", "2012-07-03"),
        relabel(p2[1, ], "non", 50, "2012-07-17", "2012-07-31"),
        relabel(p2[2, ], "kal", 5, "2012-07-17", "2012-07-31")),
    "assignment 2" = bind_rows(
        relabel(p1[2, ], "hag", 5, "2012-06-05", "2012-06-19"),
        relabel(p1[1, ], "btl", 50, "2012-06-19", "2012-07-03"),
        relabel(p2[2, ], "non", 50, "2012-07-17", "2012-07-31"),
        relabel(p2[1, ], "kal", 5, "2012-07-17", "2012-07-31")),
    "drop both rows" = p1[0, ],
    "current cleaning (both kept at original label)" = bind_rows(p1, p2))
plots <- tribble(~trans, ~dist, "hag", 5, "btl", 50, "non", 50, "kal", 5)
tab <- map_dfr(names(choices), \(nm) {
    dd <- bind_rows(base, choices[[nm]])
    plots |> mutate(choice = nm,
                    rate = map2_dbl(trans, dist, \(t, s) rate(filter(dd, trans == t, dist == s))))
}) |>
    mutate(plot = paste(trans, dist, "m")) |>
    select(choice, plot, rate) |>
    pivot_wider(names_from = plot, values_from = rate)
print(tab |> mutate(across(where(is.numeric), \(x) round(x, 1))) |> as.data.frame())
cat("\nSD of log(midge rate) across all 2012 plot-years at these four sites:",
    round(sd(log(map_dbl(split(base, ~ trans + dist, drop = TRUE), rate))), 2), "\n")
