
# What is lost by dropping both rows of E5 pairs 1-2, compared with the
# best-guess assignment (pair 1: line 864 -> hag 5 m, line 863 -> btl 50 m;
# pair 2: line 942 -> non 50 m, line 944 -> kal 5 m)? Run from project root:
#   Rscript review/15-e5-drop-cost.R > review/15-e5-drop-cost-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
options(width = 160)
d <- suppressMessages(read_csv("data/raw/infalls-08-17.csv", show_col_types = FALSE,
                               guess_max = 1e5)) |>
    unique() |> mutate(line = row_number() + 1) |>
    filter(lakeid == "myv", year(coldate) == 2012) |>
    mutate(m = bgch / fract_bgch + smch / fract_smch, days = as.numeric(coldate - setdate))
# `line` above is the row position after unique(), which matches the file line
# here only if no earlier rows were dropped; identify rows by values instead:
r863 <- d |> filter(trans == "hag", dist == 5, setdate == as.Date("2012-06-05"), smch == 291)
r864 <- d |> filter(trans == "hag", dist == 5, setdate == as.Date("2012-06-05"), smch == 468)
r942 <- d |> filter(trans == "non", dist == 50, setdate == as.Date("2012-07-17"), fract_smch == 1)
r944 <- d |> filter(trans == "non", dist == 50, setdate == as.Date("2012-07-17"), fract_smch < 1)
stopifnot(nrow(r863) == 1, nrow(r864) == 1, nrow(r942) == 1, nrow(r944) == 1)
base <- d |> anti_join(bind_rows(r863, r864, r942, r944), by = names(d))
mv <- \(r, t, ds, s, c) mutate(r, trans = t, dist = ds, setdate = as.Date(s), coldate = as.Date(c))
scen <- list(
    best_guess = bind_rows(base, mv(r864, "hag", 5, "2012-06-05", "2012-06-19"),
                           mv(r863, "btl", 50, "2012-06-19", "2012-07-03"),
                           mv(r942, "non", 50, "2012-07-17", "2012-07-31"),
                           mv(r944, "kal", 5, "2012-07-17", "2012-07-31")),
    drop_both = base)
plots <- tribble(~trans, ~dist, "hag", 5, "btl", 50, "non", 50, "kal", 5)

cat("== Samples with midge counts, by plot (periods are set dates) ==\n")
for (nm in names(scen)) {
    cat("\n", nm, ":\n", sep = "")
    scen[[nm]] |> semi_join(plots, by = c("trans", "dist")) |> filter(!is.na(m)) |>
        group_by(trans, dist) |>
        summarize(n = n(), days = sum(days),
                  periods = paste(format(sort(setdate), "%m/%d"), collapse = " "),
                  .groups = "drop") |>
        as.data.frame() |> print()
}

# Midge predictor scale used in the model (scripts/02): z of log(midges/midge_days)
mz_sd <- with(suppressMessages(read_csv("data/myv_arth.csv")) |> filter(taxon != "acar"),
              sd(log(midges / midge_days)))
cat("\nSD of log midge rate across all analysed plot-years (the midges_z scale):",
    round(mz_sd, 3), "\n")

rate <- \(dd) { dd <- filter(dd, !is.na(m)); sum(dd$m) / sum(dd$days) }
tab <- plots |>
    mutate(best_guess = map2_dbl(trans, dist, \(t, s) rate(filter(scen$best_guess, trans == t, dist == s))),
           drop_both = map2_dbl(trans, dist, \(t, s) rate(filter(scen$drop_both, trans == t, dist == s))),
           log_diff = log(drop_both / best_guess),
           midges_z_diff = log_diff / mz_sd)
cat("\n== 2012 midge rate (midges/day, NA samples' days excluded) ==\n")
print(tab |> mutate(across(where(is.double), \(x) round(x, 3))) |> as.data.frame())

cat("\n== Same-period small-midge totals at the other distances (for context) ==\n")
d |> filter((trans == "kal" & setdate == as.Date("2012-07-17")) |
                (trans == "hag" & setdate == as.Date("2012-06-05"))) |>
    select(trans, dist, setdate, coldate, small_total = smch, fract_smch) |>
    mutate(small_total = small_total / fract_smch) |> select(-fract_smch) |>
    as.data.frame() |> print()
