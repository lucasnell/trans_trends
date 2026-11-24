
# Code-review check: how much do the cleaning problems found in
# review/01-check-data.R change plot-year catch rates, and how sensitive are
# simple (non-model) trends to the within-season sampling window?
# Run from project root:
#   Rscript review/02-check-corrected-cleaning.R > review/02-check-corrected-cleaning-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
options(width = 120, pillar.width = 120)
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
taxa <- c("cara", "stap", "opil", "gnap", "sheet", "lyco")

p1 <- rd("data/raw/pitfalls-08-17.csv") |>
    select(lakeid, dist, trans, setdate, coldate, aran_other, liny,
           lyco = lyco_total, gnap, opil, stap, cara)
p2 <- rd("data/raw/pitfalls-11-19.csv") |>
    mutate(liny = ifelse(is.na(liny) & !is.na(aran), 0, liny)) |>
    select(lakeid, dist, trans, setdate, coldate, aran_other = aran, liny,
           lyco = Lyco_NEW_total, gnap, opil, stap, cara)
key <- c("lakeid", "dist", "trans", "setdate", "coldate")

# Corrected: one row per physical sample (keep p1 version where both exist)
samples <- bind_rows(p1, anti_join(p2, p1, by = key)) |>
    filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate)) |>
    mutate(year = year(coldate), sheet = aran_other + liny,
           daysout = as.numeric(coldate - setdate),
           set_doy = yday(setdate), col_doy = yday(coldate)) |>
    select(trans, dist, year, setdate, coldate, set_doy, col_doy, daysout, all_of(taxa))
stopifnot(nrow(count(samples, trans, dist, setdate, coldate) |> filter(n > 1)) == 0)

# Overlapping deployments within a plot (same trap counted twice in time)?
ovl <- samples |> arrange(trans, dist, setdate) |> group_by(trans, dist, year) |>
    mutate(overlap = setdate < lag(coldate)) |> ungroup()
cat("Samples whose setdate precedes previous sample's coldate (same plot):",
    sum(ovl$overlap, na.rm = TRUE), "\n")
print(ovl |> filter(overlap) |> count(year))

long <- samples |>
    pivot_longer(all_of(taxa), names_to = "taxon", values_to = "count")

corrected <- long |>
    filter(!is.na(count)) |>          # drop NA samples' days too
    group_by(trans, dist, year, taxon) |>
    summarize(count = sum(round(count)), season_days = sum(daysout), .groups = "drop")

orig <- rd("data/myv_arth.csv") |> filter(taxon != "acar") |>
    select(trans, dist, year, taxon, count0 = count, days0 = season_days)

cmp <- inner_join(orig, corrected, by = c("trans", "dist", "year", "taxon")) |>
    mutate(lr0 = log1p(count0 / days0), lr1 = log1p(count / season_days),
           diff = lr1 - lr0)
cat("\nRows compared:", nrow(cmp), " (orig rows:", nrow(orig), ")\n")
cat("Ratio of original to corrected season_days, by year:\n")
cmp |> group_by(year) |> summarize(days_ratio = median(days0 / season_days),
                                    count_ratio = median(count0 / pmax(count, 1)),
                                    .groups = "drop") |> print(n = 20)
cat("\n|change in log1p(rate)| > 0.05:", sum(abs(cmp$diff) > 0.05), "of", nrow(cmp), "\n")
print(cmp |> filter(abs(diff) > 0.05) |> arrange(desc(abs(diff))) |>
          select(trans, dist, year, taxon, count0, days0, count, season_days, diff), n = 40)

# ---- Sampling-window sensitivity (descriptive, not the paper's model) ----
# Within-plot OLS slope of log1p(rate) on year, per taxon:
#   (a) all samples, (b) only samples fully inside a window sampled every year.
win <- samples |> group_by(year) |>
    summarize(f = min(set_doy), l = max(col_doy), .groups = "drop")
cat("\nPer-year earliest set DOY and latest collect DOY:\n"); print(win, n = 20)
lo <- 160; hi <- 205
cat(sprintf("Common window used: DOY %d-%d (samples set >= %d and collected <= %d)\n",
            lo, hi, lo, hi))
slope_fun <- function(dd) {
    dd |> group_by(taxon) |> group_modify(\(x, ...) {
        f <- lm(lr ~ year + factor(paste(trans, dist)), data = x)
        tibble(slope = coef(f)[["year"]], se = sqrt(vcov(f)["year", "year"]),
               n = nrow(x))
    }) |> ungroup()
}
rate_tbl <- function(ll) ll |> filter(!is.na(count)) |>
    group_by(trans, dist, year, taxon) |>
    summarize(lr = log(sum(round(count) + 0.5) / sum(daysout)), .groups = "drop")
a <- slope_fun(rate_tbl(long))
b <- slope_fun(rate_tbl(long |> filter(set_doy >= lo, col_doy <= hi)))
cat("\nWithin-plot slope of log(rate) per year (OLS; SEs ignore autocorrelation):\n")
print(left_join(a, b, by = "taxon", suffix = c("_all", "_window")) |>
          mutate(across(where(is.double), \(x) round(x, 4))))
cat("Plot-years retained in window subset:",
    nrow(distinct(filter(long, set_doy >= lo, col_doy <= hi), trans, dist, year)),
    "of", nrow(distinct(long, trans, dist, year)), "\n")
cat("Years with no samples inside the window:\n")
print(setdiff(unique(samples$year),
              unique(filter(samples, set_doy >= lo, col_doy <= hi)$year)))

# ---- Infall: NA midge counts still contributing days ----
i1 <- rd("data/raw/infalls-08-17.csv"); i2 <- rd("data/raw/infalls-13-19.csv")
ii <- bind_rows(i1 |> select(lakeid, dist, trans, setdate, coldate, fract_bgch, bgch, fract_smch, smch),
                i2 |> select(lakeid, dist, trans, setdate, coldate, fract_bgch, bgch, fract_smch, smch)) |>
    distinct(lakeid, dist, trans, setdate, coldate, .keep_all = TRUE) |>
    filter(lakeid == "myv", !is.na(setdate)) |>
    mutate(year = year(coldate), m = bgch / fract_bgch + smch / fract_smch,
           days = as.numeric(coldate - setdate))
pit_plots <- distinct(samples, trans, dist, year)
na_py <- ii |> filter(is.na(m)) |> distinct(trans, dist, year) |> semi_join(pit_plots, by = c("trans", "dist", "year"))
cat("\nInfall samples with NA midge total:", sum(is.na(ii$m)),
    "; analysed plot-years affected:", nrow(na_py), "of", nrow(pit_plots), "\n")
mc <- ii |> group_by(trans, dist, year) |>
    summarize(r_orig = sum(m, na.rm = TRUE) / sum(days),
              r_fix = sum(m, na.rm = TRUE) / sum(days[!is.na(m)]), .groups = "drop") |>
    semi_join(pit_plots, by = c("trans", "dist", "year")) |>
    mutate(dlog = log(r_fix) - log(r_orig))
cat("Change in log midge rate from dropping NA samples' days: max",
    round(max(mc$dlog), 3), "; >0.1 in", sum(mc$dlog > 0.1), "plot-years\n")
print(mc |> filter(dlog > 0.1) |> arrange(desc(dlog)), n = 30)
