# E4: can 2012 infall dates inform the overlapping 2012 pitfall dates?
# (1) How often do pitfall service dates match infall service dates at the same
#     plot in other years? (2) Side-by-side 2012 pitfall and infall dates.
# Run from project root:
#   Rscript review/20-e4-infall-dates.R > review/20-e4-infall-dates-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
key <- c("lakeid", "dist", "trans", "setdate", "coldate")

# Pitfall rows as built in scripts/00-clean-data-round1.R
p1 <- rd("data/raw/pitfalls-08-17.csv") |> select(all_of(key))
p2 <- rd("data/raw/pitfalls-11-19.csv") |> select(all_of(key))
pit <- bind_rows(p1, anti_join(p2, p1, by = key)) |>
    filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate)) |>
    distinct()

# Infall rows (raw, before E5/E7 corrections, so dates are as recorded)
i1 <- rd("data/raw/infalls-08-17.csv") |> select(all_of(key))
i2 <- rd("data/raw/infalls-13-19.csv") |> select(all_of(key))
inf <- bind_rows(i1, i2) |>
    filter(lakeid == "myv", !is.na(setdate)) |>
    distinct()

# ---- (1) Agreement of service dates in all years ----
# Service dates = every set or collection date for a plot-year.
svc <- \(d) d |>
    mutate(year = year(coldate)) |>
    pivot_longer(c(setdate, coldate), values_to = "date") |>
    distinct(trans, dist, year, date)
ps <- svc(pit); is <- svc(inf)
both <- inner_join(distinct(ps, trans, dist, year),
                   distinct(is, trans, dist, year), by = c("trans", "dist", "year"))
agree <- ps |>
    semi_join(both, by = c("trans", "dist", "year")) |>
    left_join(is |> rename(idate = date), by = c("trans", "dist", "year"),
              relationship = "many-to-many") |>
    group_by(trans, dist, year, date) |>
    summarize(gap = min(abs(as.numeric(date - idate))), .groups = "drop")
cat("Pitfall service dates at plot-years that also have infall data\n")
cat("(gap = days to nearest infall service date at the same plot)\n\n")
agree |>
    mutate(yr = if_else(year == 2012, "2012", "other years"),
           gap = cut(gap, c(-Inf, 0, 1, 3, Inf), labels = c("0", "1", "2-3", ">3"))) |>
    count(yr, gap) |>
    group_by(yr) |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    print()
cat("\nBy year (% of pitfall service dates matched exactly by an infall date):\n")
agree |>
    group_by(year) |>
    summarize(n = n(), exact_pct = round(100 * mean(gap == 0), 1)) |>
    print(n = Inf)

# ---- (2) 2012 side by side ----
fmt <- \(d) d |>
    filter(year(coldate) == 2012) |>
    arrange(setdate, coldate) |>
    group_by(trans, dist) |>
    summarize(dates = paste0(format(setdate, "%m/%d"), "-", format(coldate, "%m/%d"),
                             collapse = "  "), .groups = "drop")
side <- full_join(fmt(pit) |> rename(pitfall = dates),
                  fmt(inf) |> rename(infall = dates), by = c("trans", "dist")) |>
    arrange(trans, dist)
cat("\n\n2012 samples by plot (set-collection, month/day)\n")
for (i in seq_len(nrow(side))) {
    cat(sprintf("\n%s %d m\n  pit: %s\n  inf: %s\n", side$trans[i], side$dist[i],
                side$pitfall[i], side$infall[i]))
}

# Infall overlaps in 2012 (same clustering rule as script 06)
cat("\n\n2012 infall samples that overlap within a plot:\n")
inf |>
    filter(year(coldate) == 2012) |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           grp = cumsum(is.na(prev_max) | as.numeric(setdate) >= prev_max)) |>
    group_by(trans, dist, grp) |>
    filter(n() > 1) |>
    summarize(samples = paste0(format(setdate, "%m/%d"), "-", format(coldate, "%m/%d"),
                               collapse = "; "), .groups = "drop") |>
    print(n = Inf, width = Inf)
