
# Code-review checks on data cleaning (scripts/00-clean-data-round1.R and
# scripts/02-clean-data.R). Run from project root:
#   Rscript review/01-check-data.R > review/01-check-data-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
options(width = 120, pillar.width = 120)

# Raw files as named in scripts/00:
f <- c(pit1 = "data/raw/pitfalls-08-17.csv",
       pit2 = "data/raw/pitfalls-11-19.csv",
       inf1 = "data/raw/infalls-08-17.csv",
       inf2 = "data/raw/infalls-13-19.csv")
cat("Files named in scripts/00 exist?\n")
print(file.exists(c("data/raw/pitfalls-08-17.csv", "data/raw/pitfalls-11-19.csv",
                    "data/raw/infalls-08-17.csv", "data/raw/infalls-13-19.csv")))

rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
pit1 <- rd(f[["pit1"]]); pit2 <- rd(f[["pit2"]])
inf1 <- rd(f[["inf1"]]); inf2 <- rd(f[["inf2"]])

# ---- Replicate scripts/00 pitfall steps ----
p1 <- pit1 |>
    select(lakeid, dist, trans, setdate, coldate, daysout, aran_other, liny,
           acar_total, lyco_total, gnap, opil, stap, cara)
p2 <- pit2 |>
    mutate(acar_total = trom_adult + trom_juv + acar,
           liny = ifelse(is.na(liny) & !is.na(aran), 0, liny)) |>
    select(lakeid, dist, trans, setdate, coldate, daysout,
           aran_other = aran, liny, acar_total, lyco_total = Lyco_NEW_total,
           gnap, opil, stap, cara)

cat("\n== Overlap check that treats NA as a mismatch (scripts/00 drops NA comparisons) ==\n")
key <- c("lakeid", "dist", "trans", "setdate", "coldate")
ov <- inner_join(p1, p2, by = key, suffix = c(".1", ".2"))
cat("Overlapping samples on (lakeid, dist, trans, setdate, coldate):", nrow(ov), "\n")
for (v in c("daysout", "aran_other", "liny", "acar_total", "lyco_total",
            "gnap", "opil", "stap", "cara")) {
    a <- ov[[paste0(v, ".1")]]; b <- ov[[paste0(v, ".2")]]
    mism <- sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    cat(sprintf("  %-11s mismatches incl. NA-vs-value: %d\n", v, mism))
}

pit <- full_join(p1, p2, by = join_by(lakeid, dist, trans, setdate, coldate, daysout,
                                      aran_other, liny, acar_total, lyco_total,
                                      gnap, opil, stap, cara))
dups <- pit |> filter(lakeid == "myv") |> count(trans, dist, setdate, coldate) |> filter(n > 1)
cat("\nmyv samples appearing >1 time after full_join (would be double-counted):",
    nrow(dups), "\n")
if (nrow(dups) > 0) print(dups, n = 30)

pit_m <- pit |>
    filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate)) |>
    mutate(year = year(coldate), sheet = aran_other + liny,
           daysout = as.numeric(coldate - setdate)) |>
    rename(lyco = lyco_total)

cat("\n== Pitfall samples with NA count but whose daysout is still added to season_days ==\n")
pit_m |>
    select(trans, dist, year, daysout, lyco, sheet, gnap, opil, stap, cara) |>
    pivot_longer(lyco:cara, names_to = "taxon", values_to = "count") |>
    group_by(taxon) |>
    summarize(n_samples = n(), n_na = sum(is.na(count)),
              na_days = sum(daysout[is.na(count)]),
              plot_years_affected = n_distinct(paste(trans, dist, year)[is.na(count)]),
              .groups = "drop") |>
    print()

cat("\nNegative or zero daysout after recomputation:\n")
print(pit_m |> filter(daysout <= 0) |> select(trans, dist, setdate, coldate, daysout))

# ---- Infall ----
i1 <- inf1 |> unique() |> select(lakeid, dist, trans, setdate, coldate, daysout,
                                 contains("bgch"), contains("smch"))
i2 <- inf2 |> unique() |> select(lakeid, dist, trans, setdate, coldate, daysout,
                                 contains("bgch"), contains("smch"))
inf <- full_join(i1, i2, by = join_by(lakeid, dist, trans, setdate, coldate, daysout,
                                      fract_bgch, bgch, fract_smch, smch))
idups <- inf |> filter(lakeid == "myv") |> count(trans, dist, setdate, coldate) |> filter(n > 1)
cat("\nmyv infall samples appearing >1 time after full_join:", nrow(idups), "\n")
if (nrow(idups) > 0) print(idups, n = 30)
inf_m <- inf |> filter(lakeid == "myv", !is.na(setdate)) |>
    mutate(bgch2 = bgch / fract_bgch, smch2 = smch / fract_smch)
cat("\nInfall samples: NA big / NA small / fract==0 or NA:\n")
print(c(na_big = sum(is.na(inf_m$bgch2)), na_small = sum(is.na(inf_m$smch2)),
        bad_fract = sum(is.na(inf_m$fract_bgch) | inf_m$fract_bgch == 0 |
                            is.na(inf_m$fract_smch) | inf_m$fract_smch == 0)))
cat("Infall distances present:\n"); print(table(inf_m$dist))
cat("Pitfall distances present:\n"); print(table(pit_m$dist))

# ---- Cleaned data actually used ----
cat("\n== data/myv_arth.csv ==\n")
ma <- rd("data/myv_arth.csv") |> filter(taxon != "acar")
cat("NA midges:", sum(is.na(ma$midges)), " zero midges:", sum(ma$midges == 0, na.rm = TRUE),
    " zero season_days:", sum(ma$season_days == 0), "\n")

cat("\nObservations per time series (trans x dist x taxon):\n")
ts <- ma |> group_by(trans, dist, taxon) |>
    summarize(n = n(), yrs = paste(range(year), collapse = "-"),
              gaps = sum(diff(sort(year)) > 1), .groups = "drop")
print(table(ts$n)); cat("Series with gaps:", sum(ts$gaps > 0), "of", nrow(ts), "\n")
print(ts |> filter(gaps > 0 | n < 12) |> distinct(trans, dist, n, yrs, gaps))

cat("\n== Pitfall sampling window by year (one row per plot-year) ==\n")
sw <- ma |> distinct(trans, dist, year, first_day, last_day, season_days,
                     midge_first_day, midge_last_day, midge_days) |>
    mutate(first_doy = yday(first_day), last_doy = yday(last_day),
           mfirst_doy = yday(midge_first_day), mlast_doy = yday(midge_last_day))
sw |> group_by(year) |>
    summarize(first_doy = median(first_doy), last_doy = median(last_doy),
              season_days = median(season_days),
              range_season_days = paste(range(season_days), collapse = "-"),
              midge_first = median(mfirst_doy), midge_last = median(mlast_doy),
              midge_days = median(midge_days), .groups = "drop") |>
    print(n = 20)
cat("\nCorrelation of season_days with year:",
    round(cor(sw$season_days, sw$year), 3), "\n")
cat("Correlation of last_doy with year:", round(cor(sw$last_doy, sw$year), 3), "\n")
cat("Correlation of first_doy with year:", round(cor(sw$first_doy, sw$year), 3), "\n")

# Does a later-ending window change the catch rate? (taxon-level, descriptive)
cat("\nWithin-taxon correlation of log1p(count/season_days) with last_doy and first_doy:\n")
ma |> mutate(rate = log1p(count / season_days), last_doy = yday(last_day),
             first_doy = yday(first_day)) |>
    group_by(taxon) |>
    summarize(r_last = cor(rate, last_doy), r_first = cor(rate, first_doy),
              r_days = cor(rate, season_days), .groups = "drop") |>
    print()
