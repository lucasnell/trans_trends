
# Checks the E7 set-date corrections in scripts/00-clean-data-round1.R, and
# that no analysed infall plot has overlapping samples left. Run from project
# root:
#   Rscript review/19-check-e7-fix.R > review/19-check-e7-fix-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
src <- readLines("scripts/00-clean-data-round1.R")
src <- src[!grepl("^\\s*write_csv", src)]
out <- capture.output(suppressMessages(eval(parse(text = src))))
cat(out[grepl("i1:|i2:", out)], sep = "\n")

cat("\n== Kálfaströnd infall samples, 2008/06-07 and 2011/07-08 ==\n")
inf |> filter(lakeid == "myv", trans == "kal", dist %in% c(5, 50, 150),
              (year(coldate) == 2008 & month(coldate) %in% 6:7) |
                  (year(coldate) == 2011 & month(coldate) %in% 7:8)) |>
    mutate(days = as.numeric(coldate - setdate)) |>
    arrange(year(coldate), dist, setdate) |>
    select(dist, setdate, coldate, days, bgch, smch) |> as.data.frame() |> print()

cat("\n== Overlapping samples remaining at analysed infall plots (all years) ==\n")
ov <- inf |> filter(lakeid == "myv", !is.na(setdate),
                    trans %in% c("btl", "hag", "kal", "non", "vin")) |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist, year = year(coldate)) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           overlap = !is.na(prev_max) & as.numeric(setdate) < prev_max) |>
    ungroup() |> filter(overlap)
print(as.data.frame(select(ov, trans, dist, setdate, coldate)))
cat("(The two 2008 E5 duplicates with all-NA copies appear here if present;",
    "their NA copies are removed later by the NA-count filter.)\n")

cat("\n== Kálfaströnd 2008 and 2011 midge days in rebuilt myv_arth ==\n")
myv_arth |> distinct(trans, dist, year, midges, midge_days) |>
    filter(trans == "kal", year %in% c(2008, 2011)) |> as.data.frame() |> print()
