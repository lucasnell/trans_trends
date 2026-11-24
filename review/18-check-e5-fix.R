
# Checks the E5 pair 1-2 reassignment in scripts/00-clean-data-round1.R.
# Expected (review/15, best guess): 2012 midges/day hag 5 m 155.5,
# btl 50 m 46.9, non 50 m 40.6, kal 5 m 77.8; six samples at hag 5 m,
# btl 50 m, and non 50 m, four at kal 5 m; no remaining 2012 overlaps at
# these plots. Run from project root:
#   Rscript review/18-check-e5-fix.R > review/18-check-e5-fix-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
src <- readLines("scripts/00-clean-data-round1.R")
src <- src[!grepl("^\\s*write_csv", src)]
out <- capture.output(suppressMessages(eval(parse(text = src))))
cat(out[grepl("i1:|i2:|p1:|p2:", out)], sep = "\n")

plots <- tribble(~trans, ~dist, "hag", 5, "btl", 50, "non", 50, "kal", 5)
cat("\n== 2012 infall samples kept at the four plots ==\n")
inf |> filter(lakeid == "myv", year(coldate) == 2012, !is.na(setdate)) |>
    semi_join(plots, by = c("trans", "dist")) |>
    mutate(midges = bgch / fract_bgch + smch / fract_smch) |>
    filter(!is.na(midges)) |>
    arrange(trans, dist, setdate) |>
    group_by(trans, dist) |>
    summarize(n = n(), periods = paste(format(setdate, "%m/%d"), collapse = " "),
              overlaps = sum(setdate < lag(coldate), na.rm = TRUE), .groups = "drop") |>
    as.data.frame() |> print()

cat("\n== 2012 midge rates in rebuilt myv_arth ==\n")
myv_arth |> distinct(trans, dist, year, midges, midge_days) |>
    filter(year == 2012) |> semi_join(plots, by = c("trans", "dist")) |>
    mutate(rate = round(midges / midge_days, 1)) |> as.data.frame() |> print()
cat("myv_arth rows:", nrow(myv_arth), "| NA midges:", sum(is.na(myv_arth$midges)), "\n")
