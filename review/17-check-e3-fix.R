
# Checks the E3 fix in scripts/00-clean-data-round1.R: samples with NA counts
# should be dropped before summing, so their days don't count as effort.
# Run from project root:
#   Rscript review/17-check-e3-fix.R > review/17-check-e3-fix-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
src <- readLines("scripts/00-clean-data-round1.R")
src <- src[!grepl("^\\s*write_csv", src)]
invisible(capture.output(suppressMessages(eval(parse(text = src)))))

cat("== Pitfalls ==\n")
cat("pit_clean rows:", nrow(pit_clean), "| NA counts:", sum(is.na(pit_clean$count)), "\n")
print(pit_clean |> filter(trans == "hag", dist == 5, year == 2018, taxon == "sheet") |>
          as.data.frame())

cat("\n== Infalls ==\n")
cat("inf_clean rows:", nrow(inf_clean), "| sum of midges:", sum(inf_clean$midges), "\n")
cat("myv_arth rows:", nrow(myv_arth),
    "| rows with NA midges:", sum(is.na(myv_arth$midges)),
    "| rows with midges == 0:", sum(myv_arth$midges == 0, na.rm = TRUE), "\n")

cat("\n== Spot checks against review/02 (corrected rates) and review/11 ==\n")
myv_arth |> distinct(trans, dist, year, midges, midge_days) |>
    filter((trans == "hag" & dist == 150 & year == 2014) |
               (trans == "btl" & dist == 5 & year == 2013)) |>
    mutate(rate = midges / midge_days) |> as.data.frame() |> print()
cat("Expected: hag 150 m 2014 about 73.6/day (02); btl 5 m 2013 9117 over 70 days (11).\n")
