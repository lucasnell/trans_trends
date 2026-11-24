
# Verifies the E6 fix in scripts/00-clean-data-round1.R (filter on inf2):
# exactly one row removed, and btl 5 m 2013 midge totals corrected. Run from
# project root:
#   Rscript review/11-check-e6-fix.R > review/11-check-e6-fix-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
src <- readLines("scripts/00-clean-data-round1.R")
src <- src[!grepl("^\\s*write_csv", src)]
invisible(capture.output(suppressMessages(eval(parse(text = src)))))

raw2 <- inf2 |> unique()
removed <- anti_join(raw2 |> select(all_of(names(i2))), i2, by = names(i2))
cat("inf2 rows after unique():", nrow(raw2), "| i2 rows:", nrow(i2),
    "| removed:", nrow(raw2) - nrow(i2), "\n")
print(as.data.frame(removed))
cat("\nbtl 5 m 2013 in rebuilt myv_arth (expected 9117 midges over 70 days):\n")
print(myv_arth |> filter(trans == "btl", dist == 5, year == 2013) |>
          distinct(midges, midge_days) |> as.data.frame())
cat("myv_arth rows:", nrow(myv_arth), "\n")
