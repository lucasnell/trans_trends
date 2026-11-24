# Check the E4 corrections in scripts/00-clean-data-round1.R: 2012 pitfall
# samples recorded 07/17-08/09 or 06/26-08/09 now start 07/31, and samples
# recorded 06/19-07/03 now start 06/22 or 06/23. Run from project root:
#   Rscript review/21-check-e4-july-fix.R > review/21-check-e4-july-fix-output.txt 2>&1

suppressPackageStartupMessages(library(tidyverse))
invisible(capture.output(source("scripts/00-clean-data-round1.R")))

cat("Corrected rows (pit, after correction):\n")
pit |> filter(e4_july | e4_44d | e4_june) |> select(trans, dist, setdate, coldate) |> print(n = Inf)

cat("\nRemaining 2012 pitfall overlap groups within a plot:\n")
pit |>
    filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate),
           year(coldate) == 2012) |>
    distinct(trans, dist, setdate, coldate) |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           grp = cumsum(is.na(prev_max) | as.numeric(setdate) >= prev_max)) |>
    group_by(trans, dist, grp) |>
    filter(n() > 1) |>
    summarize(samples = paste0(format(setdate, "%m/%d"), "-", format(coldate, "%m/%d"),
                               collapse = "; "), .groups = "drop") |>
    select(-grp) |>
    print(n = Inf, width = Inf)

cat("\n2012 pitfall season_days per plot (same for every taxon):\n")
pit_clean |>
    filter(year == 2012) |>
    distinct(trans, dist, season_days, first_day, last_day) |>
    print(n = Inf)
