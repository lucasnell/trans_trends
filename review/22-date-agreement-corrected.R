# Agreement between pitfall and infall service dates after the corrections in
# scripts/00-clean-data-round1.R (E4 July and 44-day, E5, and E7). Same measure
# as script 20: for each pitfall set or collection date, days to the nearest
# infall set or collection date at the same plot and year.
# Run from project root:
#   Rscript review/22-date-agreement-corrected.R > review/22-date-agreement-corrected-output.txt 2>&1

suppressPackageStartupMessages(library(tidyverse))
invisible(capture.output(suppressMessages(suppressWarnings(
    source("scripts/00-clean-data-round1.R")))))

pit_d <- pit |> filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate))
inf_d <- bind_rows(i1, i2) |> filter(lakeid == "myv", !is.na(setdate))

svc <- \(d) d |>
    mutate(year = year(coldate)) |>
    pivot_longer(c(setdate, coldate), values_to = "date") |>
    distinct(trans, dist, year, date)
ps <- svc(pit_d); is <- svc(inf_d)
both <- inner_join(distinct(ps, trans, dist, year),
                   distinct(is, trans, dist, year), by = c("trans", "dist", "year"))
agree <- ps |>
    semi_join(both, by = c("trans", "dist", "year")) |>
    left_join(is |> rename(idate = date), by = c("trans", "dist", "year"),
              relationship = "many-to-many") |>
    group_by(trans, dist, year, date) |>
    summarize(gap = min(abs(as.numeric(date - idate))), .groups = "drop")

cat("Overall: ", sum(agree$gap == 0), " of ", nrow(agree), " pitfall service dates (",
    round(100 * mean(agree$gap == 0), 1), "%) fall on an infall service date\n\n", sep = "")
agree |>
    mutate(yr = if_else(year == 2012, "2012", "other years"),
           gap = cut(gap, c(-Inf, 0, 1, 3, Inf), labels = c("0", "1", "2-3", ">3"))) |>
    count(yr, gap) |>
    group_by(yr) |>
    mutate(pct = round(100 * n / sum(n), 1)) |>
    print()
cat("\nBy year:\n")
agree |>
    group_by(year) |>
    summarize(n = n(), exact_pct = round(100 * mean(gap == 0), 1)) |>
    print(n = Inf)
cat("\n2012 pitfall service dates with no infall match:\n")
agree |> filter(year == 2012, gap > 0) |> count(date, gap) |> print(n = Inf)
cat("\nOther-year pitfall service dates with no infall match, by year:\n")
agree |> filter(year != 2012, gap > 0) |>
    group_by(year) |>
    summarize(n = n(), dates = paste(sort(unique(format(date, "%m/%d"))), collapse = " ")) |>
    print(width = Inf)
cat("\nGap sizes for other-year mismatches:\n")
agree |> filter(year != 2012, gap > 0) |> count(year, gap) |> print()
