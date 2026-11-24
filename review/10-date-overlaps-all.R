
# Within-plot overlapping deployment periods (setdate before an earlier
# sample's coldate), all years, for the combined pitfall and infall data as
# built in scripts/00-clean-data-round1.R. Run from project root:
#   Rscript review/10-date-overlaps-all.R > review/10-date-overlaps-all-output.txt 2>&1

suppressPackageStartupMessages({ library(tidyverse); library(lubridate) })
src <- readLines("scripts/00-clean-data-round1.R")
src <- src[!grepl("^\\s*write_csv", src)]
invisible(capture.output(suppressMessages(eval(parse(text = src)))))

overlaps <- function(d) d |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist, year = year(coldate)) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           overlap = !is.na(prev_max) & as.numeric(setdate) < prev_max) |>
    ungroup()

po <- overlaps(pit_kept(pit))
io <- overlaps(inf_kept(inf))
cat("Pitfall samples starting before an earlier sample at the plot ended:\n")
print(count(filter(po, overlap), year))
cat("\nInfall samples starting before an earlier sample at the plot ended:\n")
print(count(filter(io, overlap), year))
cat("\nInfall overlaps (sample and the latest earlier coldate it overlaps):\n")
io |> filter(overlap) |>
    mutate(prev_coldate = as_date(prev_max)) |>
    select(trans, dist, setdate, coldate, prev_coldate) |>
    as.data.frame() |> print()

cat("\n== btl 5 m 2013 infall rows (as combined in `inf`) ==\n")
bt <- inf |> filter(lakeid == "myv", trans == "btl", dist == 5, year(coldate) == 2013) |>
    arrange(setdate, coldate) |>
    mutate(days = as.numeric(coldate - setdate),
           midges = bgch / fract_bgch + smch / fract_smch)
print(as.data.frame(select(bt, setdate, coldate, daysout, days, fract_smch, bgch, smch, midges)))
bad <- bt$setdate == as.Date("2013-06-05") & bt$coldate == as.Date("2013-08-14")
cat(sprintf("\nWith the 06/05-08/14 row: midges %.0f over %.0f days (%.1f/day)\n",
            sum(bt$midges), sum(bt$days), sum(bt$midges) / sum(bt$days)))
cat(sprintf("Without it:              midges %.0f over %.0f days (%.1f/day)\n",
            sum(bt$midges[!bad]), sum(bt$days[!bad]),
            sum(bt$midges[!bad]) / sum(bt$days[!bad])))
cat("data/myv_arth.csv has:\n")
print(suppressMessages(read_csv("data/myv_arth.csv")) |>
          filter(trans == "btl", dist == 5, year == 2013) |>
          distinct(midges, midge_days))

cat("\n== Kálfaströnd infall samples around the 2008 and 2011 overlaps ==\n")
inf1 |> mutate(line = row_number() + 1) |>
    filter(lakeid == "myv", trans == "kal", dist %in% c(5, 50, 150),
           (setdate >= as.Date("2008-06-10") & coldate <= as.Date("2008-07-08")) |
               (setdate >= as.Date("2011-07-05") & coldate <= as.Date("2011-08-02"))) |>
    mutate(days_from_dates = as.numeric(coldate - setdate)) |>
    select(line, dist, setdate, coldate, daysout, days_from_dates, bgch, smch, comments) |>
    arrange(setdate, dist) |> as.data.frame() |> print()
cat("\n== Kálfaströnd pitfall visit dates, 2008 and 2011 (pitfalls-08-17.csv) ==\n")
pit1 |> filter(lakeid == "myv", trans == "kal", dist %in% c(5, 50, 150),
               year(coldate) %in% c(2008, 2011)) |>
    distinct(year = year(coldate), setdate, coldate) |> arrange(setdate) |>
    as.data.frame() |> print()
