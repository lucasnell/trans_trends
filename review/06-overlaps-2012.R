
# List 2012 pitfall samples whose deployment periods overlap within a plot (E4).
# One line per group of mutually overlapping samples. Run from project root:
#   Rscript review/06-overlaps-2012.R > review/06-overlaps-2012-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))
key <- c("lakeid", "dist", "trans", "setdate", "coldate")

# Pitfall rows as built in scripts/00-clean-data-round1.R (current version)
p1 <- rd("data/raw/pitfalls-08-17.csv") |>
    select(lakeid, dist, trans, setdate, coldate, daysout, lyco_total, cara)
p2 <- rd("data/raw/pitfalls-11-19.csv") |>
    select(lakeid, dist, trans, setdate, coldate, daysout,
           lyco_total = Lyco_NEW_total, cara)
pit <- bind_rows(p1 |> mutate(file = "08-17"),
                 anti_join(p2, p1, by = key) |> mutate(file = "11-19")) |>
    filter(lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate))

# Group samples into overlap clusters within each plot: a new cluster starts
# when a sample's setdate is on/after the latest coldate seen so far.
ov <- pit |>
    filter(year(coldate) == 2012) |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           new_grp = is.na(prev_max) | as.numeric(setdate) >= prev_max,
           grp = cumsum(new_grp)) |>
    group_by(trans, dist, grp) |>
    filter(n() > 1) |>
    summarize(samples = paste0(format(setdate, "%m-%d"), "–", format(coldate, "%m-%d"),
                               " (", file, ")", collapse = "; "),
              n = n(), .groups = "drop")

cat("2012 overlap groups:", nrow(ov), "| samples involved:", sum(ov$n), "\n\n")
for (i in seq_len(nrow(ov))) {
    cat(sprintf("%-3s %3d m: %s\n", ov$trans[i], ov$dist[i], ov$samples[i]))
}

# ---- Counts of analysed taxa in each sample of each overlap group ----
taxa <- c("cara", "stap", "opil", "gnap", "sheet", "lyco")
p1c <- rd("data/raw/pitfalls-08-17.csv") |>
    mutate(sheet = aran_other + liny) |>
    select(all_of(key), lyco = lyco_total, sheet, gnap, opil, stap, cara)
grp_samples <- pit |>
    filter(year(coldate) == 2012) |>
    arrange(trans, dist, setdate, coldate) |>
    group_by(trans, dist) |>
    mutate(prev_max = lag(cummax(as.numeric(coldate))),
           grp = cumsum(is.na(prev_max) | as.numeric(setdate) >= prev_max)) |>
    group_by(trans, dist, grp) |>
    filter(n() > 1) |>
    ungroup() |>
    select(all_of(key), grp) |>
    left_join(p1c, by = key) |>
    mutate(days = as.numeric(coldate - setdate))

cat("\n\nCounts per sample (days = coldate - setdate). Taxa: ",
    paste(taxa, collapse = ", "), "\n", sep = "")
grp_samples |>
    group_by(trans, dist, grp) |>
    group_walk(\(x, g) {
        cat(sprintf("\n%s %d m\n", g$trans, g$dist))
        for (j in seq_len(nrow(x))) {
            cat(sprintf("  %s–%s  %2d d | %s\n",
                        format(x$setdate[j], "%m-%d"), format(x$coldate[j], "%m-%d"),
                        x$days[j],
                        paste(sprintf("%s %4s", taxa,
                                      ifelse(is.na(unlist(x[j, taxa])), "NA",
                                             as.character(unlist(x[j, taxa])))),
                              collapse = "  ")))
        }
    })

# Summary: how often are paired samples identical across all six taxa?
pairs <- grp_samples |>
    group_by(trans, dist, grp) |>
    summarize(n = n(),
              identical_all = n_distinct(paste(cara, stap, opil, gnap, sheet, lyco)) == 1,
              .groups = "drop")
cat("\nGroups whose samples have identical counts for all six taxa:",
    sum(pairs$identical_all), "of", nrow(pairs), "\n")
