
# Details behind review/E5-infall-duplicates.md: the conflicting duplicate
# infall samples in infalls-08-17.csv, their curation comments, and which
# plots are missing a 2012 sample. Run from project root:
#   Rscript review/07-e5-infall-duplicates.R > review/07-e5-infall-duplicates-output.txt 2>&1

suppressPackageStartupMessages(library(tidyverse))
options(width = 200)
d <- suppressMessages(read_csv("data/raw/infalls-08-17.csv", show_col_types = FALSE,
                               guess_max = 1e5))
d$line <- seq_len(nrow(d)) + 1   # line in the csv, header = line 1
key <- c("lakeid", "dist", "trans", "setdate", "coldate")

cat("Exact duplicate rows (removed by unique() in scripts/00):",
    sum(duplicated(select(d, -line))), "\n\n")

cat("== Mývatn samples with the same key but different values ==\n")
dd <- d |> distinct(across(-line), .keep_all = TRUE) |>
    filter(lakeid == "myv", !is.na(setdate)) |>
    group_by(across(all_of(key))) |> filter(n() > 1) |> ungroup() |>
    arrange(trans, dist, setdate, line)
dd |> select(line, trans, dist, setdate, coldate, daysout, fract_bgch, bgch,
             fract_smch, smch, simu, comments) |> as.data.frame() |> print()
cat("\nCuration notes (comments.2):\n")
for (i in seq_len(nrow(dd))) cat(sprintf("line %d: %s\n", dd$line[i], dd$comments.2[i]))

cat("\n== Number of 2012 samples per Mývatn plot (analysed sites) ==\n")
d |> filter(lakeid == "myv", format(coldate, "%Y") == "2012",
            trans %in% c("btl", "hag", "kal", "non", "vin")) |>
    count(trans, dist) |> as.data.frame() |> print()

cat("\n== 2012 samples at btl 50 m and kal 5 m ==\n")
d |> filter(lakeid == "myv", format(coldate, "%Y") == "2012",
            (trans == "btl" & dist == 50) | (trans == "kal" & dist == 5)) |>
    arrange(trans, dist, setdate) |>
    select(line, trans, dist, setdate, coldate, daysout, bgch, smch) |>
    as.data.frame() |> print()
