
# E4: within each 2012 plot, is any pitfall sample equal to the sum of two or
# more other samples from that plot (or twice a sample equal to such a sum,
# i.e. an average)? Checked on every count column in the raw file, not only
# the analysed taxa. Run from project root:
#   Rscript review/08-e4-sum-check.R > review/08-e4-sum-check-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
options(width = 200)

d <- suppressMessages(read_csv("data/raw/pitfalls-08-17.csv", show_col_types = FALSE,
                               guess_max = 1e5)) |>
    filter(lakeid == "myv", year(coldate) == 2012, !is.na(setdate),
           dist != 200, trans != "kal2")

count_cols <- c("aran_other", "lyco_total", "egsa", "liny", "gnap", "thom", "opil",
                "acar_other", "trom", "acar_total", "cole_other", "cryp", "stap",
                "cara", "curc", "colj", "clmb")
analysed <- c("cara", "stap", "opil", "gnap", "aran_other", "liny", "lyco_total")
# One 2012 sample (hag 5 m, 05/19–06/05) has NA for every count; it can be
# neither a sum nor a part, so drop it. `curc` is NA in 74 of 75 samples, so
# drop that column. All remaining columns are then complete.
na_rows <- d |> filter(if_any(all_of(analysed), is.na))
cat("Dropped samples with NA counts:\n")
print(as.data.frame(select(na_rows, trans, dist, setdate, coldate)))
d <- d |> filter(!if_any(all_of(analysed), is.na))
cols_ok <- setdiff(count_cols, "curc")
# `egsa` is read as character (non-numeric entries in other years); its 2012
# values are all numeric:
d <- d |> mutate(egsa = as.numeric(egsa))
stopifnot(!anyNA(d[cols_ok]))
cat("Count columns compared (no NAs in 2012):", paste(cols_ok, collapse = ", "), "\n")
cat("Samples:", nrow(d), "in", nrow(distinct(d, trans, dist)), "plots\n\n")

lab <- \(x) paste0(format(x$setdate, "%m/%d"), "–", format(x$coldate, "%m/%d"))

res <- list()
for (pl in split(d, ~ trans + dist, drop = TRUE)) {
    pl <- arrange(pl, setdate, coldate)
    n <- nrow(pl)
    M <- as.matrix(pl[cols_ok])
    for (i in seq_len(n)) {
        others <- setdiff(seq_len(n), i)
        for (k in 2:length(others)) {
            for (s in combn(others, k, simplify = FALSE)) {
                tot <- colSums(M[s, , drop = FALSE])
                for (type in c("sum", "average")) {
                    target <- if (type == "sum") M[i, ] else 2 * M[i, ]
                    if (type == "average" && k != 2) next
                    all_eq <- all(abs(target - tot) < 1e-8)
                    an_eq <- all(abs(target[analysed] - tot[analysed]) < 1e-8)
                    if (all_eq || an_eq) {
                        res[[length(res) + 1]] <- tibble(
                            plot = paste(pl$trans[1], pl$dist[1], "m"),
                            sample = lab(pl[i, ]), type = type,
                            parts = paste(lab(pl[s, ]), collapse = " + "),
                            all_columns_match = all_eq,
                            analysed_match = an_eq,
                            n_nonzero_analysed = sum(target[analysed] != 0))
                    }
                }
            }
        }
    }
}
out <- bind_rows(res)
cat("Matches (sample = sum of other samples, or sample = average of two):",
    nrow(out), "\n")
if (nrow(out) > 0) print(as.data.frame(out))

# Near misses: closest subset-sum per sample, by total absolute difference
# across analysed columns, to show how far the best candidate is.
cat("\nClosest subset sum for each sample in an overlap plot (analysed columns):\n")
near <- list()
for (pl in split(d, ~ trans + dist, drop = TRUE)) {
    pl <- arrange(pl, setdate, coldate); n <- nrow(pl)
    M <- as.matrix(pl[analysed])
    for (i in seq_len(n)) {
        others <- setdiff(seq_len(n), i); best <- Inf; bs <- NULL
        for (k in 2:length(others)) for (s in combn(others, k, simplify = FALSE)) {
            dd <- sum(abs(M[i, ] - colSums(M[s, , drop = FALSE])))
            if (dd < best) { best <- dd; bs <- s }
        }
        near[[length(near) + 1]] <- tibble(plot = paste(pl$trans[1], pl$dist[1], "m"),
                                           sample = lab(pl[i, ]),
                                           total = sum(M[i, ]),
                                           best_parts = paste(lab(pl[bs, ]), collapse = " + "),
                                           abs_diff = best)
    }
}
near <- bind_rows(near)
cat("Samples whose closest subset sum differs by <= 2 animals in total:\n")
print(as.data.frame(filter(near, abs_diff <= 2)))
cat("\nDistribution of abs_diff / sample total (all 2012 samples):\n")
print(summary(with(near, abs_diff / pmax(total, 1))))
