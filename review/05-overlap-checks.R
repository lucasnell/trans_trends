
# Proposed replacement for the overlap checks in scripts/00-clean-data-round1.R
# (pitfall loop at line 41, infall loop at line 134), run on the current raw
# files. Also tests the new pitfall join (lines 56-57). Run from project root:
#   Rscript review/05-overlap-checks.R > review/05-overlap-checks-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})
options(width = 120, pillar.width = 120)
rd <- \(x) suppressMessages(read_csv(x, show_col_types = FALSE, guess_max = 1e5))

sample_key <- c("lakeid", "dist", "trans", "setdate", "coldate")

# Stop if any physical sample (key) appears more than once in `d`.
check_unique_key <- function(d, key = sample_key, label = deparse(substitute(d))) {
    dups <- d |> count(across(all_of(key))) |> filter(n > 1)
    if (nrow(dups) > 0) {
        print(semi_join(d, dups, by = key) |> arrange(across(all_of(key))), width = Inf)
        stop(label, ": ", nrow(dups), " sample key(s) appear more than once (shown above).",
             call. = FALSE)
    }
    invisible(TRUE)
}

# Compare every shared value column for samples present in both `a` and `b`.
# NA vs. value counts as a mismatch; NA vs. NA does not.
compare_overlap <- function(a, b, key = sample_key, cols = NULL) {
    if (is.null(cols)) cols <- setdiff(intersect(names(a), names(b)), key)
    ab <- inner_join(a, b, by = key, suffix = c(".a", ".b"))
    cat("Samples in both:", nrow(ab),
        "| only in first:", nrow(anti_join(a, b, by = key)),
        "| only in second:", nrow(anti_join(b, a, by = key)), "\n")
    map_dfr(cols, \(v) {
        x <- ab[[paste0(v, ".a")]]; y <- ab[[paste0(v, ".b")]]
        bad <- xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
        ab[bad, key] |> mutate(column = v, first = as.character(x[bad]),
                               second = as.character(y[bad]))
    })
}

# ---- Pitfalls (p1, p2 built exactly as in scripts/00) ----
p1 <- rd("data/raw/pitfalls-08-17.csv") |>
    select(lakeid, dist, trans, setdate, coldate, daysout, aran_other, liny,
           acar_total, lyco_total, gnap, opil, stap, cara)
p2 <- rd("data/raw/pitfalls-11-19.csv") |>
    mutate(acar_total = trom_adult + trom_juv + acar,
           liny = ifelse(is.na(liny) & !is.na(aran), 0, liny)) |>
    select(lakeid, dist, trans, setdate, coldate, daysout, aran_other = aran, liny,
           acar_total, lyco_total = Lyco_NEW_total, gnap, opil, stap, cara)

cat("== Pitfall key uniqueness (all lakes) ==\n")
for (nm in c("p1", "p2")) {
    r <- tryCatch(check_unique_key(get(nm), label = nm), error = \(e) conditionMessage(e))
    cat(nm, ":", if (isTRUE(r)) "OK" else r, "\n")
}
cat("\nSame check restricted to rows that survive cleaning (myv, setdate present):\n")
keep <- \(d) filter(d, lakeid == "myv", dist != 200, trans != "kal2", !is.na(setdate))
for (nm in c("p1", "p2")) {
    r <- tryCatch(check_unique_key(keep(get(nm)), label = nm), error = \(e) conditionMessage(e))
    cat(nm, ":", if (isTRUE(r)) "OK" else r, "\n")
}

cat("\n== Pitfall overlap mismatches by column ==\n")
pm <- compare_overlap(p1, p2)
print(count(pm, column))
cat("Mismatches outside acar_total:\n"); print(filter(pm, column != "acar_total"))

cat("\n== New join (scripts/00 lines 56-57) ==\n")
pit <- bind_rows(p1, anti_join(p2, p1, by = sample_key))
cat("Rows: p1", nrow(p1), "+ p2-only", nrow(anti_join(p2, p1, by = sample_key)),
    "=", nrow(pit), "\n")
r <- tryCatch(check_unique_key(keep(pit), label = "pit (cleaned rows)"),
              error = \(e) conditionMessage(e))
cat("Duplicate sample keys among cleaned rows:", if (isTRUE(r)) "none" else r, "\n")

# ---- Infalls (i1, i2 built exactly as in scripts/00) ----
i1 <- rd("data/raw/infalls-08-17.csv") |> unique() |>
    select(lakeid, dist, trans, setdate, coldate, daysout, contains("bgch"), contains("smch"))
i2 <- rd("data/raw/infalls-13-19.csv") |> unique() |>
    select(lakeid, dist, trans, setdate, coldate, daysout, contains("bgch"), contains("smch"))

cat("\n== Infall key uniqueness (myv, setdate present) ==\n")
keep_i <- \(d) filter(d, lakeid == "myv", !is.na(setdate))
for (nm in c("i1", "i2")) {
    r <- tryCatch(check_unique_key(keep_i(get(nm)), label = nm), error = \(e) conditionMessage(e))
    cat(nm, ":", if (isTRUE(r)) "OK" else r, "\n")
}

cat("\n== Infall overlap mismatches by column ==\n")
im <- compare_overlap(i1, i2)
print(count(im, column))
print(im, n = 50)
