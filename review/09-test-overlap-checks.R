
# Tests that the overlap checks added to scripts/00-clean-data-round1.R catch
# the problems they're meant to catch. Run from project root:
#   Rscript review/09-test-overlap-checks.R > review/09-test-overlap-checks-output.txt 2>&1

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})

# Pull the helper definitions out of the cleaning script (everything before
# the pitfall section) so the tests use the exact code in that script.
src <- readLines("scripts/00-clean-data-round1.R")
eval(parse(text = src[1:(grep("^#====Pitfall data", src) - 1)]))

expect_error <- function(expr, label) {
    ok <- tryCatch({ capture.output(expr); FALSE }, error = \(e) TRUE)
    cat(sprintf("%-60s %s\n", label, if (ok) "PASS (stopped)" else "FAIL (did not stop)"))
}
expect_ok <- function(expr, label) {
    ok <- tryCatch({ capture.output(expr); TRUE }, error = \(e) FALSE)
    cat(sprintf("%-60s %s\n", label, if (ok) "PASS (ran)" else "FAIL (stopped)"))
}

a <- tibble(lakeid = "myv", dist = c(5, 50, 150), trans = "btl",
            setdate = as.Date("2013-06-05") + c(0, 0, 0),
            coldate = as.Date("2013-06-18") + c(0, 0, 0),
            daysout = 13, cara = c(1, 2, 3), acar_total = c(5, 6, 7))
b <- a |> mutate(acar_total = NA_real_)

# Duplicate detection
expect_ok(check_unique_key(a), "unique keys pass")
expect_error(check_unique_key(bind_rows(a, a[1, ])), "duplicated sample stops")
expect_ok(check_unique_key(bind_rows(a, a[1, ]), known = a[1, sample_key]),
          "known duplicate is allowed")
expect_error(check_unique_key(bind_rows(a, a[1:2, ]), known = a[1, sample_key]),
             "unlisted duplicate still stops when another is known")

# Mismatch detection
m <- compare_overlap(a, b)
cat(sprintf("%-60s %s\n", "NA vs value is reported (old check missed this)",
            if (nrow(m) == 3 && all(m$column == "acar_total")) "PASS" else "FAIL"))
m <- compare_overlap(a, a |> mutate(cara = c(1, 2, 4)))
cat(sprintf("%-60s %s\n", "differing count is reported",
            if (nrow(m) == 1 && m$column == "cara") "PASS" else "FAIL"))
m <- compare_overlap(a, a |> mutate(daysout = 14))
cat(sprintf("%-60s %s\n", "differing daysout is compared, not used to match",
            if (nrow(m) == 3 && all(m$column == "daysout")) "PASS" else "FAIL"))
m <- compare_overlap(mutate(a, cara = NA_real_), mutate(a, cara = NA_real_))
cat(sprintf("%-60s %s\n", "NA vs NA is not a mismatch",
            if (nrow(m) == 0) "PASS" else "FAIL"))
m <- compare_overlap(mutate(a, setdate = NA), mutate(a, setdate = NA, cara = 9))
cat(sprintf("%-60s %s\n", "rows with NA setdate are not matched",
            if (nrow(m) == 0) "PASS" else "FAIL"))
