
# This file does the first round of cleaning data to generate `data/myv_arth.csv`

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})


# Columns that identify one physical sample:
sample_key <- c("lakeid", "dist", "trans", "setdate", "coldate")

# Stop if any sample (by `key`) appears more than once in `d`, except for
# duplicates listed in `known` (a data frame of key columns).
check_unique_key <- function(d, known = NULL, key = sample_key,
                             label = deparse(substitute(d))) {
    dups <- d |> count(across(all_of(key))) |> filter(n > 1)
    if (!is.null(known)) dups <- anti_join(dups, known, by = key)
    if (nrow(dups) > 0) {
        print(semi_join(d, dups, by = key) |> arrange(across(all_of(key))),
              width = Inf)
        stop(label, ": ", nrow(dups), " sample(s) appear more than once ",
             "(shown above).", call. = FALSE)
    }
    cat(label, ": no unexpected duplicated samples\n", sep = "")
    invisible(TRUE)
}

# For samples present in both `a` and `b`, return every value that differs
# between them in the shared non-key columns. NA vs. a value counts as a
# difference; NA vs. NA does not. Samples with NA in a key column are not
# matched.
compare_overlap <- function(a, b, key = sample_key) {
    cols <- setdiff(intersect(names(a), names(b)), key)
    ab <- inner_join(a, b, by = key, suffix = c(".a", ".b"),
                     na_matches = "never")
    cat("Samples in both:", nrow(ab),
        "| only in first:", nrow(anti_join(a, b, by = key)),
        "| only in second:", nrow(anti_join(b, a, by = key)), "\n")
    map(cols, \(v) {
        x <- ab[[paste0(v, ".a")]]
        y <- ab[[paste0(v, ".b")]]
        bad <- xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
        ab[bad, key] |>
            mutate(column = v, first = as.character(x[bad]),
                   second = as.character(y[bad]))
    }) |>
        list_rbind()
}



#====Pitfall data=====
# load data
pit1 <- read_csv("data/raw/pitfalls-08-17.csv")
pit2 <- read_csv("data/raw/pitfalls-11-19.csv")

#find columns that don't match
colnames(pit2)[which(!colnames(pit2) %in% colnames(pit1))]
colnames(pit1)[which(!colnames(pit1) %in% colnames(pit2))]

# pull columns used for future analyses
p1 <- pit1 |>
    select(lakeid, dist, trans, setdate, coldate, daysout,
           aran_other, liny,
           acar_total, lyco_total,
           gnap, opil, stap, cara)

p2 <- pit2 |>
    mutate(acar_total = trom_adult + trom_juv + acar,
           # linyphiidae and other sheet weavers were combined in 2018:
           liny = ifelse(is.na(liny) & !is.na(aran), 0, liny)) |>
    select(lakeid, dist, trans, setdate, coldate, daysout,
           aran_other = aran, liny,
           acar_total, lyco_total = Lyco_NEW_total,
           gnap, opil, stap, cara)


# check that each sample appears once per file (among rows kept below)
pit_kept <- \(x) filter(x, lakeid == "myv", dist != 200, trans != "kal2",
                        !is.na(setdate))
check_unique_key(pit_kept(p1), label = "p1")
check_unique_key(pit_kept(p2), label = "p2")

# check that overlapping data match before combining
pit_mismatch <- compare_overlap(p1, p2)
count(pit_mismatch, column) |> print()
# `acar_total` is NA in pitfalls-11-19.csv for 2011-2017 because `trom_adult`
# and `trom_juv` are NA there. Mites aren't analysed, and p1's values are kept
# below, so that disagreement is expected. Stop on any other mismatch:
if (any(pit_mismatch$column != "acar_total")) {
    print(filter(pit_mismatch, column != "acar_total"), n = Inf)
    stop("Overlapping pitfall samples disagree (shown above).", call. = FALSE)
}

# join new and old data to get full extent
pit <- bind_rows(p1, anti_join(p2, p1, by = c("lakeid", "dist", "trans",
                                              "setdate", "coldate")))



# clean piftall data
pit_clean <- pit |>
    # only include data from lakeid = "myv" (the other lake has gaps) and dist = 200 (which is an error)
    # remove trans "kal2" which is an error
    filter(lakeid == "myv", dist != 200, trans != "kal2") |>
    # create variable for year
    # create variables for aggregated abundances of certain groups
    mutate(year = year(coldate),
           sheet = aran_other + liny) |>
    # rename 'acar_total" as "acar" and "lyco_total" and "lyco
    rename(acar = acar_total,
           lyco = lyco_total) |>
    # eliminate sample with missing set date
    filter(!is.na(setdate)) |>
    # define span of sampling season by year and site
    mutate(daysout = as.numeric(coldate - setdate)) |>
    # select columns to keep
    select(trans, dist, year, setdate, coldate, daysout, lyco, sheet, gnap,
           opil, acar, stap, cara) |>
    # Gather days out for mean and sd calc below:
    (\(x) {
        daysout <<- x$daysout[!is.na(x$daysout)]
        return(x)
    })() |>
    # convert to long form
    # keep specific taxa that make most sense for analysis
    pivot_longer(lyco:cara, names_to = "taxon", values_to = "count") |>
    # Make pivot_longer produce output identical to previous call to gather:
    arrange(match(taxon, unique(taxon))) |>
    # calculate total year abundance for each taxon and site
    # a few "counts" are averages, so round
    group_by(trans, dist, year, taxon) |>
    # remove rows with no counts:
    filter(!is.na(count)) |>
    summarize(count = sum(round(count), na.rm = TRUE),
              first_day = min(setdate),
              last_day = max(coldate),
              season_days = sum(daysout),
              .groups = "drop") |>
    # sort by year within site and taxon
    arrange(trans, dist, taxon, year)


# mean and sd for how often pitfall traps were sampled:
daysout |> mean()
# [1] 14.54329
daysout |> sd()
# [1] 3.930394



#=====Infall data====
# load data
inf1 <- read_csv("data/raw/infalls-08-17.csv")
inf2 <- read_csv("data/raw/infalls-13-19.csv")

#find columns that don't match
colnames(inf2)[which(!colnames(inf2) %in% colnames(inf1))]
colnames(inf1)[which(!colnames(inf1) %in% colnames(inf2))]

# pull columns used for future analyses
i1 <- inf1 |>
    unique() |>  # one row is duplicated
    select(lakeid, dist, trans, setdate, coldate, daysout,
           contains("bgch"), contains("smch"))

i2 <- inf2 |>
    unique() |> # one row is duplicated
    # This row is an duplicated version of another with a slight (erroneous)
    # change:
    filter(!(trans == "btl" & dist == 5 & setdate == as.Date("2013-06-05") &
             coldate == as.Date("2013-08-14"))) |>
    select(lakeid, dist, trans, setdate, coldate, daysout,
           contains("bgch"), contains("smch"))

# check that each sample appears once per file (among rows kept below)
inf_kept <- \(x) filter(x, lakeid == "myv", !is.na(setdate))
# Four samples in infalls-08-17.csv appear twice with different values; they
# are allowed here until they're resolved against the field sheets. Remove
# each from this list once it's fixed in the raw file:
inf_known_dups <- tribble(
    ~lakeid, ~dist, ~trans,                ~setdate,                ~coldate,
    "myv",       5,  "hag", as.Date("2012-06-05"), as.Date("2012-06-19"),
    "myv",      50,  "non", as.Date("2012-07-17"), as.Date("2012-07-31"),
    "myv",     150,  "non", as.Date("2008-06-04"), as.Date("2008-06-10"),
    "myv",       5,  "vin", as.Date("2008-05-30"), as.Date("2008-06-09"))
check_unique_key(inf_kept(i1), known = inf_known_dups, label = "i1")
check_unique_key(inf_kept(i2), label = "i2")

# check that overlapping data match before combining
inf_mismatch <- compare_overlap(i1, i2)
if (nrow(inf_mismatch) > 0) {
    print(inf_mismatch, n = Inf)
    stop("Overlapping infall samples disagree (shown above).", call. = FALSE)
}

# join new and old data to get full extent
inf <- full_join(i1, i2,
                 by = join_by(lakeid, dist, trans, setdate, coldate, daysout,
                              fract_bgch, bgch, fract_smch, smch))


# clean infall data
inf_clean <- inf |>
    # only include data from lakeid = "myv" (the other lake has gaps)
    filter(lakeid=="myv") |>
    # create variable for year
    # calculate trapping days
    # correct large (bgch) and small (smch) midges for subsampling; calculate total midge abundance
    mutate(year = year(coldate),
           daysout = as.numeric(coldate - setdate),
           bgch = bgch/fract_bgch,
           smch = smch/fract_smch,
           midges = bgch + smch) |>
    # eliminate sample with missing set date
    filter(!is.na(setdate)) |>
    # remove rows with no midge counts
    filter(is.na(midges)) |>
    # calculate total year abundance for each taxon and site
    group_by(trans, dist, year) |>
    summarize(midges = sum(midges, na.rm = TRUE),
              small_midges = sum(smch, na.rm = TRUE),
              big_midges = sum(bgch, na.rm = TRUE),
              midge_first_day = min(setdate),
              midge_last_day = max(coldate),
              midge_days = sum(daysout),
              .groups = "drop")

# merge data frames and save
myv_arth <- left_join(pit_clean, inf_clean, by = join_by(trans, dist, year)) |>
    # as.data.frame() |>
    select(trans, dist, year, contains("midges"), taxon, everything())
# write_csv(myv_arth, "data/myv_arth.csv")

