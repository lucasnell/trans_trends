
# This file does the first round of cleaning data to generate `data/myv_arth.csv`

suppressPackageStartupMessages({
    library(tidyverse)
    library(lubridate)
})



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


# check that overlapping data match before combining
p2.1 <- p2
colnames(p2.1)[7:14] <- paste0(colnames(p2.1)[7:14], "_2")

for (i in 7:14){

    val1 <- colnames(p1)[i]
    val2 <- colnames(p2.1)[i]

    out <- left_join(p1, p2.1, by = join_by(lakeid, dist, trans, setdate, coldate, daysout)) |>
        (\(x) filter(x, !is.na(x[[7]])))() |>
        select(coldate, all_of(c(val1, val2))) |>
        (\(x) filter(x, x[[2]] != x[[3]]))()


    cat(val1, val2, " -- ", "nrows =", nrow(out), "\n")
}; rm(i, val1, val2, out) # no columns means no discrepancies in counts.

# join new and old data to get full extent
pit <- full_join(p1, p2,
                 by = join_by(lakeid, dist, trans, setdate, coldate, daysout,
                              aran_other, liny, acar_total, lyco_total,
                              gnap, opil, stap, cara))


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
    select(lakeid, dist, trans, setdate, coldate, daysout,
           contains("bgch"), contains("smch"))

# check that overlapping data match before combining
i2.1 <- i2
colnames(i2.1)[7:10] <- paste0(colnames(i2.1)[7:10], "_2")

for(i in 7:10){

    val1 <- colnames(i1)[i]
    val2 <- colnames(i2.1)[i]

    out <- left_join(i1, i2.1, by = join_by(lakeid, dist, trans, setdate,
                                            coldate, daysout)) |>
        (\(x) filter(x, !is.na(x[[7]])))() |>
        select(coldate, all_of(c(val1, val2))) |>
        (\(x) filter(x, x[[2]] != x[[3]]))()

    cat(val1, val2, " -- ", "nrows =", nrow(out), "\n")
}; rm(i, val1, val2, out) # no columns means no discrepancies in counts.

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

