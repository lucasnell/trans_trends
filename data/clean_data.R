# load packges
library(tidyverse)
library(lubridate)

# ── Attaching packages ───────────────────────────────────────────────────────── tidyverse 1.2.1 ──
# ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
# ✔ tibble  1.4.2     ✔ dplyr   0.7.6
# ✔ tidyr   0.8.1     ✔ stringr 1.3.1
# ✔ readr   1.1.1     ✔ forcats 0.3.0
# ── Conflicts ──────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# ✖ dplyr::filter() masks stats::filter()
# ✖ dplyr::lag()    masks stats::lag()

# load data
pit = read_csv("data/myvatn_pitfalls.csv")
inf = read_csv("data/myvatn_infalls.csv")

# clean piftall data
pit_clean = pit %>%
    # only include data from lakeid = "myv" (the other lake has gaps) and dist = 200 (which is an error)
    # remove trans "kal2" which is an error
    filter(lakeid == "myv", dist != 200, trans != "kal2") %>%
    # create variable for year
    # create variables for aggregated abundances of certain groups
    mutate(year = year(coldate),
           sheet = aran_other + liny,
           days = as.numeric(coldate - setdate)) %>%
    # rename 'acar_total" as "acar" and "lyco_total" and "lyco
    rename(lyco = lyco_total) %>%
    # select columns to keep
    select(trans, dist,setdate, coldate, year, days, lyco, sheet, gnap, opil, stap, cara) %>%
    na.omit() %>%
    pivot_longer(cols = c(lyco, sheet, gnap, opil, stap, cara)) %>%
    rename(taxon = name,
           count = value) %>%
    group_by(trans, dist, year, taxon) %>%
    summarize(days = sum(days),
              count = round(sum(count))) %>%
    ungroup() %>%
    arrange(trans, dist, taxon, year)


# clean infall data
inf_clean = inf %>%
    # only include data from lakeid = "myv" (the other lake has gaps)
    filter(lakeid=="myv") %>%
    # create variable for year
    # correct large (bgch) and small (smch) midges for subsampling; calculate total midge abundance
    mutate(year = year(coldate),
           bgch = bgch/fract_bgch,
           smch = smch/fract_smch,
           midges = bgch + smch,
           m_days  = as.numeric(coldate - setdate)) %>%
    select(trans, dist, year, m_days, midges) %>%
    na.omit() %>%
    # calculate total year abundance for each taxon and site
    group_by(trans, dist, year) %>%
    summarize(m_days = sum(m_days),
              midges = sum(midges))

# merge data frames and save
myv_arth = left_join(pit_clean, inf_clean)
# write_csv(myv_arth, "data/myv_arth.csv")
