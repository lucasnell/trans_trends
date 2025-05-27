
#'
#' This file cleans the raw data and preps it for model fitting and plotting.
#'


# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")

z_trans <- function(x) (x - mean(x)) / sd(x)


# ------------*
# Load raw data ----
# ------------*

myv_arth <- read_csv("data/myv_arth.csv", col_types = cols())


# ------------*
# Prep data ----
# ------------*


data_fit <- myv_arth |>
    filter(taxon != "acar") |>
    rename(distance = dist) |>
    # group_by(taxon) |>
    mutate(y_noz = log1p(count / season_days),
           y = z_trans(y_noz)) |>
    ungroup() |>
    mutate(across(contains("midges"), \(x) z_trans(log(x / midge_days)),
                  .names = "{.col}_z"),
           time = factor(year, levels = c(1:max(year))) |> as.numeric(),
           time = time - min(time),
           time_z = (time)/(sd(time)),
           dist  = log(distance),
           dist_z = z_trans(dist),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon, levels = taxa_lvls),
           plot = factor(paste0(trans, dist)),
           taxon_plot = factor(paste0(taxon, plot)),
           taxon_trans = factor(paste0(taxon, trans)),
           taxon_pretty = factor(paste(taxon), levels = taxa_lvls, labels = taxa_labs)) |>
    arrange(trans, distance, taxon, year)




# ------------*
# Write data to RDS file ----
# ------------*


write_rds(data_fit, data_rds)
