
#'
#' This file cleans the raw data and preps it for model fitting and plotting.
#' It creates two versions of the dataset: one that includes lagged midge
#' abundances as covariates and another that doesn't include these effects.
#' Lagged midge abundances do not have a statistically significant effect,
#' so we don't include them in the final model.
#' Using the non-lagged version also allows us to use the entire dataset,
#' since lagging causes the first element in a time series to be `NA`.
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
# Prep data with lag ----
# ------------*

data_fit_lag <- myv_arth |>
    filter(taxon != "acar") |>
    rename(distance = dist) |>
    # group_by(taxon) |>
    mutate(y_noz = log1p(count / season_days),
           y = z_trans(y_noz)) |>
    ungroup() |>
    group_by(trans, distance, taxon) |>
    arrange(trans, distance, taxon, year) |>
    mutate(midges_lag = lag(midges),
           small_midges_lag = lag(small_midges),
           big_midges_lag = lag(big_midges)) |>
    ungroup() |>
    na.omit() |>
    mutate(across(contains("midges"), \(x) z_trans(x / midge_days),
                  .names = "{.col}_z"),
           time = factor(year, levels = c(1:max(year))) |>
               as.numeric() |>
               (\(x) (x - min(x)))(),
           time_z = (time)/(sd(time)),
           dist = log(distance),
           dist_z = z_trans(dist),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon),
           plot = factor(paste0(trans, dist)),
           taxon_plot = factor(paste0(taxon, plot)),
           taxon_trans = factor(paste0(taxon, trans)),
           taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs)) |>
    arrange(trans, distance, taxon, year)



# ------------*
# Prep data without lag ----
# (this retains first year of data)
# ------------*


data_fit_nolag <- myv_arth |>
    filter(taxon != "acar") |>
    rename(distance = dist) |>
    # group_by(taxon) |>
    mutate(y_noz = log1p(count / season_days),
           y = z_trans(y_noz)) |>
    ungroup() |>
    mutate(across(contains("midges"), \(x) z_trans(x / midge_days),
                  .names = "{.col}_z"),
           time = factor(year, levels = c(1:max(year))) |> as.numeric(),
           time = time - min(time),
           time_z = (time)/(sd(time)),
           dist  = log(distance),
           dist_z = z_trans(dist),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon),
           plot = factor(paste0(trans, dist)),
           taxon_plot = factor(paste0(taxon, plot)),
           taxon_trans = factor(paste0(taxon, trans)),
           taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs)) |>
    arrange(trans, distance, taxon, year)


# # graph to check
# # use full data
# data_fit_nolag |>
#     ggplot(aes(x = year,
#                y = y))+
#     facet_wrap(~taxon)+
#     geom_line(aes(group = plot),
#               alpha =0.3,
#               size = 0.3)+
#     geom_line(data = data_fit_nolag |>
#                   group_by(year, taxon) |>
#                   summarize(y = mean(small_midges_z), .groups = "drop"),
#               size = 0.7, color = "firebrick")+
#     geom_line(data = data_fit_nolag |>
#                   group_by(year, taxon) |>
#                   summarize(y = mean(big_midges_z), .groups = "drop"),
#               size = 0.7, color = "dodgerblue")+
#     geom_line(data = data_fit_nolag |>
#                   group_by(year, taxon) |>
#                   summarize(y = mean(y), .groups = "drop"),
#               size = 0.7, color = "black")+
#     theme_classic()



# ------------*
# Write data to RDS files ----
# ------------*


write_rds(data_fit_nolag, data_rds$nolag)
write_rds(data_fit_lag, data_rds$lag)
