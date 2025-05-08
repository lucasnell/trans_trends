
#'
#' This file fits the model using the `TransTrendsPkg` package,
#' and saves these objects to `*.rds` files.
#' Below, change the `.REFIT_MODEL` object to `TRUE` if you want to refit
#' these models even if the `rds` files exist.
#'


.REFIT_MODEL <- FALSE

# =============================================================================*
# Preliminaries ----
# =============================================================================*

library(tidyverse)
library(TransTrendsPkg)

# Set up parallel processing (for brms)
options(mc.cores = max(1, parallel::detectCores()-2))


# Names of RDS files created here:
rds_files <- list(lag = "model_fits/lag-model.rds",
                  nolag = "model_fits/no_lag-model.rds",
                  nolag_reds = "model_fits/no_lag-reduced-models.rds")


z_trans <- function(x) (x - mean(x)) / sd(x)

# ------------*
# Standardize taxa order ----
# ------------*

taxa_order <- c(5:6, 4, 1, 3, 2)
taxa_lvls <- c("gnap","lyco","sheet","opil","cara", "stap")[taxa_order]
taxa_labs <- c("Ground spiders","Wolf spiders","Sheet weavers",
               "Harvestmen","Ground beetles", "Rove beetles")[taxa_order]


# ------------*
# Load data ----
# ------------*

myv_arth <- read_csv("data/myv_arth.csv", col_types = cols())






# ------------*
# Prep data with lag ----
# ------------*




# examine midge catch
data_fit_lag <- myv_arth |>
    filter(taxon != "acar") |>
    rename(distance = dist) |>
    # group_by(taxon) |>
    mutate(y = log1p(count / season_days),
           y = z_trans(y)) |>
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


# examine midge catch
data_fit_nolag <- myv_arth |>
    filter(taxon != "acar") |>
    rename(distance = dist) |>
    # group_by(taxon) |>
    mutate(y = log1p(count / season_days),
           y = (y - mean(y))/(sd(y))) |>
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




# =============================================================================*
# Fit model with lag ----
# =============================================================================*

if (! file.exists(rds_files$lag) || .REFIT_MODEL) {

    # Takes ~11 min on my machine with >= 4 threads:
    fit_lag <- armm(formula = y ~ small_midges_z + small_midges_lag_z +
                              big_midges_z + big_midges_lag_z +
                              time_z + dist_z +
                    (1 | taxon + taxon_plot + taxon_trans) +
                    (small_midges_z + small_midges_lag_z +
                     big_midges_z + big_midges_lag_z +
                     time_z + dist_z | taxon),
                time_form = ~ time | trans + distf + taxon,
                ar_form = ~ taxon,
                obs_error = TRUE,
                distr = "normal",
                data = data_fit_lag,
                x_scale = FALSE,
                y_scale = NULL,
                hmc = TRUE,
                change = TRUE,
                rstan_control = list(iter = 3000, chains = 4, seed = 1895472999,
                                     control = list(adapt_delta = 0.99)))
    write_rds(fit_lag, rds_files$lag)

} else {
    fit_lag <- read_rds(rds_files$lag)
}


# examine fit
summary(fit_lag)
coef(fit_lag)







# =============================================================================*
# Fit model without lag ----
# =============================================================================*

# This will be manipulated for reduced models below:
nolag_formula <- y ~ small_midges_z + big_midges_z + time_z + dist_z +
    (1 | taxon + taxon_plot + taxon_trans) +
    (small_midges_z + big_midges_z + time_z + dist_z | taxon)

if (! file.exists(rds_files$nolag) || .REFIT_MODEL) {

    # Takes ~8 min on my machine with >= 4 threads:
    fit_nolag <- armm(formula = nolag_formula,
                      time_form = ~ time | trans + distf + taxon,
                      ar_form = ~ taxon,
                      obs_error = TRUE,
                      distr = "normal",
                      data = data_fit_nolag,
                      x_scale = FALSE,
                      y_scale = NULL,
                      hmc = TRUE,
                      change = TRUE,
                      rstan_control = list(iter = 3000, chains = 4, seed = 130748637,
                                           control = list(adapt_delta = 0.99)))
    write_rds(fit_nolag, rds_files$nolag)

} else {
    fit_nolag <- read_rds(rds_files$nolag)
}


# examine fit
summary(fit_nolag)
coef(fit_nolag)



# using "slim" data missing first year (prepared for lag model even though lags aren't included)
# fit_nolag_slim <- armm(formula = y ~ small_midges_z + big_midges_z + time_z + dist_z +
#                       (1 | taxon + taxon_plot + taxon_trans) +
#                       (small_midges_z + big_midges_z + time_z + dist_z | taxon),
#                   time_form = ~ time | trans + distf + taxon,
#                   ar_form = ~ taxon,
#                   obs_error = TRUE,
#                   distr = "normal",
#                   data = data_fit_nolag,
#                   x_scale = FALSE,
#                   y_scale = NULL,
#                   hmc = TRUE,
#                   change = TRUE,
#                   rstan_control = list(iter = 3000, chains = 4, seed = 254498293,
#                                        control = list(adapt_delta = 0.99)))

# export
# write_rds(fit_nolag_slim, "fit_nolag_slim.rds")

# examine fit
# summary(fit_nolag_slim)
# coef(fit_nolag_slim)



#=======#=======#=========================================================================================*




# =============================================================================*
# Fit reduced models without lag ----
# =============================================================================*

# Extract random effects variables by taxon, then reconstruct formulas with
# just that random effect removed.
reduced_forms <- TransTrendsPkg:::findbars(nolag_formula) |>
    map_chr(deparse) |>
    keep(\(x) str_ends(x,"\\| taxon")) |>
    str_remove(" \\| taxon") |>
    str_split(" \\+ ") |>
    getElement(1) |>
    set_names() |>
    map(\(r) {
        ff <- deparse(TransTrendsPkg:::nobars(nolag_formula))
        rf1 <- TransTrendsPkg:::findbars(nolag_formula) |>
            map_chr(deparse) |>
            discard(\(x) str_ends(x,"\\| taxon")) |>
            map_chr(\(x) sprintf(fmt = "(%s)", x))
        rf0 <- TransTrendsPkg:::findbars(nolag_formula) |>
            map_chr(deparse) |>
            keep(\(x) str_ends(x,"\\| taxon")) |>
            str_remove(" \\| taxon") |>
            str_split(" \\+ ") |>
            getElement(1) |>
            discard(\(x) x == r) |>
            paste(collapse = " + ") |>
            sprintf(fmt = "(%s | taxon)")
        as.formula(paste(ff, rf1, rf0, sep = "+"), env = globalenv())
    })

# Random seeds for each reduced model:
reduced_seeds <- list(small_midges_z = 406770207, big_midges_z = 9178281,
                      time_z = 354495043, dist_z = 700478870)
stopifnot(identical(sort(names(reduced_seeds)), sort(names(reduced_forms))))


if (! file.exists(rds_files$nolag_reds) || .REFIT_MODEL) {

    # Takes ~31 min on my machine with >= 4 threads:
    fit_nolag_reds <- map2(reduced_forms, reduced_seeds,
                           \(red_form, seed) {
                               call_ <- fit_nolag$call
                               call_$formula <- red_form
                               call_$rstan_control$seed <- seed
                               fit <- eval(call_)
                               return(fit)
                           })
    write_rds(fit_nolag_reds, rds_files$nolag_reds)

} else {
    fit_nolag_reds <- read_rds(rds_files$nolag_reds)
}


# Takes ~12 sec
loo_full <- loo(fit_nolag)
loo_reds <- map(fit_nolag_reds, loo)
map_dbl(loo_reds, \(x) x$estimates["looic","Estimate"]) -
    loo_full$estimates["looic","Estimate"]

# Takes ~ 1 sec
waic_full <- waic(fit_nolag)
waic_reds <- map(fit_nolag_reds, waic)
map_dbl(waic_reds, \(x) x$estimates["waic","Estimate"]) -
    waic_full$estimates["waic","Estimate"]
