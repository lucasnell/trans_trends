
#'
#' This file fits the model using the `TransTrendsPkg` package,
#' and saves these objects to `*.rds` files.
#' Below, change the `.REFIT_MODELS` object to `TRUE` if you want to refit
#' these models even if the `rds` files exist.
#'


.REFIT_MODELS <- FALSE

# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")

# data_fit_lag <- read_rds(data_rds$lag)
data_fit_nolag <- read_rds(data_rds$nolag)


# # =============================================================================*
# # Fit model with lag ----
# # =============================================================================*
#
# if (! file.exists(model_rds$lag) || .REFIT_MODELS) {
#
#     # Takes ~6 min on my machine with >= 4 threads:
#     fit_lag <- armm(formula = count ~ small_midges_z + small_midges_lag_z +
#                               big_midges_z + big_midges_lag_z +
#                               time_z + dist_z +
#                     (1 | taxon + taxon_plot + taxon_trans) +
#                     (small_midges_z + small_midges_lag_z +
#                      big_midges_z + big_midges_lag_z +
#                      time_z + dist_z | taxon) +
#                         offset(log(season_days)),
#                 time_form = ~ time | trans + distf + taxon,
#                 ar_form = ~ taxon,
#                 obs_error = TRUE,
#                 distr = "lnorm_poisson",
#                 data = data_fit_lag,
#                 x_scale = FALSE,
#                 y_scale = NULL,
#                 hmc = TRUE,
#                 change = TRUE,
#                 rstan_control = list(iter = 3000, chains = 4, seed = 1895472999,
#                                      control = list(adapt_delta = 0.99)))
#     write_rds(fit_lag, model_rds$lag)
#
# } else {
#     fit_lag <- read_rds(model_rds$lag)
# }
#
#
# # # examine fit
# summary(fit_lag)
# coef(fit_lag)







# =============================================================================*
# Fit model without lag ----
# =============================================================================*

# This will be manipulated for reduced models below:
nolag_formula <- count ~ small_midges_z + big_midges_z + time_z + dist_z +
    (1 | taxon + taxon_plot + taxon_trans) +
    (small_midges_z + big_midges_z + time_z + dist_z | taxon) +
    offset(log(season_days))

if (! file.exists(model_rds$nolag) || .REFIT_MODELS) {

    # Takes ~6 min on my machine with >= 4 threads:
    fit_nolag <- armm(formula = nolag_formula,
                      time_form = ~ time | trans + distf + taxon,
                      ar_form = ~ taxon,
                      obs_error = TRUE,
                      distr = "lnorm_poisson",
                      data = data_fit_nolag,
                      x_scale = FALSE,
                      y_scale = NULL,
                      hmc = TRUE,
                      change = TRUE,
                      rstan_control = list(iter = 3000, chains = 4, seed = 130748637,
                                           control = list(adapt_delta = 0.99)))
    write_rds(fit_nolag, model_rds$nolag)

} else {
    fit_nolag <- read_rds(model_rds$nolag)
}


# # examine fit
summary(fit_nolag)
coef(fit_nolag)






# # =============================================================================*
# # Fit reduced models without lag ----
# # =============================================================================*
#
# # Extract random effects variables by taxon, then reconstruct formulas with
# # just that random effect removed.
# reduced_forms <- TransTrendsPkg:::findbars(nolag_formula) |>
#     map_chr(deparse) |>
#     keep(\(x) str_ends(x,"\\| taxon")) |>
#     str_remove(" \\| taxon") |>
#     str_split(" \\+ ") |>
#     getElement(1) |>
#     set_names() |>
#     map(\(r) {
#         ff <- deparse(TransTrendsPkg:::nobars(nolag_formula))
#         rf1 <- TransTrendsPkg:::findbars(nolag_formula) |>
#             map_chr(deparse) |>
#             discard(\(x) str_ends(x,"\\| taxon")) |>
#             map_chr(\(x) sprintf(fmt = "(%s)", x))
#         rf0 <- TransTrendsPkg:::findbars(nolag_formula) |>
#             map_chr(deparse) |>
#             keep(\(x) str_ends(x,"\\| taxon")) |>
#             str_remove(" \\| taxon") |>
#             str_split(" \\+ ") |>
#             getElement(1) |>
#             discard(\(x) x == r) |>
#             paste(collapse = " + ") |>
#             sprintf(fmt = "(%s | taxon)")
#         as.formula(paste(ff, rf1, rf0, sep = "+"), env = globalenv())
#     })
#
# # Random seeds for each reduced model:
# reduced_seeds <- list(small_midges_z = 406770207, big_midges_z = 9178281,
#                       time_z = 354495043, dist_z = 700478870)
# stopifnot(identical(sort(names(reduced_seeds)), sort(names(reduced_forms))))
#
#
# if (! file.exists(model_rds$nolag_reds) || .REFIT_MODELS) {
#
#     # Takes ~22 min on my machine with >= 4 threads:
#     fit_nolag_reds <- map2(reduced_forms, reduced_seeds,
#                            \(red_form, seed) {
#                                call_ <- fit_nolag$call
#                                call_$formula <- red_form
#                                call_$rstan_control$seed <- seed
#                                # turns off printing for individual model fits:
#                                call_$rstan_control$refresh <- 0L
#                                # Do the model fitting:
#                                fit <- eval(call_)
#                                return(fit)
#                            }, .progress = TRUE)
#     write_rds(fit_nolag_reds, model_rds$nolag_reds)
#
#     # Because the file with all the fits is quite large, I'm also going to
#     # save just the loo and waic objects
#     red_mod_comps <- list(loo = map(c(list(full = fit_nolag), fit_nolag_reds), loo),
#                           waic = map(c(list(full = fit_nolag), fit_nolag_reds), waic))
#     write_rds(red_mod_comps, model_rds$nolag_red_comps)
#
# } else {
#     fit_nolag_reds <- read_rds(model_rds$nolag_reds)
#     red_mod_comps <- read_rds(model_rds$nolag_red_comps)
# }
#
# # Each ∆LOOIC or ∆WAIC compared to the full model:
# map_dbl(red_mod_comps$loo, \(x) x$estimates["looic","Estimate"]) |>
#     (\(x) x[-1] - x[1])()
# map_dbl(red_mod_comps$waic, \(x) x$estimates["waic","Estimate"]) |>
#     (\(x) x[-1] - x[1])()
#
