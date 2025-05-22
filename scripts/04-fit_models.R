
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

data_fit_nolag <- read_rds(data_rds$nolag)




# =============================================================================*
# Fit model without lag ----
# =============================================================================*

# This will be manipulated for reduced models below:
nolag_formula <- count ~ midges_z + time_z + dist_z +
    (1 | taxon + taxon_plot + taxon_trans) +
    (midges_z + time_z + dist_z | taxon) +
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

