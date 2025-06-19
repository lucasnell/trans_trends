
#'
#' This file fits the model using the `TransTrendsPkg` package,
#' and saves the resulting `armmMod` object to a `*.rds` file.
#'
#'
#' Below, change the `.REFIT_MODELS` object to `TRUE` if you want to refit
#' these models even if the `rds` files exist.
#'


.REFIT_MODELS <- FALSE

# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")

data_df <- read_rds(data_rds)




# =============================================================================*
# Fit model ----
# =============================================================================*


if (! file.exists(model_rds) || .REFIT_MODELS) {

    # Takes ~4.2 min on my machine with >= 4 threads:
    model_fit <- armm(formula = count ~ midges_z + time_z + dist_z +
                          (1 | taxon + taxon_plot + taxon_trans) +
                          (midges_z + time_z + dist_z | taxon) +
                          offset(log(season_days)),
                      time_form = ~ time | trans + distf + taxon,
                      ar_form = ~ taxon,
                      obs_error = TRUE,
                      distr = "lnorm_poisson",
                      data = data_df,
                      x_scale = FALSE,
                      y_scale = NULL,
                      hmc = TRUE,
                      change = TRUE,
                      rstan_control = list(iter = 3000, chains = 4, seed = 130748637,
                                           control = list(adapt_delta = 0.99)))
    write_rds(model_fit, model_rds)

} else {
    model_fit <- read_rds(model_rds)
}


# # examine fit
summary(model_fit)
coef(model_fit)

