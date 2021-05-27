
#'
#' This file cleans the data, fits the model using the `TransTrendsPkg` package,
#' and summarizes the model output to some temporary files.
#' Below, change the `.REFIT_MODEL` object to `TRUE` if you want to refit
#' the full model.
#' If you already have the `fit.rds` file that contains the model fit,
#' then it's much faster to just use that.
#' That's why `.REFIT_MODEL` is currently set to `FALSE`.
#'
#' All writing to `CSV` or `RDS` files is commented out here, so that
#' unnecessary files aren't written to disk.
#' Uncomment those lines if you want to create new versions.
#'



# .REFIT_MODEL <- F




#==========
#========== Preliminaries
#==========

# load packages
library(TransTrendsPkg)
library(tidyverse)
options(mc.cores = max(1, parallel::detectCores()-4))

# load raw data
myv_arth <- read_csv("data/myv_arth.csv")






#==========
#========== Prepare data
#==========

data_fit <- myv_arth %>%
    filter(taxon != "acar") %>%
    rename(distance = dist) %>%
    mutate(catch_rate = log ((count + 1) / days),
           midge_rate = midges / m_days,
           midges_z = (midge_rate - mean(midge_rate)) / sd(midge_rate),
           dist_z = (distance - mean(distance)) / sd(distance),
           time_z = (year - min(year)) / sd(year),
           plot = factor(paste(distance, trans)),
           taxon = factor(taxon)) %>%
    group_by(plot, taxon) %>%
    arrange(plot, taxon, year) %>%
    mutate(midges_lag = lag(midges_z)) %>%
    filter(!is.na(midges_lag)) %>%
    group_by(taxon) %>%
    mutate(y = (catch_rate - mean(catch_rate)) / sd(catch_rate)) %>%
    ungroup() %>%
    arrange(trans, distance, taxon, year)

# write_csv(data_fit, "analysis/data_fit.csv")





#==========
#==========  Model formulas
#==========

# specify model names
model_names <- c("mltd_mltd",
                 "mltd_ltd",
                 "mltd_mtd",
                 "mltd_mld",
                 "mltd_mlt")

# specify model forumlas
model_list <- list(as.formula(y ~ midges_z + midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + midges_lag + time_z + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_lag + time_z + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + time_z + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + midges_lag + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + midges_lag + time_z | taxon)))

# set model names
names(model_list) <- model_names

# fit models
# lapply(model_names,
#        function(x_){
#            formula_ <- model_list[[x_]]
#            fit_ <- armm(formula = formula_,
#                        time_form = ~ year | plot + taxon,
#                        ar_form = ~ taxon,
#                        obs_error = TRUE,
#                        distr = "normal",
#                        data = data_fit,
#                        x_scale = FALSE,
#                        y_scale = NULL,
#                        hmc = TRUE,
#                        change = TRUE,
#                        rstan_control = list(iter = 4000, chains = 4, seed = 3e3,
#                                             control = list(adapt_delta = 0.97)))
#            write_rds(fit_,
#                      paste0("analysis/output/model_fits/",x_,".rds"))
#
#            rm(fit_)
#        })
