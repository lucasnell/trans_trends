
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
                 "mltd_mlt",
                 "ltd_ltd",
                 "mtd_mtd",
                 "mld_mld",
                 "mlt_mlt")

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
                                  (0 + midges_z + midges_lag + time_z | taxon)),
                   as.formula(y ~ midges_lag + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_lag + time_z + dist_z | taxon)),
                   as.formula(y ~ midges_z + time_z + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + time_z + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + dist_z +
                                  (1 | taxon + plot) +
                                  (0 + midges_z + midges_lag + dist_z | taxon)),
                   as.formula(y ~ midges_z + midges_lag + time_z+
                                  (1 | taxon + plot) +
                                  (0 + midges_z + midges_lag + time_z | taxon)))

# set model names
names(model_list) <- model_names

# fit models (only the first four; slimmed down for testing)
lapply(model_names[1:4],
       function(x_){
           formula_ <- model_list[[x_]]
           fit_ <- armm(formula = formula_,
                       time_form = ~ year | plot + taxon,
                       ar_form = ~ taxon,
                       obs_error = TRUE,
                       distr = "normal",
                       data = data_fit,
                       x_scale = FALSE,
                       y_scale = NULL,
                       hmc = TRUE,
                       change = TRUE,
                       rstan_control = list(iter = 1, chains = 1, seed = 3e3,
                                            control = list(adapt_delta = 0.95)))
           write_rds(fit_,
                     paste0("analysis/output/model_fits/",x_,".rds"))

           rm(fit_)
       })







#==========
#========== Full model
#==========


# if (.REFIT_MODEL) {
#
#     fit <- armm(formula = model_list$mltd_mltd,
#                 time_form = ~ year | plot + taxon,
#                 ar_form = ~ taxon,
#                 obs_error = TRUE,
#                 distr = "normal",
#                 data = data_fit,
#                 x_scale = FALSE,
#                 y_scale = NULL,
#                 hmc = TRUE,
#                 change = TRUE,
#                 rstan_control = list(iter = 100, chains = 4, seed = 3e3,
#                                      control = list(adapt_delta = 0.97)))
#
#     write_rds(fit, "analysis/output/fit.rds")
#
#     # Gets uncertainty intervals via quantiles, and estimates via median
#     get_fit_info <- function(.var) {
#         z <- do.call(c, lapply(1:fit$stan@sim$chains,
#                                function(i) {
#                                    fit$stan@sim$samples[[i]][[.var]][
#                                        -(1:fit$stan@sim$warmup2[i])]
#                                }))
#         return(tibble(var = .var,
#                       lo = unname(quantile(z, 0.16)),
#                       mi = median(z),
#                       hi = unname(quantile(z, 0.84))))
#     }
#
#     # summarize
#     fit_sum <- rstan::summary(fit$stan, probs = c()) %>%
#         .[["summary"]] %>%
#         as.data.frame() %>%
#         rownames_to_column() %>%
#         as_tibble() %>%
#         rename(var = rowname) %>%
#         left_join(map_dfr(.$var, get_fit_info), by = "var") %>%
#         select(var, lo, mi, hi, n_eff, Rhat)
#
#     write_csv(fit_sum, "analysis/output/fit_sum.csv")
#
# } else {
#
#     fit <- read_rds("analysis/output/fit.rds")
#     fit_sum <- read_csv("analysis/output/fit_sum.csv")
#
# }










#==========
#========== Extract estimates
#==========

# create taxon lists
# taxa_short <- tibble(taxon = data_fit$taxon %>% unique()) %>%
#     mutate(id = row_number())
# taxa_long <- data_fit %>%
#     expand(nesting(plot, taxon)) %>%
#     select(plot, taxon) %>%
#     mutate(id = row_number())
#
# # ar coefficients
# ar <- fit_sum %>%
#     filter(str_detect(fit_sum$var, "phi")) %>%
#     mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
#     full_join(taxa_short)
#
# # taxon-specific slopes
# beta <- fit_sum %>%
#     filter(str_detect(var, "beta"), !str_detect(var, "sig")) %>%
#     mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#            coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
#            coef = factor(coef, levels = 1:4,
#                          labels = c("int","midges","time","dist"))) %>%
#     full_join(taxa_long) %>%
#     filter(coef != "int") %>%
#     group_by(coef, taxon) %>%
#     summarize(lo = unique(lo),
#               mi = unique(mi),
#               hi = unique(hi)) %>%
#     ungroup()
#
# # intercepts
# int_full <- fit_sum %>%
#     filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) %>%
#     mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#            coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
#            coef = factor(coef, levels = c(1:4),
#                          labels = c("int","midges","time","dist"))) %>%
#     full_join(taxa_long) %>%
#     filter(coef == "int")
#
# # taxon-specific intercepts
# # average over plot and transect variation
# # standardize to mean for each predictor (only nonzero for time_z)
# beta_pars <- fit_sum %>%
#     filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) %>%
#     .[["var"]]
# int_taxon <- rstan::extract(fit$stan, beta_pars) %>%
#     bind_cols() %>%
#     mutate(step = row_number()) %>%
#     gather(var, val, -step) %>%
#     mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#            coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
#            coef = factor(coef, levels = c(1:4),
#                          labels = c("int","midges","time","dist"))) %>%
#     full_join(taxa_long) %>%
#     filter(coef %in% c("int","time"))  %>%
#     group_by(step, taxon, coef) %>%
#     summarize(val = mean(val)) %>%
#     ungroup() %>%
#     spread(coef, val) %>%
#     mutate(val = int + time * mean(data_fit$time_z)) %>%
#     split(.$taxon) %>%
#     map_dfr(function(.x) {
#         .dd <- tibble(taxon = .x$taxon[1],
#                       lo = unname(quantile(.x$val, 0.16)),
#                       mi = median(.x$val),
#                       hi = unname(quantile(.x$val, 0.84)))
#         return(.dd)
#     })
#
#
#
#
# # mean slopes
# alpha <- fit_sum %>%
#     filter(str_detect(fit_sum$var, "alpha")) %>%
#     mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#            coef = factor(coef, levels = 1:4,
#                          labels = c("int","midges","time","dist")))
#
# # sigmas
# sig_beta <- fit_sum %>%
#     filter(str_detect(fit_sum$var, "sig_beta")) %>%
#     mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#            coef = factor(coef, levels = 1:6,
#                          labels = c("int_tax","int_tax_plot","int_tax_trans",
#                                     "midges","time","dist")))
#
# coef_sum <- list(ar = ar,
#                  beta = beta,
#                  int_full = int_full,
#                  int_taxon = int_taxon,
#                  alpha = alpha,
#                  sig_beta = sig_beta)

# write_rds(coef_sum, "analysis/output/coef_sum.rds")
