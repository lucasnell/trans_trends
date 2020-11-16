#==========
#========== Preliminaries
#==========

# load packages
library(armmr)
library(loo)
library(tidyverse)
options(mc.cores = parallel::detectCores()-4)

# load data
myv_arth <- read_csv("data/myv_arth.csv")

# set theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0,0,0,0),
                    strip.text = element_text(size=12),
                    legend.text = element_text(size=12),
                    legend.title = element_text(size=14),
                    axis.text = element_text(size=12, color="black"),
                    axis.title.y = element_text(size=14,angle = 90,
                                                margin=margin(0,15,0,0)),
                    axis.title.x = element_text(size=14,margin=margin(15,0,0,0)),
                    strip.text.x = element_text(margin=margin(0,0,10,0)),
                    strip.text.y = element_text(margin=margin(0,0,0,10), angle=270),
                    axis.title = element_text(size=14)))





#==========
#========== Prepare data
#==========

# prepare data
data_fit <- myv_arth %>%
    filter(taxon != "acar") %>%
    rename(distance = dist) %>%
    # group_by(taxon) %>%
    mutate(y = log1p(count),
           y = (y - mean(y))/(sd(y))) %>%
    ungroup() %>%
    mutate(midges_z = log1p(midges),
           midges_z = (midges_z - mean(midges_z))/(sd(midges_z)),
           time = factor(year, levels = c(1:max(year))) %>% as.numeric(),
           time = time - min(time),
           time_z = (time)/(sd(time)),
           dist  = log(distance),
           dist_z = (dist - mean(dist))/(sd(dist)),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon),
           plot = factor(paste0(trans, dist)),
           taxon_plot = factor(paste0(taxon, plot)),
           taxon_trans = factor(paste0(taxon, trans))) %>%
    arrange(trans, distance, taxon, year)

# write_csv(data_fit, "analysis/data_fit.csv")






#==========
#========== Full model
#==========

# fit model
# start_time <- Sys.time()
# fit <- armm(formula = y ~ midges_z + time_z + dist_z +
#                 (1 | taxon + taxon_plot + taxon_trans) +
#                 (midges_z + time_z + dist_z | taxon),
#            time_form = ~  time | trans + distf + taxon,
#            ar_form = ~ taxon,
#            obs_error = T,
#            distr = "normal",
#            data = data_fit,
#            x_scale = FALSE,
#            y_scale = NULL,
#            hmc = T,
#            change = T,
#            rstan_control = list(iter = 4000, chains = 4, seed = 3e3,
#                                 control = list(adapt_delta = 0.97)))
# end_time <- Sys.time()
# end_time - start_time

# export fit
# write_rds(fit, "analysis/output/fit.rds")

# import fit
# fit <- read_rds("analysis/output/fit.rds")

# summarize
# fit_sum <- rstan::summary(fit$stan, probs = c(0.16,0.50,0.84))$summary %>%
#     as.data.frame() %>%
#     rownames_to_column() %>%
#     as_tibble() %>%
#     rename(var = rowname,
#            lo = `16%`,
#            mi = `50%`,
#            hi = `84%`) %>%
#     select(var, lo, mi, hi, n_eff, Rhat)

# write_csv(fit_sum, "analysis/output/fit_sum.csv")

# import fit summary
# fit_sum <- read_csv("analysis/output/fit_sum.csv")

# apply loo and export
# fit_loo <- loo(fit$stan, cores = 10)
# write_rds(fit_loo, "analysis/output/fit_loo.rds")
# fit_loo <- read_rds("analysis/output/fit_loo.rds")



#==========
#========== Extract estimates
#==========

# create taxon lists
taxa_short <- tibble(taxon = data_fit$taxon %>% unique()) %>%
    mutate(id = row_number())
taxa_long <- data_fit %>%
    expand(nesting(plot, taxon)) %>%
    select(plot, taxon) %>%
    mutate(id = row_number())

# ar coefficients
ar <- fit_sum %>%
    filter(str_detect(fit_sum$var, "phi")) %>%
    mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
    full_join(taxa_short)

# taxon-specific slopes
beta <- fit_sum %>%
    filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) %>%
    mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = c(1:4),
                         labels = c("int","midges","time","dist"))) %>%
    full_join(taxa_long) %>%
    filter(coef != "int") %>%
    group_by(coef, taxon) %>%
    summarize(lo = unique(lo),
              mi = unique(mi),
              hi = unique(hi))

# intercepts
int_full <- fit_sum %>%
    filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) %>%
    mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = c(1:4),
                         labels = c("int","midges","time","dist"))) %>%
    full_join(taxa_long) %>%
    filter(coef == "int")

# taxon-specific intercepts
# average over plot and transect variation
# standardize to mean for each predictor (only nonzero for time_z)
beta_pars <- {fit_sum %>%
        filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig"))}$var
int_taxon <- rstan::extract(fit$stan, beta_pars) %>%
    bind_cols() %>%
    mutate(step = row_number()) %>%
    gather(var, val, -step) %>%
    mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = c(1:4),
                         labels = c("int","midges","time","dist"))) %>%
    full_join(taxa_long) %>%
    filter(coef %in% c("int","time"))  %>%
    group_by(step, taxon, coef) %>%
    summarize(val = mean(val)) %>%
    spread(coef, val) %>%
    group_by(taxon) %>%
    mutate(val = int + time * mean(data_fit$time_z)) %>%
    summarize(lo = quantile(val, probs = 0.16),
              mi = median(val),
              hi = quantile(val, probs = 0.84))


# mean slopes
alpha <- fit_sum %>%
    filter(str_detect(fit_sum$var, "alpha")) %>%
    mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = factor(coef, levels = c(1:4),
                         labels = c("int","midges","time","dist")))

# sigmas
sig_beta <- fit_sum %>%
    filter(str_detect(fit_sum$var, "sig_beta")) %>%
    mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = factor(coef, levels = c(1:6),
                         labels = c("int_tax","int_tax_plot","int_tax_trans",
                                    "midges","time","dist")))

coef_sum <- list(ar = ar, beta = beta, int_full = int_full, int_taxon = int_taxon,
                      alpha = alpha, sig_beta = sig_beta)

# write_rds(coef_sum, "analysis/output/coef_sum.rds")





#==========
#========== Comparison to reduced models (remove random slopes by taxon)
#==========

# set model names
model_names <- c("midges_z","time_z","dist_z")

# reduced models, removing random slopes by taxon
red_re <- c("y ~ midges_z + time_z + dist_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (time_z + dist_z | taxon)",
             "y ~ midges_z + time_z + dist_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (midges_z + dist_z | taxon)",
             "y ~ midges_z + time_z + dist_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (midges_z + time_z | taxon)")

# fit models
# red_re_fits <- lapply(red_re, function(x){m0 =
#     armm(formula = as.formula(x),
#          time_form = ~  time | trans + distf + taxon,
#          ar_form = ~ taxon,
#          obs_error = T,
#          distr = "normal",
#          data = data_fit,
#          x_scale = FALSE,
#          y_scale = NULL,
#          hmc = T,
#          change = T,
#          rstan_control = list(iter = 4000, chains = 4, seed = 3e3,
#                               control = list(adapt_delta = 0.97)))
# })

# extract stan output
red_re_stan <- lapply(red_re_fits, function(x){x$stan}) %>%
    set_names(model_names)

# write_rds(red_re_stan, "analysis/output/reduced_fits_re.rds")
# red_re_stan <- read_rds("analysis/output/reduced_fits_re.rds")

# apply loo and export
# loo_re <- c("midges_z","time_z","dist_z") %>%
#     lapply(function(x){
#         loo(red_re_stan[[x]], cores = 10)}
#     ) %>%
#     set_names(model_names)

# write_rds(loo_re, "analysis/output/loo_re.rds")
# loo_re <- read_rds("analysis/output/loo_re.rds")


# calculate deviances

dev_re <- lapply(model_names, function(x){
    loo_compare(fit_loo, loo_re[[x]]) %>%
        {as_tibble(.) %>%
                mutate(model = rownames(.))} %>%
        filter(model == 2) %>%
        mutate(model = x,
               delt_looic = -2 * elpd_diff) %>%
        select(model, delt_looic)
}) %>%
    bind_rows()

# write_csv(dev_re, "analysis/output/dev_re.csv")
# dev_re <- read_csv("analysis/output/dev_re.csv")




#==========
#========== Comparison to reduced models (remove fixed + random slopes by taxon)
#==========

# set model names
model_names <- c("midges_z","time_z","dist_z")

# reduced models, removing fixed and random slopes andby taxon
red_fere <- c("y ~ time_z + dist_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (time_z + dist_z | taxon)",
              "y ~ midges_z + dist_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (midges_z + dist_z | taxon)",
              "y ~ midges_z + time_z +
                (1 | taxon + taxon_plot + taxon_trans) +
                (midges_z + time_z | taxon)")

# fit models
# red_fere_fits <- lapply(red_fere, function(x){m0 =
#     armm(formula = as.formula(x),
#          time_form = ~  time | trans + distf + taxon,
#          ar_form = ~ taxon,
#          obs_error = T,
#          distr = "normal",
#          data = data_fit,
#          x_scale = FALSE,
#          y_scale = NULL,
#          hmc = T,
#          change = T,
#          rstan_control = list(iter = 4000, chains = 4, seed = 3e3,
#                               control = list(adapt_delta = 0.97)))
# })

red_fere_names <- c("midges_z","time_z","dist_z","full")

# extract stan output
red_fere_stan = lapply(red_fere_fits, function(x){x$stan}) %>%
    set_names(model_names)

# write_rds(red_fere_stan, "analysis/output/reduced_fits_fere.rds")
# red_fere_stan <- read_rds("analysis/output/reduced_fits_fere.rds")

# apply loo and export
# loo_fere <- c("midges_z","time_z","dist_z") %>%
#     lapply(function(x){
#         loo(red_fere_stan[[x]], cores = 10)}
#     ) %>%
#     set_names(model_names)

# write_rds(loo_fere, "analysis/output/loo_fere.rds")
# loo_fere <- read_rds("analysis/output/loo_fere.rds")

# calculate deviances

dev_fere <- lapply(model_names, function(x){
    loo_compare(fit_loo, loo_fere[[x]]) %>%
        {as_tibble(.) %>%
                mutate(model = rownames(.))} %>%
        filter(model == 2) %>%
        mutate(model = x,
               delt_looic = -2 * elpd_diff) %>%
        select(model, delt_looic)
}) %>%
    bind_rows()

# write_csv(dev_fere, "analysis/output/dev_fere.csv")
# dev_fere <- read_csv("analysis/output/dev_fere.csv")










