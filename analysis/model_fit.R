#==========
#========== Preliminaries
#==========

# load packages
library(TransTrendsPkg)
library(loo)
library(tidyverse)
options(mc.cores = parallel::detectCores()-4)

# load data
myv_arth <- read_csv("data/myv_arth.csv")






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



# Gets uncertainty intervals via HPDI, and estimates via mode
get_fit_hpdi <- function(.var) {
    z <- do.call(c, lapply(1:fit$stan@sim$chains,
                           function(i) {
                               fit$stan@sim$samples[[i]][[.var]][
                                   -(1:fit$stan@sim$warmup2[i])]
                           }))
    z_hpdi <- hpdi(z, 0.68)
    .mi <- posterior_mode(z)
    return(tibble(var = .var,
                  lo = z_hpdi[["lower"]],
                  mi = .mi,
                  hi = z_hpdi[["upper"]]))
}


# # summarize
# fit_sum <- rstan::summary(fit$stan, probs = c()) %>%
#     .[["summary"]] %>%
#     as.data.frame() %>%
#     rownames_to_column() %>%
#     as_tibble() %>%
#     rename(var = rowname) %>%
#     left_join(map_dfr(.$var, get_fit_hpdi), by = "var") %>%
#     select(var, lo, mi, hi, n_eff, Rhat)



# write_csv(fit_sum, "analysis/output/fit_sum.csv")

# import fit summary
# fit_sum <- read_csv("analysis/output/fit_sum.csv")




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
    filter(str_detect(var, "beta"), !str_detect(var, "sig")) %>%
    mutate(id = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = 1:4,
                         labels = c("int","midges","time","dist"))) %>%
    full_join(taxa_long) %>%
    filter(coef != "int") %>%
    group_by(coef, taxon) %>%
    summarize(lo = unique(lo),
              mi = unique(mi),
              hi = unique(hi)) %>%
    ungroup()

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
beta_pars <- fit_sum %>%
    filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) %>%
    .[["var"]]
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
    ungroup() %>%
    spread(coef, val) %>%
    mutate(val = int + time * mean(data_fit$time_z)) %>%
    split(.$taxon) %>%
    map_dfr(function(.x) {
        ci <- hpdi(.x$val, 0.68)
        .dd <- tibble(taxon = .x$taxon[1],
                      lo = ci[["lower"]],
                      mi = posterior_mode(.x$val),
                      hi = ci[["upper"]])
        return(.dd)
    })




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
