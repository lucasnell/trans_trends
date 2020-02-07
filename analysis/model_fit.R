#==========
#========== Preliminaries
#==========

# load packages
library(lizard)
library(loo)
library(tidyverse)
options(mc.cores = parallel::detectCores()-2)

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
    group_by(taxon) %>%
    mutate(y = log1p(count),
           y = (y - mean(y))/(sd(y))) %>%
    ungroup() %>%
    mutate(midges_z = log1p(midges),
           midges_z = (midges_z - mean(midges_z))/(sd(midges_z)),
           time = factor(year, levels = c(1:max(year))) %>% as.numeric(),
           time = time - min(time),
           time_z = (time - mean(time))/(sd(time)),
           dist  = log(distance),
           dist_z = (dist - mean(dist))/(sd(dist)),
           trans = factor(trans),
           distf = factor(dist),
           taxon = factor(taxon),
           plot = factor(paste0(trans, dist))) %>%
    arrange(trans, distance, taxon, year)

# write_csv(data_fit, "analysis/data_fit.csv")






#==========
#========== Full model
#==========

# fit model
# fit <- liz_fit(formula = y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
#               (midges_z + time_z + dist_z | taxon),
#            time_form = ~  time | trans + distf + taxon,
#            ar_form = ~ taxon,
#            data = data_fit,
#            x_scale = FALSE,
#            y_scale = NULL,
#            hmc = T,
#            change = T,
#            rstan_control = list(iter = 2000, chains = 4, seed = 3e3,
#                                 control = list(adapt_delta = 0.95)))

# export fit
# saveRDS(fit, "analysis/output/fit.rds")

# import fit
# fit <- readRDS("analysis/output/fit.rds")

# summarize
fit_sum <- rstan::summary(fit$stan, probs = c(0.16,0.50,0.84))$summary %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(var = rowname,
           lo = `16%`,
           mi = `50%`,
           hi = `84%`) %>%
    select(var, lo, mi, hi)

# write_csv(fit_sum, "analysis/output/fit_sum.csv")






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
int_taxon <- int_full %>%
    group_by(taxon) %>%
    summarize(lo = mean(lo),
              mi = mean(mi),
              hi = mean(hi))


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
                         labels = c("int_tax","int_plot","int_trans",
                                    "midges","time","dist")))

coef_sum <- list(ar = ar, beta = beta, int_full = int_full, int_taxon = int_taxon,
                      alpha = alpha, sig_beta = sig_beta)

# saveRDS(coef_sum, "analysis/output/coef_sum.rds")








#==========
#========== Comparison to reduced models (remove random slopes by taxon)
#==========

# function for LOOIC
loo_fn <- function(x){

    log_lik <- extract_log_lik(x, merge_chains = FALSE)
    r_eff <- relative_eff(exp(log_lik))

    loo <- loo(log_lik, r_eff = r_eff, cores = 10)

    return(loo)
}

# reduced models, removing random slopes by taxon
red_re <- c("y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
              (time_z + dist_z | taxon)",
             "y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
              (midges_z + dist_z | taxon)",
             "y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
              (midges_z + time_z | taxon)"
             )

# fit models
# red_re_fits <- lapply(red_re, function(x){m0 =
#   liz_fit(formula = as.formula(x),
#           time_form = ~  time | trans + distf + taxon,
#           ar_form = ~ taxon,
#           data = data_fit,
#           x_scale = FALSE,
#           y_scale = NULL,
#           hmc = T,
#           change = T,
#           rstan_control = list(iter = 2000, chains = 4, seed = 3e3,
#                                control = list(adapt_delta = 0.95)))
# })

# extract stan output
red_re_stan <- lapply(red_re_fits, function(x){x$stan})
names(red_re_stan) <- c("midges_z","time_z","dist_z")

# saveRDS(red_re_stan, "analysis/output/reduced_fits_re.rds")
# red_re_stan <- readRDS("analysis/output/reduced_fits_re.rds")

# calculate deviances
dev_re <- c("midges_z","time_z","dist_z") %>%
    parallel::mclapply(function(x){
        as_tibble(loo_compare(loo_fn(red_re_stan[[x]]), loo_fn(fit$stan))) %>%
            mutate(var = c("full",x))}
    ) %>%
    bind_rows() %>%
    unique()

# write_csv(dev_re, "analysis/output/dev_re.csv")





#==========
#========== Comparison to reduced models (remove fixed + random slopes by taxon)
#==========

# function for LOOIC
loo_fn <- function(x){

    log_lik <- extract_log_lik(x, merge_chains = FALSE)
    r_eff <- relative_eff(exp(log_lik))

    loo <- loo(log_lik, r_eff = r_eff, cores = 10)

    return(loo)
}

# reduced models, removing random slopes by taxon
red_fere <- c("y ~ time_z + dist_z + (1 | taxon + plot + trans) +
            (time_z + dist_z | taxon)",
            "y ~ midges_z + dist_z + (1 | taxon + plot + trans) +
            (midges_z + dist_z | taxon)",
            "y ~ midges_z + time_z + (1 | taxon + plot + trans) +
            (midges_z + time_z | taxon)"
)

# fit models
# red_fere_fits <- lapply(red_fere, function(x){m0 =
#     liz_fit(formula = as.formula(x),
#             time_form = ~  time | trans + distf + taxon,
#             ar_form = ~ taxon,
#             data = data_fit,
#             x_scale = FALSE,
#             y_scale = NULL,
#             hmc = T,
#             change = T,
#             rstan_control = list(iter = 2000, chains = 4, seed = 3e3,
#                                 control = list(adapt_delta = 0.95)))
# })

red_fere_names <- c("midges_z","time_z","dist_z","full")

# extract stan output
red_fere_stan = lapply(red_fere_fits, function(x){x$stan})
names(red_fere_stan) <- c("midges_z","time_z","dist_z")

# saveRDS(red_fere_stan, "analysis/output/reduced_fits_fere.rds")
# red_fere_stan <- readRDS("analysis/output/reduced_fits_fere.rds")

# calculate deviances
dev_fere <- c("midges_z","time_z","dist_z") %>%
    parallel::mclapply(function(x){
        as_tibble(loo_compare(loo_fn(red_fere_stan[[x]]), loo_fn(fit$stan))) %>%
            mutate(var = c("full",x))}
    ) %>%
    bind_rows() %>%
    unique()

# write_csv(dev_fere, "analysis/output/dev_fere.csv")





#==========
#========== Comparison to reduced models (combined)
#==========

# Import for Table I:
# dev_re <- read_csv("analysis/output/dev_re.csv")
# dev_fere <- read_csv("analysis/output/dev_fere.csv")


# Table I (exclude standard errors)
dev <- dev_re %>%
    mutate(dev = 2*abs(elpd_diff),
           dev_se = 2*se_diff) %>%
    select(var, dev, dev_se) %>%
    rename( dev_re = dev, dev_re_se = dev_se) %>%
    full_join(dev_fere %>%
                  mutate(dev = 2*abs(elpd_diff),
                         dev_se = 2*se_diff) %>%
                  select(var,dev, dev_se) %>%
                  rename(dev_fere = dev, dev_fere_se = dev_se)) %>%
    filter(var != "full")





