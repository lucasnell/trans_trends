#==========
#========== Preliminaries
#==========

# load packages
library(lizard)
library(tidyverse)
options(mc.cores = parallel::detectCores()-2)

# load data
myv_arth = read_csv("data/myv_arth.csv")

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))





#==========
#========== Prepare data
#==========

# prepare data
myv_arth2 = myv_arth %>%
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
         plot = factor(paste0(trans, dist)),
         Taxon = factor(taxon,
                        levels = c("gnap","lyco","sheet","opil","acar","cara","stap"),
                        labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                   "Harvestman","Mites","Ground beetles","Rove beetles"))) %>%
    arrange(trans, distance, taxon, year)






#==========
#========== Plot data
#==========

# time
myv_arth2 %>%
  ggplot(aes(year, y))+
  facet_wrap(~Taxon, nrow = 4)+
  geom_hline(yintercept = 0, color = "gray50")+
  geom_line(aes(group = plot), size = 0.2, color = "gray70")+
  geom_point(aes(group = plot), size = 1.5, alpha = 0.5)+
  geom_line(data = myv_arth2 %>%
              group_by(Taxon, year) %>%
              summarize(y = mean(midges_z)) %>%
              mutate(type  = "Midge catch") %>%
              bind_rows(myv_arth2 %>%
                          group_by(Taxon, year) %>%
                          summarize(y = mean(y)) %>%
                          mutate(type  = "Mean activity-density")) %>%
              mutate(type = factor(type, levels = c("Midge catch", "Mean activity-density"))),
            aes(color = type), size = 1)+
  scale_color_manual("",values = c("red","blue"))+
  scale_x_continuous("Year", breaks = c(2010,2013,2016), limits = c(2007,2018))+
  scale_y_continuous("Standardized activity-density", limits = c(-4.5,4.5), breaks = c(-3,0,3))+
  theme(legend.position = c(0.75, 0.1))

# distance
myv_arth2 %>%
  ggplot(aes(distance, y))+
  facet_wrap(~Taxon, nrow = 4)+
  geom_hline(yintercept = 0, color = "gray50")+
  geom_jitter(aes(group = plot), size = 1.5, alpha = 0.5, width = 0.1)+
  geom_line(data = myv_arth2 %>%
              group_by(Taxon, distance) %>%
              summarize(y = mean(midges_z)) %>%
              mutate(type  = "Midge catch") %>%
              bind_rows(myv_arth2 %>%
                          group_by(Taxon, distance) %>%
                          summarize(y = mean(y)) %>%
                          mutate(type  = "Mean activity-density")) %>%
              mutate(type = factor(type, levels = c("Midge catch", "Mean activity-density"))),
            aes(color = type), size = 1)+
  scale_color_manual("",values = c("red","blue"))+
  scale_x_continuous("Distance", trans = "log", breaks = c(5,50,500))+
  scale_y_continuous("Standardized activity-density", limits = c(-4.5,4.5), breaks = c(-3,0,3))+
  theme(legend.position = c(0.75, 0.1))





#==========
#========== Full model
#==========

# fit model
# fit <- liz_fit(formula = y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
#               (midges_z + time_z + dist_z | taxon),
#            time_form = ~  time | trans + distf + taxon,
#            ar_form = ~ taxon,
#            data = myv_arth2,
#            x_scale = FALSE,
#            y_scale = NULL,
#            hmc = T,
#            rstan_control = list(iter = 2000, chains = 4, control = list(adapt_delta = 0.9)))

# export fit
# saveRDS(fit, "analysis/fit.rds")

# import fit
# fit <- readRDS("analysis/fit.rds")

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








#==========
#========== Extract estiamtes
#==========

# create taxon lists
taxa_short <- tibble(taxon = myv_arth2$taxon %>% unique()) %>%
  mutate(id = row_number())
taxa_long <- myv_arth2 %>%
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
         coef = factor(coef, levels = c(1:4), labels = c("int","midges","time","dist"))) %>%
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
         coef = factor(coef, levels = c(1:4), labels = c("int","midges","time","dist"))) %>%
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
         coef = factor(coef, levels = c(1:4), labels = c("int","midges","time","dist")))

# sigmas
sig_beta <- fit_sum %>%
  filter(str_detect(fit_sum$var, "sig_beta")) %>%
  mutate(coef = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         coef = factor(coef, levels = c(1:6),
                       labels = c("int_tax","int_plot","int_trans","midges","time","dist")))




#==========
#========== Plot coefficients
#==========

# taxon-specific slopes
beta %>%
  filter(coef != "int") %>%
  mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil","acar","cara","stap")))) %>%
  ggplot(aes(tx, mi))+
  facet_wrap(~coef)+
  geom_rect(data = alpha%>%
              filter(coef != "int"),
            aes(xmin = 0.5, xmax = 7.5, ymin = lo, ymax = hi), inherit.aes = F,
            alpha = 0.5, fill = "gray50", linetype =  0)+
  geom_hline(yintercept = 0, color = "gray50")+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0)+
  scale_y_continuous("Effect size (response SD)", breaks = c(-0.5, 0, 0.5))+
  scale_x_continuous("Taxon",breaks = 1:7,
                     labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                              "Harvestman","Mites","Ground beetles","Rove beetles"),
                     limits = c(0.5, 7.5))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))

# sigmas
{fit$stan %>% rstan::extract(pars = "sig_beta")}[[1]] %>%
  as.matrix() %>%
  as_tibble() %>%
  gather() %>%
  mutate(coef = factor(key, levels = c("V1","V2","V3","V4","V5","V6"),
                       labels = c("int_tax","int_plot","int_trans","midges","time","dist"))) %>%
  filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
  ggplot(aes(value, fill = coef))+
  geom_vline(xintercept = 0,  color = "gray50")+
  geom_rect(data  = sig_beta %>%
              filter(!(coef %in% c("int_tax","int_plot","int_trans"))),
            aes(xmin = lo, xmax = hi, ymin = 0, ymax = 8, fill = coef), inherit.aes = F,
            alpha = 0.3, linetype =  0)+
  geom_density(linetype = 0, alpha = 0.7)+
  scale_fill_manual("",values = c("firebrick","gray30","dodgerblue"))+
  scale_y_continuous("Posterior Density")+
  scale_x_continuous("Effect SD among taxa",  limits = c(0, 1.5))+
  theme(legend.position = c(0.8,0.8))


# ar
ar %>%
  mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil","acar","cara","stap")))) %>%
  ggplot(aes(tx, mi))+
  geom_hline(yintercept = median(ar$mi),color = "gray50")+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0)+
  scale_y_continuous("Autoregressive parameter",  limits = c(0, 1), breaks = c(0, 0.5, 1))+
  scale_x_continuous("Taxon",breaks = 1:7,
                     labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                "Harvestman","Mites","Ground beetles","Rove beetles"),
                     limits = c(0.5, 7.5))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))





#==========
#========== Comparison to reduced models
#==========

# reduced models, removing random slopes by taxon
# red_re <- c("y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
#               (time_z + dist_z | taxon)",
#              "y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
#               (midges_z + dist_z | taxon)",
#              "y ~ midges_z + time_z + dist_z + (1 | taxon + plot + trans) +
#               (midges_z + time_z | taxon)"
#              )

# fit models
# red_re_fits <- parallel::mclapply(red_re, function(x){m0 =
#   liz_fit(formula = as.formula(x),
#           time_form = ~  time | trans + distf + taxon,
#           ar_form = ~ taxon,
#           data = myv_arth2,
#           x_scale = FALSE,
#           y_scale = NULL,
#           hmc = T,
#           rstan_control = list(iter = 2000, chains = 1, control = list(adapt_delta = 0.8)))
# })

# red_re_names <- c("midges_z","time_z","dist_z","full")
#
#
# loo_fn <- function(x){
#
#   log_lik <- loo::extract_log_lik(x, merge_chains = FALSE)
#   r_eff <- loo::relative_eff(exp(log_lik))
#
#   loo <- loo::loo(log_lik, r_eff = r_eff, cores = 10)
#
#   return(loo)
# }
#
# red_re_fits2 = lapply(red_re_fits, function(x){x$stan}) %>% append(fit$stan)
#
# test <- parallel::mclapply(1:4, function(x){
#     y <- loo_fn(red_re_fits2[[x]])$estimates["looic",]
#     return(tibble(model = red_re_names[x], looic = y["Estimate"]))
#   }) %>%
#   bind_rows() %>%
#   mutate(deviance = looic - {filter(., model == "full")$looic}) %>%
#   filter(model != "full")
#
#
# test
#
# loo_fn(red_re_fits2[[4]])




#==========
#========== PCA
#==========

source("analysis/pca_funs.R")

# pca on predicted values
pred_pca <- pred_pca_fn(myv_arth2, beta, int_taxon)

# variance explained in predicted values (first three axes must explain everything)
summary(pred_pca$pca)

# plot midge effect
pred_pca$obs_rot %>%
    ggplot(aes(PC1, PC2))+
    geom_point(aes(color = midges_z), alpha = 0.6, size = 3)+
    geom_segment(data = pred_pca$taxon_vec,
                 aes(x = 0, xend = 3.5*PC1, y = 0, yend = 3.5*PC2, group = taxon),
                 arrow = arrow(length = unit(0.5, "cm")),
                 size = 0.8, color = "black")+
    geom_text(data = pred_pca$taxon_vec,
              aes(label = taxon, group = taxon, x = 3.8*PC1, y = 3.8*PC2),
              size = 4, color = "black")+
    scale_colour_gradient2("Midges",low =  "firebrick", mid = "gray70", high = "royalblue",
                           midpoint = -0.5, limits  = c(-3,2), breaks = c(-3,-0.5,2))+
    coord_equal()

# plot time effect
pred_pca$obs_rot %>%
    ggplot(aes(PC1, PC2))+
    geom_point(aes(color = time_z), alpha = 0.6, size = 3)+
    geom_segment(data = pred_pca$taxon_vec,
                 aes(x = 0, xend = 3.5*PC1, y = 0, yend = 3.5*PC2, group = taxon),
                 arrow = arrow(length = unit(0.5, "cm")),
                 size = 0.8, color = "black")+
    geom_text(data = pred_pca$taxon_vec,
              aes(label = taxon, group = taxon, x = 3.8*PC1, y = 3.8*PC2),
              size = 4, color = "black")+
    scale_colour_gradient2("Time",low =  "firebrick", mid = "gray70", high = "royalblue",
                           midpoint = -0.5, limits  = c(-3,2), breaks = c(-3,-0.5,2))+
    coord_equal()

# plot distance effect
pred_pca$obs_rot %>%
    ggplot(aes(PC1, PC2))+
    geom_point(aes(color = dist_z), alpha = 0.6, size = 3)+
    geom_segment(data = pred_pca$taxon_vec,
                 aes(x = 0, xend = 3.5*PC1, y = 0, yend = 3.5*PC2, group = taxon),
                 arrow = arrow(length = unit(0.5, "cm")),
                 size = 0.8, color = "black")+
    geom_text(data = pred_pca$taxon_vec,
              aes(label = taxon, group = taxon, x = 3.8*PC1, y = 3.8*PC2),
              size = 4, color = "black")+
    scale_colour_gradient2("Distance",low =  "firebrick", mid = "gray70", high = "royalblue",
                           midpoint = -0.5, limits  = c(-3,2), breaks = c(-3,-0.5,2))+
    coord_equal()



