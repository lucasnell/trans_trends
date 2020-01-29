#==========
#========== Preliminaries
#==========

# load packages
library(lizard)
library(tidyverse)
options(mc.cores = parallel::detectCores()-2)

# load data
data_fit <- read_csv("analysis/data_fit.csv")
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")


# set theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0,0,0,0),
                    strip.text = element_text(size=12),
                    legend.text = element_text(size=12),
                    legend.title = element_text(size=14),
                    axis.text = element_text(size=12, color="black"),
                    axis.title.y = element_text(size=14,angle = 90 ,margin=margin(0,15,0,0)),
                    axis.title.x = element_text(size=14,margin=margin(15,0,0,0)),
                    strip.text.x = element_text(margin=margin(0,0,10,0)),
                    strip.text.y = element_text(margin=margin(0,0,0,10), angle=270),
                    axis.title = element_text(size=14)))





#==========
#========== Observed data
#==========

# time
data_fit %>%
    ggplot(aes(year, y))+
    facet_wrap(~Taxon, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    geom_point(aes(group = plot), size = 1.5, alpha = 0.5)+
    geom_line(data = data_fit %>%
                  group_by(Taxon, year) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge catch") %>%
                  bind_rows(data_fit %>%
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
data_fit %>%
    ggplot(aes(distance, y))+
    facet_wrap(~Taxon, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_jitter(aes(group = plot), size = 1.5, alpha = 0.5, width = 0.1)+
    geom_line(data = data_fit %>%
                  group_by(Taxon, distance) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge catch") %>%
                  bind_rows(data_fit %>%
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
#========== Plot coefficients
#==========

# taxon-specific slopes
coef_sum$beta %>%
    filter(coef != "int") %>%
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil","acar","cara","stap")))) %>%
    ggplot(aes(tx, mi))+
    facet_wrap(~coef, nrow = 3)+
    geom_rect(data = coef_sum$alpha%>%
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
    coord_flip()

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
    geom_rect(data  = coef_sum$sig_beta %>%
                  filter(!(coef %in% c("int_tax","int_plot","int_trans"))),
              aes(xmin = lo, xmax = hi, ymin = 0, ymax = 8, fill = coef), inherit.aes = F,
              alpha = 0.3, linetype =  0)+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual("",values = c("firebrick","gray30","dodgerblue"))+
    scale_y_continuous("Posterior Density")+
    scale_x_continuous("Effect SD among taxa",  limits = c(0, 1.5))+
    theme(legend.position = c(0.8,0.8))


# ar
coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil","acar","cara","stap")))) %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = median(coef_sum$ar$mi),color = "gray50")+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0)+
    scale_y_continuous("Autoregressive parameter",  limits = c(0, 1), breaks = c(0, 0.5, 1))+
    scale_x_continuous("Taxon",breaks = 1:7,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Mites","Ground beetles","Rove beetles"),
                       limits = c(0.5, 7.5))+
    coord_flip()






#==========
#========== PCA
#==========

source("analysis/pca_funs.R")

# pca on predicted values
pred_pca <- pred_pca_fn(data_fit, coef_sum$beta, coef_sum$int_taxon)

# taxon vectors
pred_pca$taxon_vec %>%
    arrange(-abs(PC1))

# variance explained in predicted values (first three axes must explain everything)
summary(pred_pca$pca)

# variance explained in observed values
pred_pca$obs_exp

# plot midge effect
biplot_fn(pred_pca, "Midges")

# plot time effect
biplot_fn(pred_pca, "Time")

# plot distance effect
biplot_fn(pred_pca, "Distance")


# correlation between predictors and pc axes
lapply(c("midges_z","time_z","dist_z"), function(x){
    pred_pca$axes %>%
        select(id, midges_z, time_z, dist_z, PC1, PC2, PC3) %>%
        gather(pc, val, PC1, PC2, PC3) %>%
        select(x, pc, val) %>%
        split(.$pc) %>%
        lapply(function(y){
            unique(y) %>% {tibble(var = x, pc = unique(.$pc), cor = cor(.[,1], .[,3]) %>% round(2))}
        }) %>%
        bind_rows()
}) %>%
    bind_rows() %>%
    spread(pc, cor)

# variance partition of PC axes by predictors
lapply(c("PC1","PC2","PC3"), function(x){
    y = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(var = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(y, data = pred_pca$axes))[,2] %>% {.[1:3]/sum(.[1:3])} %>% round(3)
    )
}) %>%
    bind_rows() %>%
    spread(pc, cont)
