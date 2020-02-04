#==========
#========== Preliminaries
#==========

# load packages
library(lizard)
library(tidyverse)
options(mc.cores = parallel::detectCores()-2)

if (Sys.info()[["sysname"]] == "Darwin") {
    options("device" = "quartz")
    grDevices::graphics.off()
}


# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(Taxon = factor(taxon,
                          levels = c("gnap","lyco","sheet","opil","cara","stap"),
                          labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                     "Harvestman","Ground beetles","Rove beetles")))
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")


# set theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0, 0, 0, 0),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    axis.text = element_text(size = 12, color = "black"),
                    axis.title.y = element_text(size = 14, angle = 90,
                                                margin = margin(0,15,0,0)),
                    axis.title.x = element_text(size = 14, margin = margin(15,0,0,0)),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6), vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10), angle = 270),
                    axis.title = element_text(size = 14)))





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
                  mutate(type = factor(type,
                                       levels = c("Midge catch", "Mean activity-density"))),
              aes(color = type), size = 1)+
    scale_color_manual("",values = c("red","blue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016), limits = c(2007,2018))+
    scale_y_continuous("Standardized activity-density",
                       limits = c(-4.5,4.5), breaks = c(-3,0,3))+
    theme(legend.position = "top")

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
                  mutate(type = factor(type,
                                       levels = c("Midge catch",
                                                  "Mean activity-density"))),
              aes(color = type), size = 1)+
    scale_color_manual("",values = c("red","blue"))+
    scale_x_continuous("Distance", trans = "log", breaks = c(5,50,500))+
    scale_y_continuous("Standardized activity-density", limits = c(-4.5,4.5),
                       breaks = c(-3,0,3))+
    theme(legend.position = "top")





#==========
#========== Plot coefficients
#==========

# taxon-specific slopes
coef_sum$beta %>%
    filter(coef != "int") %>%
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet",
                                                    "opil","cara","stap")))) %>%
    ggplot(aes(tx, mi))+
    facet_wrap(~coef, nrow = 3)+
    geom_rect(data = coef_sum$alpha%>%
                  filter(coef != "int"),
              aes(xmin = 0.5, xmax = 6.5, ymin = lo, ymax = hi), inherit.aes = F,
              alpha = 0.5, fill = "gray50", linetype =  0)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0)+
    scale_y_continuous("Effect size (response SD)", breaks = c(-0.5, 0, 0.5),
                       limits = c(-0.6, 0.6))+
    scale_x_continuous("Taxon",breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       limits = c(0.5, 6.5))+
    coord_flip()

# sigmas
{fit$stan %>% rstan::extract(pars = "sig_beta")}[[1]] %>%
    as.matrix() %>%
    as_tibble() %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("V1","V2","V3","V4","V5","V6"),
                         labels = c("int_tax","int_plot","int_trans","midges"
                                    ,"time","dist"))) %>%
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
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil",
                                                    "cara","stap")))) %>%
    ggplot(aes(tx, mi))+
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0)+
    scale_y_continuous("Autoregressive parameter",  limits = c(0, 1),
                       breaks = c(0, 0.5, 1))+
    scale_x_continuous("Taxon",breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       limits = c(0.5, 6.5))+
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

# biplot function
biplot_fn <- function(pred_pca_, var_) {
    var_fn <- function(var__) {
        if(var_ == "Midges") return(pred_pca_$obs_rot$midges_z)
        if(var_ == "Time") return(pred_pca_$obs_rot$time_z)
        if(var_ == "Distance") return(pred_pca_$obs_rot$dist_z)
    }
    pred_pca_$obs_rot %>%
        ggplot(aes(PC1, -PC2))+
        geom_point(aes(color = var_fn(var_)), alpha = 0.6, size = 3)+
        geom_segment(data = pred_pca_$taxon_vec,
                     aes(x = 0, xend = 3.5*PC1, y = 0, yend = -3.5*PC2, group = taxon),
                     arrow = arrow(length = unit(0.5, "cm")),
                     size = 0.8, color = "black")+
        geom_text(data = pred_pca_$taxon_vec,
                  aes(label = taxon, group = taxon, x = 3.8*PC1, y = -3.8*PC2),
                  size = 4, color = "black")+
        scale_colour_gradient2(var_,low =  "firebrick", mid = "gray70", high = "royalblue",
                               midpoint = -0.5, limits  = c(-3,2), breaks = c(-3,-0.5,2))+
        coord_equal()
}

# plot midge effect
biplot_fn(pred_pca, "Midges")

# plot time effect
biplot_fn(pred_pca, "Time")

# plot distance effect
biplot_fn(pred_pca, "Distance")


# variance partition of PC axes by predictors
var_part <- lapply(c("PC1","PC2","PC3"), function(x){
    y = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(var = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(y, data = pred_pca$axes))[,2] %>% {.[1:3]/sum(.[1:3])} %>%
               round(3)
    )
}) %>%
    bind_rows() %>%
    spread(pc, cont)

var_part2 <- as.matrix(var_part[1:3,2:4]) %*% t(as.matrix(pred_pca$obs_exp[1,2:4]))
row.names(var_part2) <- var_part$var

