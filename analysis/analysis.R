#==========
#========== Preliminaries
#==========

# load packages
library(lizard)
library(tidyverse)
library(cowplot)
options(mc.cores = parallel::detectCores()-2)

# In RStudio on macOS, this makes a separate plotting window that's independent of
# the Plots pane
if (Sys.info()[["sysname"]] == "Darwin" && .Platform$GUI == "RStudio") {
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


coef_sum$beta <- coef_sum$beta %>%
    ungroup() %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "midges", "dist")),
           tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet",
                                                    "opil","cara","stap"))))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "midges", "dist")))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "midges", "dist")))
coef_sum$ar <- coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil",
                                                    "cara","stap"))))



# ------------*
# Define theme and colorblind-friendly palettes
# ------------*
coef_palette <- viridisLite::plasma(3, end = 0.8)
pca_palette <- viridisLite::viridis(3)


# set theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 12),
                    axis.text = element_text(size = 10, color = "black"),
                    axis.title.y = element_text(size = 12, angle = 90,
                                                margin = margin(0,0,0,r=6)),
                    axis.title.x = element_text(size = 12, margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 10),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6), vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10), angle = 270)))





#==========
#========== Observed data (Combine 'time' and 'distance' into Figure 1)
#==========

# time
time_p <- data_fit %>%
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
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016), limits = c(2007,2018))+
    scale_y_continuous("Standardized activity-density",
                       limits = c(-4.5,4.5), breaks = c(-3,0,3))+
    # theme(legend.position = "none") +
    NULL

# distance
dist_p <- data_fit %>%
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
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    # scale_y_continuous("Standardized activity-density", limits = c(-4.5,4.5),
    scale_y_continuous(NULL, limits = c(-4.5,4.5),
                       breaks = c(-3,0,3)) +
    NULL


time_dist_legend <- get_legend(
    time_p +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)


prow <- plot_grid(time_p + theme(legend.position = "none"),
                  dist_p + theme(legend.position = "none"),
                  labels = c("A", "B"),
                  align = "vh")

# cairo_pdf("analysis/output/fig1.pdf", width = 6, height = 6)
# plot_grid(time_dist_legend, prow, ncol = 1, rel_heights = c(0.1, 1))
# dev.off()



#==========
#========== Plot coefficients (Use `taxon-specicific-slopes` for Figure 2)
#==========



# taxon-specific slopes
slope_p <- coef_sum$beta %>%
    ggplot()+
    facet_wrap(~coef, nrow = 1)+
    geom_rect(data = coef_sum$alpha,
              aes(xmin = 0.5, xmax = 6.5, ymin = lo, ymax = hi, fill = coef),
              alpha = 0.5, linetype =  0)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(aes(tx, mi), size = 1.5)+
    geom_linerange(aes(tx, ymin = lo, ymax = hi))+
    # geom_text(data = tibble(tx = 0.5, mi = -0.6,
    #                         coef = factor(c("time", "midges", "dist"))),
    #           aes(tx, mi, label = coef), size = 10/2.83465,
    #           hjust = 0, vjust = 1, fontface = "bold") +
    scale_y_continuous(expression("Effect size" ~ ({}^{k} * beta[taxon])),
                       breaks = c(-0.5, 0, 0.5), limits = c(-0.6, 0.6))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       trans = "reverse", limits = c(6.5, 0.5)) +
    scale_fill_manual(NULL, values = coef_palette, guide = FALSE)+
    scale_color_manual(NULL, values = coef_palette, guide = FALSE)+
    coord_flip() +
    theme(#strip.text.x = element_blank(),
          axis.title.x = element_text(margin=margin(0,0,0,0))) +
    NULL

small_labels <- theme(axis.title.y = element_text(size = 10, angle = 90,
                                                  margin = margin(0,0,0,r=2)),
                      axis.title.x = element_text(size = 10, margin = margin(0,0,0,t=2)),
                      axis.text = element_text(size = 8, color = "black"))


# ar (extra)
ar_p <- coef_sum$ar %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 1.5)+
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_y_continuous(expression("AR parameter" ~ (phi)),  limits = c(0, 1),
                       breaks = c(0, 0.5, 1))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       trans = "reverse", limits = c(6.5, 0.5))+
    coord_flip() +
    theme(axis.title.x = element_text(size = 10, margin = margin(0,0,0,t=2)),
          axis.text.x = element_text(size = 8, color = "black")) +
    NULL



# sigmas (extra)
effect_sigmas_p <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    .[[1]] %>%
    as.matrix() %>%
    as_tibble() %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("V1","V2","V3","V4","V5","V6"),
                         labels = c("int_tax","int_plot","int_trans","midges",
                                    "time","dist"))) %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "midges", "dist"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_rect(data  = coef_sum$sig_beta,
              aes(xmin = lo, xmax = hi, ymin = 0, ymax = 8, fill = coef), inherit.aes = F,
              alpha = 0.3, linetype =  0)+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual(NULL, values = coef_palette)+
    guides(fill = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
    scale_y_continuous("Density", limits = c(0, 8), breaks = 0:4 * 2)+
    scale_x_continuous(expression("Response SD" ~ ({}^{k} * sigma[g]))) +
    small_labels +
    NULL


effect_mean_p <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    .[[1]] %>%
    as.matrix() %>%
    as_tibble() %>%
    gather() %>%
    filter(key != "V1") %>%  # intercept not necessary
    mutate(coef = factor(key, levels = c("V3","V2","V4"),
                         labels = c("time", "midges", "dist"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_rect(data  = coef_sum$alpha,
              aes(xmin = lo, xmax = hi, ymin = 0, ymax = 8, fill = coef), inherit.aes = F,
              alpha = 0.3, linetype =  0)+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual("", values = coef_palette)+
    scale_y_continuous("Density", limits = c(0, 8), breaks = 0:4 * 2)+
    scale_x_continuous(expression("Mean response (" * {}^{k} * alpha * ")")) +
    theme(legend.position = "none") +
    small_labels +
    NULL




# fig2 <- plot_grid(NULL,
#                   slope_p,
#                   plot_grid(NULL, NULL, NULL, labels = c("B", "C", "D"),
#                             rel_heights = c(1.2, 1, 1), ncol = 1),
#                   plot_grid(ar_p + theme(axis.text.y = element_blank()),
#                             effect_sigmas_p +
#                                 theme(legend.position = c(0.7,0.5),
#                                       plot.margin = margin(t=0,r=4,b=10,l=0)),
#                             effect_mean_p +
#                                 theme(plot.margin = margin(t=0,r=4,b=10,l=0)),
#                             rel_heights = c(1.2, 1, 1),
#                             ncol = 1, align = "h", axis = "t"),
#                   labels = c("A", "", "", ""), rel_widths = c(0.08, 1, 0.08, 0.6),
#                   nrow = 1)
#
#
# cairo_pdf("analysis/output/fig2-A.pdf", width = 6, height = 5)
# fig2
# dev.off()


fig2 <- plot_grid(NULL, slope_p + theme(plot.margin = margin(t=0,r=4,b=10,l=0)),
                  NULL, plot_grid(ar_p +
                                      theme(plot.margin = margin(t=0,r=4,b=10,l=0)),
                                  NULL,
                                  effect_sigmas_p +
                                      theme(legend.position = c(0.7,0.5),
                                            plot.margin = margin(t=0,r=4,b=10,l=0)),
                                  NULL,
                                  effect_mean_p +
                                      theme(plot.margin = margin(t=0,r=4,b=10,l=0),
                                            axis.title.y = element_blank(),
                                            axis.text.y = element_blank()),
                                  rel_widths = c(1.2, 0.15, 1.1, 0.15, 1),
                                  labels = c("", "C", "", "D", ""),
                                  nrow = 1),
                  labels = c("A", "", "B", ""), rel_widths = c(0.05, 1),
                  rel_heights = c(1.1, 1),
                  ncol = 2)


cairo_pdf("analysis/output/fig2-B.pdf", width = 6, height = 4)
fig2
dev.off()



#==========
#========== PCA (Combine biplots for `midge effect`, `time effect`, and `distance` for Fig 3
#==========

source("analysis/pca_funs.R")

# pca on predicted values
pred_pca <- pred_pca_fn(data_fit, coef_sum$beta, coef_sum$int_taxon)

# taxon vectors
pred_pca$taxon_vec %>%
    arrange(-abs(PC1))

# variance explained in predicted values (first three axes must explain everything)
summary(pred_pca$pca)

# variance explained in observed values (For Results)
pred_pca$obs_exp


# variance partition of PC axes by predictors (for Table I)
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

# overall variance accounted for by predictors (for Table I)
overall_part <- as.matrix(var_part[1:3,2:4]) %*% t(as.matrix(pred_pca$obs_exp[1,2:4]))
row.names(overall_part) <- var_part$var

# predictor vectors (example for midges)
pred_pca$axes %>%
    filter(midges_z == min(midges_z)|midges_z == max(midges_z)) %>%
    group_by(midges_z) %>%
    summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3)) %>%
    rename(val = midges_z) %>%
    arrange(val) %>%
    mutate(var = "midges_z")  %>%
    mutate(PC1 = PC1 - PC1[1],
           PC2 = PC2 - PC2[1],
           PC3 = PC3 - PC3[1])

# biplot function
biplot_fn <- function(pred_pca_, var_) {

    var_ <- match.arg(var_, c("midges", "time", "dist"))

    pred_pca_$obs_rot %>%
        ggplot(aes(PC1, -PC2))+
        geom_point(aes(color = pred_pca_$obs_rot[[paste0(var_,"_z")]]),
                   alpha = 0.6, size = 3)+
        ggtitle(var_) +
        # geom_segment(data = pred_pca_$taxon_vec,
        #              aes(x = 0, xend = 3.5*PC1, y = 0, yend = -3.5*PC2, group = taxon),
        #              arrow = arrow(length = unit(0.5, "cm")),
        #              size = 0.8, color = "black")+
        # geom_text(data = pred_pca_$taxon_vec,
        #           aes(label = taxon, group = taxon, x = 3.8*PC1, y = -3.8*PC2),
        #           size = 4, color = "black")+
        scale_colour_gradient2("value", low =  pca_palette[1],
                               mid = pca_palette[2],
                               high = pca_palette[3],
                               midpoint = -0.5, limits = c(-3,2),
                               breaks = c(-3,-0.5,2))+
        coord_equal(xlim = c(-4.1, 4.1), ylim = c(-4.1, 4.1)) +
        theme(plot.title = element_text(size = 11, hjust = 0.5, vjust = 0,
                                        face = "plain"),
              plot.margin = margin(0,0,0,0))
}



fig3a <- pred_pca$taxon_vec %>%
    mutate(PC1_lab = 5*PC1, PC2_lab = -5*PC2,
           PC2_lab = ifelse(taxon == "gnap", PC2_lab * 0.8, PC2_lab),
           PC2_lab = ifelse(taxon == "opil", 0.5, PC2_lab),
           PC1_lab = ifelse(taxon == "opil", 1, PC1_lab)) %>%
    ggplot()+
    geom_segment(aes(x = 0, xend = 3.5*PC1, y = 0, yend = -3.5*PC2, group = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 0.5, color = "black")+
    geom_text(aes(label = taxon, group = taxon, x = PC1_lab, y = PC2_lab),
              size = 9 / 2.83465, color = "black")+
    xlab("PC1") +
    ylab("-PC2") +
    coord_equal(xlim = c(-4.1, 4.1), ylim = c(-4.1, 4.1)) +
    theme(plot.margin = margin(0,0,0,0))

# plot midge effect
fig3b <- biplot_fn(pred_pca, "midges")
# plot time effect
fig3c <- biplot_fn(pred_pca, "time")
# plot distance effect
fig3d <- biplot_fn(pred_pca, "dist")


pca_legend <- get_legend(fig3b + theme(legend.title = element_blank()))



fig3 <- plot_grid(plot_grid(fig3a,
                            NULL,
                            fig3b + theme(legend.position = "none"),
                            NULL, NULL, NULL,
                            fig3c + theme(legend.position = "none"),
                            NULL,
                            fig3d + theme(legend.position = "none"),
                            labels = c("A", "", "B",
                                       "", "", "",
                                       "C", "", "D"),
                            nrow = 3, align = "vh",
                            rel_widths = c(1, 0.1, 1), rel_heights = c(1, 0.1, 1)),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))



# cairo_pdf("analysis/output/fig3.pdf", width = 6, height = 5)
# fig3
# dev.off()


