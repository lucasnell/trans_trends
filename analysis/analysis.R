#==========*
#========== Preliminaries ----
#==========*

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
                         levels = c("time", "dist", "midges")),
           tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet",
                                                    "opil","cara","stap"))))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges")))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges")))
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
                    axis.text = element_text(size = 9, color = "black"),
                    axis.title.y = element_text(size = 12, angle = 90,
                                                margin = margin(0,0,0,r=6)),
                    axis.title.x = element_text(size = 12, margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 10),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6), vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10), angle = 270)))



# cairo_pdf embeds fonts by default
save_file <- function(x, fn, ...) {
    cairo_pdf(paste0("analysis/output/", gsub(".pdf$", "", fn), ".pdf"), ...)
    print(x)
    dev.off()
}

no_x <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
}
no_y <- function(p) {
    p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
}
no_xy <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.title.y = element_blank(), axis.text.y = element_blank())
}
no_leg <- function(p) {
    p + theme(legend.position = "none")
}



#==========*
#========== Observed data (Combine 'time' and 'distance' into Figure 1) ----
#==========*

# time
time_p <- data_fit %>%
    ggplot(aes(year, y))+
    facet_wrap(~Taxon, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    # geom_point(aes(group = plot), size = 1.5, alpha = 0.5)+
    geom_line(data = data_fit %>%
                  group_by(Taxon, year) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge catch") %>%
                  bind_rows(data_fit %>%
                                group_by(Taxon, year) %>%
                                summarize(y = mean(y)) %>%
                                mutate(type  = "Mean activity-density")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge catch",
                                                  "Mean activity-density"))),
              aes(color = type), size = 0.75)+
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
    geom_jitter(aes(group = plot), size = 1, alpha = 0.5, width = 0.1, shape = 1)+
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
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick", "dodgerblue"))+
    scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    scale_y_continuous(NULL, limits = c(-4.5,4.5),
                       breaks = c(-3,0,3)) +
    NULL


time_dist_legend <- get_legend(
    time_p +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)


prow <- plot_grid(time_p %>% no_leg(),
                  dist_p  %>% no_leg(),
                  labels = c("A", "B"),
                  align = "vh")

fig1 <- plot_grid(time_dist_legend, prow, ncol = 1, rel_heights = c(0.1, 1))


save_file(fig1, "fig1", width = 6, height = 4)


#==========*
#========== Plot coefficients (Use `taxon-specicific-slopes` for Figure 2) ----
#==========*


small_labels <- theme(axis.text = element_text(size = 8, color = "black"),
                      axis.title.y = element_text(size = 10, angle = 90,
                                                  margin = margin(0,0,0,r=2)),
                      axis.title.x = element_text(size = 10, margin = margin(0,0,0,t=2)),
                      strip.text = element_text(size = 9))

# taxon-specific slopes
slope_p <- coef_sum$beta %>%
    ggplot()+
    geom_rect(data = coef_sum$alpha,
              aes(xmin = 0.5, xmax = 6.5, ymin = lo, ymax = hi, fill = coef),
              alpha = 0.5, linetype =  0)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(aes(tx, mi), size = 1.5)+
    geom_linerange(aes(tx, ymin = lo, ymax = hi))+
    scale_y_continuous(expression("Response" ~ ({}^{k} * beta[taxon])),
                       breaks = c(-0.5, 0, 0.5), limits = c(-0.6, 0.6))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       trans = "reverse", limits = c(6.5, 0.5)) +
    facet_wrap(~ coef, ncol = 1) +
    coord_flip() +
    scale_fill_manual(values = coef_palette, guide = FALSE) +
    small_labels


# ar (extra)
ar_p <- coef_sum$ar %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 1.5)+
    ggtitle(" ") +
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_y_continuous(expression("AR parameter" ~ (phi)),  limits = c(0, 1),
                       breaks = c(0, 0.5, 1))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles"),
                       trans = "reverse", limits = c(6.5, 0.5))+
    coord_flip() +
    small_labels


density_labels <- small_labels +
    theme(axis.text.x = element_text(size = 8, color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

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
                         levels = c("time", "dist", "midges"))) %>%
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
    density_labels +
    NULL


effect_mean_p <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    .[[1]] %>%
    as.matrix() %>%
    as_tibble() %>%
    gather() %>%
    filter(key != "V1") %>%  # intercept not necessary
    mutate(coef = factor(key, levels = c("V3","V2","V4"),
                         labels = c("time", "dist", "midges"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_rect(data  = coef_sum$alpha,
              aes(xmin = lo, xmax = hi, ymin = 0, ymax = 8, fill = coef), inherit.aes = F,
              alpha = 0.3, linetype =  0)+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual("", values = coef_palette)+
    scale_y_continuous("Density", limits = c(0, 8), breaks = 0:4 * 2)+
    scale_x_continuous(expression("Response mean (" * {}^{k} * alpha * ")")) +
    theme(legend.position = "none") +
    density_labels +
    NULL




# fig2 <-
plot_grid(slope_p,
          plot_grid(ar_p %>% no_y(),
                    plot_grid(),
                    effect_sigmas_p + theme(legend.position = c(0.7,0.5)),
                    effect_mean_p,
                    labels = c(LETTERS[2], ""), ncol = 1,
                    align = "h", axis = "t", rel_heights = c(1.35, 1, 1)),
          labels = c(LETTERS[1], "", ""), ncol = 2, rel_widths = c(1, 0.6))


# save_file(fig2, "fig2-A", width = 6, height = 5)




# fig2 <- plot_grid(NULL, slope_p + theme(plot.margin = margin(t=0,r=4,b=10,l=0)),
#                   NULL, plot_grid(ar_p +
#                                       theme(plot.margin = margin(t=0,r=4,b=10,l=0)),
#                                   NULL,
#                                   effect_sigmas_p +
#                                       theme(legend.position = c(0.7,0.5),
#                                             plot.margin = margin(t=0,r=4,b=10,l=0)),
#                                   NULL,
#                                   effect_mean_p +
#                                       theme(plot.margin = margin(t=0,r=4,b=10,l=0),
#                                             axis.title.y = element_blank(),
#                                             axis.text.y = element_blank()),
#                                   rel_widths = c(1.2, 0.15, 1.1, 0.15, 1),
#                                   labels = c("", "C", "", "D", ""),
#                                   nrow = 1),
#                   labels = c("A", "", "B", ""), rel_widths = c(0.05, 1),
#                   rel_heights = c(1.1, 1),
#                   ncol = 2)
#
#
# cairo_pdf("analysis/output/fig2-B.pdf", width = 6, height = 4)
# fig2
# dev.off()



#==========*
#========== PCA (Combine biplots for `midge effect`, `time effect`, and `distance` for Fig 3 ----
#==========*

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
overall_part <- as.matrix(var_part[,2:4]) %*% t(as.matrix(pred_pca$obs_exp[1,2:4]))
row.names(overall_part) <- var_part$var

# predictor vectors
pred_vec_fun <- function(pred_) {

    sym_ <- sym(paste0(pred_, "_z"))

    pred_pca$axes %>%
        filter(!!sym_ %in% range(!!sym_)) %>%
        group_by(!!sym_) %>%
        summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3)) %>%
        rename(val = !!sym_) %>%
        arrange(val) %>%
        mutate(var = pred_)  %>%
        mutate(PC1 = PC1 - PC1[1],
               PC2 = PC2 - PC2[1],
               PC3 = PC3 - PC3[1]) %>%
        filter(PC1 != 0 | PC2 != 0 | PC3 != 0) %>%
        select(var, starts_with("PC"))
}
pred_vec <- map_dfr(c("time", "dist", "midges"), pred_vec_fun)


# Which PC axes are most associated with each variable?
pc_vars <- var_part[,2:4] %>%
    as.matrix() %>%
    `*`(matrix(as.numeric(pred_pca$obs_exp[1,2:4]), 3, 3, byrow=TRUE)) %>%
    as.data.frame() %>%
    as_tibble() %>%
    set_names(paste0("PC", 1:3)) %>%
    mutate(var = var_part$var) %>%
    select(var, everything())



# biplot function
biplot_fn <- function(var_, pred_pca_ = pred_pca, pc_vars_ = pc_vars, .mult = 2) {

    var_ <- match.arg(var_, c("midges", "time", "dist"))

    var_ <- paste0(var_, "_z")

    # pc_axes <- pc_vars_ %>%
    #     filter(var == var_) %>%
    #     .[,2:4] %>%
    #     as.numeric() %>%
    #     rank() %>%
    #     {paste0("-PC", c(which(. == 3), which(. == 2)))}
    # Only doing PC 1 and 2 for main text
    pc_axes <- paste0("-PC", 1:2)

    pred_pca_$obs_rot %>%
        # Uncomment below if you want darker spots in front
        # I currently don't think we should do it
        ## arrange(desc(!!sym(var_))) %>%
        ggplot(aes_string(x = pc_axes[1], y = pc_axes[2]))+
        geom_point(aes_string(color = var_),
                   alpha = 0.25, size = 2)+
        ggtitle(gsub("_z$", "", var_)) +
        geom_segment(data = pred_vec %>% filter(var == gsub("_z$", "", var_)),
                     aes(x = 0, xend = -.mult*PC1, y = 0, yend = -.mult*PC2),
                     arrow = arrow(length = unit(6, "pt")),
                     size = 0.8, color = "black")+
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
    mutate(taxon = factor(taxon, levels = c("gnap","lyco","sheet","opil","cara","stap"),
                          labels = c("ground\nspiders","wolf\nspiders","sheet\nweavers",
                                     "harvestman","ground\nbeetles","rove\nbeetles"))) %>%
    arrange(taxon) %>%
    mutate(PC1_lab = -5*PC1 + c(1.2, 0.2, 1, 0.5, -0.8, 0.4),
           PC2_lab = -5*PC2 + c(-0.4, -0.8, -0.5, 1.7, 0.7, -0.9),
           ang = c(rep(0,3), 73.5, 0, 0)) %>%
    arrange(desc(taxon)) %>%
    ggplot()+
    geom_segment(aes(x = 0, xend = -5*PC1, y = 0, yend = -5*PC2, group = taxon,
                     color = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 1)+
    geom_text(aes(label = taxon, x = PC1_lab, y = PC2_lab,
                  color = taxon, angle = ang),
              size = 9 / 2.83465, lineheight = 0.75)+
    xlab("-PC1") +
    ylab("-PC2") +
    coord_equal(xlim = c(-4.1, 4.1), ylim = c(-4.1, 4.1)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(1:4, 6:5)],
                       guide = FALSE) +
    theme(plot.margin = margin(0,0,0,0))

# plot time effect
fig3b <- biplot_fn("time")
# plot distance effect
fig3c <- biplot_fn("dist")
# plot midge effect
fig3d <- biplot_fn("midges")


pca_legend <- get_legend(fig3b + theme(legend.title = element_blank()))



fig3 <- plot_grid(plot_grid(fig3a %>% no_x(),
                            fig3b %>% no_leg() %>% no_xy(),
                            fig3c %>% no_leg(),
                            fig3d %>% no_leg() %>% no_y(),
                            labels = LETTERS[1:4],
                            nrow = 2, align = "vh"),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))
# fig3


save_file(fig3, "fig3", width = 6, height = 5)





