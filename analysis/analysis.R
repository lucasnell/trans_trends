
#'
#' This file creates figures 1-3 and table 1.
#'




# ==================================================*
# Preliminaries ----
# ==================================================*

# load packages
library(TransTrendsPkg)
library(tidyverse)
library(cowplot)
library(egg)
library(grid)
library(viridisLite)




# ------------*
# Map abbreviated taxa names to full names, and order them properly for plots
# ------------*

taxa_order <- c(5:6, 4, 1, 3, 2)
taxa_lvls = c("gnap","lyco","sheet","opil","cara", "stap")[taxa_order]
taxa_labs = c("Ground spiders","Wolf spiders","Sheet weavers",
              "Harvestmen","Ground beetles", "Rove beetles")[taxa_order]



# ------------*
# Load and clean data from model fit
# ------------*

# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs))
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")

coef_sum$int_taxon <- coef_sum$int_taxon %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)),
           coef = "int")

# Order levels of coefficient factor:
make_coef_fct <- function(.coef) {
    factor(.coef %>% paste(),
           levels = c("time", "dist", "midges"),
           labels = c("time", "distance", "midges"))
}

coef_sum$beta <- coef_sum$beta %>%
    mutate(coef = make_coef_fct(coef),
           tx = as.numeric(factor(taxon, levels = taxa_lvls)))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = make_coef_fct(coef))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_tax_plot","int_tax_trans"))) %>%
    mutate(coef = make_coef_fct(coef))
coef_sum$ar <- coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)))



# ------------*
# Define theme, a colorblind-friendly palette, and
# helper functions/objects for plots
# ------------*
pca_palette <- viridis(3)


# set theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 12),
                    legend.background = element_blank(),
                    axis.text = element_text(size = 9, color = "black"),
                    axis.title.y = element_text(size = 12, angle = 90,
                                                margin = margin(0,0,0,r=6)),
                    axis.title.x = element_text(size = 12,
                                                margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 10),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6),
                                                vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10),
                                                angle = 270)))

# empty plot
EMPTY <- ggplot() + geom_blank() + theme_void()


# use cairo_pdf because it embeds fonts by default
save_file <- function(x, fn, ...) {
    cairo_pdf(paste0("analysis/output/", gsub(".pdf$", "", fn), ".pdf"), ...)
    print(x)
    dev.off()
}
# remove x axis from a plot
no_x <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
}
# remove y axis from a plot
no_y <- function(p) {
    p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
}
# remove x and y axes from a plot
no_xy <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.title.y = element_blank(), axis.text.y = element_blank())
}
# remove legend from a plot
no_leg <- function(p) {
    p + theme(legend.position = "none")
}



# ==================================================*
# ==================================================*

# Fig 1 - Observed data ----

# ==================================================*
# ==================================================*

# by time
time_p <- data_fit %>%
    ggplot(aes(year, y))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, year) %>%
                  summarize(y = mean(midges_z), .groups = "drop") %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, year) %>%
                                summarize(y = mean(y), .groups = "drop") %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016),
                       limits = c(2007,2018))+
    scale_y_continuous("Transformed abundance",
                       limits = c(-2, 3.1), breaks = c(-2,0,2))+
    NULL

# by distance
dist_p <- data_fit %>%
    ggplot(aes(distance, y))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_jitter(aes(group = plot), size = 1, alpha = 0.5,
                width = 0.1, shape = 1)+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, distance) %>%
                  summarize(y = mean(midges_z), .groups = "drop") %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, distance) %>%
                                summarize(y = mean(y), .groups = "drop") %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick", "dodgerblue"))+
    scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    scale_y_continuous(NULL, limits = c(-2, 3.1), breaks = c(-2,0,2)) +
    NULL

# separate legend
time_dist_legend <- get_legend(
    time_p +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)


prow <- plot_grid(time_p %>% no_leg(),
                  dist_p  %>% no_leg(),
                  labels = c("a", "b"),
                  align = "vh", label_fontface = "plain", label_size = 16)

fig1 <- plot_grid(time_dist_legend, prow, ncol = 1, rel_heights = c(0.1, 1))

# fig1

# save_file(fig1, "fig1", width = 6, height = 4)




# ==================================================*
# ==================================================*

# Fig 2 - Coefficients ----

# ==================================================*
# ==================================================*


#' Extract posterior densities for the fixed effects
#' and for the random effect SDs
slope_density_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int", "midges", "time", "distance"); .} %>%
    as_tibble() %>%
    select(-int) %>%  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 2048)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)
slope_sd_density_df <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int_tax","int_plot","int_trans","midges",
                      "time","distance"); .} %>%
    as_tibble() %>%
    select(all_of(c("midges", "time", "distance"))) %>%
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 2048)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)


#' Do the same as above, but restrict the densities to being between
#' the 68% uncertainty intervals.
#' These are used to create the darker shaded regions in fig 2.
slope_density_UI_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int", "midges", "time", "distance"); .} %>%
    as_tibble() %>%
    select(-int) %>%  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ci <- quantile(z[["value"]], c(0.16, 0.84))
        X <- density(z[["value"]], from = .ci[["16%"]], to = .ci[["84%"]])
        tibble(coef = z[["coef"]][1],
               x = c(X$x[1], X$x, tail(X$x, 1)),
               y = c(0, X$y, 0))
    }) %>%
    arrange(coef, x)
slope_sd_density_UI_df <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int_tax","int_plot","int_trans","midges",
                      "time","distance"); .} %>%
    as_tibble() %>%
    select(all_of(c("midges", "time", "distance"))) %>%
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ci <- quantile(z[["value"]], c(0.16, 0.84))
        X <- density(z[["value"]], from = .ci[["16%"]], to = .ci[["84%"]])
        tibble(coef = z[["coef"]][1],
               x = c(X$x[1], X$x, tail(X$x, 1)),
               y = c(0, X$y, 0))
    }) %>%
    arrange(coef, x)


# Creates subpanels for each predictor (fig 2 b--d)
slope_p_fun <- function(.coef) {

    .coef <- match.arg(tolower(.coef), levels(coef_sum$beta$coef))

    .ylab <- gsub("s$", "", .coef) %>%
        str_to_title() %>%
        paste("response")

    .plot <- coef_sum$beta %>%
        filter(coef == .coef) %>%
        ggplot()+
        # main density curves:
        geom_polygon(data = slope_density_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray80") +
        geom_polygon(data = slope_sd_density_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray80") +
        # uncertainty interval density curves:
        geom_polygon(data = slope_density_UI_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray60") +
        geom_polygon(data = slope_sd_density_UI_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray60") +
        geom_hline(yintercept = 0, color = "gray50")+
        geom_point(aes(tx, mi), size = 1.5, color = "black")+
        geom_linerange(aes(tx, ymin = lo, ymax = hi), color = "black")+
        scale_y_continuous(.ylab,  breaks = c(-0.5, 0, 0.5))+
        scale_x_continuous(NULL,
                           breaks = 1:6, labels = taxa_labs,
                           sec.axis =
                               sec_axis(~ .,
                                        breaks = 1:6, labels = taxa_labs)) +
        coord_cartesian(ylim = c(-0.68, 0.68), xlim = c(0.5, 7),
                        expand = FALSE) +
        theme(axis.text.y.left = element_text(size = 8,
                                              margin = margin(0,0,0,r=2)),
              axis.text.y.right = element_blank(),
              axis.ticks.y.right = element_blank(),
              axis.title.y.left = element_text(size = 10,
                                               margin = margin(0,0,0,r=6)),
              axis.title.y.right = element_text(size = 10,
                                                margin = margin(0,0,0,l=6)),
              axis.title.x = element_blank(),
              axis.text.x.top = element_text(size = 8, angle = 45,
                                             vjust = 0.1, hjust = 0.1,
                                             margin = margin(0,0,0,b=2)),
              axis.text.x.bottom = element_text(size = 8, angle = -45,
                                                vjust = 0, hjust = 0.1,
                                                margin = margin(0,0,0,t=6)),
              plot.margin = margin(t=8, r=9, b=8, l=8)) +
        NULL

    return(.plot)
}


slope_p <- lapply(levels(coef_sum$beta$coef), slope_p_fun)

#' Add labels for density curves to Fig 2b
slope_p[[1]] <- slope_p[[1]] +
    geom_text(data = tibble(x = c(0.7,    4.5),
                            y = c(-0.28, 0.45),
                            lab = c("among taxa mean", "among taxa SD")),
              aes(x, y, label = lab), size = 8 / 2.83465,
              hjust = 0, vjust = 0.5, color = "gray50")


# AR coefficient plot:
ar_p <- coef_sum$ar %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 1.5)+
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_x_continuous(NULL,
                       breaks = 1:6, labels = taxa_labs,
                       sec.axis =
                           sec_axis(~ .,
                                    breaks = 1:6, labels = taxa_labs)) +
    scale_y_continuous("AR parameter", breaks = c(0, 0.5, 1)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0.5, 7),
                    expand = FALSE) +
    theme(axis.text.y = element_text(size = 8, margin = margin(0,0,0,r=2)),
          axis.title.y = element_text(size = 10, margin = margin(0,0,0,r=6)),
          axis.title.x = element_blank(),
          axis.text.x.top = element_text(size = 8, angle = 45,
                                         vjust = 0.1, hjust = 0.1,
                                         margin = margin(0,0,0,b=2)),
          axis.text.x.bottom = element_blank(),
          plot.margin = margin(t=8, r=8, b=8, l=8)) +
    NULL





#' Helper plot functions.
#' No taxa names, top:
ntt <- function(.x) .x + theme(axis.text.x.top = element_blank())
#' No taxa names, bottom:
ntb <- function(.x) .x <- .x + theme(axis.text.x.bottom = element_blank())


fig2 <- ggarrange(ar_p,
                  slope_p[[1]] %>% ntt() %>% ntb(),
                  slope_p[[2]] %>% ntt() %>% ntb(),
                  slope_p[[3]] %>% ntt(),
                  labels = letters[1:4],
                  label.args = list(gp = gpar(font = 1, fontsize = 16),
                                    x = unit(0, "line"), hjust = 0),
                  ncol = 1, draw = FALSE)



# fig2


# save_file(fig2, "fig2", width = 3, height = 8)





# ==================================================*
# ==================================================*

# Fig 3 - PCA ----

# ==================================================*
# ==================================================*

# This file contains a bunch of functions for the PCA:
source("analysis/pca_funs.R")

# PCA on predicted values
pred_pca <- pred_pca_fn(data_fit, coef_sum$beta, coef_sum$int_taxon)


#
# This function is to switch PC2 and PC.
#
# NOTE: We're moving PC3 to PC2 bc it explains more of the observed variance!
# ----------------------*
switch_pcs <- function(.x) {
    if (inherits(.x, "data.frame")) {
        if (!"PC2" %in% colnames(.x)) return(.x)
        z <- .x[["PC2"]]
        .x[["PC2"]] <- .x[["PC3"]]
        .x[["PC3"]] <- z
    } else if (inherits(.x, "prcomp")) {
        z <- .x$sdev[2]
        .x$sdev[2] <- .x$sdev[3]
        .x$sdev[3] <- z
        z <- .x$rotation[,"PC2"]
        .x$rotation[,"PC2"] <- .x$rotation[,"PC3"]
        .x$rotation[,"PC3"] <- z
        z <- .x$x[,"PC2"]
        .x$x[,"PC2"] <- .x$x[,"PC3"]
        .x$x[,"PC3"] <- z
    } else stop("\nUnknown type in `switch_pcs`")
    return(.x)
}

pred_pca <- map(pred_pca, switch_pcs)
# Also adjust cumulative proportion for `obs_exp`:
pred_pca$obs_exp[pred_pca$obs_exp$type == "cumulative",-1] <-
    pred_pca$obs_exp[pred_pca$obs_exp$type == "individual",-1] %>%
    unlist() %>%
    cumsum() %>%
    as.list()



# taxon vectors
pred_pca$taxon_vec %>%
    arrange(-abs(PC1))

# variance explained in predicted values (first three axes must
# explain everything)
summary(pred_pca$pca)


# variance explained in observed values (For Results)
pred_pca$obs_exp



# predictor vectors
pred_vec <- map_dfr(c("time", "distance", "midges"),
                    function(pred_) {
                        if (pred_ == "distance") pred_ <- "dist"
                        sym_ <- sym(paste0(pred_, "_z"))
                        pred_pca$axes %>%
                            filter(!!sym_ %in% range(!!sym_)) %>%
                            group_by(!!sym_) %>%
                            summarize(PC1 = mean(PC1),
                                      PC2 = mean(PC2),
                                      PC3 = mean(PC3),
                                      .groups = "drop") %>%
                            rename(val = !!sym_) %>%
                            arrange(val) %>%
                            mutate(var = pred_)  %>%
                            mutate(PC1 = PC1 - PC1[1],
                                   PC2 = PC2 - PC2[1],
                                   PC3 = PC3 - PC3[1]) %>%
                            filter(PC1 != 0 | PC2 != 0 | PC3 != 0) %>%
                            select(var, starts_with("PC"))
                    })



pc_axis_lim <- 4.1
#' We're scaling some of the projections. This is to keep that consistent.
pc_mults <- list(pred = c(-2, 2, 2), taxon = c(-5, 5, 5))


pca_theme <- theme(plot.margin = margin(t=8,r=8,b=0,l=8),
                   axis.text.y = element_text(size = 8,
                                               margin = margin(0,0,0,r=2)),
                   axis.title.y = element_text(size = 10,
                                               margin = margin(0,0,0,0)),
                   axis.text.x = element_text(size = 8,
                                              margin = margin(0,0,0,t=2)),
                   axis.title.x = element_text(size = 10,
                                               margin = margin(0,0,0,b=8)))


# ----------------*
# Plots of taxon response vectors
# ----------------*


taxon_pca_fun <- function(.xPC, .yPC,
                          .nudge_x = NULL,
                          .nudge_y = NULL,
                          .segment_df = NULL) {

    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    .PCs <- c(.xPC, .yPC)

    .dd <- pred_pca$taxon_vec %>%
        mutate(taxon = factor(taxon, levels = taxa_lvls %>% rev(),
                              labels = c("ground\nspiders","wolf\nspiders",
                                         "sheet\nweavers", "harvestmen",
                                         "ground\nbeetles","rove\nbeetles") %>%
                                  .[rev(taxa_order)])) %>%
        arrange(taxon) %>%
        mutate(PC1 = pc_mults$taxon[1] * PC1,
               PC2 = pc_mults$taxon[2] * PC2,
               PC3 = pc_mults$taxon[3] * PC3)

    if (!is.null(.nudge_x) && !is.null(.nudge_y)) {
        .labs <- paste0("PC", .PCs, "_lab")
        .rto <- rev(taxa_order)
        .dd[[.labs[1]]] <- .dd[[paste0("PC", .PCs[1])]] + .nudge_x[.rto]
        .dd[[.labs[2]]] <- .dd[[paste0("PC", .PCs[2])]] + .nudge_y[.rto]
    }

    .dd <- .dd %>%
        arrange(desc(taxon))

    .p <- .dd %>%
        ggplot()+
        geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
        geom_segment(aes_string(x = 0, xend = paste0("PC", .PCs[1]),
                                y = 0, yend = paste0("PC", .PCs[2]),
                                group = "taxon", color = "taxon"),
                     arrow = arrow(length = unit(6, "pt")),
                     size = 1) +
        scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
        scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        scale_color_manual(values = RColorBrewer::brewer.pal(9, "RdYlBu") %>%
                               .[c(1:3, 7:9)],
                           guide = FALSE) +
        pca_theme

    if (!is.null(.nudge_x) && !is.null(.nudge_y)) {
        .p <- .p +
            geom_text(aes_string(label = "taxon", x = .labs[1], y = .labs[2],
                                 color = "taxon"),
                      size = 7 / 2.83465, lineheight = 0.75)
    }

    if (!is.null(.segment_df)) {
        .p <- .p +
            geom_segment(data = .segment_df,
                         aes(x = x, y = y, xend = xend, yend = yend),
                         color = "black", size = 0.25)
    }

    return(.p)
}





taxon_pca_p <- list(
    taxon_pca_fun(1, 2,
                  .nudge_x = c(-0.4,  0.5, -1.5,  2.7, -1.0,  0.5),
                  .nudge_y = c( 1.0, -1.4, -3.0, -2.6,  0.8, -1.3),
                  .segment_df = tibble(x = c(1.25,  -1.3),
                                       y = c(-0.6,   -1.3),
                                       xend = c(-0.4, -0.4),
                                       yend = c(1.2,   0.6))),
    taxon_pca_fun(1, 3),
    taxon_pca_fun(2, 3))



# ----------------*
# Plots of model predictors and observed data
# ----------------*

predictor_pca_fun <- function(.xPC, .yPC, .label_df = NULL, ...) {

    # .xPC = 1; .yPC = 2; .label_df = NULL

    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    .PCs <- c(.xPC, .yPC)

    .obs_rot <- pred_pca$obs_rot %>%
        mutate(PC1 = sign(pc_mults$pred[1]) * PC1,
               PC2 = sign(pc_mults$pred[2]) * PC2,
               PC3 = sign(pc_mults$pred[3]) * PC3)

    vars <- c("time", "distance", "midges")

    .pred_vec <- pred_vec %>%
        mutate(PC1 = pc_mults$pred[1] * PC1,
               PC2 = pc_mults$pred[2] * PC2,
               PC3 = pc_mults$pred[3] * PC3) %>%
        mutate(var = factor(var, levels = gsub("^distance$", "dist", vars),
                            labels = vars))


    pc_axes <- paste0("PC", .PCs)

    .p <- .obs_rot %>%
        ggplot(aes_string(x = pc_axes[1], y = pc_axes[2]))+
        geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
        geom_point(alpha = 0.25, size = 1, shape = 1,
                   # color = "gray60") +
                   color = viridis::viridis(1, begin = 0.9, end = 0.9)) +
        geom_segment(data = .pred_vec,
                     aes_string(x = 0, xend = pc_axes[1],
                                y = 0, yend = pc_axes[2],
                                color = "var"),
                     arrow = arrow(length = unit(6, "pt")),
                     size = 1) +
        scale_x_continuous(pc_axes[1], breaks = 3*-1:1) +
        scale_y_continuous(pc_axes[2], breaks = 3*-1:1) +
        scale_color_viridis_d(end = 0.75, guide = FALSE) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        pca_theme +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank())

    if (!is.null(.label_df)) {
        stopifnot(is.data.frame(.label_df) && nrow(.label_df) == 3)
        stopifnot(identical(c("lab", "x", "y"), sort(colnames(.label_df))))
        stopifnot(all(c("time", "midges") %in% .label_df$lab))
        stopifnot(any(grepl("^dist", .label_df$lab)))

        .p <- .p +
            geom_text(data = .label_df %>%
                          mutate(var = ifelse(grepl("^dist", lab),
                                              "distance", lab) %>%
                                     factor(levels = vars)),
                      aes(x = x, y = y, label = lab, color = var),
                      size = 8 / 2.83465,
                      ...)
    }

    return(.p)

}



pred_pca_p <- list(
    predictor_pca_fun(1, 2,
                      .label_df = tibble(x = c(-0.8, -1, pc_axis_lim),
                                         y = c(2.5, -1, 1),
                                         lab = c("midges", "time", "dist.")),
                      hjust = 1, vjust = 0.5),
    predictor_pca_fun(1, 3),
    predictor_pca_fun(2, 3))





# --------------*
# Variance partitioning
# --------------*


# variance partition of PC axes by predictors
var_part <- lapply(c("PC1","PC2","PC3"), function(x) {
    y = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(var = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(y, data = pred_pca$axes))[,2] %>%
               {.[1:3]/sum(.[1:3])})
    }) %>%
    bind_rows() %>%
    spread(pc, cont) %>%
    mutate(var = gsub("_z$", "", var))


# variance accounted for by predictors and PCs:
var_df <- var_part %>%
    mutate(PC1 = PC1 * pred_pca$obs_exp$PC1[1],
           PC2 = PC2 * pred_pca$obs_exp$PC2[1],
           PC3 = PC3 * pred_pca$obs_exp$PC3[1],
           var = make_coef_fct(var)) %>%
    pivot_longer(starts_with("PC"), names_to = "pc") %>%
    mutate(pc = gsub("PC", "", pc) %>% as.integer()) %>%
    arrange(pc, var) %>%
    mutate(cumvalue = cumsum(value),
           lag_cumvalue = lag(cumvalue, default = 0))




var_part_p <- var_df %>%
    mutate(y = pc - 1,
           ymax = 0 - y * 0.5,
           ymin = ymax - 1) %>%
    ggplot() +
    geom_rect(aes(xmin = lag_cumvalue, xmax = cumvalue,
                  ymax = ymax, ymin = ymin, fill = var)) +
    geom_text(data = var_df %>%
                  group_by(pc) %>%
                  summarize(value = (min(lag_cumvalue) + max(cumvalue)) / 2,
                            .groups = "drop") %>%
                  add_row(pc = 4, value = (1 + sum(var_df$value)) / 2) %>%
                  mutate(y = 0 - (pc-1) * 0.5,
                         lab = ifelse(pc < 4, paste0("PC", pc), "")),
              aes(value, y, label = lab),
              nudge_y = 0.3, size = 8 / 2.83465) +
    geom_point(data = tibble(var = var_df$var %>% unique() %>% sort(),
                             value = 0.79,
                             y = 0.25 - 0:2 * 0.5),
              aes(value, y, color = var),
              size = 2, shape = 15) +
    geom_text(data = tibble(var = var_df$var %>% unique() %>% sort(),
                            value = 0.82,
                            y = 0.25 - 0:2 * 0.5),
              aes(value, y, label = var, color = var),
              size = 8 / 2.83465, hjust = 0, vjust = 0.5) +
    scale_color_viridis_d(end = 0.75, aesthetics = c("fill", "color"),
                          guide = FALSE) +
    scale_x_continuous("Proportion of variance", breaks = seq(0, 1, 0.2)) +
    coord_cartesian(xlim = c(0, 1),
                    ylim = c(-2.1, 0.7),
                    expand = FALSE) +
    pca_theme +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(margin = margin(0,0,0,b=6)),
          # panel.border = element_blank(),
          plot.margin = margin(0,0,r=8,l=27)) +
    NULL





# --------------*
# Combine them all into fig 3
# --------------*
fig3 <- plot_grid(ggarrange(plots = c(taxon_pca_p, pred_pca_p),
                            nrow = 3, labels = letters[c(1,3,5,2,4,6)],
                            draw = FALSE, byrow = FALSE,
                            label.args = list(gp = gpar(font = 1, fontsize = 16),
                                              x = unit(0,"line"), hjust = 0)),
                  var_part_p,
                  labels = c("", "g"),
                  label_fontface = "plain", label_size = 16,
                  label_x = 0, label_y = 1.3,
                  ncol = 1, rel_heights = c(5.5, 1))

# fig3


save_file(fig3, "fig3", width = 3, height = 6.5)





# ==================================================*
# ==================================================*

# Fig S1 - colored by predictor ----

# ==================================================*
# ==================================================*


pred_color_pca_fun <- function(.xPC, .yPC) {

    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    .PCs <- c(.xPC, .yPC)

    vars <- c("time", "distance", "midges")

    var_z_ <- paste0(gsub("^distance$", "dist", vars), "_z")

    pc_axes <- paste0("PC", .PCs)

    .obs_rot <- pred_pca$obs_rot %>%
        mutate(PC1 = sign(pc_mults$pred[1]) * PC1,
               PC2 = sign(pc_mults$pred[2]) * PC2,
               PC3 = sign(pc_mults$pred[3]) * PC3) %>%
        select(all_of(c(var_z_, pc_axes))) %>%
        pivot_longer(all_of(var_z_), names_to = "var", values_to = "value") %>%
        mutate(var = str_replace(var, "_z$", "") %>%
                   factor(levels = gsub("^distance$", "dist", vars),
                          labels = vars))

    .pred_vec <- pred_vec %>%
        mutate(PC1 = pc_mults$pred[1] * PC1,
               PC2 = pc_mults$pred[2] * PC2,
               PC3 = pc_mults$pred[3] * PC3) %>%
        mutate(var = factor(var, levels = gsub("^distance$", "dist", vars),
                            labels = vars))


    .obs_rot %>%
        ggplot(aes_string(x = pc_axes[1], y = pc_axes[2]))+
        geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
        geom_point(aes(color = value), alpha = 0.25, size = 1)+
        geom_segment(data = .pred_vec,
                     aes_string(x = 0, xend = pc_axes[1],
                                y = 0, yend = pc_axes[2]),
                     arrow = arrow(length = unit(6, "pt")),
                     size = 0.8, color = "black") +
        facet_grid(~ var) +
        scale_colour_gradient2("Predictor\nvalue", low =  pca_palette[1],
                               mid = pca_palette[2],
                               high = pca_palette[3],
                               midpoint = -0.5, limits = c(-3,2),
                               breaks = c(-3,-0.5,2))+
        scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
        scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        pca_theme +
        theme(plot.margin = margin(0,0,0,t=8)) +
        guides(color = guide_colorbar(
            title.hjust = 0,
            title.theme = element_text(size = 10,
                                       margin = margin(0,4,4,0),
                                       lineheight = 0.75),
            label.theme = element_text(size = 8),
            barwidth = unit(0.15, "inches"),
            barheight = unit(1, "inches")))
}


figS1 <- ggarrange(pred_color_pca_fun(1, 2) %>% no_leg(),
                   pred_color_pca_fun(1, 3),
                   pred_color_pca_fun(2, 3) %>% no_leg(),
                   labels = letters[1:3],
                   draw = FALSE,
                   label.args = list(gp = gpar(font = 1, fontsize = 16),
                                     x = unit(0,"line"), hjust = 0))



# save_file(figS1, "figS1", width = 6, height = 7)





