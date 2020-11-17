#==========*
#========== Preliminaries ----
#==========*

# load packages
library(armmr)
library(tidyverse)
library(cowplot)
options(mc.cores = parallel::detectCores()-2)

# In RStudio on macOS, this makes a separate plotting window that's independent of
# the Plots pane
if (Sys.info()[["sysname"]] == "Darwin" && .Platform$GUI == "RStudio") {
    options("device" = "quartz")
    grDevices::graphics.off()
}


options(dplyr.summarise.inform = FALSE)



taxa_order <- c(5:6, 4, 1, 3, 2)
taxa_lvls = c("gnap","lyco","sheet","opil","cara", "stap")[taxa_order]
taxa_labs = c("Ground spiders","Wolf spiders","Sheet weavers",
              "Harvestman","Ground beetles", "Rove beetles")[taxa_order]

# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs))
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")



coef_sum$int_taxon <- coef_sum$int_taxon %>%
    ungroup() %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)),
           coef = "int")

coef_sum$beta <- coef_sum$beta %>%
    ungroup() %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")),
           tx = as.numeric(factor(taxon, levels = taxa_lvls)))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_tax_plot","int_tax_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")))
coef_sum$ar <- coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)))



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
                    legend.background = element_blank(),
                    axis.text = element_text(size = 9, color = "black"),
                    axis.title.y = element_text(size = 12, angle = 90,
                                                margin = margin(0,0,0,r=6)),
                    axis.title.x = element_text(size = 12, margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 10),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6), vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10), angle = 270)))

# empty plot
EMPTY <- ggplot() + geom_blank() + theme_void()


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
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    # geom_point(aes(group = plot), size = 1.5, alpha = 0.5)+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, year) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, year) %>%
                                summarize(y = mean(y)) %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016), limits = c(2007,2018))+
    scale_y_continuous("Transformed abundance",
                       limits = c(-2, 3.1), breaks = c(-2,0,2))+
    # theme(legend.position = "none") +
    NULL

# distance
dist_p <- data_fit %>%
    ggplot(aes(distance, y))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_jitter(aes(group = plot), size = 1, alpha = 0.5, width = 0.1, shape = 1)+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, distance) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, distance) %>%
                                summarize(y = mean(y)) %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick", "dodgerblue"))+
    scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    scale_y_continuous(NULL, limits = c(-2, 3.1), breaks = c(-2,0,2)) +
    NULL


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

fig1

# save_file(fig1, "fig1", width = 6, height = 4)




#==========*
#========== Plot coefficients (Use `taxon-specific-slopes` for Figure 2) ----
#==========*



slope_density_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    .[[1]] %>%
    as.matrix() %>%
    {colnames(.) <- paste0("...", 1:ncol(.)); .} %>%
    as_tibble() %>%
    gather() %>%
    filter(key != "...1") %>%  # intercept not necessary
    mutate(coef = factor(key, levels = c("...3","...4","...2"),
                         labels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 1024)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)


slope_sd_density_df <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    .[[1]] %>%
    as.matrix() %>%
    {colnames(.) <- paste0("...", 1:ncol(.)); .} %>%
    as_tibble() %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("...1","...2","...3","...4","...5","...6"),
                         labels = c("int_tax","int_plot","int_trans","midges",
                                    "time","distance"))) %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 1024)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)


# slope_density_CI_df <- slope_density_df %>%
#     split(.$coef) %>%
#     map_dfr(function(.x) {
#         .lo <- filter(coef_sum$alpha, coef == .x$coef[1])[["lo"]]
#         .hi <- filter(coef_sum$alpha, coef == .x$coef[1])[["hi"]]
#         mutate(.x, lo = abs(x - .lo), hi = abs(x - .hi)) %>%
#             filter(lo == min(lo) | hi == min(hi) |
#                        (x > .lo & x < .hi)) %>%
#             select(-lo, -hi)
#     })


slope_density_CI_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    .[[1]] %>%
    as.matrix() %>%
    {colnames(.) <- paste0("...", 1:ncol(.)); .} %>%
    as_tibble() %>%
    gather() %>%
    filter(key != "...1") %>%  # intercept not necessary
    mutate(coef = factor(key, levels = c("...3","...4","...2"),
                         labels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .lo <- filter(coef_sum$alpha, coef == z$coef[1])[["lo"]]
        .hi <- filter(coef_sum$alpha, coef == z$coef[1])[["hi"]]
        X <- density(z[["value"]], from = .lo, to = .hi)
        tibble(coef = z[["coef"]][1],
               x = c(X$x[1], X$x, tail(X$x, 1)),
               y = c(0, X$y, 0))
    }) %>%
    arrange(coef, x)



# slope_density_df %>%
#     ggplot(aes(x, y, color = coef)) +
#     geom_path() +
#     geom_polygon(data = slope_density_CI_df, aes(fill = coef),
#                  alpha = 0.5, color = NA)

annotation_custom2 <- function (grob,
                                xmin = -Inf, xmax = Inf,
                                ymin = -Inf, ymax = Inf,
                                data) {

    layer(data = data, stat = StatIdentity, position = PositionIdentity,
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob,
                                            xmin = xmin, xmax = xmax,
                                            ymin = ymin, ymax = ymax))

}

coef_p_labs <- function(.x) {
    mutate(.x,
           coef = factor(paste(coef), levels = c("time", "distance", "midges"),
                         labels = c("Time response", "Distance response",
                                    "Midge response")))
}


insets <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    .[[1]] %>%
    as.matrix() %>%
    {colnames(.) <- paste0("...", 1:ncol(.)); .} %>%
    as_tibble() %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("...1","...2","...3","...4","...5","...6"),
                         labels = c("int_tax","int_plot","int_trans","midges",
                                    "time","distance"))) %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map(function(.x) {
        .coef <- .x$coef[1]
        .color <- coef_palette[as.integer(.coef)]
        if (.coef == "midges") {
            .x_axis <- element_text(size = 8, margin = margin(0,0,0,0))
            .plot.margin <- margin(0,0,0,0)
        } else {
            .x_axis <- element_blank()
            .plot.margin <- margin(0,0,0,b=8)
        }
        p <- ggplot(.x, aes(value, fill = coef))+
            geom_vline(xintercept = 0, size = 0.5, color = "black")+
            geom_hline(yintercept = 0, size = 0.5, color = "black")+
            geom_density(linetype = 0, alpha = 0.7, fill = .color)+
            scale_y_continuous("Posterior density", limits = c(0, 10))+
            scale_x_continuous("SD among taxa", breaks = c(0, 0.5, 1)) +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = 6),
                  axis.ticks.x = element_line(size = 0.25),
                  axis.title.x = .x_axis,
                  plot.margin = .plot.margin,
                  panel.border = element_blank()) +
            coord_cartesian(xlim = c(0, 1.0), ylim = c(0, 10), expand = FALSE) +
            NULL
        annotation_custom2(
            grob = ggplotGrob(p),
            data = data.frame(coef = .coef) %>% coef_p_labs(),
            ymin = 0.25, ymax=0.6, xmin=0, xmax=-3)
    })



# taxon-specific slopes
slope_p <- coef_sum$beta %>%
    coef_p_labs() %>%
    ggplot()+
    facet_wrap(~coef, ncol = 1, scales = "free_x", strip.position = "bottom")+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_polygon(data = slope_density_df %>%
                  mutate(y = 6.75 - y / max(y) * 6.75 - 0.05) %>%
                     coef_p_labs(),
                aes(x = y, y = x, fill = coef), alpha = 0.25) +
    geom_polygon(data = slope_density_CI_df %>%
                  mutate(y = 6.75 - y / max(y) * 6.75 - 0.05) %>%
                     coef_p_labs(),
                aes(x = y, y = x, fill = coef), alpha = 0.3) +
    geom_point(aes(tx, mi, color = coef), size = 1.5)+
    geom_linerange(aes(tx, ymin = lo, ymax = hi, color = coef))+
    scale_y_continuous("Response",
                       breaks = c(-0.5, 0, 0.5))+
    scale_x_continuous(breaks = 1:6, labels = taxa_labs,
                       trans = "reverse",
                       sec.axis = sec_axis(~ .,
                                           breaks = c(0, mean(c(6.75, -0.25))),
                                           labels = c("", "Posterior density"))) +
    scale_fill_manual(NULL, values = coef_palette, guide = FALSE)+
    scale_color_manual(NULL, values = coef_palette, guide = FALSE)+
    coord_flip(ylim = c(-0.68, 0.68), xlim = c(6.75, -0.25), expand = FALSE) +
    theme(plot.margin = margin(0,0,b=6.5,t=4),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 12, margin = margin(0,0,0,b=3)),
          panel.spacing = unit(1, "lines"),
          # strip.text.x = element_blank(),
          # strip.text.x = element_text(size = 12, margin = margin(0,0,t=0,b=3)),
          axis.title.x = element_blank(),
          axis.title.y.left = element_blank(),
          axis.title.y.right = element_blank(),
          # axis.text.y.right = element_blank(),
          axis.text.y.right = element_text(size = 10, angle = -90,
                                           hjust = 0.5, vjust = 0.5),
          axis.ticks.y.right = element_blank()) +
    insets +
    NULL





int_ar_p <- coef_sum$int_taxon %>%
    # Un-scale values:
    mutate(coef = "Mean trans. abundance") %>%
    bind_rows(mutate(coef_sum$ar, coef = "AR parameter")) %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_blank(data = tibble(coef = "AR parameter", mi = 1, tx = 1)) +
    geom_point(size = 1.5)+
    facet_wrap(~ coef, ncol = 1, scales = "free", strip.position = "bottom") +
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_x_continuous(NULL, breaks = 1:6, labels = taxa_labs,
                       trans = "reverse", limits = c(6.5, 0.5)) +
    coord_flip() +
    theme(axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.spacing = unit(1, "lines"),
          strip.text.x = element_text(size = 12, margin = margin(0,0,0,b=3)),
          plot.margin = margin(t=4, r=4, b=0, l=0)) +
    NULL




# obs and process error
error_p <- fit$stan %>%
    rstan::extract(pars = c("sig_obs", "sig_proc")) %>%
    map(as.matrix) %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("sig_obs", "sig_proc"); .} %>%
    as_tibble() %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("sig_obs", "sig_proc"),
                         labels = c("sigma[obs]", "sigma[proc]"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_density(linetype = 0, alpha = 0.7)+
    geom_text(data = tibble(coef = factor(c("sigma[proc]", "sigma[obs]")),
                            value = c(0.05, 0.45),
                            y = c(6, 10)),
              aes(y = y, label = coef, color = coef), hjust = 0, parse = TRUE) +
    scale_fill_manual(NULL, values = c("firebrick", "dodgerblue"))+
    scale_y_continuous("Posterior density")+
    scale_x_continuous("Error SD") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 12, margin = margin(0,0,t=6,b=6)),
          axis.title.y = element_text(size = 10, margin = margin(0,0,0,0)),
          plot.margin = margin(t=14, r=4, b=1, l=0),
          legend.position = "none") +
    NULL





fig2 <- plot_grid(plot_grid(EMPTY, EMPTY, EMPTY,
                            labels = letters[1:3],
                            rel_heights = c(1.845, 1.845, 2),
                            ncol = 1, label_x = 0,
                            label_y = c(1, 0.965, 0.93),
                            label_fontface = "plain", label_size = 16),
                  slope_p,
                  plot_grid(EMPTY, EMPTY, EMPTY,
                            labels = letters[4:6],
                            rel_heights = c(1.845, 1.845, 2),
                            ncol = 1, label_x = 0.75, hjust = 1,
                            label_y = c(1, 0.965, 0.93),
                            label_fontface = "plain", label_size = 16),
                  plot_grid(int_ar_p,
                            error_p,
                            rel_heights = c(1.845, 1),
                            ncol = 1),
                  rel_widths = c(0.08, 1, 0.12, 0.6), nrow = 1)


fig2


# save_file(fig2, "fig2", width = 6, height = 6)





#==========*
#========== PCA (Combine biplots for `midge effect`, `time effect`, and `distance` for Fig 3 ----
#==========*

source("analysis/pca_funs.R")

# pca on predicted values
pred_pca <- pred_pca_fn(data_fit, coef_sum$beta, coef_sum$int_taxon)


#
# This function is to switch PC2 and PC3 on the fly.
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



# taxon vectors
pred_pca$taxon_vec %>%
    arrange(-abs(PC1))

# variance explained in predicted values (first three axes must explain everything)
summary(pred_pca$pca)


# variance explained in observed values (For Results)
pred_pca$obs_exp



# predictor vectors
pred_vec_fun <- function(pred_) {

    if (pred_ == "distance") pred_ <- "dist"

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
pred_vec <- map_dfr(c("time", "distance", "midges"), pred_vec_fun)





pc_axis_lim <- 4.1




# biplot function
# (Only doing PC 1 and 2 for main text)
biplot_fn <- function(var_, .mult = c(-2, 2), .PCs = 1:2) {

    # var_ = "time"
    # .mult = c(-2, 2)
    # .PCs = 2:3

    if (length(.mult) == 1) .mult <- rep(.mult, 2)
    .signs <- ifelse(.mult < 0, "-", "")
    .sign_nums <- as.numeric(paste0(.signs, "1"))

    var_ <- match.arg(var_, c("midges", "time", "distance"))

    var_z_ <- paste0(gsub("^distance$", "dist", var_), "_z")

    pc_axes <- paste0(.signs, "PC", .PCs)

    pred_pca$obs_rot %>%
        # Uncomment below if you want darker spots in front
        # I currently don't think we should do it
        ## arrange(desc(!!sym(var_z_))) %>%
        ggplot(aes_string(x = pc_axes[1], y = pc_axes[2]))+
        geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
        geom_point(aes_string(color = var_z_),
                   alpha = 0.25, size = 2)+
        geom_segment(data = pred_vec %>% filter(var == gsub("_z$", "", var_z_)),
                     aes_string(x = 0,
                                xend = paste0(".mult[1]*PC", .PCs[1]),
                                y = 0,
                                yend = paste0(".mult[2]*PC", .PCs[2])),
                     arrow = arrow(length = unit(6, "pt")),
                     size = 0.8, color = "black") +
        geom_text(data = tibble(A = .sign_nums[1] * (pc_axis_lim + 0.2),
                                B = .sign_nums[2] * (pc_axis_lim + 0.2)) %>%
                      set_names(paste0("PC", .PCs)),
                  label = var_, hjust = 1, vjust = 1, size = 12 / 2.83465) +
        scale_colour_gradient2("Predictor\nvalue", low =  pca_palette[1],
                               mid = pca_palette[2],
                               high = pca_palette[3],
                               midpoint = -0.5, limits = c(-3,2),
                               breaks = c(-3,-0.5,2))+
        scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
        scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        theme(plot.title = element_text(size = 11, hjust = 0.5, vjust = 0,
                                        face = "plain"),
              plot.margin = margin(0,0,0,0))
}





fig3a <- pred_pca$taxon_vec %>%
    mutate(taxon = factor(taxon, levels = taxa_lvls %>% rev(),
                          labels = c("ground\nspiders","wolf\nspiders","sheet\nweavers",
                                     "harvestman","ground\nbeetles","rove\nbeetles")[rev(taxa_order)])) %>%
    arrange(taxon) %>%
    mutate(PC1 = -5*PC1,
           PC2 = 5*PC2,
           PC1_lab = PC1 + c(0.0, 0.8, -0.8, 3.3, -1.1, 0.4)[rev(taxa_order)],
           PC2_lab = PC2 + c(0.7, -0.6, -2.3, -0.7, 0.4, -0.9)[rev(taxa_order)]) %>%
    arrange(desc(taxon)) %>%
    ggplot()+
    geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
    geom_segment(aes(x = 0, xend = PC1, y = 0, yend = PC2, group = taxon,
                     color = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 1)+
    geom_segment(data = tibble(PC1 = c(1.25,  -0.8),
                               PC2 = c(0.8,   -0.8),
                               PC1e = c(-0.3, -0.4),
                               PC2e = c(1,    0.6)),
                 aes(x = PC1, y = PC2, xend = PC1e, yend = PC2e),
                 color = "black", size = 0.25) +
    geom_text(aes(label = taxon, x = PC1_lab, y = PC2_lab,
                  color = taxon),
              size = 9 / 2.83465, lineheight = 0.75)+
    scale_x_continuous("PC1", breaks = 3*-1:1) +
    scale_y_continuous("PC2", breaks = 3*-1:1) +
    coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                ylim = c(-pc_axis_lim, pc_axis_lim)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(1,4,3,5,2,6)],
                       guide = FALSE) +
    theme(plot.margin = margin(0,0,0,0))


# plot time effect
fig3b <- biplot_fn("time")
# plot distance effect
fig3c <- biplot_fn("distance")
# plot midge effect
fig3d <- biplot_fn("midges")


pca_legend <- get_legend(fig3b + theme(legend.title.align = 0,
                                       legend.title = element_text(size = 11)))



fig3 <- plot_grid(plot_grid(fig3a %>% no_x(),
                            EMPTY,
                            fig3b %>% no_leg() %>% no_xy(),
                            fig3c %>% no_leg(),
                            EMPTY,
                            fig3d %>% no_leg() %>% no_y(),
                            labels = c("a", "", "b", "c", "", "d"),
                            nrow = 2, align = "vh", rel_widths = c(1, 0.1, 1),
                            label_fontface = "plain", label_size = 16),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))
# fig3


# save_file(fig3, "fig3", width = 6, height = 5)



# -----*
# With PCs 1 and 3 ----
# (perhaps for the supplement)
# -----*

figS1a <- pred_pca$taxon_vec %>%
    mutate(taxon = factor(taxon, levels = taxa_lvls %>% rev(),
                          labels = c("ground\nspiders","wolf\nspiders","sheet\nweavers",
                                     "harvestman","ground\nbeetles","rove\nbeetles")[rev(taxa_order)])) %>%
    arrange(taxon) %>%
    mutate(PC1 = -5*PC1,
           PC3 = 5*PC3,
           PC1_lab = PC1 + c(1.4, 0.6, 0,
                                -1.3, -1.2, 0.2)[rev(taxa_order)],
           PC3_lab = PC3 + c(0.8, -0.9, 0.9,
                                -2.5, -0.3, 0.8)[rev(taxa_order)]) %>%
    arrange(desc(taxon)) %>%
    ggplot()+
    geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
    geom_segment(aes(x = 0, xend = PC1, y = 0, yend = PC3, group = taxon,
                     color = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 1)+
    geom_text(aes(label = taxon, x = PC1_lab, y = PC3_lab,
                  color = taxon),
              size = 9 / 2.83465, lineheight = 0.75)+
    geom_segment(data = tibble(PC1 = -0.5,
                               PC2 = -2.5,
                               PC1e = -0.2,
                               PC2e = -0.1),
                 aes(x = PC1, y = PC2, xend = PC1e, yend = PC2e),
                 color = "black", size = 0.25) +
    scale_x_continuous("PC1") +  ## , breaks = 3*-1:1) +
    scale_y_continuous("PC3") +  ## , breaks = 3*-1:1) +
    coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                ylim = c(-pc_axis_lim, pc_axis_lim)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(1,4,3,5,2,6)],
                       guide = FALSE) +
    theme(plot.margin = margin(0,0,0,0))

# plot time effect
figS1b <- biplot_fn("time", .PCs = c(1,3), .mult = c(-2, 2))
# plot distance effect
figS1c <- biplot_fn("distance", .PCs = c(1,3), .mult = c(-2, 2))
# plot midge effect
figS1d <- biplot_fn("midges", .PCs = c(1,3), .mult = c(-2, 2))


S1_pca_legend <- get_legend(figS1b + theme(legend.title.align = 0,
                                           legend.title = element_text(size = 11)))



figS1 <- plot_grid(plot_grid(figS1a %>% no_x(),
                            EMPTY,
                            figS1b %>% no_leg() %>% no_xy(),
                            figS1c %>% no_leg(),
                            EMPTY,
                            figS1d %>% no_leg() %>% no_y(),
                            labels = c("a", "", "b", "c", "", "d"),
                            label_fontface = "plain", label_size = 16,
                            nrow = 2, align = "vh", rel_widths = c(1, 0.1, 1)),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))
# figS1


# save_file(figS1, "figS1", width = 6, height = 5)









# ===============*
# Table I ----
# ===============*


#---------*
# * LOO deviance ----
#---------*

# Table I (exclude standard errors)
loo_dev <- read_csv("analysis/output/dev_re.csv", col_types = cols()) %>%
    rename(dev_re = delt_looic) %>%
    full_join(read_csv("analysis/output/dev_fere.csv", col_types = cols()) %>%
                  rename(dev_fere = delt_looic),
              by = "model") %>%
    mutate(model = gsub("_z$", "", model)) %>%
    rename(var = model)




#---------*
# * Variance partitioning ----
#---------*


# variance partition of PC axes by predictors (for Table I)
var_part <- lapply(c("PC1","PC2","PC3"), function(x){
    y = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(var = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(y, data = pred_pca$axes))[,2] %>%
               {.[1:3]/sum(.[1:3])} %>%
               round(3))
    }) %>%
    bind_rows() %>%
    spread(pc, cont) %>%
    mutate(var = gsub("_z$", "", var))

# overall variance accounted for by predictors (for Table I)
overall_part <- as.matrix(var_part[,2:4]) %*% t(as.matrix(pred_pca$obs_exp[1,2:4]))
row.names(overall_part) <- var_part$var


# Which PC axes are most associated with each variable?
pc_vars <- var_part[,2:4] %>%
    as.matrix() %>%
    `*`(matrix(as.numeric(pred_pca$obs_exp[1,2:4]), 3, 3, byrow=TRUE)) %>%
    as.data.frame() %>%
    as_tibble() %>%
    set_names(paste0("PC", 1:3)) %>%
    mutate(var = var_part$var) %>%
    select(var, everything())


coef_order <- c("time", "dist", "midges")

fmt <- function(x, .f = "%.2f") sprintf(.f,x)

tbl1_order <- function(x, .col = NULL, .coef_order = c("time", "dist", "midges")) {
    if (inherits(x, "matrix")) {
        stopifnot(!is.null(rownames(x)))
        if (is.null(.col)) .col <- 1
        x <- as.numeric(x[.coef_order,1])
    } else if (inherits(x, "data.frame")) {
        stopifnot(!is.null(.col) && (is.numeric(.col) || is.character(.col)) &&
                      length(.col) == 1)
        stopifnot("var" %in% colnames(x))
        x <- as.numeric(x[match(.coef_order, x$var),][[.col]])
    } else stop("x must be data frame or matrix")
    return(x)
}

tibble(`coef` = c("", "time", "distance", "midges"),
       `Taxon-variation` = c("(random)", tbl1_order(loo_dev,"dev_re") %>% fmt("%.1f")),
       `Overall` = c("(fixed + random)", tbl1_order(loo_dev,"dev_fere") %>% fmt("%.1f")),
       `EXTRA` = rep("", 4),
       PC1 = c(pred_pca$obs_exp[["PC1"]][1] %>% fmt("(%.2f)"),
               tbl1_order(var_part, "PC1") %>% fmt()),
       PC2 = c(pred_pca$obs_exp[["PC2"]][1] %>% fmt("(%.2f)"),
               tbl1_order(var_part, "PC2") %>% fmt()),
       PC3 = c(pred_pca$obs_exp[["PC3"]][1] %>% fmt("(%.2f)"),
               tbl1_order(var_part, "PC3") %>% fmt()),
       Total = c(pred_pca$obs_exp[1,paste0("PC", 1:3)] %>% sum() %>% fmt("(%.2f)"),
                 fmt(overall_part[coef_order,]))) %>%
    knitr::kable(format = "latex")



