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




taxa_order <- c(5:6, 4, 1, 3, 2)


# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(Taxon = factor(taxon,
                          levels = c("gnap","lyco","sheet","opil","cara","stap")[taxa_order],
                          labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                     "Harvestman","Ground beetles","Rove beetles")[taxa_order]))
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")


coef_sum$beta <- coef_sum$beta %>%
    ungroup() %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")),
           tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet",
                                                    "opil","cara","stap")[taxa_order])))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges")))
coef_sum$ar <- coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = c("gnap","lyco","sheet","opil",
                                                    "cara","stap")[taxa_order])))



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
    facet_wrap(~Taxon, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    # geom_point(aes(group = plot), size = 1.5, alpha = 0.5)+
    geom_line(data = data_fit %>%
                  group_by(Taxon, year) %>%
                  summarize(y = mean(midges_z)) %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(Taxon, year) %>%
                                summarize(y = mean(y)) %>%
                                mutate(type  = "Mean by date/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by date/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016), limits = c(2007,2018))+
    scale_y_continuous("Transformed abundance",
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
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(Taxon, distance) %>%
                                summarize(y = mean(y)) %>%
                                mutate(type  = "Mean by date/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by date/distance"))),
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


# save_file(fig1, "fig1", width = 6, height = 4)




#==========*
#========== Plot coefficients (Use `taxon-specicific-slopes` for Figure 2) ----
#==========*



# taxon-specific slopes
slope_p <- coef_sum$beta %>%
    ggplot()+
    facet_wrap(~coef, ncol= 1)+
    geom_rect(data = coef_sum$alpha,
              aes(xmin = -1, xmax = 10, ymin = lo, ymax = hi, fill = coef),
              alpha = 0.25, linetype =  0)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(aes(tx, mi, color = coef), size = 1.5)+
    geom_linerange(aes(tx, ymin = lo, ymax = hi, color = coef))+
    geom_text(data = tibble(coef = sort(unique(coef_sum$beta[["coef"]])),
                            x = 0, y = -0.6),
              aes(x, y, label = coef),
              size = 12 / 2.83465, hjust = 0, vjust = 1, color = "black") +
    scale_y_continuous("Response",
                       breaks = c(-0.5, 0, 0.5))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles",
                                  "Rove beetles")[taxa_order],
                       trans = "reverse") +
    scale_fill_manual(NULL, values = coef_palette, guide = FALSE)+
    scale_color_manual(NULL, values = coef_palette, guide = FALSE)+
    # theme(plot.margin = margin(r=0,l=0)) +
    theme(plot.margin = margin(0,0,0,t=4), strip.text.x = element_blank()) +
    coord_flip(ylim = c(-0.6, 0.6), xlim = c(6.5, 0)) +
    NULL



small_labels <- theme(axis.title.y = element_text(size = 10, angle = 90,
                                                  margin = margin(0,0,0,0)),
                      axis.title.x = element_text(size = 10,
                                                  margin = margin(0,0,0,0)),
                      plot.margin = margin(4,4,0,0))


# ar (extra)
ar_p <- coef_sum$ar %>%
    ggplot(aes(tx, mi))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 1.5)+
    # ggtitle(" ") +
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_y_continuous("AR parameter",  limits = c(0, 1),
                       breaks = c(0, 0.5, 1))+
    scale_x_continuous(NULL, breaks = 1:6,
                       labels = c("Ground spiders","Wolf spiders","Sheet weavers",
                                  "Harvestman","Ground beetles","Rove beetles")[taxa_order],
                       trans = "reverse", limits = c(6.5, 0.5))+
    small_labels +
    coord_flip() +
    theme(axis.text.y = element_blank()) +
    NULL



effect_p_theme <- small_labels +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")

# alphas (extra)
effect_mean_p <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    .[[1]] %>%
    as.matrix() %>%
    as_tibble(.name_repair = "universal") %>%
    gather() %>%
    filter(key != "...1") %>%  # intercept not necessary
    mutate(coef = factor(key, levels = c("...3","...4","...2"),
                         labels = c("time", "distance", "midges"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual("", values = coef_palette)+
    scale_y_continuous("Posterior density", limits = c(0, 7), breaks = 0:4 * 2)+
    scale_x_continuous("Response mean") +
    effect_p_theme +
    NULL



# sigmas (extra)
effect_sigmas_p <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    .[[1]] %>%
    as.matrix() %>%
    as_tibble(.name_repair = "universal") %>%
    gather() %>%
    mutate(coef = factor(key, levels = c("...1","...2","...3","...4","...5","...6"),
                         labels = c("int_tax","int_plot","int_trans","midges",
                                    "time","distance"))) %>%
    filter(!(coef %in% c("int_tax","int_plot","int_trans"))) %>%
    mutate(coef = factor(coef %>% paste(),
                         levels = c("time", "distance", "midges"))) %>%
    ggplot(aes(value, fill = coef))+
    geom_vline(xintercept = 0,  color = "gray50")+
    geom_density(linetype = 0, alpha = 0.7)+
    scale_fill_manual(NULL, values = coef_palette)+
    guides(fill = guide_legend(keywidth = 0.75, keyheight = 0.75)) +
    scale_y_continuous("Posterior density", limits = c(0, 7), breaks = 0:4 * 2)+
    scale_x_continuous("Response SD") +
    effect_p_theme +
    NULL


fig2 <- plot_grid(EMPTY,
                  slope_p,
                  plot_grid(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
                            labels = c("B", "", "C", "", "D"),
                            rel_heights = c(1.265, 0.1, 1, 0.1, 1),
                            ncol = 1, label_x = 0.25),
                  plot_grid(ar_p,
                            EMPTY,
                            effect_mean_p,
                            EMPTY,
                            effect_sigmas_p,
                            rel_heights = c(1.265, 0.1, 1, 0.1, 1),
                            ncol = 1),
                  labels = c("A", "", "", ""), rel_widths = c(0.08, 1, 0.12, 0.6),
                  nrow = 1)

# fig2


# save_file(fig2, "fig2", width = 6, height = 6)





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


# # biplot function
# # (Only doing PC 1 and 2 for main text)
# biplot_fn <- function(var_, .mult = -2, .PCs = 1:2) {
#
#     # var_ = "time"
#     # pred_pca = pred_pca
#     # .mult = c(-2, 2)
#     # .PCs = 2:3
#
#     if (length(.mult) == 1) .mult <- rep(.mult, 2)
#     .signs <- ifelse(.mult < 0, "-", "")
#     .sign_nums <- as.numeric(paste0(.signs, "1"))
#
#     var_ <- match.arg(var_, c("midges", "time", "distance"))
#
#     var_z_ <- paste0(gsub("^distance$", "dist", var_), "_z")
#
#     pc_axes <- paste0(.signs, "PC", .PCs)
#
#     pred_pca$obs_rot %>%
#         # Uncomment below if you want darker spots in front
#         # I currently don't think we should do it
#         ## arrange(desc(!!sym(var_z_))) %>%
#         ggplot(aes_string(x = pc_axes[1], y = pc_axes[2]))+
#         geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
#         geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
#         geom_point(aes_string(color = var_z_),
#                    alpha = 0.25, size = 2)+
#         geom_segment(data = pred_vec %>% filter(var == gsub("_z$", "", var_z_)),
#                      aes_string(x = 0,
#                                 xend = paste0(".mult[1]*PC", .PCs[1]),
#                                 y = 0,
#                                 yend = paste0(".mult[2]*PC", .PCs[2])),
#                      arrow = arrow(length = unit(6, "pt")),
#                      size = 0.8, color = "black") +
#         geom_text(data = tibble(A = .sign_nums[1] * (pc_axis_lim + 0.2),
#                                 B = .sign_nums[2] * (pc_axis_lim + 0.2)) %>%
#                       set_names(paste0("PC", .PCs)),
#                   label = var_, hjust = 1, vjust = 1, size = 12 / 2.83465) +
#         scale_colour_gradient2("Predictor\nvalue", low =  pca_palette[1],
#                                mid = pca_palette[2],
#                                high = pca_palette[3],
#                                midpoint = -0.5, limits = c(-3,2),
#                                breaks = c(-3,-0.5,2))+
#         scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
#         scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
#         coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
#                     ylim = c(-pc_axis_lim, pc_axis_lim)) +
#         theme(plot.title = element_text(size = 11, hjust = 0.5, vjust = 0,
#                                         face = "plain"),
#               plot.margin = margin(0,0,0,0))
# }


size_scale <- function(.guide = "legend") {
    scale_size(limits = range(pred_pca$obs_rot$PC3 * -2.01),
               range = c(0.5, 3), breaks = c(-4,0,4), guide = .guide)
}

set.seed(2041119884)
.inds <- sample.int(nrow(pred_pca$obs_rot),
                    round(nrow(pred_pca$obs_rot) * 0.25))


# All three axes:
triplot_fn <- function(var_, .mult = c(-2, -2, -2)) {

    # var_ = "distance"
    # .mult = c(-2, -2, -2)
    .PCs <- 1:3

    stopifnot(length(.mult) %in% c(1,3))
    if (length(.mult) == 1) .mult <- rep(.mult, 3)
    .signs <- ifelse(.mult < 0, "-", "")
    .sign_nums <- as.numeric(paste0(.signs, "1"))

    var_ <- match.arg(var_, c("midges", "time", "distance"))

    var_z_ <- paste0(gsub("^distance$", "dist", var_), "_z")

    pc_axes <- paste0(.signs, "PC", .PCs)

    .shape <- ifelse(pred_vec %>%
                         filter(var == gsub("_z$", "", var_z_)) %>%
                         .[["PC2"]] %>%
                         `*`(.mult[2]) %>%
                         `>=`(0), 24, 25)

    .pca_df <- pred_pca$obs_rot %>%
        # .[.inds,] %>%x
        mutate(PC1 = .sign_nums[1] * PC1,
               PC2 = .sign_nums[2] * PC2,
               PC3 = .sign_nums[3] * PC3) %>%
        arrange(PC3)

    .path_df <- pred_vec %>%
        filter(var == gsub("_z$", "", var_z_)) %>%
        mutate(PC1 = .mult[1] * PC1,
               PC2 = .mult[2] * PC2,
               PC3 = .mult[3] * PC3) %>%
        mutate(PC1 = list(seq(0, PC1, length.out = 101)),
               PC2 = list(seq(0, PC2, length.out = 101)),
               PC3 = list(seq(min(pred_pca$obs_rot$PC3 * .mult[3]),
                              PC3, length.out = 101))) %>%
        unnest(cols = c(PC1, PC2, PC3))

    # .arrow_size <- (.path_df[nrow(.path_df),][["PC3"]] - min(.pca_df$PC3)) /
    #     (diff(range(.pca_df$PC3)))

    .pca_df %>%
        # Uncomment below if you want darker spots in front
        # I currently don't think we should do it
        ## arrange(desc(!!sym(var_z_))) %>%
        ggplot(aes(PC1, PC2, size = PC3))+
        geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
        geom_point(aes_string(color = var_z_, fill = var_z_), shape = 21) +
        geom_path(data = .path_df, color = "black", alpha = 0.5) +
        geom_point(data = .path_df[nrow(.path_df),],
                   color = "black", fill = "black", alpha = 0.75,
                   shape = .shape) +
        geom_text(data = tibble(A = (pc_axis_lim + 0.2),
                                B = (pc_axis_lim + 0.2)) %>%
                      set_names(paste0("PC", .PCs[1:2])),
                  label = var_, hjust = 1, vjust = 1, size = 12 / 2.83465) +
        scale_color_viridis_c("Predictor\nvalue", alpha = 0.5,
                              limits = c(-3,2), breaks = c(-3,-0.5,2)) +
        scale_fill_viridis_c(limits = c(-3,2), breaks = c(-3,-0.5,2),
                             alpha = 0.05, guide = FALSE) +
        size_scale() +
        scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
        scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        theme(plot.title = element_text(size = 11, hjust = 0.5, vjust = 0,
                                        face = "plain"),
              plot.margin = margin(0,0,0,0)) +
        guides(color = guide_colorbar(override.aes = list(alpha = 1)),
               size = guide_legend(override.aes = list(shape = 1)))
}



fig3a_df <- pred_pca$taxon_vec %>%
    mutate(taxon = factor(taxon, levels = c("gnap","lyco","sheet","opil",
                                            "cara","stap") %>% rev(),
                          labels = c("ground\nspiders","wolf\nspiders",
                                     "sheet\nweavers",
                                     "harvestman","ground\nbeetles",
                                     "rove\nbeetles") %>% rev())) %>%
    arrange(desc(taxon)) %>%
    mutate(PC1 = -5 * PC1,
           PC2 = -5 * PC2,
           PC3 = -5 * PC3,
           PC1_lab = PC1 + c(1.2, 0.2, 1, 0.5, -0.8, 0.4),
           PC2_lab = PC2 + c(-0.4, -0.8, -0.5, 1.7, 0.7, -0.9),
           ang = c(rep(0,3), 73.5, 0, 0),
           PC1 = map(PC1, ~ seq(0, .x, length.out = 101)),
           PC2 = map(PC2, ~ seq(0, .x, length.out = 101)),
           PC3 = map(PC3, ~ seq(min(pred_pca$obs_rot$PC3 * -2),
                                .x, length.out = 101))) %>%
    unnest(cols = c(PC1, PC2, PC3)) %>%
    arrange(taxon, PC3)

fig3a <- fig3a_df %>%
    ggplot()+
    geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
    # geom_segment(aes(x = 0, xend = -5*PC1, y = 0, yend = -5*PC2, group = taxon,
    #                  color = taxon),
    #              arrow = arrow(length = unit(6, "pt")),
    #              size = 1)+
    geom_line(aes(x = PC1, y = PC2, color = taxon, size = PC3), alpha = 0.5)+
    geom_point(data = fig3a_df %>%
                   group_by(taxon) %>%
                   filter(PC3 == max(PC3)) %>%
                   ungroup() %>%
                   mutate(shape = factor(ifelse(PC2 >= 0, 0, 1))),
               aes(x = PC1, y = PC2, color = taxon, fill = taxon, size = PC3,
                   shape = shape))+
    geom_text(data = fig3a_df %>% distinct(taxon, .keep_all = TRUE),
              aes(label = taxon, x = PC1_lab, y = PC2_lab,
                  color = taxon, angle = ang),
              size = 9 / 2.83465, lineheight = 0.75)+
    scale_x_continuous("PC1", breaks = 3*-1:1) +
    scale_y_continuous("PC2", breaks = 3*-1:1) +
    coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                ylim = c(-pc_axis_lim, pc_axis_lim)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(5:6, 4:1)],
                       guide = FALSE, aesthetics = c("color", "fill")) +
    scale_shape_manual(values = c(24, 25), guide = FALSE) +
    size_scale(FALSE) +
    theme(plot.margin = margin(0,0,0,0))

# plot time effect
fig3b <- triplot_fn("time")
# plot distance effect
fig3c <- triplot_fn("distance")
# plot midge effect
fig3d <- triplot_fn("midges")


pca_legend <- get_legend(fig3b +
                             scale_color_viridis_c("Predictor\nvalue",
                                                   alpha = 1,
                                                   limits = c(-3,2),
                                                   breaks = c(-3,-0.5,2)) +
                             theme(legend.title.align = 0,
                                       legend.title = element_text(size = 11)))



fig3 <- plot_grid(plot_grid(fig3a %>% no_x(),
                            EMPTY,
                            fig3b %>% no_leg() %>% no_xy(),
                            fig3c %>% no_leg(),
                            EMPTY,
                            fig3d %>% no_leg() %>% no_y(),
                            labels = c("A", "", "B", "C", "", "D"),
                            nrow = 2, align = "vh", rel_widths = c(1, 0.1, 1)),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))
fig3


# save_file(fig3, "fig3-wPC3", width = 6, height = 5)



# -----*
# With PCs 2 and 3 ----
# (perhaps for the supplement)
# -----*

figS1a <- pred_pca$taxon_vec %>%
    mutate(taxon = factor(taxon, levels = c("gnap","lyco","sheet","opil","cara","stap"),
                          labels = c("ground\nspiders","wolf\nspiders","sheet\nweavers",
                                     "harvestman","ground\nbeetles","rove\nbeetles"))) %>%
    arrange(taxon) %>%
    mutate(PC2_lab = -5*PC2 + c(-0.8, -1.4, -1.4,
                                1.0, 1.2, 1.2),
           PC3_lab = 5*PC3 + c(-0.8, -0.5, 0,
                                1.5, 0.5, -0.4),
           ang = c(rep(0,3), 52.18151, 0, 0)) %>%
    arrange(desc(taxon)) %>%
    ggplot()+
    geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
    geom_segment(aes(x = 0, xend = -5*PC2, y = 0, yend = 5*PC3, group = taxon,
                     color = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 1)+
    geom_text(aes(label = taxon, x = PC2_lab, y = PC3_lab,
                  color = taxon, angle = ang),
              size = 9 / 2.83465, lineheight = 0.75)+
    scale_x_continuous("PC2") +  ## , breaks = 3*-1:1) +
    scale_y_continuous("PC3") +  ## , breaks = 3*-1:1) +
    coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                ylim = c(-pc_axis_lim, pc_axis_lim)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(1:4, 6:5)],
                       guide = FALSE) +
    theme(plot.margin = margin(0,0,0,0))

# plot time effect
figS1b <- biplot_fn("time", .PCs = 2:3, .mult = c(-2, 2))
# plot distance effect
figS1c <- biplot_fn("distance", .PCs = 2:3, .mult = c(-2, 2))
# plot midge effect
figS1d <- biplot_fn("midges", .PCs = 2:3, .mult = c(-2, 2))


S1_pca_legend <- get_legend(figS1b + theme(legend.title.align = 0,
                                           legend.title = element_text(size = 11)))



figS1 <- plot_grid(plot_grid(figS1a %>% no_x(),
                            EMPTY,
                            figS1b %>% no_leg() %>% no_xy(),
                            figS1c %>% no_leg(),
                            EMPTY,
                            figS1d %>% no_leg() %>% no_y(),
                            labels = c("A", "", "B", "C", "", "D"),
                            nrow = 2, align = "vh", rel_widths = c(1, 0.1, 1)),
                  pca_legend, nrow = 1, rel_widths = c(1, 0.2))
# figS1


# save_file(figS1, "figS1", width = 6, height = 5)





# -----*
# With PCs 1 and 3 ----
# (perhaps for the supplement)
# -----*

figS2a <- pred_pca$taxon_vec %>%
    mutate(taxon = factor(taxon, levels = c("gnap","lyco","sheet","opil","cara","stap"),
                          labels = c("ground\nspiders","wolf\nspiders","sheet\nweavers",
                                     "harvestman","ground\nbeetles","rove\nbeetles"))) %>%
    arrange(taxon) %>%
    mutate(PC1_lab = -5*PC1 + c(0, 1.0, 1.0,
                                -1.3, -1.2, 0.4),
           PC3_lab = 5*PC3 + c(-0.8, -0.4, 0.8,
                               1.2, 0.4, -0.8),
           ang = c(rep(0, 3), -53.27706, rep(0, 2))) %>%
    arrange(desc(taxon)) %>%
    ggplot()+
    geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
    geom_segment(aes(x = 0, xend = -5*PC1, y = 0, yend = 5*PC3, group = taxon,
                     color = taxon),
                 arrow = arrow(length = unit(6, "pt")),
                 size = 1)+
    geom_text(aes(label = taxon, x = PC1_lab, y = PC3_lab,
                  color = taxon, angle = ang),
              size = 9 / 2.83465, lineheight = 0.75)+
    scale_x_continuous("PC1") +  ## , breaks = 3*-1:1) +
    scale_y_continuous("PC3") +  ## , breaks = 3*-1:1) +
    coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                ylim = c(-pc_axis_lim, pc_axis_lim)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(6, "Dark2")[c(1:4, 6:5)],
                       guide = FALSE) +
    theme(plot.margin = margin(0,0,0,0))

# plot time effect
figS2b <- biplot_fn("time", .PCs = c(1,3), .mult = c(-2, 2))
# plot distance effect
figS2c <- biplot_fn("distance", .PCs = c(1,3), .mult = c(-2, 2))
# plot midge effect
figS2d <- biplot_fn("midges", .PCs = c(1,3), .mult = c(-2, 2))


S2_pca_legend <- get_legend(figS2b + theme(legend.title.align = 0,
                                           legend.title = element_text(size = 11)))



figS2 <- plot_grid(plot_grid(figS2a %>% no_x(),
                             EMPTY,
                             figS2b %>% no_leg() %>% no_xy(),
                             figS2c %>% no_leg(),
                             EMPTY,
                             figS2d %>% no_leg() %>% no_y(),
                             labels = c("A", "", "B", "C", "", "D"),
                             nrow = 2, align = "vh", rel_widths = c(1, 0.1, 1)),
                   pca_legend, nrow = 1, rel_widths = c(1, 0.2))
# figS2


# save_file(figS2, "figS2", width = 6, height = 5)








# ===============*
# Table I ----
# ===============*


#---------*
# * LOO deviance ----
#---------*

# Table I (exclude standard errors)
loo_dev <- read_csv("analysis/output/dev_re.csv", col_types = cols()) %>%
    mutate(dev = 2*abs(elpd_diff),
           dev_se = 2*se_diff) %>%
    select(var, dev, dev_se) %>%
    rename( dev_re = dev, dev_re_se = dev_se) %>%
    full_join(read_csv("analysis/output/dev_fere.csv", col_types = cols()) %>%
                  mutate(dev = 2*abs(elpd_diff),
                         dev_se = 2*se_diff) %>%
                  select(var,dev, dev_se) %>%
                  rename(dev_fere = dev, dev_fere_se = dev_se),
              by = "var") %>%
    filter(var != "full") %>%
    mutate(var = gsub("_z$", "", var))



#---------*
# * Variance partitioning ----
#---------*


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
       `Taxon-variation` = c("(random)", tbl1_order(loo_dev,"dev_re") %>% fmt("%.0f")),
       `Overall` = c("(fixed + random)", tbl1_order(loo_dev,"dev_fere") %>% fmt("%.0f")),
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



