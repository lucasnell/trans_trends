#'
#' This file contains PC analyses.
#'


# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")

data_df <- read_rds(data_rds)
model_fit <- read_rds(model_rds)




# =============================================================================*
# Extract info from model fit ----
# =============================================================================*


# Gets uncertainty intervals via quantiles, and estimates via median
get_fit_info <- function(.var) {
    z <- do.call(c, lapply(1:model_fit$stan@sim$chains,
                           function(i) {
                               model_fit$stan@sim$samples[[i]][[.var]][
                                   -(1:model_fit$stan@sim$warmup2[i])]
                           }))
    return(tibble(var = .var,
                  lo = unname(quantile(z, 0.16)),
                  mi = median(z),
                  hi = unname(quantile(z, 0.84))))
}

fit_sum <- rstan::summary(model_fit$stan, probs = c()) |>
    getElement("summary") |>
    as.data.frame() |>
    rownames_to_column() |>
    as_tibble() |>
    rename(var = rowname) |>
    (\(x) left_join(x, map_dfr(x$var, get_fit_info), by = "var"))() |>
    select(var, lo, mi, hi, n_eff, Rhat)

# taxon-specific slopes
beta <- fit_sum |>
    filter(str_detect(var, "beta"), !str_detect(var, "sig")) |>
    mutate(id = strsplit(var, "\\[|\\]|,") |> map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") |> map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = 1:4,
                         labels = c("int","midges","time","dist"))) |>
    full_join(data_df |>
                  distinct(plot, taxon) |>
                  mutate(id = row_number()),
              by = "id") |>
    filter(coef != "int") |>
    group_by(coef, taxon) |>
    summarize(lo = unique(lo),
              mi = unique(mi),
              hi = unique(hi),
              .groups = "drop")


# taxon-specific intercepts
# average over plot and transect variation
# standardize to mean for each predictor (only nonzero for time_z)
int_taxon <- fit_sum |>
    filter(str_detect(fit_sum$var, "beta"), !str_detect(fit_sum$var, "sig")) |>
    getElement("var") |>
    rstan::extract(object = model_fit$stan) |>
    bind_cols() |>
    mutate(step = row_number()) |>
    pivot_longer(-step, names_to = "var", values_to = "val") |>
    mutate(id = strsplit(var, "\\[|\\]|,") |> map_int(~as.integer(.x[2])),
           coef = strsplit(var, "\\[|\\]|,") |> map_int(~as.integer(.x[3])),
           coef = factor(coef, levels = c(1:4),
                         labels = c("int","midges","time","dist"))) |>
    full_join(data_df |>
                  distinct(plot, taxon) |>
                  mutate(id = row_number()),
              by = "id") |>
    filter(coef %in% c("int","time"))  |>
    group_by(step, taxon, coef) |>
    summarize(val = mean(val), .groups = "drop") |>
    pivot_wider(names_from = coef, values_from = val) |>
    mutate(val = int + time * mean(data_df$time_z)) |>
    split(~ taxon) |>
    map(function(.x) {
        .dd <- tibble(taxon = .x$taxon[1],
                      lo = unname(quantile(.x$val, 0.16)),
                      mi = median(.x$val),
                      hi = unname(quantile(.x$val, 0.84)))
        return(.dd)
    }) |>
    list_rbind()



# =============================================================================*
# PCA ----
# =============================================================================*



#' Run pca and package into a list.
#'
#' I did it this way to avoid cluttering environment with variables used
#' multiple times for just the PCA.
#'
pred_pca <- with(list(d_ = data_df, beta_ = beta, int_ = int_taxon), {


    # predicted values from model
    pred_ <- d_ |>
        expand(taxon,
               midges_z = seq(min(midges_z), max(midges_z), 1),
               time_z = seq(min(time_z), max(time_z), 1),
               dist_z = seq(min(dist_z), max(dist_z), 1)) |>
        full_join(beta_ |>
                      select(taxon, coef, mi) |>
                      pivot_wider(names_from = coef, values_from = mi),
                  by = "taxon") |>
        full_join(int_taxon |> select(taxon, mi) |> rename(int = mi), by = "taxon") |>
        # (assume the offset, season_day, is 1 below)
        mutate(count = exp(int + midges*midges_z + time*time_z + dist*dist_z)) |>
        group_by(taxon) |>
        mutate(pred_id = row_number()) |>
        ungroup()

    # pca on predicted values
    pca_ <- pred_ |>
        select(pred_id, taxon, count) |>
        pivot_wider(names_from = taxon, values_from = count) |>
        select(-pred_id) |>
        as.matrix() |>
        prcomp(center = TRUE, scale = FALSE)

    # axes from pca
    axes_ <- pred_ |>
        select(pred_id, midges_z, time_z, dist_z) |>
        full_join(as_tibble(pca_$x) |>
                      mutate(pred_id = row_number()),
                  by = "pred_id")

    # taxon vectors from pca
    taxon_vec_ <- pca_$rotation |>
        as_tibble() |>
        mutate(taxon_id = row_number()) |>
        left_join(tibble(taxon = sort(unique(d_$taxon))) |>
                      mutate(taxon_id = row_number()), by = "taxon_id") |>
        select(taxon, matches("PC"))

    # rotation of observed data
    obs_ <- d_ |>
        group_by(taxon) |>
        mutate(obs_id = row_number()) |>
        ungroup() |>
        mutate(count = count / season_days)
    obs_rot_ <- obs_ |>
        select(obs_id, plot, trans, midges_z, time_z, dist_z) |>
        full_join(as_tibble(predict(pca_, obs_ |>
                                        select(obs_id, taxon, count) |>
                                        pivot_wider(names_from = taxon,
                                                    values_from = count) |>
                                        select(-obs_id) |>
                                        as.matrix())) |>
                      mutate(obs_id = row_number()),
                  by = "obs_id")

    # variance explained in observed data by rotation
    obs_exp_ <- obs_rot_ |>
        select(matches("PC")) |>
        apply(2, sd) |>
        (\(x) {x^2/sum(x^2)})() |>
        (\(x) rbind(x, cumsum(x)))() |>
        as_tibble() |>
        mutate(type = c("individual","cumulative")) |>
        select(type, matches("PC"))

    list(pred = pred_, pca = pca_, axes = axes_, taxon_vec = taxon_vec_,
                obs_rot = obs_rot_, obs_exp = obs_exp_)
})



# taxon vectors
pred_pca$taxon_vec |>
    arrange(-abs(PC1))

# # A tibble: 6 × 7
#   taxon       PC1       PC2       PC3      PC4       PC5          PC6
#   <fct>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>        <dbl>
# 1 lyco  -0.81730   0.52913  -0.15485  -0.16344  0.036818 -0.000050249
# 2 sheet -0.50396  -0.70621   0.41178  -0.18792 -0.14038  -0.15070
# 3 opil  -0.23467  -0.18448   0.041898  0.64197  0.46857   0.52674
# 4 gnap  -0.12164   0.044956 -0.085717  0.69502 -0.62203  -0.32537
# 5 stap  -0.070271 -0.35967  -0.74559  -0.18039 -0.32696   0.41274
# 6 cara  -0.056926 -0.23638  -0.49138   0.10124  0.51533  -0.65086



# variance explained in predicted values (first three axes must
# explain everything)
summary(pred_pca$pca)

# Importance of components:
#                           PC1    PC2     PC3     PC4      PC5      PC6
# Standard deviation     0.2836 0.1824 0.03840 0.02098 0.006776 0.003983
# Proportion of Variance 0.6952 0.2877 0.01275 0.00381 0.000400 0.000140
# Cumulative Proportion  0.6952 0.9829 0.99566 0.99947 0.999860 1.000000



# variance explained in observed values (For Results)
pred_pca$obs_exp

# # A tibble: 2 × 7
#   type           PC1     PC2      PC3      PC4      PC5      PC6
#   <chr>        <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
# 1 individual 0.55654 0.26052 0.067779 0.057361 0.024922 0.032871
# 2 cumulative 0.55654 0.81707 0.88485  0.94221  0.96713  1



# predictor vectors
pred_vec <- map_dfr(c("time", "dist", "midges"),
                    function(pred_) {
                        sym_ <- sym(paste0(pred_, "_z"))
                        pred_pca$axes |>
                            filter(!!sym_ %in% range(!!sym_)) |>
                            group_by(!!sym_) |>
                            summarize(PC1 = mean(PC1),
                                      PC2 = mean(PC2),
                                      PC3 = mean(PC3),
                                      .groups = "drop") |>
                            rename(val = !!sym_) |>
                            arrange(val) |>
                            mutate(coef = pred_)  |>
                            mutate(PC1 = PC1 - PC1[1],
                                   PC2 = PC2 - PC2[1],
                                   PC3 = PC3 - PC3[1]) |>
                            filter(PC1 != 0 | PC2 != 0 | PC3 != 0) |>
                            select(coef, starts_with("PC"))
                    })



# =============================================================================*
# Plot prelims ----
# =============================================================================*

pc_axis_lim <- 4.1
#' #' We're scaling some of the projections. This is to keep that consistent.
pc_mults <- list(pred = 6 * c(-1, 1, 1), taxon = 5 * c(-1, 1, 1))

# pc_axis_lim <- 6
# pc_mults <- list(pred = c(1, 1, 1), taxon = c(1, 1, 1))

stopifnot(identical(sign(pc_mults$pred), sign(pc_mults$taxon)))


pca_theme <- theme(plot.margin = margin(t=8,r=8,b=0,l=8),
                   axis.text.y = element_text(size = 8,
                                              margin = margin(0,0,0,r=2)),
                   axis.title.y = element_text(size = 10,
                                               margin = margin(0,0,0,0)),
                   axis.text.x = element_text(size = 8,
                                              margin = margin(0,0,0,t=2)),
                   axis.title.x = element_text(size = 10,
                                               margin = margin(0,0,0,b=8)))





# =============================================================================*
# Plots of taxon response vectors ----
# =============================================================================*


taxon_pca_fun <- function(.xPC, .yPC,
                          .nudge_x = NULL,
                          .nudge_y = NULL,
                          .segment_df = NULL) {

    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    .PCs <- c(.xPC, .yPC)

    .dd <- pred_pca$taxon_vec |>
        mutate(taxon = factor(taxon, levels = taxa_lvls |> rev(),
                              labels = taxa_labs |> rev() |> tolower() |>
                                  str_replace_all(" ", "\n"))) |>
        arrange(taxon) |>
        mutate(PC1 = pc_mults$taxon[1] * PC1,
               PC2 = pc_mults$taxon[2] * PC2,
               PC3 = pc_mults$taxon[3] * PC3)

    if (!is.null(.nudge_x) && !is.null(.nudge_y)) {
        .labs <- paste0("PC", .PCs, "_lab")
        .rto <- length(taxa_labs):1
        .dd[[.labs[1]]] <- .dd[[paste0("PC", .PCs[1])]] + .nudge_x[.rto]
        .dd[[.labs[2]]] <- .dd[[paste0("PC", .PCs[2])]] + .nudge_y[.rto]
    }

    .dd <- .dd |>
        arrange(desc(taxon))

    .p <- .dd |>
        ggplot()+
        geom_hline(yintercept = 0, color = "gray90", linewidth = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", linewidth = 0.5) +
        geom_segment(aes(x = 0, xend = .data[[paste0("PC", .PCs[1])]],
                         y = 0, yend = .data[[paste0("PC", .PCs[2])]],
                         group = taxon, color = taxon),
                     arrow = arrow(length = unit(6, "pt")),
                     linewidth = 1) +
        scale_x_continuous(paste0("PC", .PCs[1]), breaks = 3*-1:1) +
        scale_y_continuous(paste0("PC", .PCs[2]), breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        scale_color_manual(values = RColorBrewer::brewer.pal(9, "RdYlBu") |>
                               base::`[`(c(1:3, 7:9)),
                           guide = "none") +
        pca_theme

    if (!is.null(.nudge_x) && !is.null(.nudge_y)) {
        .p <- .p +
            geom_text(aes(label = taxon, x = .data[[.labs[1]]], y = .data[[.labs[2]]],
                          color = taxon),
                      size = 7 / 2.83465, lineheight = 0.75)
    }

    if (!is.null(.segment_df)) {
        .p <- .p +
            geom_segment(data = .segment_df,
                         aes(x = x, y = y, xend = xend, yend = yend),
                         color = "black", linewidth = 0.25)
    }

    return(.p)
}





taxon_pca_p <- list(
    taxon_pca_fun(1, 2,
                  .nudge_x = rep(0, 6),  # c(-0.4,  0.5, -1.5,  2.7, -1.0,  0.5),
                  .nudge_y = rep(0, 6),  # c( 1.0, -1.4, -3.0, -2.6,  0.8, -1.3),
                  # .segment_df = tibble(x = c(1.25,  -1.3),
                  #                      y = c(-0.6,   -1.3),
                  #                      xend = c(-0.4, -0.4),
                  #                      yend = c(1.2,   0.6))),
                  .segment_df = NULL),
    taxon_pca_fun(1, 3),
    taxon_pca_fun(2, 3))


# wrap_plots(taxon_pca_p)



# =============================================================================*
# Plots of model predictors and observed data ----
# =============================================================================*

predictor_pca_fun <- function(.xPC, .yPC, .label_df = NULL, ...) {

    # .xPC = 1; .yPC = 2; .label_df = NULL

    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    .PCs <- c(.xPC, .yPC)

    .obs_rot <- pred_pca$obs_rot |>
        mutate(PC1 = sign(pc_mults$pred[1]) * PC1,
               PC2 = sign(pc_mults$pred[2]) * PC2,
               PC3 = sign(pc_mults$pred[3]) * PC3)

    coefs <- c("time", "distance", "midges")

    .pred_vec <- pred_vec |>
        mutate(PC1 = pc_mults$pred[1] * PC1,
               PC2 = pc_mults$pred[2] * PC2,
               PC3 = pc_mults$pred[3] * PC3) |>
        mutate(coef = factor(coef, levels = gsub("^distance$", "dist", coefs),
                            labels = coefs))


    pc_axes <- paste0("PC", .PCs)

    .p <- .obs_rot |>
        ggplot(aes(x = .data[[pc_axes[1]]], y = .data[[pc_axes[2]]]))+
        geom_hline(yintercept = 0, color = "gray90", linewidth = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", linewidth = 0.5) +
        geom_point(alpha = 0.25, size = 1, shape = 1,
                   # color = "gray60") +
                   color = viridis::viridis(1, begin = 0.9, end = 0.9)) +
        geom_segment(data = .pred_vec,
                     aes(x = 0, xend = .data[[pc_axes[1]]],
                         y = 0, yend = .data[[pc_axes[2]]], color = var),
                     arrow = arrow(length = unit(6, "pt")),
                     linewidth = 1) +
        scale_x_continuous(pc_axes[1], breaks = 3*-1:1) +
        scale_y_continuous(pc_axes[2], breaks = 3*-1:1) +
        scale_color_viridis_d(end = 0.75, guide = "none") +
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
            geom_text(data = .label_df |>
                          mutate(var = ifelse(grepl("^dist", lab),
                                              "distance", lab) |>
                                     factor(levels = vars)),
                      aes(x = x, y = y, label = lab, color = var),
                      size = 8 / 2.83465,
                      ...)
    }

    return(.p)

}



pred_pca_p <- list(
    predictor_pca_fun(1, 2,
                      # .label_df = tibble(x = c(-0.8, -1, pc_axis_lim),
                      #                    y = c(2.5, -1, 1),
                      #                    lab = c("midges", "time", "dist.")),
                      hjust = 1, vjust = 0.5),
    predictor_pca_fun(1, 3),
    predictor_pca_fun(2, 3))

wrap_plots(pred_pca_p)




# --------------*
# NEW Plots of taxon response vectors and model predictors ----
# --------------*

pca_plotter <- function(.xPC, .yPC, incl_data = FALSE) {


    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    xPC_str <- paste0("PC", .xPC)
    yPC_str <- paste0("PC", .yPC)

    taxon_p_df <- pred_pca$taxon_vec |>
        mutate(taxon2 = factor(taxon, levels = taxa_lvls |> rev(),
                               labels = taxa_labs |> rev() |> tolower() |>
                                   str_replace_all(" ", "\n"))) |>
        arrange(taxon) |>
        mutate(PC1 = pc_mults$taxon[1] * PC1,
               PC2 = pc_mults$taxon[2] * PC2,
               PC3 = pc_mults$taxon[3] * PC3) |>
        arrange(desc(taxon))
    pred_p_df <- pred_vec |>
        mutate(PC1 = pc_mults$pred[1] * PC1,
               PC2 = pc_mults$pred[2] * PC2,
               PC3 = pc_mults$pred[3] * PC3) |>
        mutate(coef = factor(coef, levels = c("time", "dist", "midges"),
                             labels = c("time", "distance", "midges")))

    p <- taxon_p_df |>
        ggplot() +
        geom_hline(yintercept = 0, color = "gray90", linewidth = 0.5) +
        geom_vline(xintercept = 0, color = "gray90", linewidth = 0.5)

    if (incl_data) {
        p <- p +
            geom_point(data = pred_pca$obs_rot |>
                           mutate(PC1 = sign(pc_mults$pred[1]) * PC1,
                                  PC2 = sign(pc_mults$pred[2]) * PC2,
                                  PC3 = sign(pc_mults$pred[3]) * PC3),
                       aes(x =.data[[xPC_str]],
                           y = .data[[yPC_str]]),
                       alpha = 0.25, size = 1, shape = 1,
                       color = "gray80")
    }

    p <- p +
        geom_point(aes(x =.data[[xPC_str]],
                       y = .data[[yPC_str]],
                       fill = taxon),
                   size = 3, shape = 21) +
        geom_segment(data = pred_p_df,
                     aes(x = 0, xend = .data[[xPC_str]],
                         y = 0, yend = .data[[yPC_str]],
                         group = coef, color = coef),
                     arrow = arrow(length = unit(6, "pt")),
                     linewidth = 1) +
        scale_x_continuous(xPC_str, breaks = 3*-1:1) +
        scale_y_continuous(yPC_str, breaks = 3*-1:1) +
        coord_equal(xlim = c(-pc_axis_lim, pc_axis_lim),
                    ylim = c(-pc_axis_lim, pc_axis_lim)) +
        scale_color_manual(values = coef_pal, guide = "none") +
        scale_fill_manual(values = taxon_pal, guide = "none") +
        pca_theme

    return(p)
}





pca_plots <- list(pca_plotter(1, 2, TRUE),
                  pca_plotter(1, 3, TRUE),
                  pca_plotter(2, 3, TRUE))

wrap_plots(pca_plots)







# --------------*
# Variance partitioning ----
# --------------*


# variance partition of PC axes by predictors
var_part <- lapply(c("PC1","PC2","PC3"), function(x) {
    form = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(var = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(form, data = pred_pca$axes))[,2] |>
               (\(x) x[1:3]/sum(x[1:3]))())
}) |>
    bind_rows() |>
    spread(pc, cont) |>
    mutate(var = gsub("_z$", "", var))


# variance accounted for by predictors and PCs:
var_df <- var_part |>
    mutate(PC1 = PC1 * pred_pca$obs_exp$PC1[1],
           PC2 = PC2 * pred_pca$obs_exp$PC2[1],
           PC3 = PC3 * pred_pca$obs_exp$PC3[1],
           var = make_coef_fct(var)) |>
    pivot_longer(starts_with("PC"), names_to = "pc") |>
    mutate(pc = gsub("PC", "", pc) |> as.integer()) |>
    arrange(pc, var) |>
    mutate(cumvalue = cumsum(value),
           lag_cumvalue = lag(cumvalue, default = 0))




var_part_p <- var_df |>
    mutate(y = pc - 1,
           ymax = 0 - y * 0.5,
           ymin = ymax - 1) |>
    ggplot() +
    geom_rect(aes(xmin = lag_cumvalue, xmax = cumvalue,
                  ymax = ymax, ymin = ymin, fill = var)) +
    geom_text(data = var_df |>
                  group_by(pc) |>
                  summarize(value = (min(lag_cumvalue) + max(cumvalue)) / 2,
                            .groups = "drop") |>
                  add_row(pc = 4, value = (1 + sum(var_df$value)) / 2) |>
                  mutate(y = 0 - (pc-1) * 0.5,
                         lab = ifelse(pc < 4, paste0("PC", pc), "")),
              aes(value, y, label = lab),
              nudge_y = 0.3, size = 8 / 2.83465) +
    geom_point(data = tibble(var = var_df$var |> unique() |> sort(),
                             value = 0.95,
                             y = 0.25 - 0:2 * 0.5),
               aes(value, y, color = var),
               size = 2, shape = 15) +
    geom_text(data = tibble(var = var_df$var |> unique() |> sort(),
                            value = 0.98,
                            y = 0.25 - 0:2 * 0.5),
              aes(value, y, label = var, color = var),
              size = 8 / 2.83465, hjust = 0, vjust = 0.5) +
    scale_color_manual(values = coef_pal, guide = "none",
                       aesthetics = c("color", "fill")) +
    scale_x_continuous("Proportion of variance", breaks = seq(0, 1, 0.2)) +
    coord_cartesian(xlim = c(-0.1, 1.1),
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
# Combine them all into single figure
# --------------*

pca_p <- wrap_plots(c(taxon_pca_p, pred_pca_p, list(var_part_p))[c(1,4,2,5,3,6,7)]) +
    plot_layout(design = "AB\nCD\nEF\nGG") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16),
          plot.tag.position = c(-0.1, 1))
# This prevents tag for left panel from being further from plot than the others
# due to the large y-axis that's on the left:
for (i in c(1,3,5,7)) pca_p[[i]] <- pca_p[[i]]  & theme(plot.tag.position = c(0.05, 1))
pca_p[[7]] <- pca_p[[7]]  & theme(plot.tag.position = c(0.025, 1))

pca_p

# save_file(pca_p, "pca", width = 3, height = 6.5)


