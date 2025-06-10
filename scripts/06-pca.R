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
    (\(x) left_join(x, map(x$var, get_fit_info) |> list_rbind(), by = "var"))() |>
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
pred_vec <- map(c("time", "dist", "midges"),
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
                }) |>
    list_rbind()


# =============================================================================*
# Plot prelims ----
# =============================================================================*

pc_axis_lim <- 4.1
#' #' We're scaling some of the projections. This is to keep that consistent.
pc_mults <- list(pred = 6 * c(-1, 1, 1), taxon = 5 * c(-1, 1, 1))

stopifnot(identical(sign(pc_mults$pred), sign(pc_mults$taxon)))





# --------------*
# Plots of taxon response vectors and model predictors ----
# --------------*

pca_plotter <- function(.xPC, .yPC, incl_data = FALSE) {


    stopifnot(isTRUE(length(.xPC) == 1) && isTRUE(length(.yPC) == 1))

    xPC_str <- paste0("PC", .xPC)
    yPC_str <- paste0("PC", .yPC)

    taxon_p_df <- pred_pca$taxon_vec |>
        arrange(taxon) |>
        mutate(PC1 = pc_mults$taxon[1] * PC1,
               PC2 = pc_mults$taxon[2] * PC2,
               PC3 = pc_mults$taxon[3] * PC3) |>
        arrange(desc(taxon)) |>
        mutate(taxon2 = factor(taxon, levels = taxa_lvls |> rev(),
                               labels = taxa_labs |> rev() |> tolower() |>
                                   str_replace_all(" ", "\n")),
               taxon_pretty = factor(paste(taxon), levels = taxa_lvls,
                                     labels = taxa_labs))
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
                       fill = taxon_pretty,
                       shape = taxon_pretty),
                   size = 3) +
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
        scale_fill_manual(NULL, values = taxa_pal) +
        scale_shape_manual(NULL, values = rep(c(21, 22, 24), 2)) +
        theme(axis.text.y = element_text(size = 8,
                                         margin = margin(0,0,0,r=2)),
              axis.title.y = element_text(size = 10,
                                          margin = margin(0,0,0,0)),
              axis.text.x = element_text(size = 8,
                                         margin = margin(0,0,0,t=2)),
              axis.title.x = element_text(size = 10,
                                          margin = margin(0,0,0,b=8)),
              legend.position = "none",
              plot.margin = margin(0,0,0,0))

    return(p)
}





pca_plots <- list(pca_plotter(1, 2),
                  pca_plotter(1, 3),
                  pca_plotter(2, 3))

# wrap_plots(pca_plots) + plot_layout(guides = "collect") & theme(legend.position = "right")







# --------------*
# Variance partitioning ----
# --------------*


# variance accounted for by predictors and PCs:
var_df <- lapply(c("PC1","PC2","PC3"), function(x) {
    form = as.formula(paste(x, "~ midges_z + time_z + dist_z"))
    tibble(coef = c("midges_z", "time_z", "dist_z"),
           pc = x,
           cont = anova(lm(form, data = pred_pca$axes))[,2] |>
               (\(x) x[1:3]/sum(x[1:3]))())
}) |>
    bind_rows() |>
    pivot_wider(names_from = pc, values_from = cont) |>
    mutate(coef = gsub("_z$", "", coef)) |>
    mutate(PC1 = PC1 * pred_pca$obs_exp$PC1[1],
           PC2 = PC2 * pred_pca$obs_exp$PC2[1],
           PC3 = PC3 * pred_pca$obs_exp$PC3[1],
           coef = factor(coef, levels = c("time", "dist", "midges"),
                         labels = c("time", "distance", "midges"))) |>
    pivot_longer(starts_with("PC"), names_to = "pc") |>
    mutate(pc = gsub("PC", "", pc) |> as.integer()) |>
    arrange(pc, coef) |>
    mutate(cumvalue = cumsum(value),
           lag_cumvalue = lag(cumvalue, default = 0))




var_part_p <- var_df |>
    mutate(y = pc - 1,
           ymax = 0 - y * 0.5,
           ymin = ymax - 1) |>
    ggplot() +
    geom_rect(aes(xmin = lag_cumvalue, xmax = cumvalue,
                  ymax = ymax, ymin = ymin, fill = coef)) +
    geom_text(data = var_df |>
                  group_by(pc) |>
                  summarize(value = (min(lag_cumvalue) + max(cumvalue)) / 2,
                            .groups = "drop") |>
                  mutate(y = 0 - (pc-1) * 0.5,
                         lab = paste0("PC", pc)),
              aes(value, y, label = lab),
              nudge_y = 0.2, size = 10 / 2.83465) +
    scale_color_manual(NULL, values = coef_pal, aesthetics = c("color", "fill")) +
    scale_x_continuous("Proportion of variance", breaks = seq(0, 1, 0.2)) +
    coord_cartesian(xlim = c(-0.02, 1.02),
                    ylim = c(-2.1, 0.7),
                    expand = FALSE) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8, margin = margin(0,0,0,t=2)),
          axis.title.x = element_text(margin = margin(0,0,0,b=6)),
          legend.position = "none",
          plot.margin = margin(0,0,0,0))





# --------------*
# Combine them all into single figure
# --------------*

# # This version has a legend instead of labels in the plot itself.
# # I added labels in Illustrator, but this one shows all the relevant info:
# simple_pca_p <- wrap_plots(c(pca_plots, list(var_part_p))) +
#     plot_layout(design = "ABC\nDDD", heights = c(1, 0.5), guides = "collect") +
#     plot_annotation(tag_levels = "a") &
#     theme(plot.tag = element_text(size = 16),
#           plot.tag.position = c(0.015, 1),
#           legend.position = "right")
# # This prevents tag for bottom panel from being further from plot
# # than the others due to its lack of a y-axis:
# simple_pca_p[[4]] <- simple_pca_p[[4]] & theme(plot.tag.position = c(0.005, 1))
# simple_pca_p



#' For the version in the paper, I saved each panel separately, then stitched
#' them together and added labels in Illustrator:
for (i in seq_along(pca_plots)) {
    save_plot(sprintf("pca-%s", letters[i]), pca_plots[[i]], w = 2, h = 2)
}; rm(i)

save_plot("pca-var", var_part_p, w = 6, h = 2)
