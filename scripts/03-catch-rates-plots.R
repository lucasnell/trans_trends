
#'
#' This file creates figure showing catch rates through time and across
#' distance from lake.
#'

# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")


data_df <- read_rds(data_rds)

type_pal <- c("black", viridis(100)[80])





# Make summary data for plotting means of midge abundance and catch rates
make_summ_df <- function(grp_var, taxa = NA, .trans_fun = NA) {
    if (is.na(taxa)) taxa <- unique(data_df$taxon_pretty)
    else stopifnot(all(taxa %in% unique(data_df$taxon_pretty)))
    sdf <- data_df |>
        filter(taxon_pretty %in% taxa) |>
        select(taxon_pretty, {{grp_var}}, y, y_noz, midges_z) |>
        group_by(taxon_pretty, {{grp_var}}) |>
        summarize(across(y:midges_z, mean), .groups = "drop") |>
        pivot_longer(y_noz:midges_z, names_to = "type", values_to = "y_noz") |>
        # Next line makes it so that y is midge abundance when type == "midges_z"
        # but log arthropod abundance when not:
        mutate(y = ifelse(type == "y_noz", y, y_noz)) |>
        # Make type column pretty:
        mutate(type = factor(type,
                             levels = c("y_noz", "midges_z"),
                             labels = c("Mean by year/distance",
                                        "Midge abundance"))) |>
        select(taxon_pretty, {{grp_var}}, type, y, y_noz)

    if (isTRUE(is.function(.trans_fun))) sdf <- sdf |>
            mutate(y = ifelse(type == "Midge abundance", .trans_fun(y), y),
                   y_noz = ifelse(type == "Midge abundance", .trans_fun(y_noz), y_noz))

    return(sdf)
}



with_options(pillar.sigfig = 7,
             make_summ_df(year) |>
                 filter(type == "Midge abundance") |>
                 summarize(lo = min(y_noz), hi = max(y_noz), diff = hi - lo) |>
                 print())
# # A tibble: 1 × 3
#          lo        hi     diff
#       <dbl>     <dbl>    <dbl>
# 1 -1.648319 0.7983331 2.446652


with_options(pillar.sigfig = 7,
             data_df |>
                 group_by(taxon, taxon_pretty) |>
                 summarize(y_noz = max(y_noz), .groups = "drop") |>
                 arrange(desc(y_noz)) |>
                 print())
# # A tibble: 6 × 3
#   taxon taxon_pretty      y_noz
#   <fct> <fct>             <dbl>
# 1 lyco  Wolf spiders   2.568725
# 2 opil  Harvestmen     1.683546
# 3 sheet Sheet weavers  1.449808
# 4 cara  Ground beetles 1.329136
# 5 stap  Rove beetles   1.170071
# 6 gnap  Ground spiders 1.130506




# # by time
# # time_p <-
# data_df |>
#     ggplot(aes(year, y_noz))+
#     facet_wrap(~taxon_pretty, nrow = 4,
#                # scales = "free_y")+
#                scales = "fixed")+
#     geom_hline(yintercept = 0, color = "gray50")+
#     # geom_line(aes(group = plot), linewidth = 0.2, color = "gray70") +
#     geom_line(aes(group = plot), linewidth = 0.2, color = "gray30", alpha = 0.25) +
#     geom_point(size = 1, alpha = 0.1, shape = 19, color = "gray30") +
#     # geom_point(data = time_summ_df, aes(color = type), size = 2, shape = 1)+
#     # geom_line(data = time_summ_df,
#     geom_line(data = make_summ_df(year, .trans_fun = NA) |> filter(grepl("^Mean", type)),
#               aes(color = type), linewidth = 0.75) +
#     scale_color_manual(NULL, values = c("black", type_pal[1:2])) +
#     scale_x_continuous("Year", breaks = c(2008,2012,2016),
#                        limits = c(2007,2020)) +
#     coord_cartesian(ylim = c(0, 1.130)) +
#     scale_y_continuous(expression("Catch rate (days"*{}^{-1}*")"),
#                        # limits = c(0, 2.6),
#                        # breaks = log1p(c(0, 2, 8)),
#                        # labels = c(0, 2, 8),
#                        sec.axis = sec_axis(inv_trans_fun2,
#                                            "Midge abundance")) +  # ,
#                                            # breaks = log1p(20 * 40^(0:2)),
#                                            # labels = c("20", "8e2", "32e3"))) +
#     theme(legend.position = "top",
#           axis.text.y.right = element_text(color = type_pal[3]),
#           axis.ticks.y.right = element_line(color = type_pal[3]),
#           # axis.title.y.right = element_text(color = type_pal[3])
#           axis.title.y.right = element_blank()
#           # axis.text.y.right = element_blank(),
#           # axis.title.y.right = element_text(color = type_pal[3])
#           )
#


cr_plotter <- function(.df, x_var, grp_var) {
    x_str <- as_name(quo(expr = {{x_var}}))
    txn <- .df$taxon_pretty |> paste() |> unique()
    stopifnot(length(txn) == 1)
    # Index for which y limits and breaks to use for this species:
    y_idx <- list("Wolf spiders" = 1L,
                  "Harvestmen" = 2L,
                  "Sheet weavers" = 2L,
                  "Ground beetles" = 3L,
                  "Rove beetles" = 3L,
                  "Ground spiders" = 3L)[[txn]]
    # now set y limits and labels (breaks are log1p(labels)):
    y_max <- c(2.569, 1.684, 1.329)[[y_idx]]
    y_labs <- list(c(0, 2, 8), c(0, 1, 4), c(0, 1, 2.5))[[y_idx]]
    # Similar but for transforming between catch rates and midge abundances:
    b0 <- y_max / 2.45 * 0.8
    trans_fun <- function(x) (x + 1.65) * b0
    inv_trans_fun <- function(y) y / b0 - 1.65
    x_scale <- NULL
    if (x_str == "year") {
        x_scale <- scale_x_continuous("Year", breaks = 2010 + 4 * 0:2,
                                      limits = c(2007,2020))
    } else if (x_str == "dist") {
        x_scale <- scale_x_continuous("Distance (m)",
                                      breaks = log(c(5, 50, 150, 500)),
                                      labels = c(5, 50, 150, 500))
    }
    .df |>
        ggplot(aes({{ x_var }}, y_noz))+
        facet_wrap(~taxon_pretty, nrow = 4, scales = "free_y")+
        geom_hline(yintercept = 0, color = "gray50")+
        geom_line(aes(group = {{ grp_var }}), linewidth = 0.2, color = "gray30",
                  alpha = 0.25) +
        geom_point(size = 1, alpha = 0.1, shape = 19, color = "gray30") +
        geom_line(data = make_summ_df({{ x_var }}, txn, .trans_fun = trans_fun),
                  aes(color = type), linewidth = 0.75) +
        scale_color_manual(NULL, values = type_pal) +
        coord_cartesian(ylim = c(0, y_max)) +
        x_scale +
        scale_y_continuous(expression("Catch rate (days"*{}^{-1}*")"),
                           breaks = log1p(y_labs),
                           labels = y_labs,
                           sec.axis = sec_axis(inv_trans_fun,
                                               "SD change in midge abundance",
                                               breaks = -1:1)) +
        theme(axis.title.y.right = element_text(color = type_pal[2]),
              axis.text.y.right = element_text(color = type_pal[2]),
              axis.ticks.y.right = element_line(color = type_pal[2]))
}




time_plots <- data_df |>
    split(~ taxon_pretty) |>
    map(cr_plotter, x_var = year, grp_var = plot)


# do.call(wrap_plots, time_plots) +
#     plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") +
#     plot_annotation(title = "a", theme = theme(plot.title = element_text(
#         face = "bold", size = 14))) &
#     theme(legend.position = "top")


dist_plots <- data_df |>
    mutate(id = interaction(trans, year)) |>
    split(~ taxon_pretty) |>
    map(cr_plotter, x_var = dist, grp_var = id)


# do.call(wrap_plots, dist_plots) +
#     plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") +
#     plot_annotation(title = "b", theme = theme(plot.title = element_text(
#         face = "bold", size = 14))) &
#     theme(legend.position = "top")


time_p <- do.call(wrap_plots, time_plots) +
    plot_layout(ncol = 2, axis_titles = "collect") &
    theme(axis.title.y.right = element_blank())
dist_p <- do.call(wrap_plots, dist_plots) +
    plot_layout(ncol = 2, axis_titles = "collect") &
    theme(axis.title.y.left = element_blank())

obs_data_p <- guide_area() / (time_p | dist_p) +
    plot_layout(nrow = 2, guides = "collect", heights = unit(c(1, 1), c("cm", "null"))) +
    plot_annotation(tag_levels = list(c("a", rep("", 5), "b", rep("", 5)))) &
    theme(legend.position = "top")


# save_plot("catch-rates-year-dist", obs_data_p, w = 8, h = 5)


