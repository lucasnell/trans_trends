
#'
#' This file creates figure showing catch rates through time and across
#' distance from lake.
#'

# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")


data_df <- read_rds(data_rds)

type_pal <- c("black", "goldenrod2")

# Mean and standard deviation of log-transformed midge catch rate:
midges_mu <- with(read_rds(data_rds), mean(log(midges / midge_days)))
midges_sd <- with(read_rds(data_rds), sd(log(midges / midge_days)))





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
    b0 <- -midges_mu / midges_sd
    # b1 <- y_max / 2.45 * 0.8
    b1 <- y_max / 3.4
    trans_fun <- function(x) (x - b0) * b1
    inv_trans_fun <- function(y) {
        cr_lz <- y / b1 + b0 # to z_trans(log(midges / midge_days))
        cr_l <- cr_lz * midges_sd + midges_mu # to log(midges / midge_days)
        cr <- exp(cr_l) # to midges / midge_days
        return(cr)
    }
    x_scale <- NULL
    if (x_str == "year") {
        x_scale <- scale_x_continuous("Year", breaks = 2010 + 4 * 0:2,
                                      limits = c(2007,2020))
    } else if (x_str == "dist") {
        x_scale <- scale_x_continuous("Distance (m)",
                                      breaks = log(c(5, 50, 150, 500)),
                                      labels = c(5, 50, 150, 500))
    }
    p <- .df |>
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
        scale_y_continuous(expression("Catch rate (ind. day"*{}^{-1}*")"),
                           breaks = log1p(y_labs),
                           labels = y_labs,
                           sec.axis = sec_axis(inv_trans_fun,
                                               expression("Midge catch rate (ind. day"*{}^{-1}*")"),
                                               breaks = 10^(0:2))) +
        theme(axis.title.y.right = element_text(color = type_pal[2]),
              axis.text.y.right = element_text(color = type_pal[2]),
              axis.ticks.y.right = element_line(color = type_pal[2]),
              plot.margin = margin(r=3, l=3, t=-9, b=3))
    if (txn %in% levels(data_df$taxon_pretty)[c(1,3,5)]) {
        p <- p + theme(axis.text.y.right = element_blank())
    }

    return(p)
}




time_p <- data_df |>
    split(~ taxon_pretty) |>
    map(cr_plotter, x_var = year, grp_var = plot) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 2, axis_titles = "collect") &
    theme(axis.title.y.right = element_blank(),
          plot.tag.position = c(0.15, 1))

dist_p <- data_df |>
    mutate(id = interaction(trans, year)) |>
    split(~ taxon_pretty) |>
    map(cr_plotter, x_var = dist, grp_var = id) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 2, axis_titles = "collect") &
    theme(axis.title.y.left = element_blank(),
          plot.tag.position = c(0, 1))

# time_p +
#     plot_layout(guides = "collect") &
#     theme(legend.position = "top")

# dist_p +
#     plot_layout(guides = "collect") &
#     theme(legend.position = "top")


obs_data_p <- guide_area() + time_p + dist_p +
    plot_layout(guides = "collect", heights = unit(c(1, 1), c("cm", "null")),
                widths = c(1, 0.05, 1),
                design = "AAA\nB#C") +
    plot_annotation(tag_levels = list(c("a", rep("", 5), "b", rep("", 5)))) &
    theme(legend.position = "top")

# obs_data_p




# save_plot("catch-rates-year-dist", obs_data_p, w = 8, h = 5)


