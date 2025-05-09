
#'
#' This file creates figure showing data through time.
#'

# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")


data_df <- read_rds(data_rds$nolag)


# Transforming from scale of midge abundances to catch rates:
trans_fun <- function(x) (x - 3) * 0.25
# ... and back:
inv_trans_fun <- function(x) x / 0.25 + 3


time_summ_df <- data_df |>
    select(taxon_plot, year, y_noz, small_midges, big_midges) |>
    mutate(across(small_midges:big_midges, log1p)) |>
    pivot_longer(y_noz:big_midges, names_to = "type", values_to = "y_noz") |>
    mutate(type = factor(type,
                         levels = c("y_noz", "small_midges", "big_midges"),
                         labels = c("Mean by year/distance",
                                    "Small midges",
                                    "Large midges"))) |>
    group_by(taxon_plot, year, type) |>
    summarize(y_noz = mean(y_noz), .groups = "drop") |>
    # group_by(type) |> summarize(min = min(y_noz), max = max(y_noz))
    mutate(y_noz = ifelse(type == "Mean by year/distance", y_noz,
                          trans_fun(y_noz)))
dist_summ_df <- data_df |>
    select(taxon_plot, distance, y_noz, small_midges, big_midges) |>
    mutate(across(small_midges:big_midges, log1p)) |>
    pivot_longer(y_noz:big_midges, names_to = "type", values_to = "y_noz") |>
    mutate(type = factor(type,
                         levels = c("y_noz", "small_midges", "big_midges"),
                         labels = c("Mean by year/distance",
                                    "Small midges",
                                    "Large midges"))) |>
    group_by(taxon_plot, distance, type) |>
    summarize(y_noz = mean(y_noz), .groups = "drop") |>
    # group_by(type) |> summarize(min = min(y_noz), max = max(y_noz))
    mutate(y_noz = ifelse(type == "Mean by year/distance", y_noz,
                          trans_fun(y_noz)),
           dist = log(distance))



type_pal <- viridis(100)[c(70, 90, 80)]


# by time
time_p <- data_df |>
    ggplot(aes(year, y_noz))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    # geom_line(aes(group = plot), linewidth = 0.2, color = "gray70") +
    geom_line(aes(group = plot), linewidth = 0.2, color = "gray30", alpha = 0.25) +
    geom_point(size = 1, alpha = 0.1, shape = 19, color = "gray30") +
    # geom_point(data = time_summ_df, aes(color = type), size = 2, shape = 1)+
    geom_line(data = time_summ_df, aes(color = type), linewidth = 0.75) +
    scale_color_manual(NULL, values = c("black", type_pal[1:2])) +
    scale_x_continuous("Year", breaks = c(2008,2012,2016),
                       limits = c(2007,2020)) +
    scale_y_continuous(expression("Catch rate (days"*{}^{-1}*")"),
                       limits = c(0, 2.6),
                       breaks = log1p(c(0, 2, 8)),
                       labels = c(0, 2, 8),
                       sec.axis = sec_axis(inv_trans_fun,
                                           "Midge abundance",
                                           breaks = log1p(20 * 40^(0:2)),
                                           labels = c("20", "8e2", "32e3"))) +
    theme(axis.title.y.right = element_blank(),
          # axis.text.y.right = element_blank(),
          # axis.title.y.right = element_text(color = type_pal[3]),
          axis.text.y.right = element_text(color = type_pal[3]),
          axis.ticks.y.right = element_line(color = type_pal[3]),
          legend.position = "top")



# by distance
dist_p <- data_df |>
    mutate(id = interaction(trans, year)) |>
    ggplot(aes(dist, y_noz))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    # geom_point(aes(group = plot), size = 1, alpha = 0.5, shape = 1)+
    geom_line(aes(group = id), linewidth = 0.2, color = "gray30", alpha = 0.25) +
    geom_point(size = 1, alpha = 0.1, shape = 19, color = "gray30") +
    # geom_point(data = dist_summ_df, aes(color = type), size = 2, shape = 1)+
    geom_line(data = dist_summ_df, aes(color = type), linewidth = 0.75)+
    scale_color_manual(NULL, values = c("black", type_pal[1:2])) +
    # scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    scale_x_continuous("Distance (m)", breaks = log(c(5, 50, 150, 500)),
                       labels = c(5, 50, 150, 500))+
    scale_y_continuous(expression("Catch rate (days"*{}^{-1}*")"),
                       # limits = c(0, 2.6),
                       breaks = log1p(c(0, 2, 8)),
                       labels = c(0, 2, 8),
                       sec.axis = sec_axis(inv_trans_fun,
                                           "Midge abundance",
                                           breaks = log1p(20 * 40^(0:2)),
                                           labels = c("20", "8e2", "32e3"))) +
    theme(axis.title.y.right = element_text(color = type_pal[3]),
          axis.title.y.left = element_blank(),
          # axis.text.y.left = element_blank(),
          axis.text.y.right = element_text(color = type_pal[3]),
          axis.ticks.y.right = element_line(color = type_pal[3]),
          legend.position = "top")



obs_data_p <- guide_area() / (time_p + dist_p) +
    plot_layout(nrow = 2, guides = "collect", heights = unit(c(1, 1), c("cm", "null"))) +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "top")

# save_plot("data-timeseries", obs_data_p, w = 6, h = 4)


