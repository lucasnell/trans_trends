
#'
#' This file creates plots of coefficient estimates from the model fit.
#'



# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")


# Standard deviation of log-transformed midge catch rate:
midges_sd <- with(read_rds(data_rds), sd(log(midges / midge_days)))


model_fit <- read_rds(model_rds)





#' Extract posterior densities for the fixed effects
#' and for the random effect SDs
#' Column `ui` indicates whether the densities are for the 68%
#' uncertainty intervals.
extract_posteriors <- function(z) {
    .ui <- unname(quantile(z[["value"]], c(0.16, 0.84)))
    X <- density(z[["value"]], n = 2048)
    Z <- density(z[["value"]], from = .ui[1], to = .ui[2])
    tibble(coef = z[["coef"]][1],
           x = c(X$x, c(Z$x[1], Z$x, tail(Z$x, 1))),
           y = c(X$y, c(0, Z$y, 0)),
           ui = c(rep(FALSE, length(X$x)), rep(TRUE, length(Z$x)+2)))
}

slope_density_df <- model_fit$stan |>
    rstan::extract(pars = "alpha") |>
    do.call(what = cbind) |>
    as.data.frame() |>
    set_names(c("int", "midges", "time", "dist")) |>
    as_tibble() |>
    select(-int) |>  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") |>
    mutate(coef = make_coef_fct(coef)) |>
    split(~ coef) |>
    map_dfr(extract_posteriors)

slope_sd_density_df <- model_fit$stan |>
    rstan::extract(pars = "sig_beta") |>
    do.call(what = cbind) |>
    as.data.frame() |>
    set_names(c("int_tax", "int_plot", "int_trans", "midges", "time", "dist")) |>
    as_tibble() |>
    select(midges, time, dist) |>
    pivot_longer(everything(), names_to = "coef") |>
    mutate(coef = make_coef_fct(coef)) |>
    split(~ coef) |>
    map_dfr(extract_posteriors)




#'
#' Posterior draws for random effect estimates
#' (Using an inline function to avoid cluttering environment)
#'
re_draws <- model_fit |>
    (function(object) {

        stopifnot(inherits(object, "armmMod"))
        fef <- rstan::extract(object$stan, "alpha")[[1]]
        colnames(fef) <- colnames(object$stan_data$x)

        sigma_names <- names(object$stan)[grepl("^sig_beta\\[", names(object$stan))]
        z_names <- names(object$stan)[grepl("^z\\[", names(object$stan))]

        S <- rstan::extract(object$stan, sigma_names)
        S <- lapply(1:length(S),
                    function(i) {
                        matrix(as.numeric(S[[i]]), length(S[[i]]),
                               object$stan_data$lev_per_g[i])
                    })
        S <- do.call(cbind, S)
        Z <- do.call(cbind, rstan::extract(object$stan, z_names))
        R <- S * Z

        # Combining random and fixed effects:
        E <- R
        for (i in 1:nrow(object$rnd_lvl_names)) {
            .f <- fef[,object$rnd_lvl_names$Name[i]]
            E[,i] <- E[,i] + .f
        }

        out <- object$rnd_lvl_names |>
            as_tibble() |>
            mutate(iters = map(1:n(), ~ E[,.x])) |>
            filter(Name != "(Intercept)", Groups == "taxon") |>
            select(-Groups) |>
            rename(taxon = Level, coef = Name) |>
            mutate(coef = str_remove(coef, "_z$"))

        return(out)
    })()






# =============================================================================*
# Individual estimates ----
# =============================================================================*

#'
#' Creates subpanels for taxon-specific values for each predictor parameter
#'
coef_indiv_plotter <- function(.coef) {

    .coef <- match.arg(tolower(.coef), levels(slope_density_df$coef))

    .ylab <- list(time = "Time response", dist = "Distance response",
                  midges = "Midge response")[[.coef]]

    .plot <- re_draws |>
        filter(coef == .coef) |>
        mutate(taxon = factor(taxon, levels = taxa_lvls) |>
                   as.integer()) |>
        mutate(med = map_dbl(iters, ~ median(.x)),
               lo = map_dbl(iters, ~ unname(quantile(.x, 0.16))),
               hi = map_dbl(iters, ~ unname(quantile(.x, 0.84)))) |>
        select(-iters) |>
        ggplot(aes(med, taxon)) +
        geom_vline(xintercept = 0, color = "gray50")+
        geom_point(size = 1.5, color = "black")+
        geom_linerange(aes(xmin = lo, xmax = hi), color = "black") +
        scale_y_continuous(breaks = 1:6, labels = taxa_labs) +
        scale_x_continuous(.ylab, breaks = c(-1, 0, 1)) +
        coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(0.5, 7),
                        expand = FALSE) +
        theme(axis.text.x = element_text(size = 8,
                                         margin = margin(0,0,0,t=2)),
              axis.title.x = element_text(size = 10,
                                          margin = margin(0,0,0,t=0)),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 8,
                                         margin = margin(0,0,0,r=3)),
              plot.margin = margin(t=0, r=0, b=0, l=24))

    return(.plot)
}





coef_p <- levels(slope_density_df$coef) |>
    map(coef_indiv_plotter) |>
    do.call(what = patchwork::wrap_plots) +
    plot_layout(nrow = 1, axes = "collect_y") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag.position = c(-0.06, 1))
# This prevents tag for left panel from being way higher than the others
# due to the large y-axis that's on the left:
coef_p[[1]] <- coef_p[[1]] & theme(plot.tag.position = c(0.29, 1))

# coef_p

# save_plot("coeffs-indiv", coef_p, w = 6.5, h = 3)






# =============================================================================*
# Mean and SD of estimates ----
# =============================================================================*

coef_mean_p <- slope_density_df |>
    filter(!ui) |>
    ggplot(aes(x, y)) +
    geom_vline(xintercept = 0, color = "gray50")+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_path(data = slope_density_df |> filter(ui) |>
                  group_by(coef) |> filter(x %in% range(x)) |>
                  ungroup() |> distinct(coef, x) |>
                  mutate(y = -0.2 - ifelse(coef == "dist", 0.1, 0)),
              aes(color = coef), linewidth = 1) +
    geom_polygon(aes(fill = coef, color = coef), alpha = 0.5, linewidth = 0.5) +
    scale_fill_manual(values = coef_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    ylab("Density") +
    scale_x_continuous("Mean for among-taxa response", breaks = -1:1) +
    coord_cartesian(xlim = c(-1.3, 1.3))


coef_sd_p <- slope_sd_density_df |>
    filter(!ui) |>
    ggplot(aes(x, y)) +
    geom_vline(xintercept = 0, color = "gray50")+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_path(data = slope_sd_density_df |> filter(ui) |>
                  group_by(coef) |> filter(x %in% range(x)) |>
                  ungroup() |> distinct(coef, x) |>
                  mutate(y = -0.2 - ifelse(coef == "time", 0.1, 0)),
              aes(color = coef), linewidth = 1) +
    geom_polygon(aes(fill = coef, color = coef), alpha = 0.5, linewidth = 0.5) +
    geom_text(data = tibble(coef = factor(c("time", "dist", "midges")),
                            x = c(0.4, 0.8, 0.2),
                            y = c(3,   1.5, 4.2),
                            lab = c("time", "distance", "midges")),
              aes(label = lab, color = coef), hjust = 0,
              fontface = "bold", size = 10 / 2.8) +
    scale_fill_manual(values = coef_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    ylab("Density") +
    scale_x_continuous("SD for among-taxa response") +
    coord_cartesian(xlim = c(0, 1.5))


coef_among_p <- (coef_mean_p | coef_sd_p) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag.position = c(0.06, 1))

# save_plot("coefs-mean-sd", coef_among_p, w = 5, h = 3)






# Autoregressive parameters ----
# =============================================================================*

rstan::extract(model_fit$stan, "phi")[[1]] |>
    as.data.frame() |>
    set_names(gsub("^taxon", "", model_fit$ar_names)) |>
    as_tibble() |>
    pivot_longer(everything(), names_to = "taxon") |>
    mutate(taxon = factor(paste(taxon), levels = taxa_lvls, labels = taxa_labs)) |>
    group_by(taxon) |>
    summarize(med = median(value),
              lo = quantile(value, 0.16),
              hi = quantile(value, 0.84))

# # A tibble: 6 × 4
#   taxon               med       lo       hi
#   <fct>             <dbl>    <dbl>    <dbl>
# 1 Ground beetles 0.059817 0.026760 0.10946
# 2 Rove beetles   0.049684 0.021062 0.094614
# 3 Harvestmen     0.033496 0.014305 0.065400
# 4 Ground spiders 0.11583  0.051896 0.19796
# 5 Sheet weavers  0.024776 0.010810 0.046281
# 6 Wolf spiders   0.032415 0.014248 0.061300



# =============================================================================*
# Effects of 10-fold midge increase ----
# =============================================================================*

re_draws |>
    filter(coef == "midges") |>
    mutate(taxon = factor(paste(taxon), levels = taxa_lvls, labels = taxa_labs)) |>
    mutate(midge_effect = 10^(map_dbl(iters, ~ median(.x)) / midges_sd)) |>
    select(taxon, midge_effect)

# # A tibble: 6 × 2
#   taxon          midge_effect
#   <fct>                 <dbl>
# 1 Ground beetles       1.3054
# 2 Rove beetles         1.3058
# 3 Harvestmen           1.2770
# 4 Ground spiders       1.7054
# 5 Sheet weavers        1.1188
# 6 Wolf spiders         1.3630
