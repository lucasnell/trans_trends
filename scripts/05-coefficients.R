
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






# =============================================================================*
# Individual estimates ----
# =============================================================================*



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
            mutate(coef = str_remove(coef, "_z$")) |>
            mutate(coef = make_coef_fct(coef))

        return(out)
    })()


#'
#' Creates subpanels for taxon-specific values for each predictor parameter
#'
coef_indiv_plotter <- function(.coef) {

    .coef <- match.arg(tolower(.coef), levels(re_draws$coef))

    .xlab <- list(time = "Time response", dist = "Distance response",
                  midges = "Midge response")[[.coef]]

    .dd <- re_draws |>
        filter(coef == .coef) |>
        mutate(med = map_dbl(iters, ~ median(.x)),
               lo = map_dbl(iters, ~ unname(quantile(.x, 0.16))),
               hi = map_dbl(iters, ~ unname(quantile(.x, 0.84))),
               lo2 = map_dbl(iters, ~ unname(quantile(.x, 0.025))),
               hi2 = map_dbl(iters, ~ unname(quantile(.x, 0.975)))) |>
        mutate(taxon = factor(taxon, levels = taxa_lvls, labels = taxa_labs),
               taxon = fct_reorder(taxon, desc(med)),
               taxon_int = as.integer(taxon))

    .plot <- .dd |>
        ggplot(aes(med, taxon_int)) +
        geom_vline(xintercept = 0, color = "gray50")+
        geom_violin(data = .dd |> select(taxon_int, taxon, iters) |> unnest(iters),
                    aes(x = iters, fill = taxon),
                    position = "identity", alpha = 0.25, color = NA) +
        geom_segment(aes(x = lo2, xend = hi2, yend = taxon_int, color = taxon),
                     linewidth = 0.5) +
        geom_linerange(aes(xmin = lo, xmax = hi, color = taxon),
                       linewidth = 1.5) +
        geom_point(aes(fill = taxon, shape = taxon), size = 3)+
        scale_y_continuous(breaks = 1:6, labels = levels(.dd$taxon)) +
        scale_x_continuous(.xlab, breaks = c(-1, 0, 1)) +
        scale_color_manual(NULL, values = taxa_pal) +
        scale_fill_manual(NULL, values = taxa_pal) +
        scale_shape_manual(NULL, values = taxa_shapes) +
        coord_cartesian(xlim = c(-1.45, 1.45), ylim = c(0.5, 7), expand = FALSE) +
        theme(axis.text.x = element_text(size = 8,
                                         margin = margin(0,0,0,t=2)),
              axis.title.x = element_text(size = 10,
                                          margin = margin(0,0,0,t=0)),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 8,
                                         margin = margin(0,0,0,r=3)),
              plot.margin = margin(t=12, r=0, b=0, l=0),
              legend.position = "none")

    return(.plot)
}



coef_p <- levels(re_draws$coef) |>
    map(coef_indiv_plotter) |>
    do.call(what = patchwork::wrap_plots) +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag.position = c(0, 1))

coef_p

# save_plot("coeffs-indiv", coef_p, w = 3.5, h = 6)






# =============================================================================*
# Mean and SD of estimates ----
# =============================================================================*




#' Process data and use one of extract posterior densities or uncertainty
#' intervals on either fixed effects or random effect SDs:
extract_effect_info <- function(model, effect_type, output) {

    stopifnot(inherits(model, "armmMod"))

    effect_type <- match.arg(tolower(effect_type), c("fixed", "random"))
    output <- match.arg(tolower(output), c("density", "ui"))

    rnd <- effect_type == "random"

    #' Extract uncertainty intervals (68% and 95%) for the fixed effects
    #' and for the random effect SDs
    extract_post_ui <- function(z) {
        ui68 <- unname(quantile(z[["value"]], c(0.16, 0.84)))
        ui95 <- unname(quantile(z[["value"]], c(0.025, 0.975)))
        tibble(coef = z[["coef"]][1],
               x = c(ui68[1], ui95[1]),
               xend = c(ui68[2], ui95[2]),
               type = factor(c(68, 95)))
    }

    out_fn <- switch(output,
                     density = identity,
                     ui = extract_post_ui)

    if (effect_type == "fixed") {
        out <- model$stan |>
            rstan::extract(pars = "alpha") |>
            do.call(what = cbind) |>
            as.data.frame() |>
            set_names(c("int", "midges", "time", "dist")) |>
            as_tibble() |>
            select(-int) |>  # intercept not necessary
            pivot_longer(everything(), names_to = "coef") |>
            mutate(coef = make_coef_fct(coef)) |>
            split(~ coef) |>
            map(out_fn) |>
            list_rbind()
    } else {
        out <- model$stan |>
            rstan::extract(pars = "sig_beta") |>
            do.call(what = cbind) |>
            as.data.frame() |>
            set_names(c("int_tax", "int_plot", "int_trans", "midges", "time", "dist")) |>
            as_tibble() |>
            select(midges, time, dist) |>
            pivot_longer(everything(), names_to = "coef") |>
            mutate(coef = make_coef_fct(coef)) |>
            split(~ coef) |>
            map(out_fn) |>
            list_rbind()
    }

    return(out)

}





slope_density_df <- extract_effect_info(model_fit, "fixed", "density")
slope_sd_density_df <- extract_effect_info(model_fit, "random", "density")

slope_ui_df <- extract_effect_info(model_fit, "fixed", "ui") |>
    mutate(y = -0.2 - 0.2 * (as.integer(coef)-1))
slope_sd_ui_df <- extract_effect_info(model_fit, "random", "ui") |>
    mutate(y = -0.2 - 0.2 * (as.integer(coef)-1))


coef_mean_p <- slope_density_df |>
    ggplot(aes(value)) +
    geom_vline(xintercept = 0, color = "gray50")+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_segment(data = slope_ui_df,
              aes(x = x, y = y, xend = xend, yend = y, color = coef, linewidth = type)) +
    stat_density(aes(fill = coef, color = coef),
                 alpha = 0.5, linewidth = 0.5, position = "identity") +
    geom_text(data = tibble(coef = factor(c("time", "dist", "midges")),
                            value = rep(-1.3, 3),
                            y = 4.6 - (0:2 * 0.4),
                            lab = c("time", "distance", "midges")),
              aes(y = y, label = lab, color = coef), hjust = 0, vjust = 1,
              fontface = "bold", size = 10 / 2.8) +
    scale_fill_manual(values = coef_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    scale_linewidth_manual(values = c(1.5, 0.5), guide = "none") +
    ylab("Density") +
    scale_x_continuous("Mean for among-taxa response", breaks = -1:1) +
    coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(-0.6, 4.6))


coef_sd_p <- slope_sd_density_df |>
    ggplot(aes(value)) +
    geom_vline(xintercept = 0, color = "gray50") +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_segment(data = slope_sd_ui_df, aes(x = x, y = y, xend = xend, yend = y,
                                            color = coef, linewidth = type)) +
    stat_density(aes(fill = coef, color = coef), bounds = c(0, Inf),
                 alpha = 0.5, linewidth = 0.5, position = "identity") +
    scale_fill_manual(values = coef_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    scale_linewidth_manual(values = c(1.5, 0.5), guide = "none") +
    ylab("Density") +
    scale_x_continuous("SD for among-taxa response") +
    coord_cartesian(xlim = c(0, 1.5), ylim = c(-0.6, 4.6))


coef_among_p <- (coef_mean_p | coef_sd_p) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag.position = c(0.06, 1))

# coef_among_p

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
