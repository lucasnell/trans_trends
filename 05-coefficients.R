
#'
#' This file creates plots of coefficient estimates from the model fit.
#'



# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")

model_fit <- read_rds(model_rds$nolag)

# Order levels of coefficient factor:
make_coef_fct <- function(.coef) {
    factor(paste(.coef),
           levels = c("time", "dist", "small_midges", "big_midges"))
}




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
    set_names(c("int", "small_midges", "big_midges", "time", "dist")) |>
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
    set_names(c("int_tax", "int_plot", "int_trans", "small_midges",
                "big_midges", "time", "dist")) |>
    as_tibble() |>
    select(ends_with("midges"), time, dist) |>
    pivot_longer(everything(), names_to = "coef") |>
    mutate(coef = make_coef_fct(coef)) |>
    split(~ coef) |>
    map_dfr(extract_posteriors)





#' Get the posterior draws for random effect and autoregressive estimates
#'
get_re_au_draws <- function(object) {

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

    rnd_draws <- object$rnd_lvl_names |>
        as_tibble() |>
        mutate(iters = map(1:n(), ~ E[,.x]))


    autoreg_draws <- rstan::extract(model_fit$stan, "phi")[[1]] |>
        as.data.frame() |>
        set_names(gsub("^taxon", "", model_fit$ar_names)) |>
        as_tibble() |>
        pivot_longer(everything(), names_to = "Level") |>
        group_by(Level) |>
        summarize(iters = list(value)) |>
        mutate(Groups = "taxon", Name = "autoreg") |>
        select(Name, Groups, Level, iters)

    return(bind_rows(rnd_draws, autoreg_draws))
}


re_au_draws <- get_re_au_draws(model_fit) |>
    filter(Name != "(Intercept)", Groups == "taxon") |>
    select(-Groups) |>
    rename(taxon = Level, coef = Name) |>
    mutate(coef = str_remove(coef, "_z$"))








# >> LEFT OFF #2 ----
#' Add autoregressive parameter into this function below so that the axes
#' are standardized and can be combined easily

# Creates subpanels for each predictor and autoregressive parameter
slope_p_fun <- function(.coef) {

    .coef <- match.arg(tolower(.coef), c("autoreg", levels(slope_density_df$coef)))

    .ylab <- list(time = "Time response", dist = "Distance response",
         small_midges = "Small midge response",
         big_midges = "Large midge response",
         autoreg = "AR parameter")[[.coef]]

    .plot <- re_au_draws |>
        filter(coef == .coef) |>
        mutate(taxon = factor(taxon, levels = taxa_lvls) |>
                   as.integer()) |>
        mutate(med = map_dbl(iters, ~ median(.x)),
               lo = map_dbl(iters, ~ unname(quantile(.x, 0.16))),
               hi = map_dbl(iters, ~ unname(quantile(.x, 0.84)))) |>
        select(-iters) |>
        ggplot()
    if (.coef == "autoreg") {
        .ybreaks <- c(0, 0.5, 1)
        .ylims <- c(-0.1, 1.1)
    } else {
        .plot <- .plot +
            # main density curves:
            geom_polygon(data = slope_density_df |>
                             mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) |>
                             filter(coef == .coef, !ui),
                         aes(x = y, y = x), fill = "gray80") +
            geom_polygon(data = slope_sd_density_df |>
                             mutate(y = 7 - y / max(y) * (3.5 - 0.5)) |>
                             filter(coef == .coef, !ui),
                         aes(x = y, y = x), fill = "gray80") +
            # uncertainty interval density curves:
            geom_polygon(data = slope_density_df |>
                             mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) |>
                             filter(coef == .coef, ui) |>
                             arrange(x),
                         aes(x = y, y = x), fill = "gray60") +
            geom_polygon(data = slope_sd_density_df |>
                             mutate(y = 7 - y / max(y) * (3.5 - 0.5)) |>
                             filter(coef == .coef, ui),
                         aes(x = y, y = x), fill = "gray60")
        .ybreaks <- c(-0.4, 0, 0.4)
        .ylims <- c(-0.5, 0.5)
    }

    .plot <- .plot +
        geom_hline(yintercept = 0, color = "gray50")+
        geom_point(aes(taxon, med), size = 1.5, color = "black")+
        geom_linerange(aes(taxon, ymin = lo, ymax = hi), color = "black")+
        scale_y_continuous(.ylab,  breaks = .ybreaks)+
        scale_x_continuous(breaks = 1:6, labels = taxa_labs, position = "top") +
        # scale_x_continuous(breaks = 1:6, labels = taxa_labs,
        #                    sec.axis =
        #                        sec_axis(~ .,
        #                                 breaks = 1:6, labels = taxa_labs)) +
        coord_cartesian(ylim = .ylims, xlim = c(0.5, 7),
                        expand = FALSE) +
        theme(axis.text.y = element_text(size = 8,
                                         margin = margin(0,0,0,r=2)),
              axis.title.y = element_text(size = 10,
                                          margin = margin(0,0,0,r=6)),
              axis.title.x = element_blank(),
              axis.text.x.top = element_text(size = 8, angle = 45,
                                             vjust = 0.1, hjust = 0.1,
                                             margin = margin(0,0,0,b=2)),
              axis.text.x.bottom = element_text(size = 8, angle = -45,
                                                vjust = 0, hjust = 0.1,
                                                margin = margin(0,0,0,t=6)),
              plot.margin = margin(t=0, r=9, b=8, l=16))
    if (.coef == "time") {
        .plot <- .plot +
            geom_text(data = tibble(x = c(1.2,    5.5),
                                    y = c(-0.30, 0.35),
                                    lab = c("among taxa mean", "among taxa SD")),
                      aes(x, y, label = lab), size = 8 / 2.83465,
                      hjust = c(0, 1), vjust = 0.5, color = "gray50")
    }

    return(.plot)
}


slope_p <- lapply(c("autoreg", levels(slope_density_df$coef)), slope_p_fun)



coef_p <- do.call(patchwork::wrap_plots, slope_p) +
    plot_layout(ncol = 2, axes = "collect_x") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag.position = c(-0.02, 1))
# This prevents tags for top panels from being way higher than the others
# due to the large x-axis that's on top:
for (i in 1:2) coef_p[[i]] <- coef_p[[i]] & theme(plot.tag.position = c(-0.02, 0.7))

# coef_p

save_plot("coefficients", coef_p, w = 5, h = 6)
