
#'
#' This file creates figures 1-3 and table 1.
#'




# ==================================================*
# Preliminaries ----
# ==================================================*

# load packages
library(TransTrendsPkg)
library(tidyverse)
library(cowplot)
library(egg)
library(grid)
library(viridisLite)
library(rgdal)
library(ggsn)
library(broom)



# ------------*
# Map abbreviated taxa names to full names, and order them properly for plots
# ------------*

taxa_order <- c(5:6, 4, 1, 3, 2)
taxa_lvls = c("gnap","lyco","sheet","opil","cara", "stap")[taxa_order]
taxa_labs = c("Ground spiders","Wolf spiders","Sheet weavers",
              "Harvestmen","Ground beetles", "Rove beetles")[taxa_order]


# ------------*
# Load and clean data from model fits
# ------------*

# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs))


# Mapping model names with more descriptive ones for the table itself:
names_map <- list(mltd_mltd  = "Full model",
                  mltd_ltd  = "No midges RE",
                  mltd_mtd  = "No lagged midges RE",
                  mltd_mld  = "No time RE",
                  mltd_mlt  = "No distance RE")

# Full model is used most often:
full_model <- paste0("analysis/output/model_fits/",
                     names(which(names_map == "Full model")), ".rds") %>%
    readRDS()

# (`capture.output` used to suppress printing associated with `summary`)
capture.output({full_model_summ <- full_model %>% summary()}, file = "/dev/null")



# Order levels of coefficient factor:
make_coef_fct <- function(.coef) {
    factor(.coef %>% paste(),
           levels = c("time", "dist", "midges", "midges_lag"))
}


bayesian_se <- TransTrendsPkg:::bayesian_se




# ------------*
# Define theme, a colorblind-friendly palette, and
# helper functions/objects for plots
# ------------*
pca_palette <- viridis(3)


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
                    axis.title.x = element_text(size = 12,
                                                margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 10),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6),
                                                vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10),
                                                angle = 270)))

# empty plot
EMPTY <- ggplot() + geom_blank() + theme_void()


# use cairo_pdf because it embeds fonts by default
save_file <- function(x, fn, ...) {
    cairo_pdf(paste0("analysis/output/", gsub(".pdf$", "", fn), ".pdf"), ...)
    print(x)
    dev.off()
}
# remove x axis from a plot
no_x <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
}
# remove y axis from a plot
no_y <- function(p) {
    p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
}
# remove x and y axes from a plot
no_xy <- function(p) {
    p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.title.y = element_blank(), axis.text.y = element_blank())
}
# remove legend from a plot
no_leg <- function(p) {
    p + theme(legend.position = "none")
}





# ==================================================*
# ==================================================*

# Fig 2 - Observed data ----

# ==================================================*
# ==================================================*

# by time
time_p <- data_fit %>%
    ggplot(aes(year, y))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_line(aes(group = plot), size = 0.2, color = "gray70")+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, year) %>%
                  summarize(y = mean(midges_z), .groups = "drop") %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, year) %>%
                                summarize(y = mean(y), .groups = "drop") %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick","dodgerblue"))+
    scale_x_continuous("Year", breaks = c(2008,2012,2016),
                       limits = c(2007,2018))+
    scale_y_continuous("Transformed abundance",
                       # limits = c(-2, 3.1),
                       breaks = c(-2,0,2))+
    NULL

# by distance
dist_p <- data_fit %>%
    ggplot(aes(distance, y))+
    facet_wrap(~taxon_plot, nrow = 4)+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_jitter(aes(group = plot), size = 1, alpha = 0.5,
                width = 0.1, shape = 1)+
    geom_line(data = data_fit %>%
                  group_by(taxon_plot, distance) %>%
                  summarize(y = mean(midges_z), .groups = "drop") %>%
                  mutate(type  = "Midge abundance") %>%
                  bind_rows(data_fit %>%
                                group_by(taxon_plot, distance) %>%
                                summarize(y = mean(y), .groups = "drop") %>%
                                mutate(type  = "Mean by year/distance")) %>%
                  mutate(type = factor(type,
                                       levels = c("Midge abundance",
                                                  "Mean by year/distance"))),
              aes(color = type), size = 0.75)+
    scale_color_manual(NULL, values = c("firebrick", "dodgerblue"))+
    scale_x_continuous("Distance (m)", trans = "log", breaks = c(5,50,500))+
    scale_y_continuous(NULL,
                       # limits = c(-2, 3.1),
                       breaks = c(-2,0,2)) +
    NULL

# separate legend
time_dist_legend <- get_legend(
    time_p +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom")
)


prow <- plot_grid(time_p %>% no_leg(),
                  dist_p  %>% no_leg(),
                  labels = c("a", "b"),
                  align = "vh", label_fontface = "plain", label_size = 16)

obs_data_p <- plot_grid(time_dist_legend, prow, ncol = 1, rel_heights = c(0.1, 1))

# obs_data_p

# save_file(obs_data_p, "fig2", width = 6, height = 4)




# ==================================================*
# ==================================================*

# Fig 3 - Coefficients ----

# ==================================================*
# ==================================================*


#' Extract posterior densities for the fixed effects
#' and for the random effect SDs
#' Column `ui` indicates whether the densities are for the 68%
#' uncertainty intervals.
slope_density_df <- full_model$stan %>%
    rstan::extract(pars = "alpha") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int", "midges", "midges_lag", "time", "dist"); .} %>%
    as_tibble() %>%
    select(-int) %>%  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = make_coef_fct(coef)) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ui <- unname(quantile(z[["value"]], c(0.16, 0.84)))
        X <- density(z[["value"]], n = 2048)
        Z <- density(z[["value"]], from = .ui[1], to = .ui[2])
        tibble(coef = z[["coef"]][1],
               x = c(X$x, c(Z$x[1], Z$x, tail(Z$x, 1))),
               y = c(X$y, c(0, Z$y, 0)),
               ui = c(rep(FALSE, length(X$x)), rep(TRUE, length(Z$x)+2)))
    })

slope_sd_density_df <- full_model$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int_plot",
                      "int_tax",
                      "midges",
                      "midges_lag",
                      "time",
                      "dist"); .} %>%
    as_tibble() %>%
    select(all_of(c("midges", "midges_lag", "time", "dist"))) %>%
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = make_coef_fct(coef)) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ui <- unname(quantile(z[["value"]], c(0.16, 0.84)))
        X <- density(z[["value"]], n = 2048)
        Z <- density(z[["value"]], from = .ui[1], to = .ui[2])
        tibble(coef = z[["coef"]][1],
               x = c(X$x, c(Z$x[1], Z$x, tail(Z$x, 1))),
               y = c(X$y, c(0, Z$y, 0)),
               ui = c(rep(FALSE, length(X$x)), rep(TRUE, length(Z$x)+2)))
    })





#' Get the posterior draws for random effect estimates
#'
get_re_draws <- function(object) {
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

    draws <- object$rnd_lvl_names %>%
        as_tibble() %>%
        mutate(iters = map(1:n(), ~ E[,.x]))

    return(draws)
}


re_draws <- get_re_draws(full_model) %>%
    filter(Name != "(Intercept)", Groups == "taxon") %>%
    select(-Groups) %>%
    rename(taxon = Level, coef = Name)








# Creates subpanels for each predictor (fig 2 b--d)
slope_p_fun <- function(.coef) {

    # .coef = "midges_lag"

    .coef <- match.arg(tolower(.coef), levels(slope_density_df$coef))

    # Version stored in stan:
    .coef_stan <- ifelse(.coef == "midges_lag", .coef, paste0(.coef, "_z"))

    .ylab <- str_replace(.coef, "s$", "") %>%
        str_replace("dist", "distance") %>%
        str_replace("midges_lag", "lagged midge") %>%
        str_to_sentence() %>%
        paste("response")

    .plot <- re_draws %>%
        filter(coef == .coef_stan) %>%
        mutate(taxon = factor(taxon, levels = taxa_lvls) %>%
                   as.integer()) %>%
        mutate(med = map_dbl(iters, ~ median(.x)),
               lo = map_dbl(iters, ~ unname(quantile(.x, 0.16))),
               hi = map_dbl(iters, ~ unname(quantile(.x, 0.84)))) %>%
        select(-iters) %>%
        ggplot() +
        # main density curves:
        geom_polygon(data = slope_density_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef, !ui),
                     aes(x = y, y = x), fill = "gray80") +
        geom_polygon(data = slope_sd_density_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef, !ui),
                     aes(x = y, y = x), fill = "gray80") +
        # uncertainty interval density curves:
        geom_polygon(data = slope_density_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef, ui) %>%
                         arrange(x),
                     aes(x = y, y = x), fill = "gray60") +
        geom_polygon(data = slope_sd_density_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef, ui),
                     aes(x = y, y = x), fill = "gray60") +
        geom_hline(yintercept = 0, color = "gray50")+
        geom_point(aes(taxon, med), size = 1.5, color = "black")+
        geom_linerange(aes(taxon, ymin = lo, ymax = hi), color = "black")+
        scale_y_continuous(.ylab,  breaks = c(-0.4, 0, 0.4))+
        scale_x_continuous(NULL,
                           breaks = 1:6, labels = taxa_labs,
                           sec.axis =
                               sec_axis(~ .,
                                        breaks = 1:6, labels = taxa_labs)) +
        coord_cartesian(ylim = c(-0.5, 0.5), xlim = c(0.5, 7),
                        expand = FALSE) +
        theme(axis.text.y.left = element_text(size = 8,
                                              margin = margin(0,0,0,r=2)),
              axis.text.y.right = element_blank(),
              axis.ticks.y.right = element_blank(),
              axis.title.y.left = element_text(size = 10,
                                               margin = margin(0,0,0,r=6)),
              axis.title.y.right = element_text(size = 10,
                                                margin = margin(0,0,0,l=6)),
              axis.title.x = element_blank(),
              axis.text.x.top = element_text(size = 8, angle = 45,
                                             vjust = 0.1, hjust = 0.1,
                                             margin = margin(0,0,0,b=2)),
              axis.text.x.bottom = element_text(size = 8, angle = -45,
                                                vjust = 0, hjust = 0.1,
                                                margin = margin(0,0,0,t=6)),
              plot.margin = margin(t=0, r=9, b=8, l=16)) +
        NULL

    return(.plot)
}


slope_p <- lapply(unique(slope_density_df$coef), slope_p_fun)

#' Add labels for density curves to Fig 2b
slope_p[[1]] <- slope_p[[1]] +
    geom_text(data = tibble(x = c(0.7,    3.5),
                            y = c(-0.28, 0.35),
                            lab = c("among taxa mean", "among taxa SD")),
              aes(x, y, label = lab), size = 8 / 2.83465,
              hjust = 0, vjust = 0.5, color = "gray50")







# AR coefficient plot:
ar_p <- rstan::extract(full_model$stan, "phi")[[1]] %>%
    as.data.frame() %>%
    set_names(gsub("^taxon", "", full_model$ar_names)) %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = "taxon") %>%
    group_by(taxon) %>%
    summarize(med = median(value),
              lo = unname(quantile(value, 0.16)),
              hi = unname(quantile(value, 0.84))) %>%
    mutate(taxon = factor(taxon, levels = taxa_lvls) %>%
               as.integer()) %>%
    ggplot(aes(taxon, med))+
    geom_hline(yintercept = 0, color = "gray50")+
    geom_point(size = 1.5)+
    geom_linerange(aes(ymin = lo, ymax = hi))+
    scale_x_continuous(NULL,
                       breaks = 1:6, labels = taxa_labs,
                       sec.axis =
                           sec_axis(~ .,
                                    breaks = 1:6, labels = taxa_labs)) +
    scale_y_continuous("AR parameter", breaks = c(0, 0.5, 1)) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0.5, 7),
                    expand = FALSE) +
    theme(axis.text.y = element_text(size = 8, margin = margin(0,0,0,r=2)),
          axis.title.y = element_text(size = 10, margin = margin(0,0,0,r=6)),
          axis.title.x = element_blank(),
          axis.text.x.top = element_text(size = 8, angle = 45,
                                         vjust = 0.1, hjust = 0.1,
                                         margin = margin(0,0,0,b=2)),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          plot.margin = margin(t=0, r=8, b=8, l=16)) +
    NULL





#' Helper plot functions.
#' No taxa names, top:
ntt <- function(.x) .x + theme(axis.text.x.top = element_blank())
#' No taxa names, bottom:
ntb <- function(.x) .x + theme(axis.text.x.bottom = element_blank(),
                               axis.ticks.x.bottom = element_blank())

x_top <- get_x_axis(slope_p[[1]] + theme(axis.ticks.x.top = element_blank()),
                    "top") %>%
    ggdraw(clip = "on")

coef_p <- ggarrange(x_top, x_top,
                    ar_p %>% ntt(),
                    slope_p[[1]] %>% ntt() %>% ntb(),
                    slope_p[[2]] %>% ntt() %>% ntb(),
                    slope_p[[3]] %>% ntt() %>% ntb(),
                    slope_p[[4]] %>% ntt() %>% ntb(),
                    labels = c("", "", sprintf("(%s)", letters[1:5])),
                    label.args = list(gp = gpar(font = 1, fontsize = 14),
                                      x = unit(0, "line"), hjust = 0),
                    heights = c(0.5, 1, 1, 1),
                    ncol = 2, draw = FALSE)




# coef_p


save_file(coef_p, "fig3", width = 5, height = 6)







# =============================================================================*
# =============================================================================*

# Table 1 ----

# =============================================================================*
# =============================================================================*


#' The first would report main and random effects estimates (basically, the
#' stuff that is in the lower part of the `summary` method for the package)
#' for the full model.

make_summ_coef_fct <- function(.coef) {
    factor(.coef %>% paste(),
           levels = c("(Intercept)", "time_z", "dist_z", "midges_z",
                      "midges_lag"),
           labels = c("Intercept", "Time", "Distance", "Midges",
                      "Lagged midges"))
}


list(full_model_summ %>%
         .[["fixef"]] %>%
         rownames_to_column("Parameter") %>%
         mutate(Parameter = make_summ_coef_fct(Parameter)) %>%
         mutate(`Fixed effect` = sprintf("$%.3f \\pm %.3f$", Median,
                                         `Std.Error`)) %>%
         select(-Median, -`Std.Error`),
     full_model_summ %>%
         .[["ranef"]] %>%
         filter(Groups == "taxon") %>%
         mutate(Parameter = make_summ_coef_fct(Name)) %>%
         select(-Groups, -Name) %>%
         mutate(`Random effect SD` = sprintf("$%.3f \\pm %.3f$", Median,
                                             `Std.Error`)) %>%
         select(-Median, -`Std.Error`)) %>%
    c(by = "Parameter") %>%
    do.call(what = left_join) %>%
    arrange(Parameter) %>%
    knitr::kable(format = "latex", escape = FALSE, booktabs = TRUE)



# =============================================================================*
# =============================================================================*

# Table 2 ----

# =============================================================================*
# =============================================================================*


tab2 <- as_tibble(names_map) %>%
    pivot_longer(everything(), names_to = "abbrev", values_to = "latex") %>%
    mutate(LL = NA_real_, LOOIC = NA_real_, WAIC = NA_real_)

for (f in names(names_map)) {
    fit <- paste0("analysis/output/model_fits/", f, ".rds") %>%
        readRDS()
    i <- which(tab2$abbrev == f)
    l <- loo(fit)
    w <- waic(fit)
    tab2[i,"LL"] <- median(rstan::extract(fit$stan, "log_lik_sum")[[1]])
    tab2[i,"LOOIC"] <- l$estimates["looic","Estimate"]
    tab2[i,"WAIC"] <- w[["estimates"]]["waic","Estimate"]
}
rm(f, fit, i, l, w)
gc()

tab2 %>%
    .[c(1, 4, 5, 2, 3),] %>%
    mutate(`$\\Delta$ LOOIC` = LOOIC - LOOIC[latex == "Full model"],
           `$\\Delta$ WAIC` = WAIC - WAIC[latex == "Full model"]) %>%
    select(-abbrev) %>%
    rename(Model = latex, `Log likelihood` = LL) %>%
    select(Model, `Log likelihood`, LOOIC, `$\\Delta$ LOOIC`, everything()) %>%
    knitr::kable(format = "latex", escape = FALSE, digits = 1, booktabs = TRUE)
