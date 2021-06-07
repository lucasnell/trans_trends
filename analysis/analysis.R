
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
# Load and clean data from model fit
# ------------*

# load data
data_fit <- read_csv("analysis/data_fit.csv") %>%
    mutate(taxon_plot = factor(taxon, levels = taxa_lvls, labels = taxa_labs))
fit <- readRDS("analysis/output/fit.rds")
fit_sum <- read_csv("analysis/output/fit_sum.csv")
coef_sum <- readRDS("analysis/output/coef_sum.rds")

coef_sum$int_taxon <- coef_sum$int_taxon %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)),
           coef = "int")

# Order levels of coefficient factor:
make_coef_fct <- function(.coef) {
    factor(.coef %>% paste(),
           levels = c("time", "dist", "midges"),
           labels = c("time", "distance", "midges"))
}

coef_sum$beta <- coef_sum$beta %>%
    mutate(coef = make_coef_fct(coef),
           tx = as.numeric(factor(taxon, levels = taxa_lvls)))
coef_sum$alpha <- coef_sum$alpha %>%
    filter(coef != "int") %>%
    mutate(coef = make_coef_fct(coef))
coef_sum$sig_beta <- coef_sum$sig_beta %>%
    filter(!(coef %in% c("int_tax","int_tax_plot","int_tax_trans"))) %>%
    mutate(coef = make_coef_fct(coef))
coef_sum$ar <- coef_sum$ar %>%
    mutate(tx = as.numeric(factor(taxon, levels = taxa_lvls)))







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




# =============================================================================*
# =============================================================================*

# Fig 1 - Map ----

# =============================================================================*
# =============================================================================*



# From decimal degrees to UTM, assuming it's Iceland and using WGS84
to_utm <- function(.df, .lat = "lat", .lon = "long") {
    .cord.dec <- SpatialPoints(cbind(.df[[.lon]], .df[[.lat]]),
                               proj4string=CRS("+proj=longlat"))
    .cord.UTM <- spTransform(.cord.dec, CRS("+proj=utm +zone=28 ellps=WGS84"))
    .df[[.lon]] <- .cord.UTM@coords[,1]
    .df[[.lat]] <- .cord.UTM@coords[,2]
    return(.df)
}



pit_df <- "site,coord,5m,50m,150m,500m
BTL,x,407123,407104,407094,NA
BTL,y,7273389,7273347,7273251,NA
HAG,x,405379,405309,405209,404938
HAG,y,7274913,7274902,7274870,7274731
SKF,x,403834,403818,403655,403497
SKF,y,7277920,7277952,7278041,7278150
VIN,x,406346,406358,406375,NA
VIN,y,7278620,7278676,7278767,NA
FLG,x,410828,410859,410930,411208
FLG,y,7281762,7281789,7281863,7282065
NON,x,410707,410695,410690,NA
NON,y,7276185,7276145,7276040,NA
KAL,x,409441,409448,409461,NA
KAL,y,7274132,7274078,7273987,NA" %>%
    read_csv() %>%
    pivot_longer(`5m`:`500m`, names_to = "dist") %>%
    pivot_wider(names_from = coord) %>%
    mutate(site = tolower(site),
           dist = str_remove(dist, "m") %>%
               as.integer() %>%
               factor(levels = c(5, 50, 150, 500))) %>%
    filter(site %in% unique(read_csv("analysis/data_fit.csv")[["trans"]])) %>%
    mutate(site = factor(site, levels = c("vin", "kal", "hag", "non", "btl"),
                         labels = c("Vindbelgur", "Kálfaströnd", "Haganes",
                                    "Nóntangi", "Fellshóll")))


myvatn_df <- readOGR(dsn = "~/Box Sync/2020/midges/shapefiles/myvatn",
                     layer = "Myvatn_WSGUTM28") %>%
    tidy() %>%
    rename(x = long, y = lat)

# Iceland outline is from GADM data (version 3.6; https://gadm.org/)
iceland_df <- readOGR(dsn = "~/Box Sync/2020/midges/shapefiles/iceland",
                      layer = "gadm36_ISL_0") %>%
    tidy() %>%
    to_utm() %>%
    # -----------`
    # Shown below are two ways to filter this dataset, to avoid
    # plotting islands far from shore:
    # -----------`
    # 1. You can filter out islands that are very far from shore:
    # filter(!piece %in% c(90, 133, 143, 157, 215, 244, 257, 258, 260, 262))
    # 2. Filter for just the mainland:
    filter(piece == 1) %>%
    rename(x = long, y = lat)




myvatn_pit_p <- myvatn_df %>%
    ggplot(aes(x, y)) +
    geom_polygon(aes(group = group, fill = hole), color = "black", size = 0.1) +
    geom_point(data = filter(pit_df, dist == 5), aes(color = site), size = 4) +
    geom_point(data = filter(pit_df, dist == 5), size = 4, shape = 1) +
    geom_segment(data = tibble(xe = 405.5e3, x = xe - 1e3, y = 7271.5e3),
                 aes(xend = xe, yend = y), size = 1) +
    geom_text(data = tibble(x = 405.5e3-500, y = 7271.5e3 + 200),
              label = "1 km", vjust = 0, hjust = 0.5, size = 10 / 2.83465) +
    north(rename(myvatn_df, long = x, lat = y), "bottomleft", symbol = 12) +
    geom_text(data = filter(pit_df, dist == 5) %>%
                  mutate(x = case_when(site == "Kálfaströnd" ~ x + 2000,
                                       site == "Haganes" ~ x - 1400,
                                       site == "Vindbelgur" ~ x - 1200,
                                       site == "Nóntangi" ~ x + 600,
                                       site == "Fellshóll" ~ x - 400,
                                       TRUE ~ x),
                         y = case_when(site == "Nóntangi" |
                                           site == "Fellshóll" ~ y - 600,
                                       site == "Vindbelgur" ~ y + 600,
                                       TRUE ~ y)),
              aes(label = site), size = 10 / 2.83465, fontface = "plain") +
    scale_fill_manual(values = c("lightblue", "white"), guide = FALSE) +
    scale_color_viridis_d(begin = 0.1, end = 0.9, guide = FALSE) +
    coord_equal(xlim = c(NA, 411932 + 700)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL





myvatn_pit_inset_p <- iceland_df %>%
    ggplot(aes(x, y)) +
    geom_polygon(aes(group = group), color = "black", size = 0.1,
                 fill = "gray80") +
    geom_point(data = tibble(x = 403118.1, y = 7271491),
               size = 2, color = "black", shape = 16) +
    coord_equal() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    NULL





sites_p <- ggdraw() +
    draw_plot(myvatn_pit_p) +
    draw_plot(myvatn_pit_inset_p, x = 0, y = 1, hjust = 0, vjust = 1,
              width = 0.55, height = 0.3)




save_file(sites_p, "sites", width = 3.21, height = 3.6)





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

fig1 <- plot_grid(time_dist_legend, prow, ncol = 1, rel_heights = c(0.1, 1))

# fig1

# save_file(fig1, "fig1", width = 6, height = 4)




# ==================================================*
# ==================================================*

# Fig 3 - Coefficients ----

# ==================================================*
# ==================================================*


#' Extract posterior densities for the fixed effects
#' and for the random effect SDs
slope_density_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int", "midges", "time", "distance"); .} %>%
    as_tibble() %>%
    select(-int) %>%  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 2048)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)
slope_sd_density_df <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int_tax",
                      "int_plot",
                      # "int_trans",
                      "midges",
                      "time",
                      "distance"); .} %>%
    as_tibble() %>%
    select(all_of(c("midges", "time", "distance"))) %>%
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        X <- density(z[["value"]], n = 2048)
        tibble(coef = z[["coef"]][1],
               x = X$x,
               y = X$y)
    }) %>%
    arrange(coef, x)


#' Do the same as above, but restrict the densities to being between
#' the 68% uncertainty intervals.
#' These are used to create the darker shaded regions in fig 2.
slope_density_UI_df <- fit$stan %>%
    rstan::extract(pars = "alpha") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int", "midges", "time", "distance"); .} %>%
    as_tibble() %>%
    select(-int) %>%  # intercept not necessary
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ci <- quantile(z[["value"]], c(0.16, 0.84))
        X <- density(z[["value"]], from = .ci[["16%"]], to = .ci[["84%"]])
        tibble(coef = z[["coef"]][1],
               x = c(X$x[1], X$x, tail(X$x, 1)),
               y = c(0, X$y, 0))
    }) %>%
    arrange(coef, x)
slope_sd_density_UI_df <- fit$stan %>%
    rstan::extract(pars = "sig_beta") %>%
    do.call(what = cbind) %>%
    {colnames(.) <- c("int_tax",
                      "int_plot",
                      # "int_trans",
                      "midges",
                      "time",
                      "distance"); .} %>%
    as_tibble() %>%
    select(all_of(c("midges", "time", "distance"))) %>%
    pivot_longer(everything(), names_to = "coef") %>%
    mutate(coef = factor(coef, levels = c("time", "distance", "midges"))) %>%
    split(.$coef) %>%
    map_dfr(function(z) {
        .ci <- quantile(z[["value"]], c(0.16, 0.84))
        X <- density(z[["value"]], from = .ci[["16%"]], to = .ci[["84%"]])
        tibble(coef = z[["coef"]][1],
               x = c(X$x[1], X$x, tail(X$x, 1)),
               y = c(0, X$y, 0))
    }) %>%
    arrange(coef, x)


# Creates subpanels for each predictor (fig 2 b--d)
slope_p_fun <- function(.coef) {

    .coef <- match.arg(tolower(.coef), levels(coef_sum$beta$coef))

    .ylab <- gsub("s$", "", .coef) %>%
        str_to_title() %>%
        paste("response")

    .plot <- coef_sum$beta %>%
        filter(coef == .coef) %>%
        ggplot()+
        # main density curves:
        geom_polygon(data = slope_density_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray80") +
        geom_polygon(data = slope_sd_density_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray80") +
        # uncertainty interval density curves:
        geom_polygon(data = slope_density_UI_df %>%
                         mutate(y = y / max(y) * (3.5 - 0.5) + 0.5) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray60") +
        geom_polygon(data = slope_sd_density_UI_df %>%
                         mutate(y = 7 - y / max(y) * (3.5 - 0.5)) %>%
                         filter(coef == .coef),
                     aes(x = y, y = x), fill = "gray60") +
        geom_hline(yintercept = 0, color = "gray50")+
        geom_point(aes(tx, mi), size = 1.5, color = "black")+
        geom_linerange(aes(tx, ymin = lo, ymax = hi), color = "black")+
        scale_y_continuous(.ylab,  breaks = c(-0.5, 0, 0.5))+
        scale_x_continuous(NULL,
                           breaks = 1:6, labels = taxa_labs,
                           sec.axis =
                               sec_axis(~ .,
                                        breaks = 1:6, labels = taxa_labs)) +
        coord_cartesian(ylim = c(-0.68, 0.68), xlim = c(0.5, 7),
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
              plot.margin = margin(t=8, r=9, b=8, l=8)) +
        NULL

    return(.plot)
}


slope_p <- lapply(levels(coef_sum$beta$coef), slope_p_fun)

#' Add labels for density curves to Fig 2b
slope_p[[1]] <- slope_p[[1]] +
    geom_text(data = tibble(x = c(0.7,    4.5),
                            y = c(-0.28, 0.45),
                            lab = c("among taxa mean", "among taxa SD")),
              aes(x, y, label = lab), size = 8 / 2.83465,
              hjust = 0, vjust = 0.5, color = "gray50")


# AR coefficient plot:
ar_p <- coef_sum$ar %>%
    ggplot(aes(tx, mi))+
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
          plot.margin = margin(t=8, r=8, b=8, l=8)) +
    NULL





#' Helper plot functions.
#' No taxa names, top:
ntt <- function(.x) .x + theme(axis.text.x.top = element_blank())
#' No taxa names, bottom:
ntb <- function(.x) .x <- .x + theme(axis.text.x.bottom = element_blank())


fig2 <- ggarrange(ar_p,
                  slope_p[[1]] %>% ntt() %>% ntb(),
                  slope_p[[2]] %>% ntt() %>% ntb(),
                  slope_p[[3]] %>% ntt(),
                  labels = letters[1:4],
                  label.args = list(gp = gpar(font = 1, fontsize = 16),
                                    x = unit(0, "line"), hjust = 0),
                  ncol = 1, draw = FALSE)



# fig2


# save_file(fig2, "fig2", width = 3, height = 8)












