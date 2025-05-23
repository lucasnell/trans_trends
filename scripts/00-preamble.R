suppressPackageStartupMessages({
    library(TransTrendsPkg)
    library(tidyverse)
    library(patchwork)
    library(viridisLite)
})


# set ggplot2 theme
theme_set(theme_bw() %+replace%
              theme(panel.grid = element_blank(),
                    strip.background = element_blank(),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 12),
                    legend.background = element_blank(),
                    axis.text = element_text(size = 8, color = "black"),
                    axis.title.y = element_text(size = 10, angle = 90,
                                                margin = margin(0,0,0,r=6)),
                    axis.title.x = element_text(size = 10,
                                                margin = margin(0,0,0,t=6)),
                    strip.text = element_text(size = 9),
                    strip.text.x = element_text(margin = margin(b = 2, t = 6),
                                                vjust = 0),
                    strip.text.y = element_text(margin = margin(0,0,0,10),
                                                angle = 270)))


# Set up parallel processing
options(mc.cores = max(1, parallel::detectCores()-2))


# Standardize taxa order
taxa_order <- c(5:6, 4, 1, 3, 2)
taxa_lvls <- c("gnap","lyco","sheet","opil","cara", "stap")[taxa_order]
taxa_labs <- c("Ground spiders","Wolf spiders","Sheet weavers",
               "Harvestmen","Ground beetles", "Rove beetles")[taxa_order]


# Names of RDS files created in `scripts/02-clean-data.R` that contain cleaned data:
data_rds <- list(lag = "rds_files/data_fit_lag.rds",
                 nolag = "rds_files/data_fit_nolag.rds")

# Names of RDS files created in `scripts/analysis.R` that contain model fits:
model_rds <- list(lag = "rds_files/lag-model.rds",
                  nolag = "rds_files/no_lag-model.rds",
                  nolag_reds = "rds_files/no_lag-reduced-models.rds",
                  nolag_red_comps = "rds_files/no_lag-reduced-models-compares.rds")

# Save figure file into `./figures/` using cairo_pdf (that embeds fonts by default)
# Can save as PNG file in addition to or in lieu of a PDF
save_plot <- function(n, p, w, h, .pdf = TRUE, .png = FALSE,
                      pdf_args = list(), png_args = list()) {
    stopifnot(is.character(n) && length(n) == 1)
    stopifnot(is_ggplot(p) || is.function(p))
    stopifnot(is.numeric(w) && length(w) == 1)
    stopifnot(is.numeric(h) && length(h) == 1)
    stopifnot(is.logical(.pdf) && length(.pdf) == 1)
    stopifnot(is.logical(.png) && length(.png) == 1)
    stopifnot(inherits(pdf_args, "list"))
    stopifnot(inherits(png_args, "list"))
    old_warn <- getOption("warn")
    options(warn = -1)
    if (.pdf) {
        args <- list(filename = sprintf("figures/%s.pdf", n),
                     width = w, height = h, bg = NA)
        if (length(pdf_args) > 0) {
            stopifnot(!is.null(names(pdf_args)) && all(names(pdf_args) != ""))
            for (n in names(pdf_args)) args[[n]] <- pdf_args[[n]]
        }
        do.call(cairo_pdf, args)
        if (is.function(p)) {
            p()
        } else {
            plot(p)
        }
        tmp <- dev.off()
    }
    if (.png) {
        args <- list(filename = sprintf("figures/%s.png", n),
                     width = w, height = h, units = "in",
                     bg = NA, res = 300)
        if (length(png_args) > 0) {
            stopifnot(!is.null(names(png_args)) && all(names(png_args) != ""))
            for (n in names(png_args)) args[[n]] <- png_args[[n]]
        }
        do.call(png, args)
        if (is.function(p)) {
            p()
        } else {
            plot(p)
        }
        tmp <- dev.off()
    }
    options(warn = old_warn)
    invisible(NULL)
}

