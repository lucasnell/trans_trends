
#'
#' This file creates maps showing where samples were taken at Lake Myvatn.
#'


# =============================================================================*
# Preliminaries ----
# =============================================================================*

source("scripts/00-preamble.R")


library(sf)
library(ggspatial) # annotation_scale and annotation_north_arrow


# =============================================================================*
# Read data ----
# =============================================================================*

myvatn_map <- st_read("data/Myvatn_WSGUTM28.geojson")

# Iceland outline is from GADM data (version 4.1; https://gadm.org/)
iceland_map <- st_read("data/gadm41_ISL_0.geojson") |>
    st_transform(st_crs("+proj=utm +zone=28 ellps=WGS84"))
geom_nrows <- map_int(iceland_map$geometry[[1]], \(x) nrow(x[[1]]))
st_geometry(iceland_map) <- iceland_map$geometry[[1]][[which(geom_nrows == max(geom_nrows))]] |>
    st_polygon() |>
    st_sfc(crs = st_crs(iceland_map))

pit_df <- read_csv("data/pitfall_locations.csv", col_types = cols()) |>
    pivot_longer(`5m`:`500m`, names_to = "dist") |>
    pivot_wider(names_from = coord) |>
    mutate(site = tolower(site),
           dist = str_remove(dist, "m") |>
               as.integer() |>
               factor(levels = c(5, 50, 150, 500))) |>
    filter(site %in% unique(read_csv("data/myv_arth.csv",
                                     col_types = cols())[["trans"]])) |>
    mutate(site = factor(site, levels = c("vin", "kal", "hag", "non", "btl"),
                         labels = c("Vindbelgur", "Kálfaströnd", "Haganes",
                                    "Nóntangi", "Fellshóll")))


# =============================================================================*
# Create maps ----
# =============================================================================*

myvatn_p <- myvatn_map |>
    ggplot(aes(x, y)) +
    geom_sf(color = "gray20", linewidth = 0.25, fill = "lightskyblue1", inherit.aes = FALSE) +
    geom_point(data = filter(pit_df, dist == 5), color = "firebrick1", size = 4) +
    geom_point(data = filter(pit_df, dist == 5), size = 4, shape = 1) +
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
    annotation_north_arrow(which_north = "true", location = "bl",
                           height = unit(0.1, "npc"), width = unit(0.1, "npc"),
                           pad_x = unit(0.02, "npc"), pad_y = unit(0.02, "npc"),
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "bl",
                     text_pad = unit(0.02, "npc"), height = unit(0.03, "npc"),
                     pad_x = unit(0.15, "npc"), pad_y = unit(0.045, "npc")) +
    scale_fill_manual(values = c("lightblue", "white"), guide = "none") +
    scale_color_viridis_d(begin = 0.1, end = 0.9, guide = "none") +
    coord_sf(xlim = c(NA, 411932 + 700)) +
    theme_void() +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent"))



iceland_p <- iceland_map |>
    ggplot(aes(x, y)) +
    geom_sf(color = "gray20", linewidth = 0.25, fill = "gray80", inherit.aes = FALSE) +
    geom_point(data = tibble(x = 403118.1, y = 7271491),
               size = 2, color = "black", shape = 16) +
    theme_void() +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent"))


sites_p <- myvatn_p +
    inset_element(iceland_p, left = -0.09, right = -0.09 + 0.55,
                  bottom = 0.7 + 0.03, top = 1 + 0.03,
                  align_to = "full", clip = FALSE)

# sites_p

# save_plot("sample-map", sites_p, w = 3.21, h = 3.6)

