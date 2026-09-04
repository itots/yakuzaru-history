library(sf)
library(tidyverse)
library(rnaturalearth)
library(forcats)
library(ggrepel)
library(ggsflabel)
library(ggplot2)

# Input path
input_path <- list(
  info = "../list/population.csv",
  morph_pop = "../list/population_morph.csv"
)

# Color palettes
color <- list(three = c("#1784bf", "#d1990d", "#07360c"))

# Load data
tbs <- list(
  info = read_csv(input_path$info) %>%
    mutate(
      population = fct_inorder(population)
    ),
  morph_pop = read_csv(input_path$morph_pop) %>%
    mutate(
      population = fct_inorder(population)
    )
)

# Load World Map
map <- list(
  world = ne_countries(scale = "large")
)

# Create SF objects & Transform CRS
sf <- list(
  points = st_as_sf(tbs$info, coords = c("longitude", "latitude"), crs = 4326) %>%
    mutate(
      population = fct_inorder(population)
    ) %>%
    st_transform(crs = 6677),
  
  points_morph = st_as_sf(tbs$morph_pop %>%
                            filter(!is.na(longitude) & !is.na(latitude)), 
                          coords = c("longitude", "latitude"), crs = 4326) %>%
    mutate(
      population = fct_inorder(population)
    ) %>%
    st_transform(crs = 6677)
)

map$world_proj <- st_transform(map$world, crs = 6677)


lim_lon <- c(127.5, 145)
lim_lat <- c(28, 45)

bbox_geo <- st_sfc(st_point(c(lim_lon[1], lim_lat[1])), 
                   st_point(c(lim_lon[2], lim_lat[2])), 
                   crs = 4326)
bbox_proj <- st_transform(bbox_geo, crs = 6677)
coords_proj <- st_coordinates(bbox_proj)

xlim_proj <- c(coords_proj[1, "X"], coords_proj[2, "X"])
ylim_proj <- c(coords_proj[1, "Y"], coords_proj[2, "Y"])


# Plotting
plots <- list(
  map = ggplot() +
    geom_sf(data = map$world_proj, linewidth = 0.2) +
    geom_sf(data = sf$points, size = 2.5, stroke = 0.35,
            aes(shape = population, color = cluster)) +
    geom_sf(data = sf$points_morph, size = 1,
            aes(color = cluster)) +
    geom_sf_label_repel(data = sf$points, aes(label = population), 
                        box.padding = 1.2,
                        point.padding = 0.4,
                        segment.size = 0.1,
                        linewidth = 0.1,
                        min.segment.length = 0.5, max.overlaps = Inf) +
    theme_minimal() +
    coord_sf(
      xlim = xlim_proj,
      ylim = ylim_proj,
      default = TRUE) + 
    theme(legend.position = c(0.01, 0.99),
          legend.justification = c("left", "top"),
          legend.background = element_rect(fill = alpha("white", 0.6), color = NA), # ←ここを追加（枠線はNAで消すのがおすすめ）
          axis.title = element_blank()) +
    scale_shape_manual(values = 0:11, guide = "none") +
    scale_color_manual(values = color$three) +
    guides(fill = guide_legend(nrow = 2))
)

# Save output
ggsave("output/map.pdf", plots$map, width = 120, height = 120, units = "mm", device = cairo_pdf)