library(tidyverse)
library(sf)
library(USAboundaries)

source("01-map-functions.R")

wa_state <- us_states(state = "WA", resolution = "high")
wa_crs <- st_crs(32149)

cow_huc <- read_wbd(8) |>
  select(name, huc8) |>
  filter(huc8 %in% c("17080004", "17080005")) |>
  st_zm() |>
  st_transform(wa_crs)

cow_wshed <- cow_huc |>
  st_union() |>
  st_sf() |>
  st_transform(wa_crs)

cow_area <- read_nhd("Area") |>
  st_zm() |>
  st_transform(wa_crs) |>
  st_filter(cow_huc, .predicate = st_intersects)

col_area <- read_nhd("Area") |>
  st_zm() |>
  st_transform(wa_crs) |>
  st_filter(cow_area, .predicate = st_is_within_distance, 1e3)

## cow_fl <- read_nhd("Flowline") |>
##   st_zm() |>
##   st_filter(cow_huc, .predicate = st_intersects)

## cow_line <- read_nhd("Line") |>
##   st_zm() |>
##   st_filter(cow_huc, .predicate = st_intersects)

cow_wb <- read_nhd("Waterbody") |>
  st_zm() |>
  st_transform(wa_crs) |>
  st_filter(cow_huc, .predicate = st_intersects) |>
  ## Only include LakePond
  filter(grepl("Mayfield", gnis_name) |
         grepl("Riffe", gnis_name) |
         grepl("Scanewa", gnis_name))

cow_dam <- read_sf("data/Power_Plants/Power_Plants.shp") |>
  filter(State == "Washington") |>
  st_transform(wa_crs) |>
  st_filter(cow_wshed, .predicate = st_within) |>
  filter(Plant_Name != "Packwood")

## cow_pt <- read_nhd("Point") |>
##   st_zm() |>
##   st_filter(cow_huc, .predicate = st_intersects)


inset_map <- ggplot() +
  geom_sf(data = wa_state, fill = "white") +
  geom_sf(data = cow_wshed) +
  coord_sf(crs = wa_crs, datum = NA) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())

cow_labs <- tribble(
~ name, ~ geometry,
"Columbia River", st_point(c(3e5, 8.75e4), dim = "XY"),
"Cowlitz River", st_point(c(3.14e5, 12.7e4), dim = "XY"),
"Toutle River", st_point(c(3.23e5, 11.15e4), dim = "XY"),
"Cispus River", st_point(c(3.965e5, 12.25e4), dim = "XY")) |>
  st_sf(crs = wa_crs)

cow_wb_labs <- tribble(
~ gnis_name, ~ x, ~ y,
"Riffe Lake", 4e3, 2.5e3,
"Lake Scanewa", 2e3, 3e3,
"Mayfield Lake", -5.5e3, 2e3)

cow_dam_labs <- cow_dam |>
  select(Plant_Name) |>
  mutate(name = paste(Plant_Name, "Dam"))

ggplot() +
  geom_sf(data = wa_state, fill = NA) +
  geom_sf(data = cow_wshed, fill = "gray95", color = "gray50") +
  geom_sf(data = col_area, fill = "skyblue", color = "skyblue") +
  ## geom_sf(data = cow_line, color = "skyblue", fill = "skyblue") +
  geom_sf(data = cow_wb, color = "skyblue", fill = "skyblue") +
  geom_sf(data = cow_dam, size = 2) +
  geom_sf_text(data = cow_wb, aes(label = gnis_name),
               nudge_x = cow_wb_labs$x,
               nudge_y = cow_wb_labs$y) +
  geom_sf_label(data = cow_dam_labs, aes(label = name),
                nudge_x = c(0, 0, 0),
                nudge_y = c(-3e3, 3e3, -3e3)) +
  geom_sf_text(data = cow_labs, aes(label = name)) +
  coord_sf(xlim = st_bbox(cow_wshed)[c(1, 3)],
           ylim = st_bbox(cow_wshed)[c(2, 4)],
           crs = wa_crs) +
  theme_minimal() +
  theme(axis.title = element_blank()) +
  annotation_custom(ggplotGrob(inset_map),
                    xmin = 2.95e5, xmax = 3.5e5,
                    ymin = 15e4, ymax = 17.5e4)
ggsave("figs/cow-map.png", width = 12, height = 8)

### Lake Scanewa map -----------------------------------------------------------
scan_wb <- cow_wb |>
  filter(gnis_name == "Lake Scanewa")

scan_bb <- st_bbox(scan_wb)
lon_mean <- mean(scan_bb[c(1, 3)])
lon_range <- diff(scan_bb[c(1, 3)])
lat_mean <- mean(scan_bb[c(2, 4)])
lat_range <- diff(scan_bb[c(2, 4)])
expand <- 0.05
lon_expand <- c(-1, 1) *  expand * lon_range + scan_bb[c(1, 3)]
lat_expand <- c(-1, 1) * expand * lat_range + scan_bb[c(2, 4)]
scan_bb_exp <- c(lon_expand[1], lat_expand[1], lon_expand[2], lat_expand[2])

scan_dup <- c(-122.09509, 46.48144) |>
  st_point(dim = "XY") |>
  st_sfc(crs = 4326) |>
  st_sf() |>
  st_transform(crs = wa_crs) |>
  mutate(name = "Day Use Park")

cf_dam <- matrix(
  c(-122.108854, 46.467103,
    -122.108777, 46.465617),
  nrow = 2, ncol = 2, byrow = TRUE) |>
  st_linestring(dim = "XY") |>
  st_sfc(crs = 4326) |>
  st_transform(crs = wa_crs) |>
  st_sf() |>
  mutate(name = "Cowlitz Falls Dam")

cf_coll <- matrix(
  c(-122.107785, 46.467101,
    -122.108776, 46.467087),
  ncol = 2, nrow = 2, byrow = TRUE) |>
  st_linestring(dim = "XY") |>
  st_sfc(crs = 4326) |>
  st_transform(wa_crs) |>
  st_sf() |>
  mutate(name = "Fish Collector")

cf_labs <- list(
  c(-122.079298, 46.483492),
  c(-122.074018, 46.470744)) |>
  map(st_point, dim = "XY") |>
  st_sfc(crs = 4326) |>
  st_transform(wa_crs) |>
  st_sf() |>
  mutate(name = c("Cowlitz River", "Cispus River"))


ggplot() +
  geom_sf(data = cow_area, fill = "skyblue", color = "skyblue") +
  geom_sf(data = scan_wb, fill = "skyblue", color = "skyblue") +
  geom_sf_text(data = scan_wb, aes(label = gnis_name),
               nudge_x = -725, nudge_y = 150) +
  geom_sf(data = scan_dup, size = 5) +
  geom_sf_text(data = scan_dup, aes(label = name), nudge_y = 75) +
  geom_sf(data = cf_dam, linewidth = 4, fill = "black") +
  geom_sf_text(data = cf_dam, aes(label = name), nudge_y = -200) +
  geom_sf(data = cf_coll, linewidth = 3, color = "firebrick3") +
  geom_sf_text(data = cf_coll, aes(label = name), nudge_y = 50) +
  geom_sf_text(data = cf_labs, aes(label = name)) +
  coord_sf(xlim = scan_bb_exp[c(1, 3)],
           ylim = scan_bb_exp[c(2, 4)],
           expand = TRUE,
           crs = wa_crs) +
  theme_minimal() +
  theme(axis.title = element_blank())
ggsave("figs/scan-map.png", width = 12, height = 8)
