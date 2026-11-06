#### RACCOON PREP ####

library(tidyverse)
library(sf)
library(readxl)
library(rgbif)
library(data.table)
library(rnaturalearth)
library(patchwork)

#### SETUP ####

#output directory
out <- "revision_equalarea"

for (x in c("data", "tables", "figures", "diagnostics", "gbif")) {
  dir.create(file.path(out, x), recursive = TRUE, showWarnings = FALSE)
}

#establish equal-area projection centered on North America
grid_crs <- paste("+proj=laea +lat_0=40 +lon_0=-100",
                  "+datum=WGS84 +units=m +no_defs")

grid_size <- 50000 #50km grid size
set.seed(1120)


#### READ & CLEAN RACCOON DATA ####

temp_df <- read_excel(file.path("data", "raw", "rac.xlsx"), n_max = 1)
types <- ifelse(grepl("field:", names(temp_df)), "text", "guess")

df <- read_excel(file.path("data", "raw", "rac.xlsx"), col_types = types) %>%
  mutate(id = as.character(id),
         across(starts_with("field:"), ~na_if(trimws(.x), "")))

#filtering stages for sample accounting
df_coords <- df %>%
  filter(is.finite(longitude), is.finite(latitude),
         between(longitude, -180, 180),
         between(latitude, -90, 90))

df_noswim <- df_coords %>%
  filter(is.na(`field:swimming.`))

#existing behavioral rules and geographic filter
df_clean <- df_noswim %>%
  mutate(arboreal = case_when(
    !is.na(`field:inside tree`) ~ 1L,
    !is.na(`field:vertical`) ~ 1L,
    !is.na(`field:arboreal.`) ~ 1L,
    !is.na(`field:in water`) ~ 0L,
    !is.na(`field:terrestrial`) ~ 0L,
    TRUE ~ 0L)) %>%
  filter(between(longitude, -170, -30), between(latitude, -5, 85))

#sample flow
sample_flow <- tibble(stage = c("Source workbook", "Valid coordinates",
            "Swimming excluded", "Geographic filter"),
            n = c(nrow(df), nrow(df_coords), nrow(df_noswim), nrow(df_clean)))

write_csv(sample_flow, file.path(out, "tables", "Table_S1_prep_flow.csv"))

#review extreme southern records
df_clean %>%
  filter(latitude < 25) %>%
  select(id, place_guess, latitude, longitude, url) %>%
  write_csv(file.path(out, "tables", "southern_records_review.csv"))

#field values for checking behavioral definitions
df_clean %>%
  select(id, starts_with("field:")) %>%
  pivot_longer(-id, names_to = "field", values_to = "value") %>%
  filter(!is.na(value)) %>%
  count(field, value, sort = TRUE) %>%
  write_csv(file.path(out, "tables", "coding_field_values.csv"))

#coordinate-quality summary
df_clean %>%
  summarise(
    n = n(),
    obscured = sum(tolower(as.character(coordinates_obscured)) == "true", na.rm = TRUE),
    accuracy_missing = sum(is.na(positional_accuracy)),
    accuracy_gt_1km = sum(positional_accuracy > 1000, na.rm = TRUE),
    accuracy_gt_10km = sum(positional_accuracy > 10000, na.rm = TRUE),
    accuracy_gt_50km = sum(positional_accuracy > 50000, na.rm = TRUE)) %>%
  write_csv(file.path(out, "tables", "Table_S1_coordinate_quality.csv"))


#### JOIN ENVIRONMENTAL DATA ####

#preserve IDs and restrict join to environmental columns
env_cols <- c("human_modification", "canopy_height_m", "tree_cover_pct", "peak_ndvi",
              "bio01", "bio02", "bio03", "bio04", "bio07", "bio12", "bio15")

df_gee <- read_csv(file.path("data", "raw", "Raccoon_GEE_Export_7Jun26.csv"),
                   col_types = cols(id = col_character()),
                   show_col_types = FALSE) %>%
  select(id, all_of(env_cols)) %>%
  mutate(environment_match = TRUE)

df_clean <- df_clean %>%
  left_join(df_gee, by = "id") %>%
  mutate(environment_match = replace_na(environment_match, FALSE))

#unmatched observations
df_clean %>%
  filter(!environment_match) %>%
  select(id, latitude, longitude) %>%
  write_csv(file.path(out, "tables", "environment_unmatched.csv"))

#rm(df, df_coords, df_noswim, df_gee, temp_df)


#### CREATE GRID SYSTEM ####

#project to equal-area coordinates
rac_projected <- st_as_sf(df_clean, coords = c("longitude", "latitude"), 
                          crs = 4326, remove = FALSE) %>% 
  st_transform(grid_crs)

xy <- st_coordinates(rac_projected)

#fixed origin, assign each observation to one cell
obs_with_geo <- df_clean %>%
  mutate(grid_col = as.integer(floor(xy[, 1] / grid_size)),
         grid_row = as.integer(floor(xy[, 2] / grid_size)),
         grid_id = paste(grid_col, grid_row, sep = "_"))

occupied_cells <- obs_with_geo %>%
  distinct(grid_id, grid_col, grid_row) %>%
  arrange(grid_row, grid_col)

#create polygons for occupied cells
make_cell <- function(col, row) {
  x <- col * grid_size
  y <- row * grid_size
  
  st_polygon(list(rbind(
    c(x, y), c(x + grid_size, y),
    c(x + grid_size, y + grid_size),
    c(x, y + grid_size), c(x, y)
  )))
}

study_grid <- st_sf(occupied_cells,
                    geometry = st_sfc(Map(make_cell, occupied_cells$grid_col, occupied_cells$grid_row),
                    crs = grid_crs))

grid_coords <- occupied_cells %>%
  transmute(grid_id,
            X_coord = (grid_col + 0.5) * grid_size,
            Y_coord = (grid_row + 0.5) * grid_size)

#save grid and observation assignments
saveRDS(study_grid, file.path(out, "data", "study_grid_50km.rds"))
saveRDS(obs_with_geo, file.path(out, "data", "observations_prepared.rds"))
write_csv(obs_with_geo, file.path(out, "data", "rac_clean.csv"))

saveRDS(list(crs = grid_crs, cellsize = grid_size, origin = c(0, 0),
             assignment = "floor", grid_id = "column_row"),
        file.path(out, "data", "grid_definition.rds"))

rm(rac_projected, xy, df_clean)


#### GRID-LEVEL MODEL AGGREGATION ####

#return NA when all values are missing
med <- function(x) {
  if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
}

#retain all occupied cells before spatial-model restriction
grid_summary <- obs_with_geo %>%
  group_by(grid_id) %>%
  summarise(
    n_arboreal = sum(arboreal == 1L),
    n_total = n(),
    latitude = median(latitude),
    n_environment_matches = sum(environment_match),
    human_mod = med(human_modification),
    canopy_height = med(canopy_height_m),
    tree_cover = med(tree_cover_pct),
    ndvi = med(peak_ndvi),
    mean_temp = med(bio01),
    md_range = med(bio02),
    isotherm = med(bio03),
    temp_seas = med(bio04),
    ann_range = med(bio07),
    ann_precip = med(bio12),
    precip_seas = med(bio15),
    .groups = "drop") %>%
  left_join(grid_coords, by = "grid_id")

#number of environmental values contributing to each cell
environment_support <- obs_with_geo %>%
  group_by(grid_id) %>%
  summarise(across(all_of(env_cols), ~sum(is.finite(.x)), 
                   .names = "n_{.col}"), .groups = "drop") %>%
  left_join(grid_summary %>% select(grid_id, n_total, n_environment_matches), by = "grid_id")

write_csv(environment_support, file.path(out, "tables", "Table_S2_environment_support.csv"))


#### GBIF EFFORT CACHE ####
#reuse the published download's equal-area cell counts
effort_grid <- readRDS(file.path("data", "derived", "grid_definition.rds"))
stopifnot(isTRUE(st_crs(grid_crs) == st_crs(effort_grid$crs)),
          isTRUE(all.equal(as.numeric(grid_size), as.numeric(effort_grid$cellsize))),
          identical(as.numeric(effort_grid$origin), c(0, 0)),
          identical(effort_grid$assignment, "floor"))
effort_dt <- read.csv(file.path("data", "derived", "gbif_effort_grid.csv"),
  colClasses = c(grid_id = "character"), check.names = FALSE)
stopifnot(!anyDuplicated(effort_dt$grid_id),
          setequal(effort_dt$grid_id, occupied_cells$grid_id),
          all(is.finite(effort_dt$effortRG)), all(is.finite(effort_dt$effortRG_mammals)))

#### MERGE & SAVE ####

#replace missing effort
grid_summary <- grid_summary %>%
  left_join(as_tibble(effort_dt), by = "grid_id") %>%
  mutate(across(c(effortRG, effortRG_mammals), ~replace_na(.x, 0)))

#check number of cells and obs
sum(grid_summary$n_total)
summary(grid_summary$n_total)

#minimum sample size for spatial model = 3
df_grid_model <- grid_summary %>% filter(n_total >= 3)

sum(df_grid_model$n_total)

#eligibility flag for maps and sample accounting
model_vars <- c("latitude", "human_mod", "tree_cover", "ndvi",
                "ann_precip", "precip_seas", "effortRG", "effortRG_mammals",
                "X_coord", "Y_coord")

grid_summary$spatial_eligible <- grid_summary$n_total >= 3 &
  apply(as.data.frame(grid_summary[, model_vars]), 1, function(x) all(is.finite(x)))

saveRDS(grid_summary, file.path(out, "data", "all_occupied_grid_data.rds"))
saveRDS(df_grid_model, file.path(out, "data", "raccoon_grid_model.rds"))
write_csv(df_grid_model, file.path(out, "data", "raccoon_grid_model.csv"))

#distinguish observations from model cells
bind_rows(
  grid_summary %>%
    summarise(dataset = "All occupied cells", n_cells = n(),
              n_observations = sum(n_total), n_arboreal = sum(n_arboreal)),
  grid_summary %>%
    filter(n_total >= 3) %>%
    summarise(dataset = "Cells with >=3 observations", n_cells = n(),
              n_observations = sum(n_total), n_arboreal = sum(n_arboreal)),
  grid_summary %>%
    filter(spatial_eligible) %>%
    summarise(dataset = "Complete spatial-model cells", n_cells = n(),
              n_observations = sum(n_total), n_arboreal = sum(n_arboreal))) %>%
  write_csv(file.path(out, "tables", "Table_S1_grid_samples.csv"))

#calculate correlation
effort_r <- cor(df_grid_model$effortRG, df_grid_model$effortRG_mammals,
                method = "pearson")

write_csv(tibble(comparison = "All taxa vs mammals", pearson_r = effort_r),
          file.path(out, "tables", "effort_correlation.csv"))

sample_flow


#### MAPS ####

#cell summaries
map_data <- study_grid %>%
  left_join(grid_summary, by = "grid_id") %>%
  mutate(prop_arboreal = n_arboreal / n_total)

#cell centers; larger circles first
map_points <- map_data %>%
  st_drop_geometry() %>%
  arrange(desc(n_total)) %>%
  st_as_sf(coords = c("X_coord", "Y_coord"), crs = st_crs(study_grid))

#full basemap; crop only the plotting window
world_map <- ne_countries(scale = "medium", returnclass = "sf")

#display extent
crop_bb <- st_bbox(st_transform(map_data, 4326))
crop_lon <- unname(c(crop_bb["xmin"], crop_bb["xmax"])) + c(-5, 15)
crop_lat <- unname(c(crop_bb["ymin"], crop_bb["ymax"]))

p_fig1 <- ggplot() +
  geom_sf(data = world_map, fill = "gray93", color = NA) +
  geom_sf(data = map_points,
          aes(size = n_total, color = prop_arboreal),
          shape = 16, alpha = 0.5, stroke = 0) +
  coord_sf(crs = st_crs(4326), default_crs = st_crs(4326),
           xlim = crop_lon, ylim = crop_lat, expand = FALSE) +
  scale_size_area(
    name = "observations\nper cell",
    max_size = 10,
    breaks = scales::breaks_pretty(n = 4),
    labels = scales::comma) +
  scale_color_gradient(
    name = "arboreal\nproportion",
    low = "black", high = "green",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(color = "gray30", alpha = 0.5)),
    color = guide_colorbar(order = 2)) +
  theme_void(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.50),
    legend.justification = c(1, 0.5),
    legend.box = "vertical",
    legend.box.background = element_rect(
      fill = scales::alpha("white", 0.8), color = NA),
    legend.background = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(5, 5, 5, 5),
    plot.margin = margin(8, 8, 8, 8))

p_fig1

ggsave(file.path(out, "figures", "Figure_1_sampling.jpg"),
       p_fig1, width = 9, height = 7, dpi = 300, bg = "white")

#source data 
write_csv(st_drop_geometry(map_data), file.path(out, "tables", "Figure_1_source_data.csv"))


#### FINAL CHECKS ####

stopifnot(
  sum(grid_summary$n_total) == nrow(obs_with_geo),
  sum(grid_summary$n_arboreal) == sum(obs_with_geo$arboreal),
  !anyDuplicated(grid_summary$grid_id),
  all(abs(as.numeric(st_area(study_grid)) - grid_size^2) < 1))

capture.output(sessionInfo(),  file = file.path(out, "diagnostics", "session_prep.txt"))
