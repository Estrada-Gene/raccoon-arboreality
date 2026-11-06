#### REBUILD GBIF EFFORT FROM THE ORIGINAL DOWNLOAD ####
#optional; first run 01_prep.R with its supplied small effort cache
library(tidyverse)
library(sf)
library(rgbif)
library(data.table)
out <- "revision_equalarea"
definition <- readRDS(file.path("data", "derived", "grid_definition.rds"))
grid_crs <- definition$crs
grid_size <- definition$cellsize
study_grid <- readRDS(file.path(out, "data", "study_grid_50km.rds"))
occupied_cells <- st_drop_geometry(study_grid)
dir.create(file.path(out, "gbif"), recursive = TRUE, showWarnings = FALSE)

#### GBIF DOWNLOAD ####

#request bounds
request_boundary <- st_as_sfc(st_bbox(study_grid)) %>%
  st_segmentize(dfMaxLength = 10000) %>%
  st_transform(4326)

bb <- st_bbox(request_boundary)

request_bounds <- c( xmin = max(-180, unname(bb["xmin"]) - 0.1),
                     xmax = min(180, unname(bb["xmax"]) + 0.1),
                     ymin = max(-90, unname(bb["ymin"]) - 0.1),
                     ymax = min(90, unname(bb["ymax"]) + 0.1))

write_csv(as_tibble_row(request_bounds), file.path(out, "tables", "gbif_request_bounds.csv"))

#reuse the key saved by a previous run
key_file <- file.path(out, "gbif", "download_key.txt")
gbif_key <- trimws(readLines(file.path("data", "metadata", "gbif_download_key.txt"), warn = FALSE)[1])

if (file.exists(key_file)) {
  gbif_key <- trimws(paste(readLines(key_file, warn = FALSE), collapse = ""))
}

if (!nzchar(gbif_key)) {
  
  credentials <- Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))
  
  if (any(!nzchar(credentials))) {
    stop("Set GBIF_USER, GBIF_PWD and GBIF_EMAIL in .Renviron, then restart R.")
  }
  
  gbif_dl <- occ_download(
    pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
    pred("hasCoordinate", TRUE),
    pred_gte("decimalLatitude", request_bounds["ymin"]),
    pred_lte("decimalLatitude", request_bounds["ymax"]),
    pred_gte("decimalLongitude", request_bounds["xmin"]),
    pred_lte("decimalLongitude", request_bounds["xmax"]),
    user = credentials["GBIF_USER"],
    pwd = credentials["GBIF_PWD"],
    email = credentials["GBIF_EMAIL"],
    format = "SIMPLE_CSV"
  )
  
  gbif_key <- as.character(gbif_dl)
  writeLines(gbif_key, key_file)
  rm(credentials)
}

occ_download_wait(gbif_key)

gbif_file <- occ_download_get(gbif_key, path = file.path(out, "gbif"), overwrite = FALSE)

#save metadata, including DOI and query
saveRDS(occ_download_meta(gbif_key), file.path(out, "gbif", "download_metadata.rds"))

inat_raw <- occ_download_import(gbif_file, 
                                select = c("decimalLongitude", "decimalLatitude", "class")) %>%
  as.data.table()

inat_raw <- inat_raw[is.finite(decimalLongitude) & is.finite(decimalLatitude) &
                       decimalLongitude >= -180 & decimalLongitude <= 180 &
                       decimalLatitude >= -90 & decimalLatitude <= 90]


#### GBIF EFFORT AGGREGATION ####

#transform in chunks
chunk_size <- 250000L
starts <- seq_len(ceiling(nrow(inat_raw) / chunk_size))
effort_parts <- vector("list", length(starts))

if (!nrow(inat_raw)) stop("No valid GBIF coordinates.")

for (j in starts) {
  first <- (j - 1L) * chunk_size + 1L
  last <- min(j * chunk_size, nrow(inat_raw))
  part <- copy(inat_raw[first:last])
  
  pts <- st_as_sf(
    as.data.frame(part),
    coords = c("decimalLongitude", "decimalLatitude"), crs = 4326
  ) %>%
    st_transform(grid_crs)
  
  coords <- st_coordinates(pts)
  
  part[, grid_id := paste(
    floor(coords[, 1] / grid_size),
    floor(coords[, 2] / grid_size), sep = "_"
  )]
  
  part <- part[grid_id %in% occupied_cells$grid_id]
  
  effort_parts[[j]] <- part[, .(
    effortRG = .N,
    effortRG_mammals = sum(class == "Mammalia", na.rm = TRUE)
  ), by = grid_id]
  
  message("Effort chunk ", j, "/", length(starts))
}

effort_dt <- rbindlist(effort_parts)[, .(
  effortRG = sum(effortRG),
  effortRG_mammals = sum(effortRG_mammals)
  ), by = grid_id]

rm(inat_raw, effort_parts, part, pts, coords)
gc()


#write a separate result; preserve the supplied input snapshot
write_csv(as_tibble(effort_dt), file.path(out, "data", "gbif_effort_grid_rebuilt.csv"))
message("Rebuilt effort table saved under revision_equalarea/data/.")

