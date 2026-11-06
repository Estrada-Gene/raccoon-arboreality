#### SETUP ####
#open raccoon-arboreality.Rproj before running scripts
required_files <- c("data/raw/rac.xlsx", "data/raw/Raccoon_GEE_Export_7Jun26.csv",
  "data/derived/gbif_effort_grid.csv", "data/derived/grid_definition.rds",
  "data/derived/era5_temperature_cache.rds")
if (!file.exists("raccoon-arboreality.Rproj")) stop("Open the RStudio Project first.")
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = ", "))
manifest <- read.csv("data/metadata/input_manifest.csv", stringsAsFactors = FALSE)
hashes <- unname(tools::md5sum(manifest$file))
if (anyNA(hashes) || !identical(hashes, manifest$md5)) {
  stop("An input file differs from the packaged snapshot. Check data/metadata/input_manifest.csv.")
}
#install missing analysis packages; rgee/reticulate are optional for fresh extraction
packages <- c("tidyverse", "sf", "readxl", "rgbif", "data.table", "rnaturalearth",
  "rnaturalearthdata", "patchwork", "lubridate", "suncalc", "spaMM", "spdep",
  "DHARMa", "corrplot", "mgcv", "MASS", "scales")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) install.packages(missing_packages, repos = "https://cloud.r-project.org")
if (any(!vapply(packages, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Some packages did not install. Resolve their installation errors before fitting models.")
}
if (utils::packageVersion("ggplot2") < "3.5.0") stop("Update ggplot2 to version 3.5.0 or later.")
for (x in c("data", "models", "tables", "figures", "diagnostics", "gbif")) {
  dir.create(file.path("revision_equalarea", x), recursive = TRUE, showWarnings = FALSE)
}
for (x in list.files("scripts", pattern = "^[0-9].*\\.R$", full.names = TRUE)) parse(file = x)
capture.output(sessionInfo(), file = "revision_equalarea/diagnostics/session_setup.txt")
message("Inputs and script syntax checked. Run scripts 01-05 in order.")

