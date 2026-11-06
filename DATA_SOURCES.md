# Data sources

| Input | Role and provenance |
|---|---|
| `data/raw/rac.xlsx` | Original visually coded iNaturalist observation workbook; script 01 applies the documented cleaning and behavioral rules. |
| `data/raw/Raccoon_GEE_Export_7Jun26.csv` | Saved observation-ID environmental extraction; the original JavaScript is in `gee/environmental_extraction.js`. |
| `data/derived/gbif_effort_grid.csv` | All-taxa and mammal counts copied from the completed equal-area `all_occupied_grid_data.rds`, with one row per occupied cell. No refitting or new download is used to create this table. |
| `data/derived/grid_definition.rds` | Projection, 50-km cell size, fixed origin, and floor-based cell assignment used for the effort snapshot. |
| `data/derived/era5_temperature_cache.rds` | Completed extraction cache, including unavailable pairs. This contains extraction results, not fitted models. |

`data/metadata/gbif_download_key.txt` identifies the GBIF snapshot actually used.
When available, its original metadata and request bounds are included alongside
it. The manuscript identifies the download DOI as
https://doi.org/10.15468/dl.phgwkv; verify it against the saved key and metadata if
the input download changes. A fresh occurrence download can produce different
counts and is not a substitute for this recorded snapshot.

The environmental code uses WorldClim V1 BIO12/BIO15, Hansen `treecover2000`, the
2020 forest-height layer, a June-August 2023 mean of MOD13Q1 Collection 6.1 NDVI,
and circa-2016 global human modification. Some exported bands were excluded from
the final GLMM. The environmental export is a fixed spatial reference rather
than an observation-date-matched time series. Sampling requested 30 m, which does
not increase the information resolution of coarser source layers.

ERA5-Land temperatures use the hourly collection and a 1991-2020 climatology from
the monthly-by-hour collection, matched by calendar month and UTC hour. The
temperature script retains its original exact-pixel/15-km fallback pair rules.

Source catalogs:

- https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_BIO
- https://developers.google.com/earth-engine/datasets/catalog/UMD_hansen_global_forest_change_2023_v1_11
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
- https://developers.google.com/earth-engine/datasets/catalog/CSP_HM_GlobalHumanModification
- https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY
- https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_MONTHLY_BY_HOUR

`input_manifest.csv` records file sizes and MD5 hashes for the supplied data and
metadata. `script_manifest.csv` records hashes of both the source scripts and
their portable copies; source paths are recorded as basenames, not personal
working-directory paths. Checksums identify file changes and are not validation
of the scientific measurements themselves.

