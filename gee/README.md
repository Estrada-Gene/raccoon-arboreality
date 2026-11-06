# Environmental extraction

`environmental_extraction.js` preserves the environmental GEE code from the local
working directory. The R analysis uses the saved CSV export, so normal analysis
does not require uploading or running this code again.

For a fresh extraction, upload an observation-coordinate table to Earth Engine
and import its FeatureCollection with the variable name `raccoons`, as expected
by the code. Observation IDs must remain unchanged. A coordinate table can be
written after running script 01:

```r
d <- readRDS("revision_equalarea/data/observations_prepared.rds")
readr::write_csv(dplyr::select(d, id, latitude, longitude),
                "revision_equalarea/data/gee_coordinates.csv")
```

The extraction code preserves the original reference layers, date filters, and
masking behavior. A fresh export may differ if upstream assets change. Preserve
the supplied CSV for reproducing the submitted analysis. Temperature extraction
code is retained in the unmatched-record section of `scripts/05_temperature.R`.

