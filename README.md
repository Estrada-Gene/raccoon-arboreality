# Raccoon arboreality

Code and input data for the North American raccoon arboreality study. The workflow
rebuilds the equal-area grid and fits the spatial, diel, seasonal, and temperature
analyses from the coded observation workbook and saved environmental extractions.

## Start here

1. Open `raccoon-arboreality.Rproj` in RStudio. This sets the working directory to
   the repository root; do not change it to `scripts/`.
2. Run `source("scripts/00_setup.R")` once. This checks input checksums, installs
   missing analysis packages, and checks script syntax without fitting models.
3. Run the five numbered analysis scripts in the order below, using the Source
   button or `source()` from the console. Each reads the saved inputs it needs.

| Order | Script | Main saved products |
|---|---|---|
| 1 | `scripts/01_prep.R` | Equal-area grid, prepared observations, grid summaries, sampling map |
| 2 | `scripts/02_spatial.R` | Spatial GLMMs, coefficients, conditional curves, diagnostics |
| 3 | `scripts/03_diel.R` | Parsed timestamps, solar variables, descriptive diel summaries |
| 4 | `scripts/04_seasonal.R` | Seasonal GAMs, geographic-basis sensitivity, temporal curves |
| 5 | `scripts/05_temperature.R` | Temperature GAMs, odds-ratio contrasts, sensitivity comparisons |

```r
source("scripts/00_setup.R")
source("scripts/01_prep.R")
source("scripts/02_spatial.R")
source("scripts/03_diel.R")
source("scripts/04_seasonal.R")
source("scripts/05_temperature.R")
```

Alternatively, `source("scripts/run_all.R")` runs the complete workflow. Model
fitting can take substantial time and memory. This command fits models; it is not
needed merely to prepare the repository or publish the files.

## Folder layout

- `scripts/`: numbered analysis scripts, setup, and optional reconstruction code.
- `data/raw/`: coded observations and the environmental GEE export.
- `data/derived/`: frozen GBIF effort counts, their grid definition, and the
  completed ERA5-Land extraction cache used as analysis inputs.
- `data/metadata/`: input checksums, script provenance, download key, and available
  extraction metadata and reference R sessions.
- `gee/`: environmental extraction JavaScript and instructions.
- `revision_equalarea/`: generated data, models, tables, figures, diagnostics, and
  optional large downloads. The scripts create this folder; Git ignores it.

The repository does not need previously fitted model RDS files, the multi-GB GBIF
archive, legacy script versions, or the historical temperature request and
intermediate files. Those are replaced by the supplied compact input snapshots.

## Extraction snapshots and external services

The default workflow uses saved environmental extractions and effort counts.
It does not require a GBIF account or Earth Engine authentication. Internet
access may still be needed for package installation and basemap resources.

The GBIF cache contains the original download's counts in the study's fixed
50-km equal-area cells. Script 01 checks the projection, cell size, origin,
assignment rule, and cell IDs before using it. The optional
`scripts/optional/rebuild_gbif_effort.R` reconstructs counts from the recorded
download key after script 01 has created the grid. It downloads the large original
archive and writes a separate rebuilt table under `revision_equalarea/data/`;
it does not replace the supplied input snapshot.

Script 05 checks the temperature cache against observation IDs, coordinates, and
rounded UTC hours. Records without complete temperatures remain in the cache so
the original exclusions can be reproduced. It uses a working copy under the
output folder. If changed inputs require fresh extraction, the existing rgee
extraction block also requires `rgee`, `reticulate`, its configured Python
environment, Earth Engine access, and `GEE_PROJECT`. These are not required for
the supplied observation snapshot.

See `DATA_SOURCES.md` for provenance and `docs/packaging_changes.md` for the limited
portability changes made while copying the local scripts. Current local script
edits, including figure settings, are copied into this project.

## Software versions

Setup installs missing packages but does not pin their versions. Available
reference `sessionInfo()` files are retained under `data/metadata/reference_sessions/`.
New sessions are saved in the generated diagnostics folder. After confirming that
the project works in the intended R environment, a lockfile can be added with
`install.packages("renv")` followed by `renv::init()`. Commit the resulting
`renv.lock`, `.Rprofile`, and renv activation files; local package libraries stay
ignored. No lockfile is supplied with unverified package versions.

## Publication and reuse

Repository URL: `[ADD GITHUB REPOSITORY URL]`.
Select an appropriate code license before public release. Upstream data retain
their applicable source terms; a code license does not relicense third-party data.
Photographs are not included. The coded observations contain source identifiers
and links so records can be traced to iNaturalist.

CSV and figure filenames retain the analysis scripts' existing output names.
Their supplementary numbers may differ from the final manuscript numbering;
match exported content to the final captions when preparing the submission.

