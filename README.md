# Raccoon arboreality

Code and data for Greiner et al. Latitude and anthropogenic modification predict geographic variation in raccoon arboreality. 
The code here rebuilds the equal-area grid and fits the spatial, diel, seasonal, and temperature
analyses from the coded raccoon observation table and environmental variables.

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

## Software versions

Setup installs missing packages but does not pin their versions. Available
reference `sessionInfo()` files are retained under `data/metadata/reference_sessions/`.
