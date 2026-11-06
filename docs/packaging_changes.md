# Packaging changes

The setup utility copies the latest revised scripts in the local `coond` folder.
It does not use the older `coon*.R` scripts and does not run analyses.

Changes to copied files are limited to:

1. Numbered filenames under `scripts/` and removal of `setwd()` lines.
2. Relative paths to the coded workbook and saved environmental export.
3. Keeping generated results under the ignored `revision_equalarea/` folder.
4. Reading the original GBIF effort counts from a compact CSV. The original
   download and aggregation section is preserved as an optional reconstruction
   script using the saved download key.
5. Loading the completed temperature cache instead of rebuilding it from legacy
   request/intermediate RDS files. Extraction updates write to an output copy.
6. Replacing a personal Earth Engine project string with `GEE_PROJECT`, used only
   when fresh extraction is needed.
7. Disabling a leftover `compare_ML <- TRUE` setting if it is still present, in
   accordance with the manuscript's omission of that optional comparison.

Behavioral coding rules, grid geometry, eligibility filters, predictors, model
formulas, fitting settings, simulations, and locally edited figure code are
otherwise retained. Before writing project files, the utility parses every
copied R script and checks the required cache columns and grid metadata. Setup
checks hashes and script syntax without fitting models. These checks do not
constitute an end-to-end model rerun; compare the regenerated sample accounting
and key results when performing that verification in R.

