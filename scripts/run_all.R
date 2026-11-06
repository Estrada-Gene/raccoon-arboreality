#### RUN ALL ANALYSES ####
#fits every manuscript model; run individual scripts for step-by-step work
if (!file.exists("raccoon-arboreality.Rproj")) stop("Open the RStudio Project first.")
scripts <- c("00_setup.R", "01_prep.R", "02_spatial.R", "03_diel.R",
             "04_seasonal.R", "05_temperature.R")
for (script in scripts) {
  message("Running ", script)
  #isolate scripts while keeping package attachment behavior unchanged
  source(file.path("scripts", script), local = new.env(parent = globalenv()), echo = FALSE)
}

