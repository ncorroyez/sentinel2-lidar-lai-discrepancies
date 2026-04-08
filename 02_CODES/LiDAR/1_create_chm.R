# ---
# title: "0_create_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library("lidR")
library("terra")
library("future")
library("future.apply")

# --------------------------- Import useful functions --------------------------

source("../libraries/functions_chm.R")

# --------------------------------- Setup --------------------------------------
# data_path <- "../../01_DATA"
data_path <- "/media/corroyez/MyPassport/01_DATA"
results_dir <- "../../03_RESULTS"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois", "Mormal")
sites <- c("Reine", "Hayes")
sites <- "Hayes"
# sites <- "Mormal" # Aigoual Blois Mormal
resolutions <- 1 # 1
# resolutions <- c(10, 1)
leaf_state <- "leaf_off"
# -------------------------------- Process -------------------------------------
# for (site in sites){
#   cat("Processing site:", site, "\n")
#   for (resolution in resolutions){
#     cat("Processing resolution:", resolution, "\n")
#     create_chm(data_path, results_dir, site, resolution)
#   }
# }
# plan(multisession(workers = 4))
plan(sequential)
for (resolution in resolutions) {
  cat("Processing resolution:", resolution, "\n")
  for (site in sites) {
    cat("Processing site:", site, "\n")
    s2_rast <- switch(site,
                      "Aigoual" = terra::rast(file.path(results_dir, site, "L2A_T31TEJ_A031608_20210711T104217/Reflectance/res_10_m/L2A_T31TEJ_A031608_20210711T104217_Refl"))[[1]],
                      "Blois" = terra::rast(file.path(results_dir, site, "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31TCN_A031222_20210614T105443_Refl"))[[1]],
                      "Mormal" = terra::rast(file.path(results_dir, site, "L2A_T31UER_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31UER_A031222_20210614T105443_Refl"))[[1]],
                      "Hayes" = rast(file.path(data_path, site, "Sentinel-2", "2023-12-03/raster_samples", "Hayes_001_2023-12-03.tiff")),
                      # "Reine" = rast(file.path(data_path, site, "Sentinel-2", "2024-03-07/raster_samples", "Reine_001_2024-03-07.tiff"))
                      # "Hayes" = rast(file.path(data_path, site, "Sentinel-2", "2024-08-11/raster_samples", "Hayes_001_2024-08-11.tiff")),
                      "Reine" = rast(file.path(data_path, site, "Sentinel-2", "2024-07-30/raster_samples", "Reine_001_2024-07-30.tiff"))
    )
    create_chm(data_path, results_dir, site, s2_rast, leaf_state, resolution)
  }
}
