# ---
# title: "main_viz_norms.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-11-15"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("../libraries/functions_lidar.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_chm.R")

# Parellel
# plan(multisession)

# Directories
data <- "/media/corroyez/My Passport/01_DATA"
sites <- c("Mormal") # Mormal Blois Aigoual
# sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual", "Blois")
lidar <- "LiDAR"
results <- "../../03_RESULTS"

for (site in sites){

  # Directories setup
  data_site <- file.path(data, site)
  # UTM
  las_utm <- file.path(lidar,"2-las_utm")
  las_dir <- file.path(data_site, las_utm)
  las_utm_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_dir, showWarnings = FALSE, recursive = TRUE)

  # Normalized UTM by DTM
  las_norm_utm <- file.path(lidar,"3-las_normalized_utm")
  las_norm_dir <- file.path(data_site, las_norm_utm)
  las_norm_utm_files <- list.files(las_norm_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_norm_dir, showWarnings = FALSE, recursive = TRUE)

  # Normalized UTM by DSM
  las_norm_utm_dsm <- file.path(lidar,"4-las_normalized_utm_dsm")
  las_norm_dsm_dir <- file.path(data_site, las_norm_utm_dsm)
  las_norm_dsm_files <- list.files(las_norm_dsm_dir, pattern = "\\.las$", full.names = TRUE)
  dir.create(path = las_norm_dsm_dir, showWarnings = FALSE, recursive = TRUE)

  utm <- readLAS(las_utm_files[[20]])
  dtm <- readLAS(las_norm_utm_files[[20]])
  dsm <- readLAS(las_norm_dsm_files[[20]])
  print(summary(utm@data$Z))
  print(summary(dtm@data$Z))
  print(summary(dsm@data$Z))
  plot(utm)
  plot(dtm)
  plot(dsm)
}
