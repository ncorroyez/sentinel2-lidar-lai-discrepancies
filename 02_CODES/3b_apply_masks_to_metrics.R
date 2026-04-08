# ---
# title: "3b_apply_masks_to_metrics"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-14"
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
library("lmom")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Import useful functions
source("libraries/functions_create_masks.R")

# Directories 
# sites <- "Aigoual" # Mormal Blois Aigoual
sites <- c("Aigoual", "Blois", "Mormal")

for (site in sites){
  cat("Processing site:", site, "\n")
  results_path <- file.path("../03_RESULTS", site)
  masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
  
  # Output
  results_path <- file.path("../03_RESULTS", site)
  metrics_dir <- file.path(results_path, "Metrics")
  raw_dir <- file.path(metrics_dir, "Raw")
  pad_dir <- file.path(raw_dir, "PAD_Profiles")
  
  metrics <- list.files(raw_dir, pattern = ".tif$")
  pad_files <- list.files(pad_dir, pattern = ".tif$")
  
  for(metric_name in metrics){
    cat("Processing metric:", metric_name, "\n")
    metric_raster <- terra::rast(file.path(raw_dir, metric_name))
    apply_and_save_masks(metric_raster,
                         metric_name,
                         masks_dir,
                         metrics_dir)
  }
  for(pad_file in pad_files){
    cat("Processing PAD Profile:", pad_file, "\n")
    pad_raster <- terra::rast(file.path(raw_dir, 
                                        "PAD_Profiles",
                                        pad_file))
    apply_and_save_masks(pad_raster,
                         pad_file,
                         masks_dir,
                         metrics_dir,
                         pad_bool = TRUE)
  }
}
