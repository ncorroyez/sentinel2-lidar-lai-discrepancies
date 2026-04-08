# ---
# title: "1.create_masks.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
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

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_lidar.R")
source("../libraries/functions_create_masks.R")

# Pre-processing Parameters & Directories
# data_dir <- '../../01_DATA'
data_dir <- '/media/corroyez/My Passport/01_DATA'
results_dir <- '../../03_RESULTS'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Mormal" # Mormal Blois Aigoual

# Main
for (site in sites){
  shapefiles_dir <- file.path(data_dir, site, "Shapefiles")
  results_path <- file.path(results_dir, site)
  masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
  create_vegetation_forest_mask(data_dir,
                                site,
                                shapefiles_dir,
                                masks_dir,
                                results_path,
                                resolution = 10)
}
