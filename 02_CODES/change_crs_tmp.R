# ---
# title: "change_crs_tmp.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-11-29"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(terra)
library(sgsR)
library(ggplot2)
library(sf)
library(sp)
library(spdep)
library(dplyr)
library(geojsonio)
library(ranger)
library(VSURF)
library(progress)
library(corrplot)
library(patchwork)

source("libraries/functions_general_tools.R")
source("libraries/functions_plots.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual", "Blois")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/variogram"
pp_dir <- "../04_FIGURES/pp"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
set.seed(42)

# Lists
for (site in sites){
  date <- switch(site,
                 "Aigoual" = '2021-07-11',
                 "Blois" = '2021-06-14',
                 "Mormal" = '2021-06-14')
  lai_path <- paste0("/home/corroyez/Downloads/PROSAIL_Inversions/03_RESULTS/", 
                     site, "/HybridInversion_2/", site, "_",
                     date, "_lai.tiff")
  lai_raster <- terra::rast(lai_path)
  good_lai_crs_path <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/",
                              site, "/Metrics/Deciduous_Only/s2lai_summer_res_10_m.tif")
  good_lai_crs <- terra::rast(good_lai_crs_path)
  
  lai_raster_reprojected <- terra::project(lai_raster, good_lai_crs)
  lai_raster_cropped <- terra::crop(lai_raster_reprojected, good_lai_crs)
  lai_raster_masked <- terra::mask(lai_raster_cropped, good_lai_crs)
  
  output_path <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/",
                        site, "/Metrics/Deciduous_Only/s2lai_summer_v2_res_10_m.tif")
  terra::writeRaster(lai_raster_masked, output_path, overwrite = TRUE)
}