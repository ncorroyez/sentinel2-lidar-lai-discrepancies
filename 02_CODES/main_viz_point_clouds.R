# ---
# title: "main_viz_point_clouds.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-12-09"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(lidR)
library(terra)
library(sf)
library(spdep)

sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois", "Mormal")
sites <- c("Aigoual")
data_dir <- "/media/corroyez/My Passport/01_DATA"
results_dir <- "../03_RESULTS"
metrics_dir <- "Metrics/Not_Masked"

# Loop over sites
for (site in sites) {
  las_utm_dir <- file.path(data_dir, site, "LiDAR", "3-las_normalized_utm")
  ctg <- readLAScatalog(las_utm_dir)
  lcv <- rast(file.path(results_dir, site, metrics_dir, "lcv_res_10_m.tif"))
  
  # Extract centroids of LAS file geometries
  ctg_centroids <- st_centroid(st_geometry(ctg@data))
  
  # Extract Lcv values at the centroid coordinates
  coords <- st_coordinates(ctg_centroids)
  ctg@data$Lcv <- extract(lcv, coords)[, 1]
  
  # Find the LAS file with the highest and lowest Lcv
  highest_lcv_row <- which.max(ctg@data$Lcv)
  lowest_lcv_row <- which.min(ctg@data$Lcv)
  
  highest_lcv_file <- ctg@data$filename[highest_lcv_row]
  lowest_lcv_file <- ctg@data$filename[lowest_lcv_row]
  
  # Output results
  cat("Highest Lcv file:", highest_lcv_file, "\n")
  cat("Lowest Lcv file:", lowest_lcv_file, "\n")
  print(ctg@data$Lcv)
  # aa
  # Open the files
  las_high <- readLAS(highest_lcv_file)
  las_low <- readLAS(lowest_lcv_file)
  las_high_clip <- clip_circle(las_high, 
                               xcenter = (las_high@header$`Max X` + las_high@header$`Min X` ) / 2, 
                               ycenter = (las_high@header$`Max Y` + las_high@header$`Min Y` ) / 2, 
                               radius = 100)
  las_low_clip <- clip_circle(las_low, 
                              xcenter = (las_low@header$`Max X` + las_low@header$`Min X` ) / 2, 
                              ycenter = (las_low@header$`Max Y` + las_low@header$`Min Y` ) / 2, 
                              radius = 100)
  plot(las_high_clip, main = "Highest Lcv File")
  plot(las_low_clip, main = "Lowest Lcv File")
}
