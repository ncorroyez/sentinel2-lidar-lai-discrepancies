# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Load necessary libraries
library(lidR)
library(lidRviewer)
library(raster)
library(ggplot2)
library(rgl)

# Setup
figures_dir <- "../04_FIGURES"

# Load the LAS file
las <- readLAS("/home/corroyez/Documents/NC_Full/01_DATA/Aigoual/LiDAR/2-las_utm/LAS_742000_6336000.las")

# Create a DTM from the LiDAR data
dtm <- grid_terrain(las, res = 1, algorithm = tin())
# writeRaster(dtm, "path_to_save_dtm.tif", format = "GTiff"

# Plot the DTM
# png(file.path(figures_dir, "DTM.png"))
# plot(dtm, main = "Digital Terrain Model (DTM)")
# dev.off()

# Create a DSM from the LiDAR data
dsm <- grid_canopy(las, res = 1, algorithm = p2r())

# Plot the raw LiDAR data
# png(file.path(figures_dir, "Raw_LiDAR_Data.png"))
# plot(las, color = "Z", size = 4, bg = "white", main = "Raw LiDAR Data")
x <- lidR::plot(las, bg = "white", size = 3, backend = 'rgl', legend = TRUE)
add_dtm3d(x, dtm)
# lidR::plot(las) |> add_dtm3d(dtm)
# dev.off()

# Normalize the LiDAR data using the DTM
las_normalized_dtm <- normalize_height(las, dtm)

# Plot the normalized DTM LiDAR data
png(file.path(figures_dir, "Normalized_LiDAR_Data.png"))
plot(las_normalized_dtm, color = "Z", size = 1, bg = "white", main = "Normalized DTM LiDAR Data", legend = TRUE)
dev.off()

# Plot the normalized translated DTM LiDAR data
las_normalized_dtm$Z <- las_normalized_dtm$Z + 40 - max(las_normalized_dtm$Z)
plot(las_normalized_dtm, color = "Z", size = 1, bg = "white", main = "Normalized DTM Translated LiDAR Data", legend = TRUE)

# Normalize the LiDAR data using the DSM
las_normalized_dsm <- normalize_height(las, algo = dsm)

# Plot the normalized DSM LiDAR data
png(file.path(figures_dir, "Normalized_LiDAR_Data.png"))
plot(las_normalized_dsm, color = "Z", size = 1, bg = "white", main = "Normalized DSM LiDAR Data", legend = TRUE)
dev.off()

# Plot the normalized shifted DSM LiDAR data
las_normalized_dsm$Z <- las_normalized_dsm$Z + 40
plot(las_normalized_dsm, color = "Z", size = 1, bg = "white", main = "Normalized DSM Translated LiDAR Data", legend = TRUE)
