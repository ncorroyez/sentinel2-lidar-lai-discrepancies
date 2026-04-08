# ---
# title: "2.heterogeneity_lai_lidar.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-11-20"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("data.table")
library("raster")
library("rgdal")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

results_path <- "../../03_RESULTS/Mormal"
deciles_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks/Deciles/")
pai_path <- file.path(results_path, "LiDAR/PAI/lidR/pai_las_cleaned.tif")
forest_path <- file.path(results_path, 'LiDAR/Heterogeneity_Masks/L2A_T31UER_A031222_20210614T105443_forest.envi')
low_vegetation_path <- file.path(results_path, 'LiDAR/Heterogeneity_Masks/mne_summed_95_p_kept_new.tif')

forest_mask <- open_raster_file(forest_path)
low_vegetation_mask <- open_raster_file(low_vegetation_path)

pai_las <- terra::rast(pai_path)
pai_las_forest <- terra::project(pai_las, low_vegetation_mask)
pai_las_forest <- mask(pai_las_forest, low_vegetation_mask)
pai_las_low_vegetation <- terra::project(pai_las, forest_mask)
pai_las_low_vegetation <- mask(pai_las_low_vegetation, forest_mask)

# Convert raster to a matrix
pai_matrix <- terra::as.matrix(pai_las_cleaned)

# Flatten the matrix to a vector
pai_values <- as.vector(pai_matrix)

# Calculate the deciles
deciles <- quantile(pai_values, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE)

# Display the deciles
print(deciles)

inc <- 0.1
for (low_quantile in seq(0.0, 0.99, by = inc)) {
  quantiles <- terra::rast(file.path(deciles_dir, 
                                    paste0("heter_binary_mask_", 
                                           low_quantile, 
                                           "_", 
                                           low_quantile + inc,
                                           ".tif")))
  pai_resample <- resample(pai_las_low_vegetation, quantiles)
  pai_masked <- mask(pai_resample, quantiles)
  plot(pai_masked, main = paste("pai", low_quantile, "to", low_quantile + inc))
  writeRaster(pai_masked, 
              paste(deciles_dir, "pai_masked_", low_quantile, 
                    "_", low_quantile + inc, ".tif", sep = ''),
              overwrite = T)
}