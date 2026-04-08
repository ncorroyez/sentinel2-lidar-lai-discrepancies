# ---
# title: "4_create_s2_metric_time_series.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-24"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
library("readr")
library("terra")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_general_tools.R")
source("../libraries/functions_create_masks.R")

# Pre-processing Parameters
data_dir <- '/media/corroyez/My Passport/01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../../03_RESULTS"

# Read S2 Acquisition Dates (1/month)
file_path <- file.path(data_dir, "S2_Dates.csv")
data <- read_csv(file_path)

# Function to stack and save rasters
stack_and_save <- function(rasters_list, results_site, output_filename,
                           masks_dir, metrics_dir) {
  stacked_raster <- rast(rasters_list)
  writeRaster(stacked_raster, file.path(results_site, output_filename),
              overwrite = TRUE)
  cat("Stacked raster saved to:", output_filename, "\n")
  apply_and_save_masks(stacked_raster,
                       output_filename,
                       masks_dir,
                       metrics_dir)
}

for (site in sites){
  
  # Setup 
  results_dir <- '../../03_RESULTS'
  results_site <- file.path(results_dir, site)
  masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics")
  
  # Get a list of subdirectories starting with "L2A"
  subdirs <- list.dirs(results_site, recursive = FALSE)
  subdirs <- subdirs[grepl("^L2A", basename(subdirs))]
  
  # Initialize an empty list to store SpatRaster objects
  lai_rasters_list <- evi_rasters_list <- ndvi_rasters_list <- savi_rasters_list <- list()
  
  # Loop through each subdirectory
  for (subdir in subdirs) {
    # Extract the base name of the subdirectory (e.g., "L2A_T31TEJ_A022485_20210626T103707")
    dir_name <- basename(subdir)
    
    # Extract the date (first 8 characters) from the directory name after the third '_'
    date <- strsplit(dir_name, "_")[[1]][4]
    date <- substr(date, 1, 8)  # Keep only the first 8 characters (YYYYMMDD)
    
    # Construct the file path to the desired .envi file
    lai_file <- file.path(subdir, "res_10m", "PRO4SAIL_INVERSION_atbd", 
                          "addmult", "bands_3_4_8", 
                          paste0(dir_name, "_Refl_lai.envi"))
    evi_file <- file.path(subdir, "SpectralIndices", "res_10_m",
                          paste0(dir_name, "_EVI"))
    ndvi_file <- file.path(subdir, "SpectralIndices", "res_10_m",
                           paste0(dir_name, "_NDVI"))
    savi_file <- file.path(subdir, "SpectralIndices", "res_10_m",
                           paste0(dir_name, "_SAVI"))
    
    # Read the files as SpatRaster
    lai_raster <- rast(lai_file)
    evi_raster <- rast(evi_file)
    ndvi_raster <- rast(ndvi_file)
    savi_raster <- rast(savi_file)
    
    # Append them to the list of rasters
    lai_rasters_list[[date]] <- lai_raster
    evi_rasters_list[[date]] <- evi_raster
    ndvi_rasters_list[[date]] <- ndvi_raster
    savi_rasters_list[[date]] <- savi_raster
    
    # Sort the rasters list by date
    date_keys <- names(lai_rasters_list)
    dates <- as.Date(date_keys, format = "%Y%m%d")
    ordered_indices <- order(dates)
    lai_rasters_list <- lai_rasters_list[ordered_indices]
    evi_rasters_list <- evi_rasters_list[ordered_indices]
    ndvi_rasters_list <- ndvi_rasters_list[ordered_indices]
    savi_rasters_list <- savi_rasters_list[ordered_indices]
    
    # Save the stacked rasters to new files
    lai_output_name <- "s2lai_stacked_raster.tif"
    evi_output_name <- "evi_stacked_raster.tif"
    ndvi_output_name <- "ndvi_stacked_raster.tif"
    savi_output_name <- "savi_stacked_raster.tif"
    
    # Stack and save each type of raster
    stack_and_save(lai_rasters_list, results_site, lai_output_name, masks_dir, metrics_dir)
    stack_and_save(evi_rasters_list, results_site, evi_output_name, masks_dir, metrics_dir)
    stack_and_save(ndvi_rasters_list, results_site, ndvi_output_name, masks_dir, metrics_dir)
    stack_and_save(savi_rasters_list, results_site, savi_output_name, masks_dir, metrics_dir)
  }
}
  