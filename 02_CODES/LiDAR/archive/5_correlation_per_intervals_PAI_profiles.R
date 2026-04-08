# ---
# title: "3.correlation_per_intervals_PAI_profiles.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-19"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("plotly")
library("terra")
library("viridis")
library("future")
library("RColorBrewer")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../functions_plots.R")
source("functions_lidar.R")

# Pre-processing Parameters & Directories
results_dir <- '../../03_RESULTS'
site <- "Mormal"
results_path <- file.path("../../03_RESULTS", site)
lidar_dir <- "LiDAR"
masks_dir <- file.path(results_path, lidar_dir, "Heterogeneity_Masks")
deciles_dir <- file.path(masks_dir, "Deciles") 
equal_intervals_dir <- file.path(masks_dir, "Equal_Intervals") 
profiles_dir <- file.path(results_path,lidar_dir, "Profiles")

# S-2
mask_10m <- terra::rast(file.path(masks_dir, "artifacts_low_vegetation_majority_90_p_res_10_m.envi"))
lai_sentinel <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/res_10m/PRO4SAIL_INVERSION_atbd/addmult/bands_3_4_8/L2A_T31UER_A031222_20210614T105443_Refl_lai.envi")
lai_sentinel_masked <- lai_sentinel * mask_10m

# Function to plot correlation curves
plot_correlation_curves <- function(correlation_values, file_name) {
  
  # Plotting
  png(file_name, width = 3000, height = 2000, units = "px", res = 200)
  
  # Initialize the plot
  plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.75), 
       xlab = "PAI Profile Depth", ylab = "Correlation Value",
       main = "Correlation Values of Intervals PAI Profiles of Varying Depth with Sentinel-2 LAI for an Error")
  
  # Define line styles
  # line_styles <- c("solid", "dashed", "dotted", "dotdash", "longdash",
  #                  "twodash", "dashdot", "dashdotdot", "longdashdot", "longdashdotdot")
  
  # Define line styles
  # line_styles <- c("solid", "dashed", "dotted", "dotdash")
  # 
  # # Repeat line styles to ensure length is a multiple of 4
  # while (length(line_styles) < length(correlation_values)) {
  #   line_styles <- c(line_styles, line_styles)
  # }
  
  # Plot all correlation curves
  # colors <- rainbow(length(correlation_values))
  colors <- brewer.pal(length(correlation_values), "Paired")
  for (i in seq_along(correlation_values)) {
    key <- names(correlation_values)[i]
    values <- correlation_values[[key]]
    values <- rev(values)
    # Check for NAs in the data
    if (any(is.na(values))) {
      next  # Skip this interval if there are NAs in the data
    }
    # lines(seq(0, 44, by = 1), values, 
    #       type = "l", col = colors[i],
    #       lty = line_styles[i %% length(line_styles) + 1],  # Cycling through line styles
    #       lwd = 2)
    lines(seq(0, 44, by = 1), values, 
          type = "l", col = colors[i], lwd = 2)
  }
  
  # Add legend
  legend("bottomright", legend = names(correlation_values), col = colors, 
         # lty = 1, 
         lwd = 2, cex = 0.8)  # Ensure consistent line type in legend
  
  dev.off()
}

# Function to load correlation values based on error type and interval
load_correlation_values <- function(error_type, interval_dir, interval_key) {
  correlation_values <- list()
  for (i in seq(46.5, 2.5, by = -1)) {
    error <- terra::rast(file.path(interval_dir, 
                                   paste0(sprintf("%s_heter_binary_mask_res_10_m_%s.tif", error_type, interval_key))))
    # Load lidar file
    lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
    lidar <- terra::rast(lidar_file)
    lidar <- terra::project(lidar, mask_10m)
    lidar_masked <- lidar * mask_10m
    lidar_masked <- mask(lidar_masked, error)
    lidar_values <- values(lidar_masked)
    
    # Load S-2 data
    lai_sentinel_masked_error <- mask(lai_sentinel_masked, error)
    sentinel_values <- values(lai_sentinel_masked_error)
    
    # Calculate correlation
    correlation_value <- cor(lidar_values, 
                             sentinel_values,
                             use = "pairwise.complete.obs")
    
    # Store correlation value
    correlation_values[[i]] <- correlation_value
  }
  correlation_values <- correlation_values[!sapply(correlation_values, is.null)]
  return(correlation_values)
}

intervals_10 <- seq(0, 0.9, by = 0.1)
for (error_type in c("mean_heights", "cv", "std")) {
  print(error_type)
  correlation_values <- list()
  for (interval in intervals_10) {
    print(interval)
    key_start <- ifelse(interval == 0, 0, format(interval, nsmall = 1))
    key_end <- ifelse(interval + 0.1 == 1.0, 1, format(interval + 0.1, nsmall = 1))
    key <- paste0(key_start, "_", key_end)
    correlation_values[[as.character(paste0(interval, "_", interval + 0.1))]] <- load_correlation_values(error_type, deciles_dir, key)
  }
  
  # Plot correlation curves for the current error type
  file_name <- file.path(profiles_dir, paste0("correlation_lidar_s2_", error_type, "_10_deciles.png"))
  plot_correlation_curves(correlation_values, file_name)
}

# Loop over intervals for 5 intervals (0-0.2, 0.2-0.4, ..., 0.8-1)
intervals_5 <- seq(0, 0.8, by = 0.2)
for (error_type in c("mean_heights", "cv", "std")) {
  print(error_type)
  correlation_values <- list()
  for (interval in intervals_5) {
    print(interval)
    key_start <- ifelse(interval == 0, 0, format(interval, nsmall = 1))
    key_end <- ifelse(interval + 0.2 == 1.0, 1, format(interval + 0.2, nsmall = 1))
    key <- paste0(key_start, "_", key_end)
    correlation_values[[as.character(paste0(interval, "_", interval + 0.2))]] <- load_correlation_values(error_type, equal_intervals_dir, key)
  }
  
  # Plot correlation curves for the current error type
  file_name <- file.path(profiles_dir, paste0("correlation_lidar_s2_", error_type, "_5_equal_intervals.png"))
  plot_correlation_curves(correlation_values, file_name)
}



# Parametrize
correlation_values <- list()
# Loop through each raster file
for (k in seq(0, 0.9, by = 0.1)){
  error <- terra::rast(file.path(deciles_dir, 
                                 paste0(sprintf("std_heter_binary_mask_res_10_m_%s_%s.tif",
                                                k, k + 0.1))))
  correlations <- c()  # Initialize correlations vector for this interval
  for (i in seq(46.5, 2.5, by = -1)) {
    # Generate lidar file name
    lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
    lidar <- terra::rast(lidar_file)
    lidar <- terra::project(lidar, mask_10m)
    lidar_masked <- lidar * mask_10m
    lidar_masked <- mask(lidar_masked, error)
    lidar_values <- values(lidar_masked)
    
    # S-2
    lai_sentinel_masked_error <- mask(lai_sentinel_masked, error)
    sentinel_values <- values(lai_sentinel_masked_error)
    
    # Calculate correlation and store
    correlation_value <- correlation_test_function(lidar_values, sentinel_values)
    correlations <- c(correlations, correlation_value)
  }
  correlation_values[[paste0(k, "-", k + 0.1)]] <- correlations 
}

# Plotting
png(file.path(profiles_dir, "correlation_lidar_s2_for_mean_heights_deciles_profiles.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

# Initialize the plot
plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.65), xlab = "PAI Profile Depth", ylab = "Correlation Value",
     main = "Correlation Values of PAI Profiles of Varying Depth with Sentinel-2 LAI")

line_styles <- c("solid", "dashed", "dotted", "dotdash", "longdash")

# Plot all correlation curves
colors <- rainbow(length(correlation_values))
for (i in seq_along(correlation_values)) {
  key <- names(correlation_values)[i]
  values <- correlation_values[[i]]
  # Check for NAs in the data
  if (any(is.na(values))) {
    next  # Skip this interval if there are NAs in the data
  }
  
  lines(seq(0, 44, by = 1),
        values, 
        type = "l", 
        col = colors[i],
        lty = line_styles[i %% length(line_styles) + 1],
        lwd = 2)
}

# Add legend
legend("bottomright", legend = names(correlation_values), col = colors, lty = 1, lwd = 2, cex = 0.8)
dev.off()

