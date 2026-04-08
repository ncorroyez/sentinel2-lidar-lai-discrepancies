# ---
# title: "main_explore_init.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-30"
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
library("terra")
library("ggplot2")
library("dplyr")
# library("gstat")
library("sp")
library("data.table")
library("spdep")
library("cluster")
library("drf")
library("bigmemory")
library("ff")
library("ranger")
library("randomForest")
library("bench")
library("Matrix")
source("libraries/functions_general_tools.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/explore_init"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)

vars <- c("lcv", "mean", "max", "vci", "cv_lad", "cv_lad_dtm", "fCover", "lskew", "lidarlai", "s2lai", "cv", "deltaLAI", "slope")
vars <- c("dsm_sd", "dsm_cv")

IQR_outliers <- function(DistVal,weightIRQ = 1.5){
  range_IQR <- c(stats::quantile(DistVal, 0.25,na.rm=TRUE),stats::quantile(DistVal, 0.75,na.rm=TRUE))
  iqr <- diff(range_IQR)
  outlier_IQR <- c(range_IQR[1]-weightIRQ*iqr,range_IQR[2]+weightIRQ*iqr)
  return(outlier_IQR)
}

IQR_inliers <- function(DistVal, weightIRQ = 1.5) {
  range_IQR <- c(stats::quantile(DistVal, 0.25, na.rm = TRUE), stats::quantile(DistVal, 0.75, na.rm = TRUE))
  iqr <- diff(range_IQR)
  outlier_IQR <- c(range_IQR[1] - weightIRQ * iqr, range_IQR[2] + weightIRQ * iqr)
  
  # Filter values within the IQR range
  inliers <- DistVal[DistVal >= outlier_IQR[1] & DistVal <= outlier_IQR[2]]
  
  return(inliers)
}

# Initialize a list to store results
extracted_values <- list()

# Loop through each variable and site
for (var in vars) {
  for (site in sites) {
    # Construct full path to the raster file
    full_path <- file.path(results_dir, site, metrics_dir)
    raster_file <- file.path(full_path, paste0(var, "_res_10_m.tif"))
    
    # Special handling for deltaLAI
    if (var == "deltaLAI") {
      # Compute deltaLAI dynamically
      lidar_file <- file.path(full_path, "lidarlai_res_10_m.tif")
      s2lai_file <- file.path(full_path, "s2lai_res_10_m.tif")
      
      if (file.exists(lidar_file) && file.exists(s2lai_file)) {
        # Load rasters and compute deltaLAI
        lidar_raster <- terra::rast(lidar_file)
        s2lai_raster <- terra::rast(s2lai_file)
        raster_data <- lidar_raster - s2lai_raster
      } else {
        message(paste("Missing lidar or s2lai raster for deltaLAI at site:", site))
        next
      }
    } else {
      # Load the standard variable raster
      if (file.exists(raster_file)) {
        raster_data <- terra::rast(raster_file)
      } else {
        message(paste("File does not exist:", raster_file))
        next
      }
    }
    
    # Extract raster values and filter inliers
    values <- values(raster_data, na.rm = TRUE)
    values <- as.data.frame(IQR_inliers(values))
    names(values) <- "Values"
    
    # Check if there are any valid values
    if (length(values$Values) == 0) {
      message(paste("No valid values found in raster for:", var, "at site:", site))
      next
    }
    
    # Calculate binwidth based on the range of values
    binwidth <- diff(range(values$Values)) / 30
    
    # Create a data frame for this site and variable
    df <- data.frame(
      Site = site,
      Variable = var,
      Values = values$Values,
      Binwidth = binwidth
    )
    
    # Store the data frame in the list
    extracted_values[[paste(site, var, sep = "_")]] <- df
    
    # Plot and save the individual histogram
    histogram_plot <- ggplot(df, aes(x = Values)) +
      geom_histogram(binwidth = binwidth, fill = "lightblue", color = "black") +
      labs(title = paste("Histogram of", var, "for", site),
           x = paste(var, "Values"),
           y = "Frequency") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = file.path(output_dir, paste0(var, "_histogram_", site, ".png")), 
           plot = histogram_plot, width = 8, height = 6)
  }
}

# Combine all data into a single data frame for easier plotting
combined_df <- do.call(rbind, extracted_values)

# Loop through each variable to create a combined plot with 3 panels (one per site)
for (var in unique(combined_df$Variable)) {
  var_df <- combined_df[combined_df$Variable == var, ]
  
  # Plot combined histogram for the variable with individual binwidths per site
  var_plot <- ggplot(var_df, aes(x = Values)) +
    geom_histogram(data = subset(var_df, Site == "Aigoual"), 
                   binwidth = unique(var_df$Binwidth[var_df$Site == "Aigoual"]), 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_histogram(data = subset(var_df, Site == "Blois"), 
                   binwidth = unique(var_df$Binwidth[var_df$Site == "Blois"]), 
                   fill = "skyblue", color = "black", alpha = 0.7) +
    geom_histogram(data = subset(var_df, Site == "Mormal"), 
                   binwidth = unique(var_df$Binwidth[var_df$Site == "Mormal"]), 
                   fill = "steelblue", color = "black", alpha = 0.7) +
    facet_wrap(~ Site, ncol = 3, scales = "free") +  # Free scales allow xlim and ylim to adapt to each site
    labs(title = paste("Histogram of", var, "across Sites"),
         x = paste(var, "Values"),
         y = "Frequency") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the combined plot for the variable to a file
  ggsave(filename = file.path(output_dir, paste0(var, "_facet_histogram.png")), plot = var_plot, width = 12, height = 6)
}

# Loop through each variable to create an overlaid normalized histogram
for (var in unique(combined_df$Variable)) {
  var_df <- combined_df[combined_df$Variable == var, ]
  
  # Create a density plot to compare distributions across sites
  overlay_plot <- ggplot(var_df, aes(x = Values, color = Site, fill = Site)) +
    geom_density(alpha = 0.7, adjust = 1) +  # Use density for normalized comparisons
    labs(title = paste("Normalized Histogram (Density) of", var, "across Sites"),
         x = paste(var, "Values"),
         y = "Density") +
    scale_color_manual(values = c("Aigoual" = "lightblue", "Blois" = "skyblue", "Mormal" = "steelblue")) +
    scale_fill_manual(values = c("Aigoual" = "lightblue", "Blois" = "skyblue", "Mormal" = "steelblue")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the overlaid density plot for the variable to a file
  ggsave(filename = file.path(output_dir, paste0(var, "_overlaid_density.png")), 
         plot = overlay_plot, width = 8, height = 6)
}







# Set up a 1-row, 3-column layout for the histograms
# par(mfrow = c(1, 3))  # Adjust rows/columns as needed

# Loop through each site and plot the histogram
# for (var in vars){
#   for (site in sites) {
#     full_path <- file.path(results_dir, site, metrics_dir)
#     raster_file <- file.path(full_path, paste0(var, "_res_10_m.tif"))
#     raster_data <- terra::rast(raster_file)
# 
#     hist(raster_data[],
#          main = paste("Histogram of", var, "for", site),
#          xlab = paste(var, "Values"),
#          ylab = "Frequency",
#          col = "lightblue",
#          border = "black")
#   }
# }

# Reset plotting layout
# par(mfrow = c(1, 1))