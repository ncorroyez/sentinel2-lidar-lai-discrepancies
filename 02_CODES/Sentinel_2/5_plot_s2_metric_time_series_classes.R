# ---
# title: "5_plot_s2_metric_time_series_classes.R"
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
library("ggplot2")
library("dplyr")
library("zoo")
source("../libraries/functions_plot_time_series.R")
source("../libraries/functions_general_tools.R")

# Pre-processing Parameters
data_dir <- '/media/corroyez/My Passport/01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../../03_RESULTS"
output_dir <- "../../04_FIGURES/lai_time_series"

# Read S2 Acquisition Dates (1/month)
file_path <- file.path(data_dir, "S2_Dates.csv")
data <- read_csv(file_path)

# List of directories to keep
dirs_to_keep <- c(
  # "Full_Composition",
  # "Not_Masked",
  "Deciduous_Only"
  # "Deciduous_Flex",
  # "Coniferous_Only",
  # "Coniferous_Flex"
)

structural_metrics <- c(
  "lcv",
  "mean"
)

# Initialize empty data frames to store LAI metrics for plotting
lai_data <- data.frame(site = character(),
                       composition = character(),
                       quantile = integer(),
                       layer = numeric(),
                       lai_mean = numeric(),
                       lai_std = numeric(),
                       structural_metric = character()
)
lai_data_values <- data.frame(site = character(),
                              composition = character(),
                              quantile = integer(),
                              layer = numeric(),
                              lai_values = numeric(),
                              structural_metric = character()
)

# lai_data <- prepare_lai_data(sites, results_dir, dirs_to_keep, structural_metrics)
# 
# plot_lai_facet_sites(lai_data, output_dir)
# 
# plot_lai_separate(lai_data, output_dir)

# Preparation
for (site in sites){
  
  cat("Processing site: ", site, "\n")
  
  # Open Metrics Dir
  masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks/Quantiles")
  metrics_dir <- file.path(results_dir, site, "Metrics")
  
  # Keep: Full_Composition, Not_Masked, Deciduous_Only, Deciduous_Flex, Coniferous_Only, Coniferous_Flex
  composition_dirs <- list.dirs(metrics_dir, recursive = FALSE)
  composition_dirs <- composition_dirs[grepl(paste(dirs_to_keep, collapse = "|"), composition_dirs)]
  
  for (composition_dir in composition_dirs){
    
    composition <- basename(composition_dir)
    cat("Processing composition: ", composition, "\n")
    # a
    # Open Stack Raster for a given composition
    lai_rasters <- terra::rast(file.path(composition_dir, 
                                         "stack",
                                         paste0(site, "_lai_stacked_raster.tif")))
    
    # Extract layer names (dates)
    layer_names <- names(lai_rasters)
    
    for (structural_metric in structural_metrics){
      
      cat("Processing metric: ", structural_metric, "\n")
      
      # Open Quantiles Masks
      quantiles_dir <- file.path(masks_dir, composition, structural_metric,
                                 "Deciles", "3_Quantiles")
      quantiles_low <- c(0, 0.33, 0.67)
      quantiles_high <- c(0.33, 0.66, 1)
      
      for (quantile_idx in seq_along(quantiles_low)) {  # use the index here
        
        cat("Processing quantile: ", quantile_idx, "\n")
        quantile_low <- quantiles_low[quantile_idx]  # Get the quantile increment (0, 1/3, 2/3)
        quantile_high <- quantiles_high[quantile_idx]  # Get the quantile increment (0, 1/3, 2/3)
        
        # Mask Stack Raster for a given height
        qt <- terra::rast(file.path(quantiles_dir,
                                    paste0(structural_metric, 
                                           "_heter_binary_mask_res_10_m_",
                                           quantile_low,
                                           "_",
                                           quantile_high,
                                           ".tif")))
        
        # Get Layer (Date)
        for (layer_idx in 1:nlyr(lai_rasters)){
          
          # cat("Processing layer: ", layer_idx, "\n")
          lai_raster <- lai_rasters[[layer_idx]]
          
          # Mask the LAI raster using the quantile mask
          lai_masked <- mask(lai_raster, qt)
          
          # Define the window size (e.g., 3x3 pixels)
          window_size <- 3
          
          # Apply focal to calculate the mean
          lai_mean_raster <- terra::focal(lai_masked, 
                                          w = matrix(1,
                                                     nrow = window_size,
                                                     ncol = window_size),
                                          fun = mean, 
                                          na.policy = "omit",
                                          na.rm = TRUE)
          
          # Apply focal to calculate the standard deviation
          lai_sd_raster <- terra::focal(lai_masked, 
                                        w = matrix(1,
                                                   nrow = window_size,
                                                   ncol = window_size),
                                        fun = sd, 
                                        na.policy = "omit",
                                        na.rm = TRUE)
          
          # Get all LAI values for current layer
          lai_values <- as.vector(values(lai_masked, na.rm = TRUE))
          names(lai_values) <- "lai_values"
          lai_mean_values <- values(lai_mean_raster, na.rm = TRUE)
          lai_sd_values <- values(lai_sd_raster, na.rm = TRUE)
          lai_mean <- mean(lai_mean_values, na.rm = TRUE)
          lai_std <- mean(lai_sd_values, na.rm = TRUE)
          
          # Get all LAI values
          # lai_values <- values(lai_masked, na.rm = T)
          
          # Calculate mean and standard deviation of LAI
          # lai_mean <- mean(lai_values, na.rm = TRUE)
          # lai_std <- sd(lai_values, na.rm = TRUE)
          
          # Store the data in a data frame
          lai_data <- rbind(lai_data, 
                            data.frame(site = site,
                                       composition = composition,
                                       quantile = quantile_idx,
                                       layer = layer_names[layer_idx],
                                       lai_mean = lai_mean,
                                       lai_std = lai_std,
                                       structural_metric = structural_metric)
          )
          lai_data_values <- rbind(lai_data_values, 
                                   data.frame(site = site,
                                              composition = composition,
                                              quantile = quantile_idx,
                                              layer = layer_names[layer_idx],
                                              lai_values = lai_values,
                                              structural_metric = structural_metric)
          )
        }
      }
    }
  }
}

# Convert `layer` from YYYYMMDD to YYYY-MM (only year and month)
# Use `substr` to extract the year and month, then convert to yearmon format
lai_data$layer <- as.yearmon(substr(lai_data$layer, 1, 6), "%Y%m")
lai_data_values$layer <- as.yearmon(substr(lai_data_values$layer, 1, 6), "%Y%m")


# Loop through each unique composition and structural metric for plotting
for (composition_name in unique(lai_data_values$composition)) {
  for (metric_name in unique(lai_data_values$structural_metric)) {
    
    # Filter data for the current composition and structural metric
    composition_metric_data <- lai_data_values[lai_data_values$composition == composition_name & 
                                                 lai_data_values$structural_metric == metric_name, ]
    
    # Check if there's data to plot
    if (nrow(composition_metric_data) > 0) {
      # Create the plot using all raw values for the boxplot
      p <- ggplot(composition_metric_data, aes(x = as.factor(layer), y = lai_values, fill = factor(quantile))) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "Month", y = "LAI Values", fill = "Quantile Class",
             title = paste("LAI for", composition_name, "using", metric_name)) +
        facet_wrap(~ site) +  # Separate panels for each site
        theme_bw() +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14))  # Center the title
      
      # Print the plot to view in RStudio
      print(p)
      
      # Save the plot as an image file (e.g., PNG)
      ggsave(filename = paste0(output_dir, "LAI_plot_", composition_name, "_", metric_name, ".png"), 
             plot = p, width = 10, height = 8)
    }
  }
}

# Loop through each unique site, composition, and structural metric for plotting
for (site_name in unique(lai_data_values$site)) {
  for (composition_name in unique(lai_data_values$composition)) {
    for (metric_name in unique(lai_data_values$structural_metric)) {
      
      # Filter data for the current site, composition, and structural metric
      site_composition_metric_data <- lai_data_values[lai_data_values$site == site_name &
                                                        lai_data_values$composition == composition_name & 
                                                        lai_data_values$structural_metric == metric_name, ]
      
      # Check if there's data to plot
      if (nrow(site_composition_metric_data) > 0) {
        # Create the plot for the current site, composition, and structural metric
        p <- ggplot(site_composition_metric_data, aes(x = as.factor(layer), y = lai_values, fill = factor(quantile))) +
          geom_boxplot(outlier.shape = NA) +
          labs(x = "Month", y = "Mean LAI", color = "Quantile Class",
               title = paste("Mean LAI with Standard Deviation for", site_name, 
                             "\nComposition:", composition_name, "| Metric:", metric_name)) +
          # scale_x_yearmon(n.breaks = 20, format = "%b %Y") +  # Use scale_x_yearmon for year-month format
          theme_bw() +
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14))  # Center the title
        
        # Print the plot to view in RStudio
        print(p)
        
        # Save the plot as an image file (e.g., PNG)
        ggsave(filename = paste0(output_dir, "LAI_plot_", site_name, "_", composition_name, "_", metric_name, ".png"), 
               plot = p, width = 10, height = 8)
      }
    }
  }
}







# Loop through each unique composition and structural metric for plotting
for (composition_name in unique(lai_data$composition)) {
  for (metric_name in unique(lai_data$structural_metric)) {
    
    # Filter data for the current composition and structural metric
    composition_metric_data <- lai_data[lai_data$composition == composition_name & 
                                          lai_data$structural_metric == metric_name, ]
    
    # Check if there's data to plot
    if (nrow(composition_metric_data) > 0) {
      # Create the plot for the current composition and structural metric, with facets for each site
      p <- ggplot(composition_metric_data, aes(x = layer, y = lai_mean, color = factor(quantile))) +
        # geom_boxplot() +
        geom_line() +  # Line plot for LAI mean
        geom_point() +  # Add points for the mean
        geom_errorbar(aes(ymin = lai_mean - lai_std, ymax = lai_mean + lai_std), width = 0.2) +  # Error bars for SD
        labs(x = "Month", y = "Mean LAI", color = "Quantile Class",
             title = paste("Mean LAI with Standard Deviation for", composition_name, "using", metric_name)) +
        facet_wrap(~ site) +  # Separate panels for each site
        theme_bw() +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14))  # Center the title
      
      # Print the plot to view in RStudio
      print(p)
      
      # Save the plot as an image file (e.g., PNG)
      ggsave(filename = paste0(output_dir, "LAI_plot_", composition_name, "_", metric_name, ".png"), plot = p, width = 10, height = 8)
    }
  }
}

# Loop through each unique site, composition, and structural metric for plotting
for (site_name in unique(lai_data$site)) {
  for (composition_name in unique(lai_data$composition)) {
    for (metric_name in unique(lai_data$structural_metric)) {
      
      # Filter data for the current site, composition, and structural metric
      site_composition_metric_data <- lai_data[lai_data$site == site_name &
                                                 lai_data$composition == composition_name & 
                                                 lai_data$structural_metric == metric_name, ]
      
      # Check if there's data to plot
      if (nrow(site_composition_metric_data) > 0) {
        # Create the plot for the current site, composition, and structural metric
        p <- ggplot(site_composition_metric_data, aes(x = layer, y = lai_mean, color = factor(quantile))) +#group = interaction(layer, quantile), fill = factor(quantile))) +
          # geom_boxplot() +
          geom_line() +  # Line plot for LAI mean
          geom_point() +  # Add points for the mean
          geom_errorbar(aes(ymin = lai_mean - lai_std, ymax = lai_mean + lai_std), width = 0.2) +  # Error bars for SD
          labs(x = "Month", y = "Mean LAI", color = "Quantile Class",
               title = paste("Mean LAI with Standard Deviation for", site_name, 
                             "\nComposition:", composition_name, "| Metric:", metric_name)) +
          scale_x_yearmon(n.breaks = 20, format = "%b %Y") +  # Use scale_x_yearmon for year-month format
          theme_bw() +
          theme(legend.position = "top",
                axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14))  # Center the title
        
        # Print the plot to view in RStudio
        print(p)
        
        # Save the plot as an image file (e.g., PNG)
        ggsave(filename = paste0(output_dir, "LAI_plot_", site_name, "_", composition_name, "_", metric_name, ".png"), 
               plot = p, width = 10, height = 8)
      }
    }
  }
}