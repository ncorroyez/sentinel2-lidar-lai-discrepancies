# ---
# title: "functions_heterogeneity_depth_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-06-12"
# ---

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("sf")
library("dplyr")

#' Perform Heterogeneity Analysis
#'
#' This function performs a heterogeneity analysis by calculating Pearson correlation values 
#' between LiDAR and Sentinel-2 metrics across specified quantile ranges and generates a plot 
#' of these correlation values.
#'
#' @param site Character. The name of the site being analyzed.
#' @param composition_mask Character. The composition mask applied to the data.
#' @param metric Character. The metric used for the analysis (e.g., "lai", "pai").
#' @param quantiles_dir Character. The directory containing the quantile raster files.
#' @param n_interval Integer. The number of intervals for the quantile ranges.
#' @param inc Numeric. The increment value for the quantile ranges.
#' @param figures_dir Character. The directory where the output plot will be saved.
#'
#' @return A plot saved in the specified figures directory showing the correlation values 
#'         for different quantile ranges.
#'
#' @importFrom terra rast values global
#' @importFrom stats cor.test
#' @export
perform_heterogeneity_analysis <- function(site,
                                           composition_mask,
                                           metric,
                                           quantiles_dir,
                                           n_interval,
                                           inc,
                                           figures_dir){
  
  correlation_values <- quantiles_range <- c()
  
  # Quantiles range will have a size of n_intervals + 1
  # -> Add the + 1 value (can be min of 1st quantile or max of last quantile)
  # Here we choose min of 1st quantile
  quantiles_range <- c(quantiles_range,
                       round(global(terra::rast(file.path(quantiles_dir,
                                                          paste0(metric,
                                                                 "_heter_raster_res_10_m_",
                                                                 0,
                                                                 "_",
                                                                 round(inc, digits = 2),
                                                                 ".tif"))),
                                    "min", na.rm = TRUE), 2)[["min"]])
  if (inc == 0.2) {
    max_val <- 0.8
  } else if (inc == 0.1) {
    max_val <- 0.9
  } else if (inc == 1/3) {
    max_val <- 2/3
  } else {
    stop("inc is not 0.1 or 0.2 or 1/3")
  }
  for (low_value in seq(0, max_val, by = inc)) {
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    lidar_intervals <- terra::rast(file.path(quantiles_dir, 
                                             paste0(metric,
                                                    "_pai_masked_res_10_m_", 
                                                    low_value, 
                                                    "_", 
                                                    high_value,
                                                    ".tif")))
    lidar_values <- values(lidar_intervals)
    
    s2_intervals <- terra::rast(file.path(quantiles_dir, 
                                          paste0(metric,
                                                 "_lai_s2_masked_res_10_m_", 
                                                 low_value, 
                                                 "_", 
                                                 high_value,
                                                 ".tif")))
    s2_values <- values(s2_intervals)
    
    correlation_value <- cor.test(lidar_values, 
                                  s2_values, 
                                  method = "pearson")$estimate
    correlation_values <- c(correlation_values, correlation_value)
    
    quantiles_range <- c(quantiles_range, 
                         round(global(terra::rast(file.path(quantiles_dir, 
                                                            paste0(metric,
                                                                   "_heter_raster_res_10_m_", 
                                                                   low_value, 
                                                                   "_", 
                                                                   high_value,
                                                                   ".tif"))),
                                      "max", na.rm = TRUE), 2)[["max"]])
    
  }
  # print(file.path(figures_dir, 
  #                 paste0(site,
  #                        "_",
  #                        composition_mask,
  #                        "_",
  #                        metric,
  #                        "_correlation_values_for_",
  #                        n_interval,
  #                        "_quantiles_res_10_m.png")))
  # a
  png(filename = file.path(figures_dir, 
                           paste0(site,
                                  "_",
                                  composition_mask,
                                  "_",
                                  metric,
                                  "_correlation_values_for_",
                                  n_interval,
                                  "_quantiles_res_10_m.png")),
      width = 1920, height = 1080)  # Adjust width and height as needed
  cor_plot <- plot(1, type = "n", 
                   xlab = "Quantile", 
                   ylab = "Correlation Value", 
                   xlim = c(0, 1),
                   ylim = c(0, 0.8),
                   main = paste("Correlation Values for", 
                                metric, 
                                site, 
                                composition_mask,
                                n_interval,
                                "quantiles"),
                   xaxt = "n",
                   cex.main = 1.2,   # Increase title font size
                   cex.lab = 1.2,    # Normal font size for axis labels
                   cex.axis = 1.2)
  midpoints <- seq(inc/2, 1-inc/2, by = inc)
  lines(midpoints, correlation_values, type = "l")
  points(midpoints, correlation_values, pch = 16, col = "red")
  axis(side = 1, at = seq(0, 1, by = inc), labels = quantiles_range, cex.axis = 1.2)
  dev.off()
}


perform_heterogeneity_analysis_2 <- function(site,
                                             composition_mask,
                                             metric,
                                             quantiles_dir,
                                             n_interval,
                                             inc,
                                             figures_dir){
  
  correlation_values <- rmse_values <- r2_values <- slope_values <- quantiles_range <- c()
  
  # Quantiles range will have a size of n_intervals + 1
  # -> Add the + 1 value (can be min of 1st quantile or max of last quantile)
  # Here we choose min of 1st quantile
  quantiles_range <- c(quantiles_range,
                       round(global(terra::rast(file.path(quantiles_dir,
                                                          paste0(metric,
                                                                 "_heter_raster_res_10_m_",
                                                                 0,
                                                                 "_",
                                                                 round(inc, digits = 2),
                                                                 ".tif"))),
                                    "min", na.rm = TRUE), 2)[["min"]])
  
  # Set max_val based on the increment value
  if (inc == 0.2) {
    max_val <- 0.8
  } else if (inc == 0.1) {
    max_val <- 0.9
  } else if (inc == 1/3) {
    max_val <- 2/3
  } else {
    stop("inc is not 0.1 or 0.2 or 1/3")
  }
  
  for (low_value in seq(0, max_val, by = inc)) {
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    
    # Load lidar and Sentinel-2 (S2) data
    lidar_intervals <- terra::rast(file.path(quantiles_dir, 
                                             paste0(metric,
                                                    "_pai_masked_res_10_m_", 
                                                    low_value, 
                                                    "_", 
                                                    high_value,
                                                    ".tif")))
    lidar_values <- values(lidar_intervals)
    
    s2_intervals <- terra::rast(file.path(quantiles_dir, 
                                          paste0(metric,
                                                 "_lai_s2_masked_res_10_m_", 
                                                 low_value, 
                                                 "_", 
                                                 high_value,
                                                 ".tif")))
    s2_values <- values(s2_intervals)
    
    # Filter out NA values for both lidar and S2
    valid_idx <- !is.na(lidar_values) & !is.na(s2_values)
    lidar_values <- lidar_values[valid_idx]
    s2_values <- s2_values[valid_idx]
    
    # Correlation
    correlation_value <- cor.test(lidar_values, 
                                  s2_values, 
                                  method = "pearson")$estimate
    correlation_values <- c(correlation_values, correlation_value)
    
    # RMSE
    rmse_value <- sqrt(mean((lidar_values - s2_values)^2))
    rmse_values <- c(rmse_values, rmse_value)
    
    # Regression (to compute R² and slope)
    lm_model <- lm(lidar_values ~ s2_values)
    r2_value <- summary(lm_model)$r.squared
    slope_value <- coef(lm_model)[2]  # The slope is the coefficient of lidar_values
    
    r2_values <- c(r2_values, r2_value)
    slope_values <- c(slope_values, slope_value)
    
    # Quantiles range update
    quantiles_range <- c(quantiles_range, 
                         round(global(terra::rast(file.path(quantiles_dir, 
                                                            paste0(metric,
                                                                   "_heter_raster_res_10_m_", 
                                                                   low_value, 
                                                                   "_", 
                                                                   high_value,
                                                                   ".tif"))),
                                      "max", na.rm = TRUE), 2)[["max"]])
    
  }
  
  # Save the plots (adjust size and labels as needed)
  png(filename = file.path(figures_dir, 
                           paste0(site,
                                  "_",
                                  composition_mask,
                                  "_",
                                  metric,
                                  "_analysis_for_",
                                  n_interval,
                                  "_quantiles_res_10_m.png")),
      width = 1920, height = 1080)
  
  # Set up for plotting multiple metrics in one figure
  par(mfrow = c(2, 2))  # Set to 4 plots
  
  midpoints <- seq(inc/2, 1-inc/2, by = inc)
  
  # Correlation plot
  plot(1, type = "n", 
       xlab = "Quantile", 
       ylab = "Correlation Value", 
       xlim = c(0, 1),
       ylim = range(correlation_values, na.rm = TRUE),
       main = paste("Correlation Values for", metric),
       xaxt = "n")
  lines(midpoints, correlation_values, type = "l")
  points(midpoints, correlation_values, pch = 16, col = "red")
  axis(side = 1, at = seq(0, 1, by = inc), labels = quantiles_range)
  
  # RMSE plot
  plot(1, type = "n", 
       xlab = "Quantile", 
       ylab = "RMSE", 
       xlim = c(0, 1),
       ylim = range(rmse_values, na.rm = TRUE),
       main = paste("RMSE for", metric),
       xaxt = "n")
  lines(midpoints, rmse_values, type = "l")
  points(midpoints, rmse_values, pch = 16, col = "blue")
  axis(side = 1, at = seq(0, 1, by = inc), labels = quantiles_range)
  
  # R2 plot
  plot(1, type = "n", 
       xlab = "Quantile", 
       ylab = "R2", 
       xlim = c(0, 1),
       ylim = range(r2_values, na.rm = TRUE),
       main = paste("R2 for", metric),
       xaxt = "n")
  lines(midpoints, r2_values, type = "l")
  points(midpoints, r2_values, pch = 16, col = "green")
  axis(side = 1, at = seq(0, 1, by = inc), labels = quantiles_range)
  
  # Slope plot
  plot(1, type = "n", 
       xlab = "Quantile", 
       ylab = "Slope", 
       xlim = c(0, 1),
       ylim = range(slope_values, na.rm = TRUE),
       main = paste("Slope for", metric),
       xaxt = "n")
  lines(midpoints, slope_values, type = "l")
  points(midpoints, slope_values, pch = 16, col = "purple")
  axis(side = 1, at = seq(0, 1, by = inc), labels = quantiles_range)
  
  dev.off()
}




#' Perform Depth Analysis
#'
#' This function calculates the correlation between LiDAR LAI profiles of varying depth and Sentinel-2 LAI.
#'
#' @param site A character string representing the site name.
#' @param pad_profiles_dir A character string representing the directory containing PAD profile files.
#' @param figures_dir A character string representing the directory where the output figure will be saved.
#'
#' @return A PNG plot showing the correlation values between LiDAR LAI profiles of varying depth and Sentinel-2 LAI.
#' @export
perform_depth_analysis <- function(site,
                                   chm,
                                   pad_profiles_dir,
                                   figures_dir){
  correlation_values <- c()
  pad_files <- list.files(pad_profiles_dir, pattern = "\\.tif$", full.names = T)
  
  # Extract the numbers
  numbers <- as.numeric(sub(".*/PAD_(\\d+\\.\\d+)_.*", "\\1", pad_files))
  
  # Sort the numbers in descending order
  sorted_indices <- order(numbers, decreasing = TRUE)
  
  # Reorder the file list based on the sorted indices
  sorted_pad_files <- pad_files[sorted_indices]
  
  i <- 1
  for (pad_file in sorted_pad_files) {
    # Create a logical mask where CHM is less than i
    chm_mask <- chm < i
    
    s2_lai <- terra::rast(file.path(pad_profiles_dir, "../s2lai_res_10_m.tif"))
    s2_lai <- mask(s2_lai, chm_mask, maskvalues = TRUE)
    s2_lai_values <- values(s2_lai, na.rm = T)
    
    lidar_pad <- terra::rast(pad_file)
    lidar_pad <- mask(lidar_pad, chm_mask, maskvalues = TRUE)
    lidar_pad_values <- values(lidar_pad, na.rm = T)
    
    if (length(lidar_pad_values) < 100){
      correlation_value <- NA
    }
    else {
      correlation_value <- cor(lidar_pad_values,
                               s2_lai_values,
                               use = "pairwise.complete.obs")
    }
    i <- i + 1
    correlation_values <- c(correlation_values, correlation_value)
  }
  correlation_values <- Filter(Negate(is.null), correlation_values)
  
  png(file.path(figures_dir, 
                paste0("vegetation_depth_", 
                       site, 
                       "_", 
                       composition_mask, 
                       ".png")), 
      width = 1920, 
      height = 1080, 
      units = "px", res = 200)
  
  plot(seq(0, length(pad_files) - 1, by = 1), correlation_values, 
       type = "l", 
       xlim = c(0, length(pad_files) + 2), 
       ylim = c(0, 0.8), 
       xlab = "LiDAR LAI Profile Depth", 
       ylab = "Correlation Value",
       # main = "Correlation Values between LiDAR LAI Profiles of Varying Depth and Sentinel-2 LAI", 
       lwd = 2,
       xaxt = "n",
       cex.main = 1.2,   # Increase title font size
       cex.lab = 1.2,    # Normal font size for axis labels
       cex.axis = 1.2)
  title(main = paste("Pearson Correlation between LiDAR LAI Profiles\n",
                     "Integrated over Varying Depth and Sentinel-2 LAI"))
  axis(side = 1, 
       at = seq(0, length(pad_files) - 1, by = 2.5),
       labels = seq(0, length(pad_files) - 1, by = 2.5))
  points(seq(0, length(pad_files) - 1, by = 1), correlation_values, pch = 1)
  dev.off()
}

perform_depth_analysis_2 <- function(site,
                                     composition_mask,
                                     chm,
                                     pad_profiles_dir,
                                     figures_dir){
  correlation_values <- rmse_values <- r2_values <- slope_values <- c()
  pad_files <- list.files(pad_profiles_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # Extract the numbers (depth values) from the filenames
  numbers <- as.numeric(sub(".*/PAD_(\\d+\\.\\d+)_.*", "\\1", pad_files))
  
  # Sort the depth values in descending order
  sorted_indices <- order(numbers, decreasing = T)
  sorted_pad_files <- pad_files[sorted_indices]
  
  i <- 1
  for (pad_file in sorted_pad_files) {
    # Create a logical mask where CHM is less than the current depth (i)
    chm_mask <- chm < i
    
    # Load Sentinel-2 LAI raster
    s2_lai <- terra::rast(file.path(pad_profiles_dir, "../s2lai_summer_res_10_m.tif"))
    s2_lai <- mask(s2_lai, chm_mask, maskvalues = T)
    s2_lai_values <- values(s2_lai, na.rm = T)
    
    # Load LiDAR PAD raster for the current depth
    lidar_pad <- terra::rast(pad_file)
    lidar_pad <- mask(lidar_pad, chm_mask, maskvalues = T)
    lidar_pad_values <- values(lidar_pad, na.rm = T)
    
    # Filter out NA values for both LiDAR and S2
    # valid_idx <- !is.na(lidar_pad_values) & !is.na(s2_lai_values)
    # lidar_pad_values <- lidar_pad_values[valid_idx]
    # s2_lai_values <- s2_lai_values[valid_idx]
    
    # print(str(s2_lai_values))
    # print(str(lidar_pad_values))
    
    if (length(lidar_pad_values) < 100) {
      correlation_value <- NA
      rmse_value <- NA
      r2_value <- NA
      slope_value <- NA
    } else {
      # Pearson correlation
      correlation_value <- cor(lidar_pad_values,
                               s2_lai_values,
                               use = "pairwise.complete.obs")
      
      # RMSE calculation
      rmse_value <- sqrt(mean((lidar_pad_values - s2_lai_values) ^ 2))
      
      # Linear regression to get R² and slope
      lm_model <- lm(s2_lai_values ~ lidar_pad_values)
      r2_value <- summary(lm_model)$r.squared
      slope_value <- coef(lm_model)[2]  # The slope of the regression line
    }
    
    # Store the calculated values
    correlation_values <- c(correlation_values, correlation_value)
    rmse_values <- c(rmse_values, rmse_value)
    r2_values <- c(r2_values, r2_value)
    slope_values <- c(slope_values, slope_value)
    
    i <- i + 1
  }
  
  # Remove any NULL or NA values from the lists
  correlation_values <- Filter(Negate(is.null), correlation_values)
  rmse_values <- Filter(Negate(is.null), rmse_values)
  r2_values <- Filter(Negate(is.null), r2_values)
  slope_values <- Filter(Negate(is.null), slope_values)
  
  # Plot Correlation
  png(file.path(figures_dir, paste0("correlation_", site, composition_mask, "_depth_analysis.png")), 
      width = 1920, height = 1080, units = "px", res = 200)
  plot(seq(0, length(pad_files) - 1, by = 1), correlation_values, 
       type = "l", 
       xlab = "LiDAR LAI Profile Depth", 
       ylab = "Correlation Value", 
       ylim = c(0, 1), 
       lwd = 2, 
       xaxt = "n", 
       cex.lab = 1.2, 
       cex.axis = 1.2)
  title(main = "Pearson Correlation")
  axis(side = 1, at = seq(0, length(pad_files) - 1, by = 1), 
       labels = seq(0, length(pad_files) - 1, by = 1))
  points(seq(0, length(pad_files) - 1, by = 1), correlation_values, pch = 1)
  dev.off()
  
  # Plot RMSE
  png(file.path(figures_dir, paste0("rmse_", site, composition_mask, "_depth_analysis.png")), 
      width = 1920, height = 1080, units = "px", res = 200)
  plot(seq(0, length(pad_files) - 1, by = 1), rmse_values, 
       type = "l", 
       xlab = "LiDAR LAI Profile Depth", 
       ylab = "RMSE", 
       ylim = range(rmse_values, na.rm = TRUE), 
       lwd = 2, 
       xaxt = "n", 
       cex.lab = 1.2, 
       cex.axis = 1.2)
  title(main = "RMSE")
  axis(side = 1, at = seq(0, length(pad_files) - 1, by = 1), 
       labels = seq(0, length(pad_files) - 1, by = 1))
  points(seq(0, length(pad_files) - 1, by = 1), rmse_values, pch = 1)
  dev.off()
  
  # Plot R²
  png(file.path(figures_dir, paste0("r2_", site, composition_mask, "_depth_analysis.png")), 
      width = 1920, height = 1080, units = "px", res = 200)
  plot(seq(0, length(pad_files) - 1, by = 1), r2_values, 
       type = "l", 
       xlab = "LiDAR LAI Profile Depth", 
       ylab = "R²", 
       ylim = range(r2_values, na.rm = TRUE), 
       lwd = 2, 
       xaxt = "n", 
       cex.lab = 1.2, 
       cex.axis = 1.2)
  title(main = "R²")
  axis(side = 1, at = seq(0, length(pad_files) - 1, by = 1), 
       labels = seq(0, length(pad_files) - 1, by = 1))
  points(seq(0, length(pad_files) - 1, by = 1), r2_values, pch = 1)
  dev.off()
  
  # Plot Slope
  png(file.path(figures_dir, paste0("slope_", site, composition_mask, "_depth_analysis.png")), 
      width = 1920, height = 1080, units = "px", res = 200)
  plot(seq(0, length(pad_files) - 1, by = 1), slope_values, 
       type = "l", 
       xlab = "LiDAR LAI Profile Depth", 
       ylab = "Slope", 
       ylim = range(slope_values, na.rm = TRUE), 
       lwd = 2, 
       xaxt = "n", 
       cex.lab = 1.2, 
       cex.axis = 1.2)
  title(main = "Slope of Regression Line")
  axis(side = 1, at = seq(0, length(pad_files) - 1, by = 1), 
       labels = seq(0, length(pad_files) - 1, by = 1))
  points(seq(0, length(pad_files) - 1, by = 1), slope_values, pch = 1)
  dev.off()
}

perform_heterogeneity_depth_analysis_cvlad <- function(site,
                                                       composition_mask,
                                                       chm,
                                                       quantiles_dir,
                                                       n_interval,
                                                       inc,
                                                       figures_dir){
  correlation_values <- list()
  
  if (inc == 0.2) {
    max_val <- 0.8
  } else if (inc == 0.1) {
    max_val <- 0.9
  } else if (inc == 1/3) {
    max_val <- 2/3
  } else {
    stop("inc is not 0.1 or 0.2 or 1/3")
  }
  
  i_quantile <- 1
  for (low_value in seq(0, max_val, by = inc)) {
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    correlations <- numeric(37)
    i <- 1
    for (min_depth in seq(38.5, 2.5)) {
      # Create a logical mask where CHM is greater than i
      chm_mask <- chm < i + 1
      
      lidar_pad_masked <- terra::rast(file.path(quantiles_dir, 
                                                paste0("CV_LAD_",
                                                       min_depth,
                                                       "_40.tif_pai_masked_res_10_m_",
                                                       low_value,
                                                       "_",
                                                       high_value,
                                                       ".tif")))
      lidar_pad_masked <- mask(lidar_pad_masked, chm_mask, maskvalues = TRUE)
      lidar_pad_values <- values(lidar_pad_masked, na.rm = T)
      
      s2_lai_masked <- terra::rast(file.path(quantiles_dir, 
                                             paste0("CV_LAD_",
                                                    min_depth,
                                                    "_40.tif_lai_s2_masked_res_10_m_",
                                                    low_value,
                                                    "_",
                                                    high_value,
                                                    ".tif")))
      s2_lai_masked <- mask(s2_lai_masked, chm_mask, maskvalues = TRUE)
      s2_lai_values <- values(s2_lai_masked, na.rm = T)
      
      if (length(lidar_pad_values) < 100){
        correlation_value <- NA
      }
      else {
        correlation_value <- cor(lidar_pad_values,
                                 s2_lai_values,
                                 use = "pairwise.complete.obs")
      }
      correlations[i] <- correlation_value
      i <- i + 1
    }
    correlation_values[[i_quantile]] <- correlations
    i_quantile <- i_quantile + 1
  }
  correlation_values <- lapply(correlation_values, function(x) x[x != 0])
  
  png(file.path(figures_dir,
                paste0(site, 
                       "_correlation_lidar_s2_for_profiles_", 
                       composition_mask,
                       "_",
                       "CV_LAD", 
                       "_",
                       n_interval,
                       "_quantiles",
                       ".png")), 
      width = 1920, 
      height = 1080, 
      units = "px", res = 200)
  
  plot(NA, 
       type = "n",
       xlim = c(2, 38 + 3), 
       ylim = c(0, 0.8),  
       xlab = "LiDAR LAI Profile Depth",
       ylab = "Correlation Value",
       # main = "Correlation Values between LiDAR LAI Profiles of Varying Depth with Sentinel-2 LAI\n for 3 Standard Deviation Classes",
       lwd = 2,
       xaxt = "n",
       cex.main = 1.2,   # Increase title font size
       cex.lab = 1.2,    # Normal font size for axis labels
       cex.axis = 1.2)
  title(main = paste("Pearson Correlation between LiDAR LAI Profiles Integrated over Varying Depth\n",
                     "and Sentinel-2 LAI for 3 CV_LAD Classes\n",
                     "over Sentinel-2 pixels (10m resolution)"))
  
  colors <- brewer.pal(n_interval, "Paired")
  # colors <- c("#51a343", "#e6e6e6", "#9a4c00")
  for (i in seq_along(correlation_values)) {
    values <- correlation_values[[i]]
    lines(seq(2, 38, by = 1), values, 
          type = "l", col = colors[i], lwd = 2)
  }
  axis(side = 1, 
       at = seq(2, 38, by = 2.5),
       labels = seq(2, 38, by = 2.5))
  if (n_interval == 3){
    legend <- c("1: Lowest Heter", 
                "2",
                "3: Highest Heter")
  }
  else if (n_interval == 5){
    legend <- c("1: Lowest Heter", 
                "2",
                "3",
                "4",
                "5: Highest Heter")
  }
  else if (n_interval == 10){
    legend <- c("1: Lowest Heter", 
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10: Highest Heter")
  }
  else {
    stop("n_interval is not 3, 5, or 10")
  }
  legend("bottomright", legend = legend, col = colors, lwd = 2, cex = 1)
  dev.off()
}

perform_heterogeneity_depth_analysis <- function(site,
                                                 composition_mask,
                                                 metric,
                                                 chm,
                                                 quantiles_dir,
                                                 n_interval,
                                                 inc,
                                                 pad_profiles_dir,
                                                 figures_dir){
  correlation_values <- list()
  pad_files <- list.files(pad_profiles_dir, pattern = "\\.tif$", full.names = T)
  
  # Extract the numbers
  numbers <- as.numeric(sub(".*/PAD_(\\d+\\.\\d+)_.*", "\\1", pad_files))
  
  # Sort the numbers in descending order
  sorted_indices <- order(numbers, decreasing = TRUE)
  
  # Reorder the file list based on the sorted indices
  sorted_pad_files <- pad_files[sorted_indices]
  
  if (inc == 0.2) {
    max_val <- 0.8
  } else if (inc == 0.1) {
    max_val <- 0.9
  } else if (inc == 1/3) {
    max_val <- 2/3
  } else {
    stop("inc is not 0.1 or 0.2 or 1/3")
  }
  
  i_quantile <- 1
  for (low_value in seq(0, max_val, by = inc)) {
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    
    quantile <- terra::rast(file.path(quantiles_dir, 
                                      paste0(metric,
                                             "_heter_raster_res_10_m_",
                                             low_value,
                                             "_",
                                             high_value,
                                             ".tif")))
    correlations <- numeric(length(numbers))
    i <- 1
    for (pad_file in sorted_pad_files) {
      # Create a logical mask where CHM is greater than i
      chm_mask <- chm < i
      
      lidar_pad <- terra::rast(pad_file)
      lidar_pad_masked <- mask(lidar_pad, quantile)
      lidar_pad_masked <- mask(lidar_pad_masked, chm_mask, maskvalues = TRUE)
      lidar_pad_values <- values(lidar_pad_masked, na.rm = T)
      
      s2_lai <- terra::rast(file.path(pad_profiles_dir, "../s2lai_res_10_m.tif"))
      s2_lai_masked <- mask(s2_lai, quantile)
      s2_lai_masked <- mask(s2_lai_masked, chm_mask, maskvalues = TRUE)
      s2_lai_values <- values(s2_lai_masked, na.rm = T)
      
      if (grepl("cvladnorm", heterogeneity_metric)) {
        tmpdir <- paste0(dirname(pad_file), 
                         "/cvladnorm/", 
                         low_value, 
                         "_", 
                         high_value)
        
        if (!dir.exists(tmpdir)) {
          dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        }
        writeRaster(lidar_pad_masked,
                    filename = file.path(tmpdir,
                                         paste0("lidar", basename(pad_file))),
                    overwrite = T)
        writeRaster(s2_lai_masked,
                    filename = file.path(tmpdir,
                                         paste0("s2", basename(pad_file))),
                    overwrite = T)
      }
      
      if (length(lidar_pad_values) < 100){
        correlation_value <- NA
      }
      else {
        correlation_value <- cor(lidar_pad_values,
                                 s2_lai_values,
                                 use = "pairwise.complete.obs")
      }
      correlations[i] <- correlation_value
      i <- i + 1
    }
    correlation_values[[i_quantile]] <- correlations
    i_quantile <- i_quantile + 1
  }
  correlation_values <- lapply(correlation_values, function(x) x[x != 0])
  
  png(file.path(figures_dir,
                paste0(site, 
                       "_correlation_lidar_s2_for_profiles_", 
                       composition_mask,
                       "_",
                       metric, 
                       "_",
                       n_interval,
                       "_quantiles",
                       ".png")), 
      width = 1920, 
      height = 1080, 
      units = "px", res = 200)
  
  plot(NA, 
       type = "n",
       xlim = c(0, length(pad_files) + 2), 
       ylim = c(0, 0.8),  
       xlab = "LiDAR LAI Profile Depth",
       ylab = "Correlation Value",
       # main = "Correlation Values between LiDAR LAI Profiles of Varying Depth with Sentinel-2 LAI\n for 3 Standard Deviation Classes",
       lwd = 2,
       xaxt = "n",
       cex.main = 1.2,   # Increase title font size
       cex.lab = 1.2,    # Normal font size for axis labels
       cex.axis = 1.2)
  title(main = paste("Pearson Correlation between LiDAR LAI Profiles Integrated over Varying Depth\n",
                     "and Sentinel-2 LAI for 3 Local Heterogeneity Classes (Defined as the CHM\n",
                     metric, "(1m Resolution) over Sentinel-2 pixels (10m resolution)"))
  
  colors <- brewer.pal(n_interval, "Paired")
  # colors <- c("#51a343", "#e6e6e6", "#9a4c00")
  for (i in seq_along(correlation_values)) {
    values <- correlation_values[[i]]
    lines(seq(0, length(pad_files) - 1, by = 1), values, 
          type = "l", col = colors[i], lwd = 2)
  }
  axis(side = 1, 
       at = seq(0, length(pad_files) - 1, by = 2.5),
       labels = seq(0, length(pad_files) - 1, by = 2.5))
  if (n_interval == 3){
    legend <- c("1: Lowest Heter", 
                "2",
                "3: Highest Heter")
  }
  else if (n_interval == 5){
    legend <- c("1: Lowest Heter", 
                "2",
                "3",
                "4",
                "5: Highest Heter")
  }
  else if (n_interval == 10){
    legend <- c("1: Lowest Heter", 
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10: Highest Heter")
  }
  else {
    stop("n_interval is not 3, 5, or 10")
  }
  legend("bottomright", legend = legend, col = colors, lwd = 2, cex = 1)
  dev.off()
}

perform_heterogeneity_depth_analysis_2 <- function(site,
                                                   composition_mask,
                                                   metric,
                                                   chm,
                                                   quantiles_dir,
                                                   n_interval,
                                                   inc,
                                                   pad_profiles_dir,
                                                   figures_dir,
                                                   mask_chm_bool = FALSE){
  int_dir <- paste0(figures_dir, "/", n_interval, "_Quantiles")
  if (!dir.exists(int_dir)) {
    dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)
  }
  scatter_dir <- file.path(int_dir, "Scatterplots")
  if (!dir.exists(scatter_dir)) {
    dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  correlation_values <- list()
  # r2_values <- list()
  # rmse_values <- list()
  # slope_values <- list()
  pad_files <- list.files(pad_profiles_dir, pattern = "\\.tif$", full.names = T)
  
  # Extract the numbers
  numbers <- as.numeric(sub(".*/PAD_(\\d+\\.\\d+)_.*", "\\1", pad_files))
  
  # Sort the numbers in descending order
  sorted_indices <- order(numbers, decreasing = TRUE)
  
  # Reorder the file list based on the sorted indices
  sorted_pad_files <- pad_files[sorted_indices]
  
  if (inc == 0.2) {
    max_val <- 0.8
  } else if (inc == 0.1) {
    max_val <- 0.9
  } else if (inc == 1/3) {
    max_val <- 2/3
  } else {
    stop("inc is not 0.1 or 0.2 or 1/3")
  }
  
  i_quantile <- 1
  for (low_value in seq(0, max_val, by = inc)) {
    low_value <- round(low_value, digits = 2)
    high_value <- round(low_value + inc, digits = 2)
    
    quantile <- terra::rast(file.path(quantiles_dir, 
                                      paste0(metric,
                                             "_heter_raster_res_10_m_",
                                             low_value,
                                             "_",
                                             high_value,
                                             ".tif")))
    correlations <- numeric(length(numbers))
    # r2 <- numeric(length(numbers))
    # rmse <- numeric(length(numbers))
    # slopes <- numeric(length(numbers))
    
    i <- 1
    for (pad_file in sorted_pad_files) {
      # Create a logical mask where CHM is greater than i
      chm_mask <- chm < i
      
      lidar_pad <- terra::rast(pad_file)
      s2_lai <- terra::rast(file.path(pad_profiles_dir, "../s2lai_summer_res_10_m.tif"))
      
      na_index_lidar <- terra::countNA(lidar_pad)
      na_index_s2 <- terra::countNA(s2_lai)
      na_index <- na_index_lidar + na_index_s2
      lidar_pad[na_index > 0] <- NA
      s2_lai[na_index > 0] <- NA
      
      lidar_pad_masked <- mask(lidar_pad, quantile)
      s2_lai_masked <- mask(s2_lai, quantile)
      
      if (mask_chm_bool){
        lidar_pad_masked <- mask(lidar_pad_masked, chm_mask, maskvalues = TRUE)
        s2_lai_masked <- mask(s2_lai_masked, chm_mask, maskvalues = TRUE)
        lidar_pad_values <- values(lidar_pad_masked, na.rm = T)
        s2_lai_values <- values(s2_lai_masked, na.rm = T)
        mask_write <- "mask_chm_true"
        if (length(lidar_pad_values) > 100){
          # plot_density_scatterplot(s2_lai_values,
          #                          lidar_pad_values,
          #                          "Sentinel-2 LAI",
          #                          "LiDAR PAD",
          #                          paste("Density Scatterplot at", i, "meters from top of canopy"),
          #                          scatter_dir,
          #                          paste0("zscatter_", i, "_meters_from_top",
          #                                 n_interval, "_classes_from_",
          #                                 low_value, "_to_", high_value))
          tmp <- 0
        }
      } else {
        mask_write <- "mask_chm_false"
        lidar_pad_values <- values(lidar_pad_masked, na.rm = T)
        s2_lai_values <- values(s2_lai_masked, na.rm = T)
      }
      
      # if (grepl("cvladnorm", heterogeneity_metric)) {
      #   tmpdir <- paste0(dirname(pad_file), 
      #                    "/cvladnorm/", 
      #                    low_value, 
      #                    "_", 
      #                    high_value)
      #   
      #   if (!dir.exists(tmpdir)) {
      #     dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
      #   }
      #   writeRaster(lidar_pad_masked,
      #               filename = file.path(tmpdir,
      #                                    paste0("lidar", basename(pad_file))),
      #               overwrite = T)
      #   writeRaster(s2_lai_masked,
      #               filename = file.path(tmpdir,
      #                                    paste0("s2", basename(pad_file))),
      #               overwrite = T)
      # }
      
      if (length(lidar_pad_values) < 100){
        correlation_value <- NA
        # r2_value <- NA
        # rmse_value <- NA
        # slope_value <- NA
      }
      else {
        # Correlation
        correlation_value <- cor(lidar_pad_values, s2_lai_values, use = "pairwise.complete.obs")
        
        # Linear Regression for R², RMSE, and Slope
        # model <- lm(s2_lai_values ~ lidar_pad_values)
        # r2_value <- summary(model)$r.squared
        # rmse_value <- sqrt(mean((s2_lai_values - predict(model))^2))
        # slope_value <- coef(model)[2]
      }
      
      correlations[i] <- correlation_value
      # r2[i] <- r2_value
      # rmse[i] <- rmse_value
      # slopes[i] <- slope_value
      
      i <- i + 1
    }
    
    correlation_values[[i_quantile]] <- correlations
    # r2_values[[i_quantile]] <- r2
    # rmse_values[[i_quantile]] <- rmse
    # slope_values[[i_quantile]] <- slopes
    
    i_quantile <- i_quantile + 1
  }
  
  correlation_values <- lapply(correlation_values, function(x) x[x != 0])
  # r2_values <- lapply(r2_values, function(x) x[x != 0])
  # rmse_values <- lapply(rmse_values, function(x) x[x != 0])
  # slope_values <- lapply(slope_values, function(x) x[x != 0])
  
  # Colors for plotting
  colors <- brewer.pal(n_interval, "Paired")
  
  # Plot Correlation
  png(file.path(int_dir, paste0("correlation_", site, 
                                "_", composition_mask, 
                                "_depth_analysis_", mask_write, ".png")), 
      width = 1920, height = 1080, units = "px", res = 200)
  plot_metric(correlation_values, pad_files, "Correlation", colors, n_interval)
  dev.off()
  
  # Plot R²
  # png(file.path(int_dir, paste0("r2_", site, 
  #                               "_", composition_mask, 
  #                               "_depth_analysis_", mask_write, ".png")), 
  #     width = 1920, height = 1080, units = "px", res = 200)
  # plot_metric(r2_values, pad_files, "R²", colors, n_interval)
  # dev.off()
  
  # Plot RMSE
  # png(file.path(int_dir, paste0("rmse_", site, 
  #                               "_", composition_mask, 
  #                               "_depth_analysis_", mask_write, ".png")), 
  #     width = 1920, height = 1080, units = "px", res = 200)
  # plot_metric(rmse_values, pad_files, "RMSE", colors, n_interval)
  # dev.off()
  
  # Plot Slope
  # png(file.path(int_dir, paste0("slope_", site, 
  #                               "_", composition_mask, 
  #                               "_depth_analysis_", mask_write, ".png")), 
  #     width = 1920, height = 1080, units = "px", res = 200)
  # plot_metric(slope_values, pad_files, "Slope", colors, n_interval)
  # dev.off()
}

# Helper function to plot each metric
plot_metric <- function(metric_values, pad_files, metric_name, colors, n_interval) {
  plot(NA, 
       type = "n",
       xlim = c(0, length(pad_files) + 2), 
       ylim = c(0, 1),  # Adjust `ylim` as needed for each metric
       xlab = "LiDAR LAI Profile Depth",
       ylab = paste(metric_name, "Value"),
       lwd = 2,
       xaxt = "n",
       cex.lab = 1.2,
       cex.axis = 1.2)
  title(main = paste(metric_name, "between LiDAR LAI Profiles and Sentinel-2 LAI for Varying Depth"))
  
  for (i in seq_along(metric_values)) {
    values <- metric_values[[i]]
    lines(seq(1, length(pad_files), by = 1), values, 
          type = "l", col = colors[i], lwd = 2)
  }
  axis(side = 1, at = seq(1, length(pad_files), by = 1), labels = seq(1, length(pad_files), by = 1))
  
  # Add legend
  if (n_interval == 3) {
    legend_labels <- c("Class 1: min to Q1", 
                       "Class 2: Q1 to Q2", 
                       "Class 3: Q2 to max")
  } else if (n_interval == 5) {
    legend_labels <- c("1: Lowest Heterogeneity", 
                       "2: Low Heterogeneity", 
                       "3: Medium Heterogeneity",
                       "4: High Heterogeneity", 
                       "5: Highest Heterogeneity")
  } else if (n_interval == 10) {
    legend_labels <- paste(1:10, "Heterogeneity")
  } else {
    stop("n_interval is not 3, 5, or 10")
  }
  
  legend("bottomright", legend = legend_labels, col = colors, lwd = 2, cex = 0.8)
}
