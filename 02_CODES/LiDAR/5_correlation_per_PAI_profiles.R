# ---
# title: "3.correlation_per_PAI_profiles.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-18"
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

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../functions_plots.R")
source("functions_lidar.R")

# Pre-processing Parameters & Directories
results_dir <- '../../03_RESULTS'
site <- "Blois" # Mormal Blois Aigoual
results_path <- file.path("../../03_RESULTS", site)
lidar_dir <- "LiDAR"
masks_dir <- file.path(results_path, lidar_dir, "Heterogeneity_Masks")
dir.create(path = masks_dir, showWarnings = F, recursive = T)
profiles_dir <- file.path(results_path, lidar_dir, "Profiles")

# Intervals Directories
deciles_dir <- "Deciles"
equal_intervals_dir <- "Equal_Intervals"
deciles_dir <- file.path(masks_dir, "Deciles")
equal_intervals_dir <- file.path(masks_dir, "Equal_Intervals")
class_dir <- file.path(deciles_dir, "10_Quantiles")

# Open LiDAR S-2 LAI
mask_10m <- terra::rast(file.path(masks_dir,
                                  "artifacts_low_vegetation_majority_90_p_res_10_m.envi"))
lai_lidar_raster <- terra::rast(file.path(masks_dir,
                                          "lai_lidar_masked_res_10_m.envi"))
lai_lidar <- values(lai_lidar_raster)
lai_s2_raster <- terra::rast(file.path(masks_dir,
                                       "lai_s2_masked_res_10_m.envi"))
lai_s2 <- values(lai_s2_raster)

# Initial correlation
cor <- correlation_test_function(lai_lidar, lai_s2)
plot_density_scatterplot(var_x = lai_s2,
                         var_y = lai_lidar,
                         xlab = "LAI PROSAIL",
                         ylab = "LAI LiDAR",
                         dirname = profiles_dir,
                         filename = paste0("initial_correlation_lidar_s2",
                                           ".tif"),
                         xlimits = c(0, 10),
                         ylimits = c(0, 10)
)

# -------------------- Heterogeneity Correlations: mean ------------------------
correlation_values <- c()
inc <- 0.1
for (low_value in seq(0, 0.9, by = inc)) {
  lidar_intervals <- terra::rast(file.path(class_dir, 
                                           paste0("mean_heights_pai_masked_res_10_m_", 
                                                  low_value, 
                                                  "_", 
                                                  low_value + inc,
                                                  ".tif")))
  lidar_values <- values(lidar_intervals)
  
  s2_intervals <- terra::rast(file.path(class_dir, 
                                        paste0("mean_heights_lai_s2_masked_res_10_m_", 
                                               low_value, 
                                               "_", 
                                               low_value + inc,
                                               ".tif")))
  s2_values <- values(s2_intervals)
  
  correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
  
  # plot_density_scatterplot(var_x = s2_values,
  #                          var_y = lidar_values,
  #                          xlab = "LAI PROSAIL",
  #                          ylab = "LAI LiDAR",
  #                          dirname = profiles_dir,
  #                          filename = paste0("mean_heights_scatter_lai_lidar_prosail_res_10_m_", 
  #                                            low_value, 
  #                                            "_", 
  #                                            low_value + inc,
  #                                            ".tif"),
  #                          xlimits = c(0, 10),
  #                          ylimits = c(0, 10)
  # )
}

png(filename = file.path(profiles_dir, 
                         "mean_heights_correlation_values_for_10_quantiles_res_10_m.png"),
    width = 1920, height = 1080)  # Adjust width and height as needed
cor_plot <- plot(1, type = "n", 
                 xlab = "Quantile", 
                 ylab = "Correlation Value", 
                 xlim = c(0, 1), ylim = c(0, 0.8),
                 main = "Correlation Values for 10 Mean Heights Quantiles",
                 xaxt = "n",
                 cex.main = 1.2,   # Increase title font size
                 cex.lab = 1.2,    # Normal font size for axis labels
                 cex.axis = 1.2)
midpoints <- seq(0.05, 0.95, by = inc)
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2), cex.axis = 1.2)
dev.off()

# -------------------- Heterogeneity Correlations: max -------------------------
correlation_values <- c()
inc <- 0.1
for (low_value in seq(0, 0.9, by = inc)) {
  lidar_intervals <- terra::rast(file.path(class_dir, 
                                           paste0("max_heights_pai_masked_res_10_m_", 
                                                  low_value, 
                                                  "_", 
                                                  low_value + inc,
                                                  ".tif")))
  lidar_values <- values(lidar_intervals)
  
  s2_intervals <- terra::rast(file.path(class_dir, 
                                        paste0("max_heights_lai_s2_masked_res_10_m_", 
                                               low_value, 
                                               "_", 
                                               low_value + inc,
                                               ".tif")))
  s2_values <- values(s2_intervals)
  
  correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
  
  # plot_density_scatterplot(var_x = s2_values,
  #                          var_y = lidar_values,
  #                          xlab = "LAI PROSAIL",
  #                          ylab = "LAI LiDAR",
  #                          dirname = profiles_dir,
  #                          filename = paste0("mean_heights_scatter_lai_lidar_prosail_res_10_m_", 
  #                                            low_value, 
  #                                            "_", 
  #                                            low_value + inc,
  #                                            ".tif"),
  #                          xlimits = c(0, 10),
  #                          ylimits = c(0, 10)
  # )
}

png(filename = file.path(profiles_dir, 
                         "max_heights_correlation_values_for_10_quantiles_res_10_m.png"),
    width = 1920, height = 1080)  # Adjust width and height as needed
cor_plot <- plot(1, type = "n", 
                 xlab = "Quantile", 
                 ylab = "Correlation Value", 
                 xlim = c(0, 1), ylim = c(0, 0.8),
                 main = "Correlation Values for 10 Max Quantiles",
                 xaxt = "n",
                 cex.main = 1.2,   # Increase title font size
                 cex.lab = 1.2,    # Normal font size for axis labels
                 cex.axis = 1.2)
midpoints <- seq(0.05, 0.95, by = inc)
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2), cex.axis = 1.2)
dev.off()

# --------------------- Heterogeneity Correlations: cv -------------------------
correlation_values <- c()
inc <- 0.1
for (low_value in seq(0, 0.9, by = inc)) {
  lidar_intervals <- terra::rast(file.path(class_dir, 
                                           paste0("cv_pai_masked_res_10_m_", 
                                                  low_value, 
                                                  "_", 
                                                  low_value + inc,
                                                  ".tif")))
  lidar_values <- values(lidar_intervals)
  
  s2_intervals <- terra::rast(file.path(class_dir, 
                                        paste0("cv_lai_s2_masked_res_10_m_", 
                                               low_value, 
                                               "_", 
                                               low_value + inc,
                                               ".tif")))
  s2_values <- values(s2_intervals)
  
  correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
  
  plot_density_scatterplot(var_x = s2_values,
                           var_y = lidar_values,
                           xlab = "LAI PROSAIL",
                           ylab = "LAI LiDAR",
                           dirname = profiles_dir,
                           filename = paste0("cv_scatter_lai_lidar_prosail_res_10_m_", 
                                             low_value, 
                                             "_", 
                                             low_value + inc,
                                             ".tif"),
                           xlimits = c(0, 10),
                           ylimits = c(0, 10)
  )
}

png(filename = file.path(profiles_dir, 
                         "cv_correlation_values_for_10_quantiles_res_10_m.png"),
    width = 800, height = 600)  # Adjust width and height as needed
cor_plot <- plot(1, type = "n", 
                 xlab = "Quantile", 
                 ylab = "Correlation Value", 
                 xlim = c(0, 1), ylim = c(0, 0.8),
                 main = "Correlation Values for 10 CV Quantiles")
midpoints <- seq(0.05, 0.95, by = inc)
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
dev.off()

# --------------------- Heterogeneity Correlations: std ------------------------
correlation_values <- c()
inc <- 0.1
for (low_value in seq(0, 0.9, by = inc)) {
  lidar_intervals <- terra::rast(file.path(class_dir, 
                                           paste0("std_pai_masked_res_10_m_", 
                                                  low_value, 
                                                  "_", 
                                                  low_value + inc,
                                                  ".tif")))
  lidar_values <- values(lidar_intervals)
  
  s2_intervals <- terra::rast(file.path(class_dir, 
                                        paste0("std_lai_s2_masked_res_10_m_", 
                                               low_value, 
                                               "_", 
                                               low_value + inc,
                                               ".tif")))
  s2_values <- values(s2_intervals)
  
  correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
  
  plot_density_scatterplot(var_x = s2_values,
                           var_y = lidar_values,
                           xlab = "LAI PROSAIL",
                           ylab = "LAI LiDAR",
                           dirname = profiles_dir,
                           filename = paste0("std_scatter_lai_lidar_prosail_res_10_m_", 
                                             low_value, 
                                             "_", 
                                             low_value + inc,
                                             ".tif"),
                           xlimits = c(0, 10),
                           ylimits = c(0, 10)
  )
}

png(filename = file.path(profiles_dir, 
                         "std_correlation_values_for_10_quantiles_res_10_m.png"),
    width = 800, height = 600)  # Adjust width and height as needed
cor_plot <- plot(1, type = "n", 
                 xlab = "Quantile", 
                 ylab = "Correlation Value", 
                 xlim = c(0, 1), ylim = c(0, 0.8),
                 main = "Correlation Values for 10 Std Quantiles")
midpoints <- seq(0.05, 0.95, by = inc)
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
dev.off()

# -------------------- Heterogeneity Correlations: rumple ----------------------
correlation_values <- list()
classs_dir <- file.path(deciles_dir, "10_Quantiles")
inc <- 0.1

for (low_quantile in seq(0, 0.9, by = inc)){
  low_quantile <- round(low_quantile, digits = 2)
  high_quantile <- round(low_quantile + inc, digits = 2)
  lidar_intervals <- terra::rast(file.path(classs_dir, 
                                           paste0("rumple_pai_masked_res_10_m_", 
                                                  low_quantile, "_", 
                                                  high_quantile, ".tif")))
  lidar_values <- values(lidar_intervals)
  
  s2_intervals <- terra::rast(file.path(classs_dir, 
                                        paste0("rumple_lai_s2_masked_res_10_m_", 
                                               low_quantile, "_", 
                                               high_quantile, ".tif")))
  s2_values <- values(s2_intervals)
  
  correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
  correlation_values <- c(correlation_values, correlation_value)
  
  plot_density_scatterplot(var_x = s2_values,
                           var_y = lidar_values,
                           xlab = "LAI PROSAIL",
                           ylab = "LAI LiDAR",
                           dirname = profiles_dir,
                           filename = paste0("rumple_scatter_lai_lidar_prosail_res_10_m_", 
                                             low_quantile, "_", 
                                             high_quantile, ".tif"),
                           xlimits = c(0, 10),
                           ylimits = c(0, 10)
  )
}

png(filename = file.path(profiles_dir, 
                         "rumple_correlation_values_for_10_quantiles_res_10_m.png"),
    width = 800, height = 600)  # Adjust width and height as needed
cor_plot <- plot(1, type = "n", 
                 xlab = "Quantile", 
                 ylab = "Correlation Value", 
                 xlim = c(0, 1), ylim = c(0, 0.8),
                 main = "Correlation Values for 3 Rumple Quantiles")
midpoints <- seq(0.05, 0.95, by = inc)
lines(midpoints, correlation_values, type = "l")
points(midpoints, correlation_values, pch = 16, col = "red")
axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
dev.off()

# ---------------------- Vegetation Depth Correlations -------------------------
correlation_values <- list()

# for (i in seq(46.5, 2.5, by = -1)) {
for (i in seq(43.5, 2.5, by = -1)) {
  # pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
  pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_44.18", ".tif"))
  pad_lidar <- terra::rast(pad_lidar_file)
  pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
  
  lidar_masked <- pad_lidar_proj * mask_10m
  lidar_values <- values(lidar_masked)
  
  correlation_value <- cor(lidar_values, lai_s2, use = "pairwise.complete.obs")
  correlation_values[i] <- correlation_value
  print(i)
}
correlation_values <- Filter(Negate(is.null), correlation_values)

png(file.path(profiles_dir, "initial_correlation_lidar_s2_for_profiles.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

plot(seq(0, 41, by = 1), rev(correlation_values), type = "l", 
     xlim = c(0, 44), ylim = c(0, 0.65), 
     xlab = "LiDAR LAI Profile Depth", ylab = "Correlation Value",
     # main = "Correlation Values between LiDAR LAI Profiles of Varying Depth and Sentinel-2 LAI", 
     lwd = 2,
     xaxt = "n",
     cex.main = 1.2,   # Increase title font size
     cex.lab = 1.2,    # Normal font size for axis labels
     cex.axis = 1.2)
title(main = paste("Pearson Correlation between LiDAR LAI Profiles\n",
                   "Integrated over Varying Depth and Sentinel-2 LAI"))
axis(side = 1, at = seq(0, 40.5, by = 2.5), labels = seq(0, 41, by = 2.5))
points(seq(0, 41, by = 1), rev(correlation_values), pch = 1)
dev.off()

# ---------- Vegetation Depth + Heterogeneity Correlations (std 3 classes) ---------------
correlation_values <- list()
classs_dir <- file.path(deciles_dir, "3_Quantiles")
inc <- 1/3

for (low_quantile in seq(0, 2/3, by = inc)){
  low_quantile <- round(low_quantile, digits = 2)
  high_quantile <- round(low_quantile + inc, digits = 2)
  quantile <- terra::rast(file.path(classs_dir, 
                                    sprintf("std_heter_raster_res_10_m_%s_%s.tif",
                                            low_quantile,
                                            high_quantile)))
  correlations <- numeric(length(seq(43.5, 2.5, by = -1)))
  
  for (i in seq(43.5, 2.5, by = -1)) {
    pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_44.18", ".tif"))
    pad_lidar <- terra::rast(pad_lidar_file)
    pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
    
    lidar_masked <- pad_lidar_proj * mask_10m
    lidar_masked <- mask(lidar_masked, quantile)
    lidar_values <- values(lidar_masked)
    
    # Density Scatterplot
    # plot_density_scatterplot(var_x = sentinel_values,
    #                          var_y = lidar_values,
    #                          xlab = "LAI PROSAIL",
    #                          ylab = sprintf("LAI LiDAR %s to 46.69", i),
    #                          dirname = profiles_dir,
    #                          filename = sprintf("std_depth_scatter_lidar_%s_46.69_sentinel2_%s_%s",
    #                                             i, low_quantile, high_quantile),
    #                          xlimits = c(0, 10),
    #                          ylimits = c(0, 20)
    # )
    print(i)
    # residuals <- lidar_masked - lai_sentinel_masked
    # plot(residuals)
    # writeRaster(residuals, 
    #             file.path(profiles_dir,
    #                       "Residuals",
    #                       sprintf("residuals_%s_to_46.69.envi", i)))
    
    correlation_value <- cor(lidar_values, lai_s2, use = "pairwise.complete.obs")
    correlations[i] <- correlation_value
  }
  correlation_values[[length(correlation_values) + 1]] <- correlations
}
correlation_values <- lapply(correlation_values, function(x) x[x != 0])

png(file.path(profiles_dir, "correlation_lidar_s2_for_profiles_std.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.65), 
     xlab = "LiDAR LAI Profile Depth", ylab = "Correlation Value",
     # main = "Correlation Values between LiDAR LAI Profiles of Varying Depth with Sentinel-2 LAI\n for 3 Standard Deviation Classes",
     lwd = 2,
     xaxt = "n",
     cex.main = 1.2,   # Increase title font size
     cex.lab = 1.2,    # Normal font size for axis labels
     cex.axis = 1.2)
title(main = paste("Pearson Correlation between LiDAR LAI Profiles Integrated over Varying Depth\n",
                   "and Sentinel-2 LAI for 3 Local Heterogeneity Classes (Defined as the CHM\n",
                   "Standard Deviation (1m Resolution) over Sentinel-2 pixels (10m resolution)"))

# colors <- brewer.pal(length(correlation_values), "Paired")
colors <- c("#51a343", "#e6e6e6", "#9a4c00")
for (i in seq_along(correlation_values)) {
  values <- correlation_values[[i]]
  values <- rev(values)
  lines(seq(0, 41, by = 1), values, 
        type = "l", col = colors[i], lwd = 2)
}
axis(side = 1, at = seq(0, 40.5, by = 2.5), labels = seq(0, 41, by = 2.5))
legend("bottomright", legend = c("Low CHM Heterogeneity", 
                                 "Medium CHM Heterogeneity",
                                 "High CHM Heterogeneity"), col = colors, 
       lwd = 2, cex = 1)
dev.off()

# ---------- Vegetation Depth + Heterogeneity Correlations (cv 3 classes) ----------------
correlation_values <- list()
classs_dir <- file.path(deciles_dir, "3_Quantiles")
inc <- 1/3

for (low_quantile in seq(0, 2/3, by = inc)){
  low_quantile <- round(low_quantile, digits = 2)
  high_quantile <- round(low_quantile + inc, digits = 2)
  quantile <- terra::rast(file.path(classs_dir, 
                                    sprintf("cv_heter_raster_res_10_m_%s_%s.tif",
                                            low_quantile,
                                            high_quantile)))
  correlations <- numeric(length(seq(46.5, 2.5, by = -1)))
  
  for (i in seq(46.5, 2.5, by = -1)) {
    pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
    pad_lidar <- terra::rast(pad_lidar_file)
    pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
    
    lidar_masked <- pad_lidar_proj * mask_10m
    lidar_masked <- mask(lidar_masked, quantile)
    lidar_values <- values(lidar_masked)
    
    lai_sentinel_maskedd <- mask(lai_sentinel_masked, quantile)
    sentinel_values <- values(lai_sentinel_maskedd)
    
    # Density Scatterplot
    # plot_density_scatterplot(var_x = sentinel_values,
    #                          var_y = lidar_values,
    #                          xlab = "LAI PROSAIL",
    #                          ylab = sprintf("LAI LiDAR %s to 46.69", i),
    #                          dirname = profiles_dir,
    #                          filename = sprintf("cv_depth_scatter_lidar_%s_46.69_sentinel2_%s_%s",
    #                                             i, low_quantile, high_quantile),
    #                          xlimits = c(0, 10),
    #                          ylimits = c(0, 20)
    # )
    print(i)
    # residuals <- lidar_masked - lai_sentinel_masked
    # plot(residuals)
    # writeRaster(residuals, 
    #             file.path(profiles_dir,
    #                       "Residuals",
    #                       sprintf("residuals_%s_to_46.69.envi", i)))
    
    correlation_value <- cor(lidar_values, sentinel_values, use = "pairwise.complete.obs")
    correlations[i] <- correlation_value
  }
  correlation_values[[length(correlation_values) + 1]] <- correlations
}
correlation_values <- lapply(correlation_values, function(x) x[x != 0])

png(file.path(profiles_dir, "correlation_lidar_s2_for_profiles_cv.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.7), 
     xlab = "PAI Profile Depth", ylab = "Correlation Value",
     main = "Correlation Values of Intervals PAI Profiles of Varying Depth with Sentinel-2 LAI for Cv")

colors <- brewer.pal(length(correlation_values), "Paired")
for (i in seq_along(correlation_values)) {
  values <- correlation_values[[i]]
  values <- rev(values)
  lines(seq(0, 44, by = 1), values, 
        type = "l", col = colors[i], lwd = 2)
}
legend("bottomright", legend = c("Homogeneous Stand 1", 
                                 "Homogeneous Stand 2",
                                 "Homogeneous Stand 3"), col = colors, 
       lwd = 2, cex = 0.8)
dev.off()
# ---------- Vegetation Depth + Heterogeneity Correlations (std 10 classes) ---------------
correlation_values <- list()
classs_dir <- file.path(deciles_dir, "10_Quantiles")
inc <- 0.1

for (low_quantile in seq(0, 0.9, by = inc)){
  low_quantile <- round(low_quantile, digits = 2)
  high_quantile <- round(low_quantile + inc, digits = 2)
  quantile <- terra::rast(file.path(classs_dir, 
                                    sprintf("std_heter_raster_res_10_m_%s_%s.tif",
                                            low_quantile,
                                            high_quantile)))
  correlations <- numeric(length(seq(46.5, 2.5, by = -1)))
  
  for (i in seq(46.5, 2.5, by = -1)) {
    pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
    pad_lidar <- terra::rast(pad_lidar_file)
    pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
    
    lidar_masked <- pad_lidar_proj * mask_10m
    lidar_masked <- mask(lidar_masked, quantile)
    lidar_values <- values(lidar_masked)
    
    lai_sentinel_maskedd <- mask(lai_sentinel_masked, quantile)
    sentinel_values <- values(lai_sentinel_maskedd)
    
    # Density Scatterplot
    # plot_density_scatterplot(var_x = sentinel_values,
    #                          var_y = lidar_values,
    #                          xlab = "LAI PROSAIL",
    #                          ylab = sprintf("LAI LiDAR %s to 46.69", i),
    #                          dirname = profiles_dir,
    #                          filename = sprintf("std_depth_scatter_lidar_%s_46.69_sentinel2_%s_%s",
    #                                             i, low_quantile, high_quantile),
    #                          xlimits = c(0, 10),
    #                          ylimits = c(0, 20)
    # )
    print(i)
    # residuals <- lidar_masked - lai_sentinel_masked
    # plot(residuals)
    # writeRaster(residuals, 
    #             file.path(profiles_dir,
    #                       "Residuals",
    #                       sprintf("residuals_%s_to_46.69.envi", i)))
    
    correlation_value <- cor(lidar_values, sentinel_values, use = "pairwise.complete.obs")
    correlations[i] <- correlation_value
  }
  correlation_values[[length(correlation_values) + 1]] <- correlations
}
correlation_values <- lapply(correlation_values, function(x) x[x != 0])

png(file.path(profiles_dir, "correlation_lidar_s2_for_profiles_std10c.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.7), 
     xlab = "PAI Profile Depth", ylab = "Correlation Value",
     main = "Correlation Values of Intervals PAI Profiles of Varying Depth with Sentinel-2 LAI for Std10c")

colors <- brewer.pal(length(correlation_values), "Paired")
for (i in seq_along(correlation_values)) {
  values <- correlation_values[[i]]
  values <- rev(values)
  lines(seq(0, 44, by = 1), values, 
        type = "l", col = colors[i], lwd = 2)
}
legend("bottomright", legend = c("Homogeneous Stand 1", 
                                 "Homogeneous Stand 2",
                                 "Homogeneous Stand 3",
                                 "Homogeneous Stand 4", 
                                 "Homogeneous Stand 5",
                                 "Homogeneous Stand 6",
                                 "Homogeneous Stand 7", 
                                 "Homogeneous Stand 8",
                                 "Homogeneous Stand 9",
                                 "Homogeneous Stand 10"), col = colors, 
       lwd = 2, cex = 0.8)
dev.off()

# ---------- Vegetation Depth + Heterogeneity Correlations (cv 10 classes) ----------------
correlation_values <- list()
classs_dir <- file.path(deciles_dir, "10_Quantiles")
inc <- 0.1

for (low_quantile in seq(0, 0.9, by = inc)){
  low_quantile <- round(low_quantile, digits = 2)
  high_quantile <- round(low_quantile + inc, digits = 2)
  quantile <- terra::rast(file.path(classs_dir, 
                                    sprintf("cv_heter_raster_res_10_m_%s_%s.tif",
                                            low_quantile,
                                            high_quantile)))
  correlations <- numeric(length(seq(46.5, 2.5, by = -1)))
  
  for (i in seq(46.5, 2.5, by = -1)) {
    pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
    pad_lidar <- terra::rast(pad_lidar_file)
    pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
    
    lidar_masked <- pad_lidar_proj * mask_10m
    lidar_masked <- mask(lidar_masked, quantile)
    lidar_values <- values(lidar_masked)
    
    lai_sentinel_maskedd <- mask(lai_sentinel_masked, quantile)
    sentinel_values <- values(lai_sentinel_maskedd)
    
    # Density Scatterplot
    # plot_density_scatterplot(var_x = sentinel_values,
    #                          var_y = lidar_values,
    #                          xlab = "LAI PROSAIL",
    #                          ylab = sprintf("LAI LiDAR %s to 46.69", i),
    #                          dirname = profiles_dir,
    #                          filename = sprintf("cv_depth_scatter_lidar_%s_46.69_sentinel2_%s_%s",
    #                                             i, low_quantile, high_quantile),
    #                          xlimits = c(0, 10),
    #                          ylimits = c(0, 20)
    # )
    print(i)
    # residuals <- lidar_masked - lai_sentinel_masked
    # plot(residuals)
    # writeRaster(residuals, 
    #             file.path(profiles_dir,
    #                       "Residuals",
    #                       sprintf("residuals_%s_to_46.69.envi", i)))
    
    correlation_value <- cor(lidar_values, sentinel_values, use = "pairwise.complete.obs")
    correlations[i] <- correlation_value
  }
  correlation_values[[length(correlation_values) + 1]] <- correlations
}
correlation_values <- lapply(correlation_values, function(x) x[x != 0])

png(file.path(profiles_dir, "correlation_lidar_s2_for_profiles_cv10c.png"), 
    width = 1920, 
    height = 1080, 
    units = "px", res = 200)

plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.7), 
     xlab = "PAI Profile Depth", ylab = "Correlation Value",
     main = "Correlation Values of Intervals PAI Profiles of Varying Depth with Sentinel-2 LAI for Cv10c")

colors <- brewer.pal(length(correlation_values), "Paired")
for (i in seq_along(correlation_values)) {
  values <- correlation_values[[i]]
  values <- rev(values)
  lines(seq(0, 44, by = 1), values, 
        type = "l", col = colors[i], lwd = 2)
}
legend("bottomright", legend = c("Homogeneous Stand 1", 
                                 "Homogeneous Stand 2",
                                 "Homogeneous Stand 3",
                                 "Homogeneous Stand 4", 
                                 "Homogeneous Stand 5",
                                 "Homogeneous Stand 6",
                                 "Homogeneous Stand 7", 
                                 "Homogeneous Stand 8",
                                 "Homogeneous Stand 9",
                                 "Homogeneous Stand 10"), col = colors, 
       lwd = 2, cex = 0.8)
dev.off()
# ------------------------------------------------------------------------------




# 
# # Vector to store correlations
# # correlations <- numeric(0)
# # correlation_values <- c()
# correlation_values <- list()
# 
# error <- terra::rast(file.path(masks_dir, "mnc_std_res_10_m.envi"))
# 
# # Test
# # q1 <- quantile(values(error), probs = c(0, 1/3), na.rm = TRUE)
# # q2 <- quantile(values(error), probs = c(1/3, 2/3), na.rm = TRUE)
# # q3 <- quantile(values(error), probs = c(2/3, 3/3), na.rm = TRUE)
# # mask_q1 <- mask(error, mask = error >= q1[2] | error < q1[1], maskvalue = 1)
# # mask_q2 <- mask(error, mask = error >= q2[2] | error < q2[1], maskvalue = 1)
# # mask_q3 <- mask(error, mask = error >= q3[2] | error < q3[1], maskvalue = 1)
# # writeRaster(mask_q1, file.path(masks_dir, "0q1.envi"), overwrite=T)
# # writeRaster(mask_q2, file.path(masks_dir, "q1q2.envi"), overwrite=T)
# # writeRaster(mask_q3, file.path(masks_dir, "q2max.envi"), overwrite=T)
# # errors <- list(mask_q1, mask_q2, mask_q3)
# 
# for (error in errors){
#   correlations <- numeric(length(seq(46.5, 2.5, by = -1)))
#   # Loop through each raster file
#   for (i in seq(46.5, 2.5, by = -1)) {
#     # Open LiDAR PAD
#     pad_lidar_file <- file.path(profiles_dir, paste0("PAD_", i, "_46.69", ".tif"))
#     pad_lidar <- terra::rast(pad_lidar_file)
#     pad_lidar_proj <- terra::project(pad_lidar, mask_10m)
#     lidar_masked <- pad_lidar_proj * mask_10m
#     lidar_masked <- mask(lidar_masked, error)
#     lidar_values <- values(lidar_masked)
#     
#     lai_sentinel_maskedd <- mask(lai_sentinel_masked, error)
#     sentinel_values <- values(lai_sentinel_maskedd)
#     
#     # Density Scatterplot
#     # plot_density_scatterplot(var_x = sentinel_values,
#     #                          var_y = lidar_values,
#     #                          xlab = "LAI PROSAIL",
#     #                          ylab = sprintf("LAI LiDAR %s to 46.69", i),
#     #                          dirname = profiles_dir,
#     #                          filename = sprintf("scatter_lidar_%s_46.69_sentinel2", 
#     #                                             i),
#     #                          xlimits = c(0, 10),
#     #                          ylimits = c(0, 20)
#     # )
#     print(i)
#     # residuals <- lidar_masked - lai_sentinel_masked
#     # plot(residuals)
#     # writeRaster(residuals, 
#     #             file.path(profiles_dir,
#     #                       "Residuals",
#     #                       sprintf("residuals_%s_to_46.69.envi", i)))
#     
#     # Store correlation
#     # correlations <- c(correlations, cor(lidar_values, sentinel_values, 
#     #                                     use = "pairwise.complete.obs"))
#     correlation_value <- cor(lidar_values, sentinel_values, use = "pairwise.complete.obs")
#     correlations[i] <- correlation_value
#   }
#   correlation_values[[length(correlation_values) + 1]] <- correlations
#   correlation_values <- lapply(correlation_values, function(x) x[x != 0])
# }
# png(file.path(profiles_dir, "correlation_lidar_s2_for_profiles_std.png"), 
#     width = 1920, 
#     height = 1080, 
#     units = "px", res = 200)
# # Initialize the plot
# plot(NA, type = "n", xlim = c(0, 44), ylim = c(0, 0.7), 
#      xlab = "PAI Profile Depth", ylab = "Correlation Value",
#      main = "Correlation Values of Intervals PAI Profiles of Varying Depth with Sentinel-2 LAI for Std")
# 
# colors <- brewer.pal(length(correlation_values), "Paired")
# 
# for (i in seq_along(correlation_values)) {
#   values <- correlation_values[[i]]
#   values <- rev(values)
#   # Check for NAs in the data
#   if (any(is.na(values))) {
#     next  # Skip this interval if there are NAs in the data
#   }
#   # lines(seq(0, 44, by = 1), values, 
#   #       type = "l", col = colors[i],
#   #       lty = line_styles[i %% length(line_styles) + 1],  # Cycling through line styles
#   #       lwd = 2)
#   lines(seq(0, 44, by = 1), values, 
#         type = "l", col = colors[i], lwd = 2)
# }
# 
# # Add legend
# legend("bottomright", legend = c("Homogeneous Stand 1", 
#                                  "Homogeneous Stand 2",
#                                  "Homogeneous Stand 3"), col = colors, 
#        # lty = 1, 
#        lwd = 2, cex = 0.8)  # Ensure consistent line type in legend
# 
# dev.off()
# 






# 
# # Plot correlation curve
# plot(seq(0, 44, by = 1),
#      correlation_values, 
#      type = "o", 
#      xlab = "PAI Profile Depth", 
#      ylab = "Correlation Value",
#      main = "Correlation Values of PAI Profiles of Varying Depth with Sentinel-2 LAI")
# axis(side = 1, at = seq(0, 44, by = 2.5), labels = seq(0, 44, by = 2.5))
# dev.off()
# plot(seq(0, 44, by = 1),
#      correlation_values, 
#      type = "o", 
#      xlab = "PAI Profile Depth", 
#      ylab = "Correlation Value",
#      main = "Correlation Values of PAI Profiles of Varying Depth with Sentinel-2 LAI")
# axis(side = 1, at = seq(0, 44, by = 2.5), labels = seq(0, 44, by = 2.5))