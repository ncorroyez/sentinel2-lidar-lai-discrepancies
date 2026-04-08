# ---
# title: "4.prepare_heterogeneity_depth_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-03"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../libraries/functions_plots.R")
source("../libraries/functions_create_masks.R")

# Pre-processing Parameters & Directories
data_dir <- '../../01_DATA'
results_dir <- '../../03_RESULTS'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
# sites <- c("Mormal", "Blois")

choice_equal_deciles <- "deciles" #equal_intervals #deciles
# choice_equal_deciles <- c("deciles", "equal_intervals")

for (site in sites){
  cat("Processing site:", site, "\n")
  # Setup
  results_path <- file.path(results_dir, site)
  # Which mask
  composition_masks <- list(
    "Not_Masked"
    # "Full_Composition"
    # "Deciduous_Flex",
    # "Deciduous_Only",
    # "Coniferous_Flex",
    # "Coniferous_Only",
    # "Beech_Only",
    # "Oak_Only"
    )
  for (composition_mask in composition_masks){
    cat("Processing mask:", composition_mask, "\n")
    if ((composition_mask == "Oak_Only") && (site == "Aigoual")
        || (composition_mask == "Beech_Only") && (site == "Blois")){
      cat(paste(composition_mask, "is not present in", site, "\n"))
    } 
    else {
      metrics_dir <- file.path(results_path, "Metrics")
      composition_mask_dir <- file.path(metrics_dir, composition_mask)
      # cv_lad_dir <- file.path(composition_mask_dir, 
      #                         "PAD_Profiles_updated/CV_LAD")
      # cv_lad_files <- list.files(cv_lad_dir, pattern = "\\.tif$", full.names = T)
      # numbers <- as.numeric(sub(".*/CV_LAD_(\\d+\\.\\d+)_.*", "\\1", cv_lad_files))
      # sorted_indices <- order(numbers, decreasing = TRUE)
      # sorted_cv_lad_files <- cv_lad_files[sorted_indices]
      
      # Open heterogeneity metrics 
      mean <- terra::rast(file.path(composition_mask_dir, "mean_res_10_m.tif"))
      max <- terra::rast(file.path(composition_mask_dir, "max_res_10_m.tif"))
      std <- terra::rast(file.path(composition_mask_dir, "std_res_10_m.tif"))
      cv <- terra::rast(file.path(composition_mask_dir, "cv_res_10_m.tif"))
      # cvlad <- terra::rast(file.path(composition_mask_dir, 
      #                                "cvlad_res_10_m.tif"))
      # cvladnorm <- terra::rast(file.path(composition_mask_dir, 
      #                                    "cvladnorm_res_10_m.tif"))
      rumple <- terra::rast(file.path(composition_mask_dir, 
                                      "rumple_res_10_m.tif"))
      lcv <- terra::rast(file.path(composition_mask_dir, "lcv_res_10_m.tif"))
      lskew <- terra::rast(file.path(composition_mask_dir, 
                                     "lskew_res_10_m.tif"))
      vci <- terra::rast(file.path(composition_mask_dir, "vci_res_10_m.tif"))
      cv_lad <- terra::rast(file.path(composition_mask_dir, "cv_lad_res_10_m.tif"))
      vdr <- terra::rast(file.path(composition_mask_dir, 
                                   "vdr_res_10_m.tif"))
      dtm <- terra::rast(file.path(composition_mask_dir, 
                                   "dtm.tif"))
      gap_fraction <- terra::rast(file.path(composition_mask_dir, 
                                   "gap_fraction_res_10_m.tif"))
      shade <- gap_fraction <- terra::rast(file.path(composition_mask_dir, 
                                                     "shade_res_10_m.tif"))
      slope <- gap_fraction <- terra::rast(file.path(composition_mask_dir, 
                                                     "slope_res_10_m.tif"))
      # shade <- terra::rast(file.path(composition_mask_dir, 
      #                                 "shade_res_10_m.tif"))
      heterogeneity_metrics <- list(
                                    # mean, 
                                    # max, 
                                    # std, 
                                    # cv, 
                                    # cv_lad,
                                    # dtm,
                                    # gap_fraction,
                                    # shade,
                                    # slope
                                    # cvladnorm,
                                    # rumple,
                                    # lcv,
                                    # lskew,
                                    # vci,
                                    # vdr
                                    shade
      )
      
      # heterogeneity_metrics <- list(cvlad, cvladnorm)
      
      for (heterogeneity_metric_raster in heterogeneity_metrics){
        if (identical(heterogeneity_metric_raster, std)) {
          heterogeneity_metric_name <- "std"
        }
        else if (identical(heterogeneity_metric_raster, cv)) {
          heterogeneity_metric_name <- "cv"
        }
        else if (identical(heterogeneity_metric_raster, cv_lad)) {
          heterogeneity_metric_name <- "cv_lad"
        }
        else if (identical(heterogeneity_metric_raster, dtm)) {
          heterogeneity_metric_name <- "dtm"
        }
        else if (identical(heterogeneity_metric_raster, slope)) {
          heterogeneity_metric_name <- "slope"
        }
        else if (identical(heterogeneity_metric_raster, gap_fraction)) {
          heterogeneity_metric_name <- "gap_fraction"
        }
        else if (identical(heterogeneity_metric_raster, shade)) {
          heterogeneity_metric_name <- "shade"
        }
        else if (identical(heterogeneity_metric_raster, mean)) {
          heterogeneity_metric_name <- "mean"
        }
        else if (identical(heterogeneity_metric_raster, rumple)) {
          heterogeneity_metric_name <- "rumple"
        }
        else if (identical(heterogeneity_metric_raster, max)) {
          heterogeneity_metric_name <- "max"
        }
        else if (identical(heterogeneity_metric_raster, lcv)) {
          heterogeneity_metric_name <- "lcv"
        }
        else if (identical(heterogeneity_metric_raster, lskew)) {
          heterogeneity_metric_name <- "lskew"
        }
        else if (identical(heterogeneity_metric_raster, vci)) {
          heterogeneity_metric_name <- "vci"
        }
        else if (identical(heterogeneity_metric_raster, vdr)) {
          heterogeneity_metric_name <- "vdr"
        }
        else if (identical(heterogeneity_metric_raster, shade)) {
          heterogeneity_metric_name <- "shade"
        }
        else{
          stop("Error is not Standard Deviation or Coefficient of Variation
           or Mean or Max or Lcv or Lskew or Rumple or Shade")
        }
        cat("Processing metric:", heterogeneity_metric_name, "\n")
        
        masks_dir <- file.path(results_path, 
                               "LiDAR/Heterogeneity_Masks")
        if (!dir.exists(masks_dir)) {
          dir.create(masks_dir, showWarnings = FALSE, recursive = TRUE)
        }
        quantiles_dir <- file.path(masks_dir,
                                   "Quantiles",
                                   composition_mask,
                                   heterogeneity_metric_name)
        if (!dir.exists(quantiles_dir)) {
          dir.create(quantiles_dir, showWarnings = FALSE, recursive = TRUE)
        }
        
        # save_x_y_plot(xvar = mean,
        #               yvar = std,
        #               dirname = masks_dir,
        #               filename = "std_vs_mean_heights.png",
        #               xlab = "Mean Heights",
        #               ylab = "Standard Deviation",
        #               title = "Standard Deviation vs Mean Heights")
        # save_x_y_plot(xvar = mean,
        #               yvar = cv,
        #               dirname = masks_dir,
        #               filename = "cv_vs_mean_heights.png",
        #               xlab = "Mean Heights",
        #               ylab = "Coefficient of Variation",
        #               title = "Coefficient of Variation vs Mean Heights")
        # save_x_y_plot(xvar = std,
        #               yvar = cv,
        #               dirname = masks_dir,
        #               filename = "cv_vs_std.png",
        #               xlab = "Standard Deviation",
        #               ylab = "Coefficient of Variation",
        #               title = "Coefficient of Variation vs Standard Deviation")
        # save_x_y_plot(xvar = std,
        #               yvar = rumple,
        #               dirname = masks_dir,
        #               filename = "rumple_vs_std.png",
        #               xlab = "Standard Deviation",
        #               ylab = "Rumple Index",
        #               title = "Rumple Index vs Standard Deviation")
        # save_x_y_plot(xvar = mean,
        #               yvar = rumple,
        #               dirname = masks_dir,
        #               filename = "rumple_vs_mean.png",
        #               xlab = "Mean",
        #               ylab = "Rumple Index",
        #               title = "Rumple Index vs Mean")
        # save_x_y_plot(xvar = cvladnorm,
        #               yvar = std,
        #               dirname = masks_dir,
        #               filename = "std_vs_cvladnorm.png",
        #               xlab = "Standard Deviation",
        #               ylab = "Cvlad",
        #               title = "Std vs CvLAD")
        # plot_density_scatterplot(values(cvladnorm), 
        #                          values(std),
        #                          "cvladnorm",
        #                          "std")
        for (choice_equal_decile in choice_equal_deciles){
          cat("Processing quantile type:", choice_equal_decile, "\n")
          for (inc in c(0.1, 0.2, 1/3)){
            cat("Processing inc:", inc, "\n")
            create_heterogeneity_quantiles(heterogeneity_metric_raster,
                                           heterogeneity_metric_name,
                                           site,
                                           quantiles_dir,
                                           composition_mask_dir,
                                           choice_equal_decile,
                                           inc = inc,
                                           resolution = 10)
            # for (cv_lad_file in sorted_cv_lad_files){
            #   cat("Processing cv_lad_file:", cv_lad_file, "\n")
            #   quantiles_dir <- file.path(masks_dir,
            #                              "Quantiles",
            #                              composition_mask,
            #                              "CV_LAD")
            #   cv_lad_raster <- terra::rast(cv_lad_file)
            #   cv_lad_name <- basename(cv_lad_file)
            #   min_depth <- as.numeric(sub(".*CV_LAD_([0-9.]+)_40.tif.*", "\\1",
            #                               cv_lad_name)) 
            #   create_cv_lad_quantiles(cv_lad_raster,
            #                           cv_lad_name,
            #                           min_depth,
            #                           site,
            #                           quantiles_dir,
            #                           composition_mask_dir,
            #                           choice_equal_decile,
            #                           inc = 1/3,
            #                           resolution = 10)
            # }
          }
        }
      }
    }
  }
}




# Heterogeneity Mask

# 
# plot(mean_heights, standard_deviation_proj, main = "Standard Deviation vs Mean Heights")
# plot(mean_heights, coeff_variation, main = "Coefficient of Variation vs Mean Heights")
# plot(standard_deviation_proj, coeff_variation, main = "Coefficient of Variation vs Standard Deviation")
# 
# 
# 
# # Calculate deciles on heterogeneity raster
# errors <- list(standard_deviation_proj, coeff_variation, mean_heights)
# # errors <- list(mean_heights)
# 
# 
# # Test Deciles
# correlation_values <- c()
# for (low_value in seq(range[1], range[2], by = inc)) {
#   low_value <- round(low_value, digits = 2)
#   high_value <- round(low_value + inc, digits = 2)
#   lidar_intervals <- terra::rast(file.path(class_dir, 
#                                            paste0(sprintf("%s_pai_masked_res_%s_m_", 
#                                                           err, 
#                                                           resolution),
#                                                   low_value,  "_", 
#                                                   high_value, ".tif")))
#   lidar_values <- values(lidar_intervals)
#   
#   s2_intervals <- terra::rast(file.path(class_dir, 
#                                         paste0(sprintf("%s_lai_s2_masked_res_%s_m_", 
#                                                        err, 
#                                                        resolution),
#                                                low_value, "_", 
#                                                high_value, ".tif")))
#   s2_values <- values(s2_intervals)
#   
#   correlation_value <- correlation_test_function(lidar_values, s2_values, method = "pearson")
#   correlation_values <- c(correlation_values, correlation_value)
#   
#   plot_density_scatterplot(var_x = s2_values,
#                            var_y = lidar_values,
#                            xlab = "LAI PROSAIL",
#                            ylab = "LAI LiDAR",
#                            dirname = class_dir,
#                            filename = paste0(sprintf("%s_scatter_lai_lidar_prosail_res_%s_m_", 
#                                                      err, resolution),
#                                              low_value, "_", 
#                                              high_value, ".tif"),
#                            xlimits = c(0, 10),
#                            ylimits = c(0, 10)
#   )
# }
# first_range <- sapply(quantiles_list, function(x) x[1])
# quantile_labels <- c(first_range, max(values(error), na.rm=T))
# quantile_labels <- round(quantile_labels, 1)
# png(filename = file.path(class_dir, 
#                          sprintf("z_%s_correlation_values_for_%s_res_%s_m.png", 
#                                  err, choice, resolution)),
#     width = 1920, height = 1080, 
#     units = "px", res = 200
# )
# plot(1, type = "l", 
#      xlim = c(0, 1), ylim = c(0, 0.65),
#      xlab = choice_equal_deciles_for_plot, 
#      ylab = "Correlation Value", 
#      # main = "Correlation Values between LiDAR LAI and Sentinel-2 LAI for 10 Mean Heights Classes",
#      # main = paste("Correlation Values between LiDAR LAI and\n",
#      #              "Sentinel-2 LAI for 10 Mean Heights Classes"),
#      lwd = 2,
#      xaxt = "n",
#      cex.main = 1.2,   # Increase title font size
#      cex.lab = 1.2,    # Normal font size for axis labels
#      cex.axis = 1.2
# )
# title(main = paste("Pearson Correlation between LiDAR LAI and Sentinel-2 LAI\n",
#                    "for 10 Height Classes (Defined as the CHM Mean Heights\n",
#                    "(1m resolution)) over Sentinel-2 pixels (10m resolution)"))
# lines(midpoints, correlation_values, type = "l")
# axis(side = 1, at = seq(0, 1, by = 0.1), labels = quantile_labels)
# points(midpoints, correlation_values, pch = 16, col = "red")
# # axis(side = 1, at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
# dev.off()
