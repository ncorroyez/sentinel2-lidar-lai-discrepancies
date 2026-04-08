# ---
# title: "main_DRF.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-14"
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
library("gstat")
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
# sites <- "Mormal" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/RF"
forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
set.seed(42) # Ensure reproducibility

load_metrics_for_site <- function(site, results_dir, forest_composition) {
  # Define the metrics directory for the given site
  metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
  
  # Define the names of the metrics
  metric_files <- c("lidarlai_res_10_m.tif", "s2lai_res_10_m.tif", "lcv_res_10_m.tif",
                    "mean_res_10_m.tif", "max_res_10_m.tif", "lskew_res_10_m.tif",
                    "shade_res_10_m.tif", "slope_res_10_m.tif", "aspect_res_10_m.tif",
                    "cv_res_10_m.tif", "rumple_res_10_m.tif", "vci_res_10_m.tif",
                    "vdr_res_10_m.tif")
  
  # Function to load and extract raster values
  load_metric <- function(filename) {
    rast <- terra::rast(file.path(metrics_dir, filename))
    values(rast)
  }
  
  # Load all the metrics in a list
  metrics_values <- lapply(metric_files, load_metric)
  
  # Assign meaningful names to the metrics
  names(metrics_values) <- c("lidar_lai", "s2_lai", "lcv", 
                             "mean", "max", "lskew", 
                             "shade", "slope", "aspect", 
                             "cv", "rumple", "vci", 
                             "vdr")
  
  # Extract coordinates (assuming all rasters have the same structure)
  lidar_lai <- terra::rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  coordinates <- xyFromCell(lidar_lai, 1:ncell(lidar_lai))
  
  # Add longitude and latitude to the list of metrics
  metrics_values$lon <- coordinates[, 1]
  metrics_values$lat <- coordinates[, 2]
  
  # Filter out NA values across all metrics
  valid_idx <- complete.cases(do.call(cbind, metrics_values))
  
  # Filter metrics to retain only valid (non-NA) rows
  valid_metrics <- lapply(metrics_values, function(x) x[valid_idx])
  
  # Combine all predictor variables into a data frame
  predictors <- as.data.frame(valid_metrics)
  predictors$site <- site
  
  return(predictors)
}

# Initialize an empty list to store the data frames
all_predictors <- list()

# Loop through each site and load metrics
for (site in sites) {
  site_data <- load_metrics_for_site(site, results_dir, forest_composition)
  all_predictors[[site]] <- site_data  # Store the data frame in the list
}

# Combine all the data frames into one
combined_predictors <- do.call(rbind, all_predictors)

# Create dummy variables for the 'site' column
combined_predictors$site <- as.factor(combined_predictors$site)  # Convert to factor
dummy_vars <- model.matrix(~ site - 1, data = combined_predictors)  # Create dummy variables

# Rename dummy variable columns to avoid confusion
colnames(dummy_vars) <- gsub("site", "site_", colnames(dummy_vars))

# Combine dummy variables with the original predictors, removing the original 'site' column
combined_predictors <- cbind(combined_predictors[, !colnames(combined_predictors) %in% "site"], as.data.frame(dummy_vars))

str(combined_predictors)

# Sparse
# sparse_combined_predictors <- Matrix(as.matrix(combined_predictors), sparse = TRUE)
# Y <- as.vector(sparse_combined_predictors[, "lidar_lai"])
# X <- sparse_combined_predictors[, !colnames(sparse_combined_predictors) %in% "lidar_lai"]

dt <- data.table(combined_predictors)
density_estimate <- density(combined_predictors[, 1], bw = "nrd0")

# combined_predictors_drf <- combined_predictors[sample(nrow(combined_predictors), 30000), ]
# combined_predictors_drf2 <- combined_predictors[sample(nrow(combined_predictors), 30000), ]
# density_estimate_drf <- density(combined_predictors_drf[, 1], bw = "nrd0")
drf_modele <- drf(dt[, -1],
                  dt[, 1],
                  num.trees = 500,
                  num.features = 5,
                  mtry = ncol(combined_predictors) - 1,
                  splitting.rule = "FourierMMD", # FourierMMD FastMMD
                  bandwidth = density_estimate$bw,
                  num.threads = 12)
# saveRDS(drf_modele, paste0(output_dir, "/drf_model_densitybw.rds"))

# saved_drf <- readRDS(paste0(output_dir, "/drf_model_densitybw.rds"))
# test <- predict(saved_drf, num.threads=1)

w <- predict(saved_drf)$weights
pred.oob <- predict(drf_modele, 
                    functional = "mean",
                    num.threads = 12)$mean
err.oob <- apply((pred.oob-as.data.frame(combined_predictors[, 1]))^2,2, "mean")
varex <- 1 - err.oob/apply(as.data.frame(combined_predictors[, 1]), 2, "var")
VI_MMD <- variable_importance(drf_modele)
VI_MMD <- sort(VI_MMD, decreasing = T)
names(VI_MMD) <- names(combined_predictors_drf[, -1])
barplot(VI_MMD)
a
# scaled_predicted_test_lidar_lai <- predict(multiple_regression_model, newdata = scaled_test_predictors)
# predicted_test_lidar_lai <- scaled_predicted_test_lidar_lai * train_lidar_sd + train_lidar_mean

# Evaluate the model on the testing set
# r_squared_test <- cor(test_lidar_lai, predicted_test_lidar_lai)^2
# model_rmse_test <- rmse(test_lidar_lai, predicted_test_lidar_lai)
# 
# plot_density_scatterplot(test_lidar_lai,
#                          predicted_test_lidar_lai,
#                          "Observed LiDAR LAI (Test)",
#                          "Predicted LiDAR LAI (Test)",
#                          paste("Observed vs Predicted LiDAR LAI for", site, "(Test Data)"))
# 
# # Print results for each site
# cat("Multiple Linear Regression\n")
# cat("Site:", site, "\n")
# cat("Intercept:", regression_coefficients[1], "\n")
# for (i in 2:length(regression_coefficients)) {
#   cat("Coefficient for predictor", names(regression_coefficients)[i], ":", regression_coefficients[i], "\n")
# }
# cat("R-squared for training set:", r_squared_train, "\n")
# cat("RMSE for training set:", model_rmse_train, "\n")
# cat("R-squared for testing set:", r_squared_test, "\n")
# cat("RMSE for testing set:", model_rmse_test, "\n")
# cat("-------------------------\n")







