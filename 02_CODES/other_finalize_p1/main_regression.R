# ---
# title: "main_regression.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-23"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Libraries
library(terra)
library(caret)
library(randomForest)
source("libraries/functions_plots.R")

# Function to calculate RMSE
rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

# Setup
results_dir <- "../03_RESULTS"
forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Blois")
set.seed(123)

# Initial check
for (site in sites){
  metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
  lidar_lai <- terra::rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  lidar_lai_values <- values(lidar_lai)
  
  s2_lai <- terra::rast(file.path(metrics_dir, "s2lai_res_10_m.tif"))
  s2_lai_values <- values(s2_lai)
  
  # Filter out NA values for both LiDAR and S2
  valid_idx <- !is.na(lidar_lai_values) & !is.na(s2_lai_values)
  lidar_lai_values <- lidar_lai_values[valid_idx]
  s2_lai_values <- s2_lai_values[valid_idx]
  
  # Perform linear regression: LiDAR as predictor and S2 as response variable
  regression_model <- lm(s2_lai_values ~ lidar_lai_values)
  
  # Extract intercept and slope from the model
  intercept <- coef(regression_model)[1]
  slope <- coef(regression_model)[2]
  
  # Print results for each site
  cat("Linear Regression\n")
  cat("Site:", site, "\n")
  cat("Intercept:", intercept, "\n")
  cat("Slope:", slope, "\n")
  cat("-------------------------\n")
}

# Scaled
for (site in sites) {
  # Load the LiDAR LAI data (response variable)
  metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
  
  lidar_lai <- terra::rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  lidar_lai_values <- values(lidar_lai)
  
  # Load Sentinel-2 LAI data (predictor)
  s2_lai <- terra::rast(file.path(metrics_dir, "s2lai_res_10_m.tif"))
  s2_lai_values <- values(s2_lai)
  
  # Load additional metrics as predictors
  lcv <- terra::rast(file.path(metrics_dir, "lcv_res_10_m.tif"))
  mean <- terra::rast(file.path(metrics_dir, "mean_res_10_m.tif"))
  max <- terra::rast(file.path(metrics_dir, "max_res_10_m.tif"))
  lskew <- terra::rast(file.path(metrics_dir, "lskew_res_10_m.tif"))
  shade <- terra::rast(file.path(metrics_dir, "shade_res_10_m.tif"))
  slope <- terra::rast(file.path(metrics_dir, "slope_res_10_m.tif"))
  aspect <- terra::rast(file.path(metrics_dir, "aspect_res_10_m.tif"))
  aspect_cos <- terra::rast(file.path(metrics_dir, "aspect_cos_res_10_m.tif"))
  aspect_sin <- terra::rast(file.path(metrics_dir, "aspect_sin_res_10_m.tif"))
  
  lcv_values <- values(lcv)
  mean_values <- values(mean)
  max_values <- values(max)
  lskew_values <- values(lskew)
  shade_values <- values(shade)
  slope_values <- values(slope)
  aspect_values <- values(aspect)
  aspect_cos_values <- values(aspect_cos)
  aspect_sin_values <- values(aspect_sin)
  
  # Filter out NA values for all variables (lidar, S2, and additional metrics)
  valid_idx <- !is.na(lidar_lai_values) & !is.na(s2_lai_values) & 
    !is.na(lcv_values) & !is.na(mean_values) & !is.na(max_values) &
    !is.na(lskew_values) & !is.na(shade_values) & !is.na(slope_values) &
    !is.na(aspect_values) & !is.na(aspect_cos_values) & !is.na(aspect_sin_values)
  
  lidar_lai_values <- lidar_lai_values[valid_idx]
  s2_lai_values <- s2_lai_values[valid_idx]
  lcv_values <- lcv_values[valid_idx]
  mean_values <- mean_values[valid_idx]
  max_values <- max_values[valid_idx]
  lskew_values <- lskew_values[valid_idx]
  shade_values <- shade_values[valid_idx]
  slope_values <- slope_values[valid_idx]
  aspect_values <- aspect_values[valid_idx]
  aspect_cos_values <- aspect_cos_values[valid_idx]
  aspect_sin_values <- aspect_sin_values[valid_idx]
  
  # Combine all predictor variables into a data frame
  predictors <- data.frame(
                           # s2_lai_values,
                           lcv_values,
                           mean_values,
                           # max_values,
                           lskew_values,
                           # shade_values,
                           slope_values,
                           aspect_values
                           # aspect_cos_values,
                           # aspect_sin_values
                           )
  
  train_correlation_matrix <- cor(predictors)
  corrplot::corrplot(train_correlation_matrix,
                     method = "number",
                     type = "upper")
  
  # Split the data into training (80%) and testing sets (20%)
  train_index <- createDataPartition(lidar_lai_values, p = 0.8, list = FALSE)
  
  train_predictors <- predictors[train_index, ]
  test_predictors <- predictors[-train_index, ]
  
  train_lidar_lai <- lidar_lai_values[train_index] - s2_lai_values[train_index]
  test_lidar_lai <- lidar_lai_values[-train_index] - s2_lai_values[-train_index]
  
  # Manually scale the training data
  train_means <- apply(train_predictors, 2, mean)
  train_sds <- apply(train_predictors, 2, sd)
  
  scaled_train_predictors <- sweep(sweep(train_predictors, 2, train_means, "-"), 2, train_sds, "/")
  
  # Scale the response variable (LiDAR LAI) manually
  train_lidar_mean <- mean(train_lidar_lai)
  train_lidar_sd <- sd(train_lidar_lai)
  scaled_train_lidar_lai <- (train_lidar_lai - train_lidar_mean) / train_lidar_sd
  
  # rf <- randomForest(scaled_train_lidar_lai ~ ., data = scaled_train_predictors,
  #                    ntree = 30, importance = T, do.trace = T)
  # 
  # importance_df <- as.data.frame(importance(rf))
  # importance_df$Variable <- rownames(importance_df)
  
  # Fit the multiple linear regression model on the scaled training data
  multiple_regression_model <- lm(scaled_train_lidar_lai ~ ., data = scaled_train_predictors)
  
  # Extract regression coefficients (intercept and slopes)
  regression_coefficients <- coef(multiple_regression_model)
  
  # Make predictions on the scaled training data
  scaled_predicted_train_lidar_lai <- predict(multiple_regression_model, newdata = scaled_train_predictors)
  predicted_train_lidar_lai <- scaled_predicted_train_lidar_lai * train_lidar_sd + train_lidar_mean
  
  # Evaluate the model on the training set
  r_squared_train <- summary(multiple_regression_model)$r.squared
  model_rmse_train <- rmse(train_lidar_lai, predicted_train_lidar_lai)
  
  # Manually scale the test data using training mean and sd
  scaled_test_predictors <- sweep(sweep(test_predictors, 2, train_means, "-"), 2, train_sds, "/")
  
  # Predict on the test set
  scaled_predicted_test_lidar_lai <- predict(multiple_regression_model, newdata = scaled_test_predictors)
  predicted_test_lidar_lai <- scaled_predicted_test_lidar_lai * train_lidar_sd + train_lidar_mean
  
  # Evaluate the model on the testing set
  r_squared_test <- cor(test_lidar_lai, predicted_test_lidar_lai)^2
  model_rmse_test <- rmse(test_lidar_lai, predicted_test_lidar_lai)
  
  plot_density_scatterplot(test_lidar_lai,
                           predicted_test_lidar_lai,
                           "Observed LiDAR LAI (Test)",
                           "Predicted LiDAR LAI (Test)",
                           paste("Observed vs Predicted LiDAR LAI for", site, "(Test Data)"))
  
  # Print results for each site
  cat("Multiple Linear Regression\n")
  cat("Site:", site, "\n")
  cat("Intercept:", regression_coefficients[1], "\n")
  for (i in 2:length(regression_coefficients)) {
    cat("Coefficient for predictor", names(regression_coefficients)[i], ":", regression_coefficients[i], "\n")
  }
  cat("R-squared for training set:", r_squared_train, "\n")
  cat("RMSE for training set:", model_rmse_train, "\n")
  cat("R-squared for testing set:", r_squared_test, "\n")
  cat("RMSE for testing set:", model_rmse_test, "\n")
  cat("-------------------------\n")
}
