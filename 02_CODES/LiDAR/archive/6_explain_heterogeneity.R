# ------------------------------ Code Info -------------------------------------
# title: "3.explain_heterogeneity.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, 
# CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-05"

# -------------- (Optional) Clear the environment and free memory --------------
rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ------------------------------ Libraries -------------------------------------
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rgl")
library("randomForest")
library("caret")
library("e1071")
library("lme4")
library("glmnet")
library("leaps")
library("pdp")
library("ggplot2")
library("DALEX")
library("tuneRanger")
library("mlr")
library("OpenML")
library("pls")

# ------ Define working dir as the directory where the script is located -------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# --------------------------- Import useful functions --------------------------
source("../libraries/functions_lidar.R")
source("../libraries/functions_explain_heterogeneity.R")
source("../libraries/functions_plots.R")

# ----------------------------- Plotting Parameters ----------------------------
default_par <- par(no.readonly = TRUE)
par(mar=c(5, 5, 4, 2) + 0.1) # Adjust margins as per your preference
par(oma=c(0, 0, 0, 0)) # Adjust outer margins as per your preference
par(cex=1.6, cex.axis=1.6, cex.names=1.6)
# par(default_par)
# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
lidar_dir <- "LiDAR"

# Mormal
mormal_results_path <- file.path("../../03_RESULTS", "Mormal")
mormal_lidr_dir <- file.path(mormal_results_path, lidar_dir, "PAI/lidR")
mormal_masks_dir <- file.path(mormal_results_path, 
                              lidar_dir, 
                              "Heterogeneity_Masks")
# dir.create(path = mormal_lidr_dir, showWarnings = F, recursive = T)
# dir.create(path = mormal_masks_dir, showWarnings = F, recursive = T)

# Blois
blois_results_path <- file.path("../../03_RESULTS", "Blois")
blois_lidr_dir <- file.path(blois_results_path, lidar_dir, "PAI/lidR")
blois_masks_dir <- file.path(blois_results_path, 
                             lidar_dir, 
                             "Heterogeneity_Masks")
# dir.create(path = blois_lidr_dir, showWarnings = F, recursive = T)
# dir.create(path = blois_masks_dir, showWarnings = F, recursive = T)

# Aigoual
aigoual_results_path <- file.path("../../03_RESULTS", "Aigoual")
aigoual_lidr_dir <- file.path(aigoual_results_path, lidar_dir, "PAI/lidR")
aigoual_masks_dir <- file.path(aigoual_results_path, 
                               lidar_dir, 
                               "Heterogeneity_Masks")
# dir.create(path = aigoual_lidr_dir, showWarnings = F, recursive = T)
# dir.create(path = aigoual_masks_dir, showWarnings = F, recursive = T)

# Results Reproducibility
set.seed(0) 
# total_sample_size <- 205000
# training_sample_size <- 5000
# test_sample_size <- 200000

# Masks
full_composition_mask <- terra::rast(file.path(mormal_masks_dir,
                                             "artifacts_full_composition_low_vegetation_majority_90_p_res_10_m.envi"))
deciduous_only_mask <- terra::rast(file.path(mormal_masks_dir,
                                             "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"))
deciduous_flex_mask <- terra::rast(file.path(mormal_masks_dir,
                                             "artifacts_deciduous_flex_low_vegetation_majority_90_p_res_10_m.envi"))

# ------------------------ Data Preparation: Mormal ----------------------------
mormal_predictors_values <- extract_predictors_values(mormal_masks_dir, 
                                                      mormal_lidr_dir)

# Initial Correlation
plot_density_scatterplot(mormal_predictors_values$lai_s2,
                         mormal_predictors_values$lai_lidar,
                         "Sentinel-2 LAI",
                         "LiDAR LAI",
                         title = "Mormal Sentinel-2 LAI vs LiDAR LAI Values",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

# Train Test Creation
mormal_train_test_datasets <- create_train_test_datasets(
  mormal_predictors_values, 
  total_sample_size = length(mormal_predictors_values[[1]]),
  training_sample_size = 10000,
  test_sample_size = length(mormal_predictors_values[[1]]) - 10000)
mormal_training_data <- mormal_train_test_datasets$training_data
mormal_test_data <- mormal_train_test_datasets$test_data
mormal_lai_lidar_train_mean <- mormal_train_test_datasets$lai_lidar_train_mean
mormal_lai_lidar_train_sd <- mormal_train_test_datasets$lai_lidar_train_sd

# ------------------------ Data Preparation: Blois -----------------------------
blois_predictors_values <- extract_predictors_values(blois_masks_dir,
                                                     blois_lidr_dir)

# Initial Correlation
plot_density_scatterplot(blois_predictors_values$lai_s2,
                         blois_predictors_values$lai_lidar,
                         "Sentinel-2 LAI",
                         "LiDAR LAI",
                         title = "Blois Sentinel-2 LAI vs LiDAR LAI Values",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

blois_train_test_datasets <- create_train_test_datasets(
  blois_predictors_values, 
  total_sample_size = length(blois_predictors_values[[1]]),
  training_sample_size = 10000,
  test_sample_size = length(blois_predictors_values[[1]]) - 10000)
blois_training_data <- blois_train_test_datasets$training_data
blois_test_data <- blois_train_test_datasets$test_data
blois_lai_lidar_train_mean <- blois_train_test_datasets$lai_lidar_train_mean
blois_lai_lidar_train_sd <- blois_train_test_datasets$lai_lidar_train_sd

# ---------------------- Data Preparation: Aigoual  ----------------------------
aigoual_predictors_values <- extract_predictors_values(aigoual_masks_dir,
                                                       aigoual_lidr_dir)

# Initial Correlation
plot_density_scatterplot(aigoual_predictors_values$lai_s2,
                         aigoual_predictors_values$lai_lidar,
                         "Sentinel-2 LAI",
                         "LiDAR LAI",
                         title = "Aigoual Sentinel-2 LAI vs LiDAR LAI Values",
                         # xlimits = c(0,10),
                         # ylimits = c(0,10)
)
cor(aigoual_predictors_values$lai_s2, aigoual_predictors_values$lai_lidar)

aigoual_train_test_datasets <- create_train_test_datasets(
  aigoual_predictors_values, 
  total_sample_size = length(aigoual_predictors_values[[1]]),
  training_sample_size = 10000,
  test_sample_size = length(aigoual_predictors_values[[1]]) - 10000)
aigoual_training_data <- aigoual_train_test_datasets$training_data
aigoual_test_data <- aigoual_train_test_datasets$test_data
aigoual_lai_lidar_train_mean <- aigoual_train_test_datasets$lai_lidar_train_mean
aigoual_lai_lidar_train_sd <- aigoual_train_test_datasets$lai_lidar_train_sd

# Data Preparation: Mormal + Blois
mix_training_data <- rbind(mormal_training_data, blois_training_data, aigoual_training_data)
mix_training_data <- mix_training_data[sample(nrow(mix_training_data)), ]

# -------------------------------- Simple RF -----------------------------------

# Set hyperparameters
ntree <- 50
mtry <- 3

# Create formula
formula <- lai_lidar ~ 
  lai_s2 +
  mean_h +
  cv +
  rumple +
  vci +
  lcv +
  lskew

rf_model <- randomForest(formula,
                         data = mormal_training_data, # mormal_ blois_ aigoual_ mix_
                         ntree = ntree,
                         mtry = mtry,
                         importance = TRUE,
                         do.trace = TRUE
)

# Print the number of trees used
num_trees <- length(rf_model$mse)
print(paste("Number of trees:", num_trees))
print(summary(rf_model))

importance <- varImp(rf_model)
normalized_importance <- importance$Overall / sum(importance$Overall) * 100
names(normalized_importance) <- rownames(importance)
normalized_importance <- sort(normalized_importance, decreasing = TRUE)
df_importance <- data.frame(variable = names(normalized_importance), 
                            importance = normalized_importance)
barplot(df_importance$importance, 
        names.arg = df_importance$variable, 
        main = "Normalized Importance Score",
        xlab = "Variable", 
        ylab = "Importance",
        cex.lab = 2,
        cex.axis = 2,
        cex.main = 2,
        col = "#e6e6e6"
)
print(normalized_importance)

# Predict using the random forest model
# On Mormal
lai_lidar_mormal <- terra::rast(file.path(mormal_masks_dir,
                                          "lai_lidar_masked_res_10_m.envi")) 
mormal_predictions <- predict(rf_model, newdata = mormal_test_data)
mormal_comparison <- data.frame(Actual = mormal_test_data$lai_lidar,
                                Predicted = mormal_predictions)

# Compute metrics
correlation <- cor(mormal_comparison$Actual, mormal_comparison$Predicted)
rsquared <- caret::R2(mormal_comparison$Actual, mormal_comparison$Predicted)
rmse <- Metrics::rmse(mormal_comparison$Actual, mormal_comparison$Predicted)
bias <- Metrics::bias(mormal_comparison$Actual, mormal_comparison$Predicted)
mae <- Metrics::mae(mormal_comparison$Actual, mormal_comparison$Predicted)

# Print metrics
cat("Mormal Metrics:",
    "\nCorrelation:", correlation, 
    "\nR-squared:", rsquared, 
    "\nRoot Mean Squared Error:", rmse,
    "\nBias:", bias, 
    "\nMAE:", mae)

# On Blois
blois_predictions <- predict(rf_model, newdata = blois_test_data)
blois_comparison <- data.frame(Actual = blois_test_data$lai_lidar,
                               Predicted = blois_predictions)

# Compute metrics
correlation <- cor(blois_comparison$Actual, blois_comparison$Predicted)
rsquared <- caret::R2(blois_comparison$Actual, blois_comparison$Predicted)
rmse <- Metrics::rmse(blois_comparison$Actual, blois_comparison$Predicted)
bias <- Metrics::bias(blois_comparison$Actual, blois_comparison$Predicted)
mae <- Metrics::mae(blois_comparison$Actual, blois_comparison$Predicted)

# Print metrics
cat("Blois Metrics:",
    "\nCorrelation:", correlation, 
    "\nR-squared:", rsquared, 
    "\nRoot Mean Squared Error:", rmse,
    "\nBias:", bias, 
    "\nMAE:", mae)

# On Aigoual
aigoual_predictions <- predict(rf_model, newdata = aigoual_test_data)
aigoual_comparison <- data.frame(Actual = aigoual_test_data$lai_lidar,
                                 Predicted = aigoual_predictions)

# Compute metrics
correlation <- cor(aigoual_comparison$Actual, aigoual_comparison$Predicted)
rsquared <- caret::R2(aigoual_comparison$Actual, aigoual_comparison$Predicted)
rmse <- Metrics::rmse(aigoual_comparison$Actual, aigoual_comparison$Predicted)
bias <- Metrics::bias(aigoual_comparison$Actual, aigoual_comparison$Predicted)
mae <- Metrics::mae(aigoual_comparison$Actual, aigoual_comparison$Predicted)

# Print metrics
cat("Aigoual Metrics:",
    "\nCorrelation:", correlation, 
    "\nR-squared:", rsquared, 
    "\nRoot Mean Squared Error:", rmse,
    "\nBias:", bias, 
    "\nMAE:", mae)

# Unstandardize actual and predictions using mean and standard deviation 
# of original lai_lidar variable

# Mormal
mormal_visualization <- data.frame(Actual = unstandardize_variable(
  mormal_test_data$lai_lidar,
  mormal_lai_lidar_train_mean, 
  mormal_lai_lidar_train_sd),
  Predicted = unstandardize_variable(
    mormal_predictions,
    mormal_lai_lidar_train_mean, 
    mormal_lai_lidar_train_sd))
plot_density_scatterplot(mormal_visualization$Actual,
                         mormal_visualization$Predicted,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10)
)

# Blois
blois_visualization <- data.frame(Actual = unstandardize_variable(
  blois_test_data$lai_lidar,
  blois_lai_lidar_train_mean, 
  blois_lai_lidar_train_sd),
  Predicted = unstandardize_variable(
    blois_predictions,
    blois_lai_lidar_train_mean, 
    blois_lai_lidar_train_sd))

plot_density_scatterplot(blois_visualization$Actual,
                         blois_visualization$Predicted,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10)
)

# Aigoual
aigoual_visualization <- data.frame(Actual = unstandardize_variable(
  aigoual_test_data$lai_lidar,
  aigoual_lai_lidar_train_mean, 
  aigoual_lai_lidar_train_sd),
  Predicted = unstandardize_variable(
    aigoual_predictions,
    aigoual_lai_lidar_train_mean, 
    aigoual_lai_lidar_train_sd))

plot_density_scatterplot(aigoual_visualization$Actual,
                         aigoual_visualization$Predicted,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         # xlimits = c(0,10),
                         # ylimits = c(0,10)
)

# Visualize predicted raster
# new_raster[non_na_indices] <- predictions
# plot_histogram(list(predictions, predictions_ranger$data$response, test_data$lai_lidar, test_data$lai_s2),
#                var_labs = c("Predicted", "Predicted max", "Actual", "S2"))
# # plot_histogram(list(predictions, predictions_ranger$data$response),
# #                var_labs = c("Predicted", "Predicted max"))
# 
# # 1. Calculate RMSE of the trained model
# train_predictions <- predict(rf_model, newdata = training_data)
# train_rmse <- sqrt(mean((training_data$lai_lidar - train_predictions)^2))
# print(paste("RMSE on training data:", train_rmse))
# 
# # 2. Make predictions on the validation data
# predictions <- predict(rf_model, newdata = test_data)
# 
# # 3. Calculate RMSE on the predicted data
# validation_rmse <- sqrt(mean((test_data$lai_lidar - predictions)^2))
# print(paste("RMSE on validation data:", validation_rmse))


# ----------------------------------- ES ---------------------------------------
# Set hyperparameters
ntree <- 100  # Maximum number of trees
mtry <- 3
early_stopping_rounds <- 10  # Number of rounds to wait for improvement

# Create formula
formula <- lai_lidar ~ 
  lai_s2 +
  mean_h +
  cv +
  rumple +
  vci +
  lcv +
  lskew

# Fit the random forest model with monitoring for early stopping
rf_model <- randomForest(formula,
                         data = mormal_training_data,  # mormal_ blois_
                         ntree = ntree,
                         mtry = mtry,
                         importance = TRUE,
                         do.trace = TRUE,
                         keep.inbag = TRUE)

# Extract OOB error for each tree
oob_error <- rf_model$mse

# Implement early stopping
best_error <- min(oob_error)
best_iter <- which.min(oob_error)
rounds_without_improvement <- 0

for (i in seq_along(oob_error)) {
  if (oob_error[i] < best_error) {
    best_error <- oob_error[i]
    best_iter <- i
    rounds_without_improvement <- 0
  } else {
    rounds_without_improvement <- rounds_without_improvement + 1
  }
  if (rounds_without_improvement >= early_stopping_rounds) {
    cat("Early stopping at iteration:", i, "\n")
    break
  }
}

# Refit the model with the optimal number of trees
rf_model_optimal <- randomForest(formula,
                                 data = mormal_training_data,  # mormal_ blois_
                                 ntree = best_iter,
                                 mtry = mtry,
                                 importance = TRUE,
                                 do.trace = TRUE)

# Print the number of trees used
num_trees <- length(rf_model_optimal$mse)
print(paste("Number of trees:", num_trees))
print(summary(rf_model_optimal))

importance <- varImp(rf_model_optimal)
normalized_importance <- importance$Overall / sum(importance$Overall) * 100
names(normalized_importance) <- rownames(importance)
normalized_importance <- sort(normalized_importance, decreasing = TRUE)
df_importance <- data.frame(variable = names(normalized_importance), 
                            importance = normalized_importance)
barplot(df_importance$importance, 
        names.arg = df_importance$variable, 
        main = "Normalized Importance Score",
        xlab = "Variable", 
        ylab = "Importance",
        cex.lab = 2,
        cex.axis = 2,
        cex.main = 2,
        col = "#e6e6e6"
)
print(normalized_importance)

# Predict using the random forest model
# On Mormal
mormal_predictions <- predict(rf_model_optimal, newdata = mormal_test_data)
mormal_comparison <- data.frame(Actual = mormal_test_data$lai_lidar,
                                Predicted = mormal_predictions)

# Compute metrics
correlation <- cor(mormal_comparison$Actual, mormal_comparison$Predicted)
rmse <- Metrics::rmse(mormal_comparison$Actual, mormal_comparison$Predicted)
rsquared <- caret::R2(mormal_comparison$Actual, mormal_comparison$Predicted)
bias <- Metrics::bias(mormal_comparison$Actual, mormal_comparison$Predicted)
mae <- Metrics::mae(mormal_comparison$Actual, mormal_comparison$Predicted)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
    "\nMAE:", mae)
# ---------------------------- RF Train Test  ----------------------------------
# Define the range of number of trees to evaluate
num_trees <- seq(1, 1000, by = 1)

predictors = training_data[, -which(names(training_data) == "lai_lidar")]
outcome = training_data$lai_lidar

# Initialize vectors to store RMSE values
train_rmse_values <- numeric(length(num_trees))
val_rmse_values <- numeric(length(num_trees))
bias_values <- numeric(length(num_trees))
variance_values <- numeric(length(num_trees))

# CV
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 5,
                     search = "random")
mtry <- sqrt(ncol(predictors))
tunegrid <- expand.grid(.mtry=c(1:15))
metric <- "RMSE"
# rf_default <- train(x = predictors,
#                     y = outcome, 
#                     method="rf", 
#                     ntree = 15,
#                     metric=metric, 
#                     tuneGrid=tunegrid, 
#                     trControl=ctrl)
# print(rf_default)
# plot(rf_default)

# Train models and store RMSE values
ctrl <- trainControl(method = "LOOCV",
                     number = 3)
# Define random forest model
# rf_model2 <- train(
#   x = predictors,
#   y = outcome,
#   method = "rf",
#   ntree = num_trees[i], 
#   trControl = ctrl,
#   tuneGrid = data.frame(mtry = 4),
#   do.trace = TRUE
# )

rf_model <- randomForest(formula,
                         data = training_data,
                         ntree = 1000,  # or any other number of trees
                         mtry = mtry,
                         importance = TRUE,
                         do.trace = TRUE)

# Initialize an empty list to store subset models
subset_models <- list()

# Iterate through different numbers of trees
for (i in seq_along(num_trees)) {
  print(i)
  # Extract a subset of trees from the forest
  rf_model_sub <- randomForest::combine(rf_model, k = num_trees[i])
  subset_models[[i]] <- rf_model_sub
  
  # Make predictions using the subset of trees
  predictions_train <- predict(rf_model_sub, newdata = training_data)
  predictions_val <- predict(rf_model_sub, newdata = test_data)
  
  # Calculate performance metrics
  train_rmse_values[i] <- RMSE(predictions_train, training_data$lai_lidar)
  val_rmse_values[i] <- RMSE(predictions_val, test_data$lai_lidar)
  bias_values[i] <- Metrics::bias(predictions_train, training_data$lai_lidar)
  variance_values[i] <- var(predictions_train)
}

# Print the performance metrics for each tree
for (i in seq_along(num_trees)) {
  cat("Number of trees:", num_trees[i], "\n")
  cat("Train RMSE:", train_rmse_values[i], "\n")
  cat("Validation RMSE:", val_rmse_values[i], "\n")
  cat("Bias:", bias_values[i], "\n")
  cat("Variance:", variance_values[i], "\n\n")
}

# Plot RMSE values against number of trees
plot(num_trees, train_rmse_values, type = "l", col = "#51a343", 
     ylim = c(0,1),
     xlab = "Number of Trees", ylab = "RMSE",
     main = paste("RMSE vs. Number of Trees\n",
                  "K-fold Repeated Cross Validation\n",
                  "5 folds / 5 repeats"),
     cex.lab = 1.5,    # Normal font size for axis labels
     cex.axis = 1.5,
     cex.main = 1.5, 
     lwd = 2
     # main = "RMSE vs. Number of Trees"
)
lines(num_trees, val_rmse_values, type = "l", col = "#9a4c00", lwd = 2)
abline(v = 15, col = "red", lwd = 2)
legend("topright", 
       legend = c("Train", "Test", "NTree Choice"), 
       col = c("#51a343", "#9a4c00", "red"), 
       lty = c(1, 1, 1),
       lwd = c(2, 2, 2),
       cex = 1.4
)
# title(main = paste("RMSE vs. Number of Trees\n",
#                    "K-fold Repeated Cross Validation\n",
#                    "5 folds / 5 repeats"))
# legend("topright", 
#        legend = c("Train", "Test"), 
#        col = c("#51a343", "#9a4c00"), 
#        lty = 1)


# ------------------------------- Regression -----------------------------------
# Create formula
formula <- lai_lidar ~ 
  lai_s2 +
  mean_h + 
  # max_h +
  # std +
  cv +
  # variance +
  rumple +
  vci +
  lcv +
  lskew
# vdr

# Fit the multivariate regression model
lm_model <- lm(formula, data = mormal_training_data)

# Summary of the model
summary(lm_model)

# Make predictions on the test dataset
test_predictions <- predict(lm_model, newdata = mormal_test_data)

# Calculate RMSE (Root Mean Squared Error) to evaluate model performance
rmse <- Metrics::rmse(mormal_test_data$lai_lidar, test_predictions)
print(paste("RMSE:", rmse))

# Plot actual vs predicted values
plot(mormal_test_data$lai_lidar, test_predictions, 
     xlab = "Actual LAI", ylab = "Predicted LAI",
     main = "Actual vs Predicted LAI")

# Add a diagonal line to represent perfect prediction
abline(0, 1, col = "red")

# Perform stepwise regression
stepwise_model <- step(lm_model, direction = "both")

# Summary of the stepwise model (AIC)
summary(stepwise_model)

# Fit the PLS regression model
pls_model <- plsr(formula, data = mix_training_data, ncomp = 7) 

# Summary of the PLS model
summary(pls_model)

# Make predictions on the test dataset
test_predictions <- predict(pls_model, newdata = mormal_test_data, ncomp=7)

comparison <- data.frame(Actual = mormal_test_data$lai_lidar,
                         Predicted = as.double(unname(test_predictions[, , 1])))

# Compute metrics
correlation <- cor(comparison$Actual, comparison$Predicted)
rmse <- Metrics::rmse(comparison$Actual, comparison$Predicted)
rsquared <- caret::R2(comparison$Actual, comparison$Predicted)
bias <- Metrics::bias(comparison$Actual, comparison$Predicted)
mae <- Metrics::mae(comparison$Actual, comparison$Predicted)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
    "\nMAE:", mae)

visualization <- data.frame(Actual = unstandardize_variable(
  mormal_test_data$lai_lidar,
  mormal_lai_lidar_train_mean, 
  mormal_lai_lidar_train_sd),
  Predicted = unstandardize_variable(
    unname(test_predictions[, , 1]),
    mormal_lai_lidar_train_mean, 
    mormal_lai_lidar_train_sd))
plot_density_scatterplot(visualization$Actual,
                         visualization$Predicted,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10)
)