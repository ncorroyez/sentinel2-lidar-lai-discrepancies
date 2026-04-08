# ------------------------------ Code Info -------------------------------------
# title: "ml_archive.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, 
# CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-17"

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
source("functions_lidar.R")
source("functions_explain_heterogeneity.R")
source("../functions_plots.R")

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
dir.create(path = mormal_lidr_dir, showWarnings = F, recursive = T)
dir.create(path = mormal_masks_dir, showWarnings = F, recursive = T)

# Blois
blois_results_path <- file.path("../../03_RESULTS", "Blois")
blois_lidr_dir <- file.path(blois_results_path, lidar_dir, "PAI/lidR")
blois_masks_dir <- file.path(blois_results_path, 
                             lidar_dir, 
                             "Heterogeneity_Masks")
dir.create(path = blois_lidr_dir, showWarnings = F, recursive = T)
dir.create(path = blois_masks_dir, showWarnings = F, recursive = T)

# Results Reproducibility
set.seed(0) 
total_sample_size <- 205000
training_sample_size <- 5000
test_sample_size <- 200000

# ------------------------ Data Preparation: Mormal ----------------------------
mormal_predictors_values <- extract_predictors_values(mormal_masks_dir, 
                                                      mormal_lidr_dir)
print(length(mormal_predictors_values[[1]]))

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
  training_sample_size = 5000,
  test_sample_size = length(mormal_predictors_values[[1]]) - 5000)
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
  total_sample_size = total_sample_size,
  training_sample_size = training_sample_size,
  test_sample_size = test_sample_size)
blois_training_data <- blois_train_test_datasets$training_data
blois_test_data <- blois_train_test_datasets$test_data
blois_lai_lidar_train_mean <- blois_train_test_datasets$lai_lidar_train_mean
blois_lai_lidar_train_sd <- blois_train_test_datasets$lai_lidar_train_sd

# --------------------------- RF (Cross Validation) ----------------------------

# Initialization
ctrl <- trainControl(method = "LOOCV",
                     number = 5)
# Create formula
formula <- lai_lidar ~ 
  lai_s2 + 
  mean_h + 
  std + 
  cv +
  variance + 
  rumple +
  vci + 
  lcv +
  lskew +
  vdr

# Train the Random Forest model with cross-validation
mtry <- sqrt(ncol(predictors))
tunegrid <- expand.grid(.mtry=c(2:8))
metric <- "RMSE"
predictors = training_data[, -which(names(training_data) == "lai_lidar")]
outcome = training_data$lai_lidar
rf_default <- train(x = predictors,
                    y = outcome, 
                    method="rf", 
                    ntree = 15,
                    metric=metric, 
                    tuneGrid=tunegrid, 
                    trControl=ctrl)
# Check results
print(rf_model)
importance <- varImp(rf_model)
print(importance)
plot(importance)

# ------------------------------- RF tuneRanger --------------------------------
# Define your regression task
regression_task <- makeRegrTask(data = training_data, target = "lai_lidar")

# Estimate run time
estimateTimeTuneRanger(regression_task)

# Tuning parameters
tuned_ranger <- tuneRanger(
  regression_task,
  measure = list(mse),                # Specify the evaluation measure (e.g., RMSE)
  num.trees = 1000,                   # Number of trees
  num.threads = 4,                    # Number of threads for parallelization
  iters = 100,                         # Number of tuning iterations
  iters.warmup = 30                   # Number of warm-up iterations
)

# Display the results
tuned_ranger
# Print the recommended parameter settings
tuned_ranger$results

# Get the tuned Ranger model
tuned_ranger_model <- tuned_ranger$model

# Make predictions on the test data
predictions_ranger <- predict(tuned_ranger_model, newdata = test_data)
comparison_ranger <- data.frame(Actual = test_data$lai_lidar,
                                Predicted = predictions_ranger)

# Compute metrics
correlation <- cor(comparison_ranger$Actual, comparison_ranger$Predicted.response)
rmse <- Metrics::rmse(comparison_ranger$Actual, comparison_ranger$Predicted.response)
rsquared <- R2(comparison_ranger$Actual, comparison_ranger$Predicted.response)
bias <- Metrics::bias(comparison_ranger$Actual, comparison_ranger$Predicted.response)
mae <- Metrics::mae(comparison_ranger$Actual, comparison_ranger$Predicted.response)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
    "\nMAE:", mae)

plot_density_scatterplot(comparison_ranger$Actual,
                         comparison_ranger$Predicted.response,
                         "Actual LiDAR LAI",
                         "Predicted LiDAR LAI",
                         title = "Predicted LiDAR LAI vs Actual LiDAR LAI Values on the Test Set",
                         xlimits = c(0,10),
                         ylimits = c(0,10))

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
                         data = mormal_training_data, # mormal_ blois_
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
mormal_predictions <- predict(rf_model, newdata = mormal_test_data)
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

# On Blois
blois_predictions <- predict(rf_model, newdata = blois_test_data)
blois_comparison <- data.frame(Actual = blois_test_data$lai_lidar,
                               Predicted = blois_predictions)

# Compute metrics
correlation <- cor(blois_comparison$Actual, blois_comparison$Predicted)
rmse <- Metrics::rmse(blois_comparison$Actual, blois_comparison$Predicted)
rsquared <- R2(blois_comparison$Actual, blois_comparison$Predicted)
bias <- Metrics::bias(blois_comparison$Actual, blois_comparison$Predicted)
mae <- Metrics::mae(blois_comparison$Actual, blois_comparison$Predicted)

# Print metrics
cat("Root Mean Squared Error:", rmse, 
    "\nR-squared:", rsquared, 
    "\nBias:", bias, 
    "\nCorrelation:", correlation, 
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

# Visualize predicted raster
new_raster[non_na_indices] <- predictions
plot_histogram(list(predictions, predictions_ranger$data$response, test_data$lai_lidar, test_data$lai_s2),
               var_labs = c("Predicted", "Predicted max", "Actual", "S2"))
# plot_histogram(list(predictions, predictions_ranger$data$response),
#                var_labs = c("Predicted", "Predicted max"))

# 1. Calculate RMSE of the trained model
train_predictions <- predict(rf_model, newdata = training_data)
train_rmse <- sqrt(mean((training_data$lai_lidar - train_predictions)^2))
print(paste("RMSE on training data:", train_rmse))

# 2. Make predictions on the validation data
predictions <- predict(rf_model, newdata = test_data)

# 3. Calculate RMSE on the predicted data
validation_rmse <- sqrt(mean((test_data$lai_lidar - predictions)^2))
print(paste("RMSE on validation data:", validation_rmse))

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

# ------------------------------- Simple RF (ES) -------------------------------

max_trees <- 100
max_iter_without_improvement <- 10
min_delta_error <- 0.0001

# Initialize variables
best_error <- Inf
iter_without_improvement <- 0

# Training loop
for (i in 1:max_trees) {
  # Train random forest model
  rf_model <- randomForest(
    lai_lidar_val ~ ., 
    data = training_data, 
    ntree = i,  
    mtry = sqrt(ncol(training_data) - 1),  
    importance = TRUE,
    do.trace = TRUE
  )
  
  # Compute out-of-bag error rate
  oob_error <- rf_model$mse[length(rf_model$mse)]
  
  # Check if the change in error rate is below threshold
  delta_error <- best_error - oob_error
  
  if (delta_error < min_delta_error) {
    iter_without_improvement <- iter_without_improvement + 1
  } else {
    best_error <- oob_error
    iter_without_improvement <- 0
  }
  
  # Check early stopping condition
  if (iter_without_improvement >= max_iter_without_improvement) {
    print("Early stopping due to lack of improvement.")
    break
  }
}

# Print the number of trees used
num_trees <- length(rf_model$mse)
print(paste("Number of trees:", num_trees))
print(summary(rf_model))

# Normalize importance scores
importance <- varImp(rf_model)
normalized_importance <- importance$Overall / sum(importance$Overall) * 100
names(normalized_importance) <- rownames(importance)
df_importance <- data.frame(variable = names(normalized_importance), 
                            importance = normalized_importance)
barplot(df_importance$importance, 
        names.arg = df_importance$variable, 
        main = "Normalized Importance Scores", 
        xlab = "Variable", 
        ylab = "Importance"
)
print(normalized_importance)

# Extract MSE and Variance for each tree
mse_values <- rf_model$mse
# variance_values <- 100 - (rf_model$rsq * 100)  # Correct calculation for Variance

# Plot MSE and Variance for each tree
tree_numbers <- 1:length(mse_values)

# Plot MSE values according to tree number
plot(tree_numbers, mse_values, type = "l", col = "blue", xlab = "Tree", ylab = "MSE", main = "MSE by Tree")

# -------------------------------- RF Tuning -----------------------------------
# Define grid of hyperparameters
param_grid <- expand.grid(
  ntree = c(10, 15, 20, 30, 100, 200, 500, 1000),
  mtry = c(3, 4, 5, 6, 7, 8),
  nodesize_range <- c(1, 5, 10),
  sampsize_range <- c(50, 100, 150),
  splitrule_range <- c("gini", "extratrees"),
  stringsAsFactors = FALSE
)

# Perform grid search
best_rmse <- Inf
best_model <- NULL
train_control <- trainControl(method = "cv", number = 5)

# Formula
formula <- as.formula(lai_lidar ~ 
                        lai_s2 +
                        mean_h +
                        max_h +
                        std +
                        cv +
                        # variance +
                        rumple +
                        vci +
                        lcv +
                        lskew +
                        vdr)

# Perform grid search
for (i in 1:nrow(param_grid)) {
  # Train random forest model
  rf_model <- train(
    formula,
    data = training_data,
    method = "rf",
    trControl = train_control,
    tuneGrid = data.frame(
      ntree = param_grid$ntree[i],
      mtry = param_grid$mtry[i],
      nodesize = param_grid$nodesize[i],
      sampsize = param_grid$sampsize[i],
      splitrule = param_grid$splitrule[i]
    )
  )
  
  # Evaluate model
  predictions <- predict(rf_model, newdata = test_data)
  rmse <- Metrics::rmse(test_data$lai_lidar, predictions)
  
  # Check if current model is the best
  if (rmse < best_rmse) {
    best_rmse <- rmse
    best_model <- rf_model
  }
}

# Evaluate performance of the best model
cat("Best mtry:", best_mtry, 
    "\nBest ntree:", best_ntree,
    "\nBest RMSE:", best_rmse, "\n")
best_model_performance <- predict(best_model, test_data)

# ---------------------------- RF Tuning Custom --------------------------------
customRF_regression <- list(type = "Regression", 
                            library = "randomForest", 
                            loop = NULL)
customRF_regression$parameters <- data.frame(parameter = c("mtry", 
                                                           "ntree",
                                                           "nodesize",
                                                           "sampsize",
                                                           "splitrule"), 
                                             class = rep("numeric", 5), 
                                             label = c("mtry", 
                                                       "ntree",
                                                       "nodesize",
                                                       "sampsize",
                                                       "splitrule"
                                             )
)
customRF_regression$grid <- function(x, 
                                     y,
                                     len = NULL, 
                                     search = "grid") {}
customRF_regression$fit <- function(x, 
                                    y, 
                                    wts,
                                    param, 
                                    lev, 
                                    last,
                                    weights, 
                                    classProbs, 
                                    ...) {
  randomForest(x, 
               y, 
               mtry = param$mtry,
               ntree=param$ntree, 
               ...)
}
customRF_regression$predict <- function(modelFit, 
                                        newdata, 
                                        preProc = NULL, 
                                        submodels = NULL)
  predict(modelFit, newdata)
customRF_regression$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) NULL
customRF_regression$sort <- function(x) x[order(x[,1]),]
customRF_regression$levels <- function(x) x$classes

# Train model
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=c(3, 4, 5, 6, 7, 8), 
                        .ntree=c(10, 15, 20, 30, 100, 200, 500, 1000),
                        .nodesize=c(1, 5, 10),
                        .sampsize=c(50, 100, 150),
                        .splitrule=c("gini", "extratrees")
)
start_time <- Sys.time()
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time
custom_regression <- caret::train(x = predictors,
                                  y = outcome,
                                  method = customRF_regression,
                                  metric = "RMSE",
                                  trControl = control,
                                  tuneGrid = tunegrid,
                                  do.trace = TRUE
)
summary(custom_regression)
plot(custom_regression)

# ---------------------------- RF Tuning ranger --------------------------------
# Train a GT model
n_features <- 11
num.trees <- n_features * 10
ref <- ranger(
  lai_lidar ~ ., 
  data = training_data,
  num.trees = num.trees,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order",
  seed = 0
)

# get OOB RMSE
default_rmse <- sqrt(ref$prediction.error)

# create hyperparameter grid
hyper_grid <- expand.grid(
  # num.trees = c(10, 15, 20, 30, 50, 100, 200),œ
  mtry = c(2, 3, 4, 5, 6, 7, 8),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)

# Execute full cartesian grid search
for (i in seq_len(nrow(hyper_grid))) {
  # Fit model for ith hyperparameter combination
  fit <- ranger(
    formula         = lai_lidar ~ ., 
    data            = training_data, 
    num.trees       = num.trees, # hyper_grid$num.trees[i],
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = TRUE,
    seed            = 0,
    respect.unordered.factors = 'order'
  )
  
  # Export OOB error 
  hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
  
  # Make predictions on the test dataset
  predictions <- predict(fit, data = test_data)$predictions
  
  # Calculate RMSE on test dataset
  test_rmse <- sqrt(mean((test_data$lai_lidar - predictions)^2))
  
  # Store test RMSE in hyper_grid
  hyper_grid$test_rmse[i] <- test_rmse
}

# Assess top 10 models
top_models <- hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)

# Print top models
print(top_models)
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
pls_model <- plsr(formula, data = mormal_training_data, ncomp = 4)  # Adjust ncomp as needed

# Summary of the PLS model
summary(pls_model)

# Make predictions on the test dataset
test_predictions <- predict(pls_model, newdata = mormal_test_data, ncomp=4)

comparison <- data.frame(Actual = mormal_test_data$lai_lidar,
                         Predicted = as.double(unname(test_predictions[, , 1])))

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

visualization <- data.frame(Actual = unstandardize_variable(
  mormal_test_data$lai_lidar,
  mormal_lai_lidar_train_mean, 
  mormal_lai_lidar_train_sd),
  Predicted = unstandardize_variable(
    as.double(unname(test_predictions[, , 1])),
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

# -------------------------------- Bonus ---------------------------------------
# Calculate permutation importance without variable names
calculate_permutation_importance <- function(model, data, target_column_index) {
  perm_importance <- numeric(ncol(data))
  baseline_rmse <- Metrics::rmse(predict(model, data), data[[target_column_index]])
  
  for (i in seq_along(data)) {
    perm_data <- data
    perm_data[[i]] <- sample(perm_data[[i]])
    perm_rmse <- Metrics::rmse(predict(model, perm_data), data[[target_column_index]])
    perm_importance[i] <- baseline_rmse - perm_rmse
  }
  names(perm_importance) <- names(validation_data)
  perm_importance <- sort(perm_importance, decreasing = TRUE)
  return(data.frame(variable = names(perm_importance), 
                    importance = perm_importance))
}

# Calculate permutation importance
perm_importance <- calculate_permutation_importance(rf_model, validation_data, 1)
barplot(perm_importance$importance[-1], 
        names.arg = perm_importance$variable[-1],
        main = "Permutation Importance",
        xlab = "Predictor Variable Index",
        ylab = "Importance")
print(perm_importance)

# Combine rankings into a single dataframe
combined_rankings <- data.frame(
  variable = names(normalized_importance),
  norm_rank = 1:length(normalized_importance),
  perm_rank = match(names(normalized_importance), perm_importance$variable)
)
combined_rankings$perm_rank <- nrow(combined_rankings) +2 - combined_rankings$perm_rank

# Set up colors for the bars
colors <- c("blue", "red")

# Plot the rankings for each variable
barplot(
  rbind(combined_rankings$norm_rank, combined_rankings$perm_rank),
  beside = TRUE,
  names.arg = combined_rankings$variable,
  col = colors,
  main = "Variable Rankings",
  xlab = "Variable",
  ylab = "Ranking"
)
legend("topright", 
       legend = c("Normalized Importance", "Permutation Importance"), 
       fill = colors)

# Bonus Plots
explainer <- DALEX::explain(model = rf_model,
                            data = as.matrix(training_data[-1]),
                            y = training_data$lai_lidar_val,
                            verbose = FALSE)
pdp_plot <- predict_profile(explainer, pred.var = predictor_vars, new_observation = training_data[1, ])
plot(pdp_plot)
complexity_plot <- model_profile(explainer)
plot(complexity_plot)
