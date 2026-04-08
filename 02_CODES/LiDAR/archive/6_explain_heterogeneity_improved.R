# ------------------------------ Code Info -------------------------------------
# title: "6_explain_heterogeneity.R"
# authors: Nathan CORROYEZ & Jean-Baptiste FERET, 
# UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, 
# F-34196, Montpellier, France
# output: html_document
# last_update: "2024-06-03"

# -------------- (Optional) Clear the environment and free memory --------------
rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# ------------------------------ Libraries -------------------------------------
library(foreach)
# library(doFuture)
# library(caret)
library(future)
# library(parallel)

# ------ Define working dir as the directory where the script is located -------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# --------------------------- Import useful functions --------------------------
source("../libraries/functions_JBF.R")
source("../libraries/functions_plots.R")

# ----------------------------- RF hyperparameters ----------------------------
ntree <- 50
mtry <- 3

# ----------------------------- Plotting Parameters ----------------------------
default_par <- par(no.readonly = TRUE)
par(mar=c(5, 5, 4, 2) + 0.1) # Adjust margins as per your preference
par(oma=c(0, 0, 0, 0)) # Adjust outer margins as per your preference
par(cex=1.6, cex.axis=1.6, cex.names=1.6)

# --------------------------- Results Reproducibility --------------------------
training_sample_size <- 5000
test_sample_size <- 50000

# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
lidar_dir <- "LiDAR"
# sites <- c('Mormal', 'Blois', 'Aigoual')
sites <- c('Aigoual')
composition_masks <- c("Full_Composition", "Deciduous_Flex", "Deciduous_Only")

# ---------------------------------- Process  ---------------------------------- 
# nb workers in parallel for SFS
nbWorkers <- 4

# create site directories
results_path <- list()

for (site in sites){
  results_path[[site]] <- file.path('../../03_RESULTS', site)
}

# --------------------------- Data processing ----------------------------------
regression_data <- train_test_data <- list()
SelectedVars <- EvolCorr <- SelectedVars_RF <- EvolCorr_RF <- list()

# process each site independently
for (site in sites){
  figure_path_site <- file.path(results_path, "ML_Figures")
  if (!dir.exists(figure_path_site)) {
    dir.create(figure_path_site, showWarnings = FALSE, recursive = TRUE)
  }
  
  for (composition_mask in composition_masks) {
    
    # Metrics Path for training
    metrics_path_train <- file.path(results_path[[site]],
                                    'Metrics',
                                    composition_mask)
    metrics_train <- list.files(metrics_path_train, pattern = ".tif$")
    
    # Create a DataFrame for predicted value's path: LiDAR LAI (for training)
    lidar_lai_file <- grep("lidarlai", metrics_train, value = TRUE)
    # lidar_lai_name <- unlist(lapply(
    #   stringr::str_split(string = lidar_lai_file, pattern = '_'), '[[', 1))
    # predicted_path_train <- data.frame('name' = lidar_lai_name,
    #                                    'file' = file.path(metrics_path_train,
    #                                                       lidar_lai_file))
    predicted_path <- data.frame('lidar_LAI' = file.path(metrics_path_train,
                                                         lidar_lai_file))
    
    # Remove PAD_Profiles directory
    pad_profiles_dir <- grep("PAD_Profles", metrics_train, value = TRUE)
    metrics_train <- metrics_train[!metrics_train %in% lidar_lai_file]
    
    # Create a DataFrame for predictors' path (for training)
    metrics_names_train <- unlist(lapply(
      stringr::str_split(string = metrics_train, pattern = '_'), '[[', 1))
    predictors_path <- data.frame('name' = metrics_names_train, 
                                        'file' = file.path(metrics_path_train, 
                                                           metrics_train))
    
    # Extract values to be used during regression (for training)
    regression_data_train <- extract_raster_info(
      predicted_path = predicted_path, 
      predictors_path = predictors_path)
    
    # Scatterplot for direct comparison between S2 and LiDAR LAI (for training)
    filename <- paste0('S2_vs_LiDAR_LAI_', site, '_', composition_mask)
    plot_density_scatterplot(regression_data_train$predictors_val$s2lai,
                             regression_data_train$predicted_val$lidar_LAI,
                             "Sentinel-2 LAI", 
                             "LiDAR LAI",
                             title = paste("S2 LAI vs LiDAR LAI:", 
                                           site, 
                                           composition_mask),
                             xlimits = c(0, 15),
                             ylimits = c(0, 15), 
                             dirname = figure_path_site,
                             filename = filename)
  
  # ----------- Apply RF regression model Using S2 & LiDAR data ----------------
  # Create Train and Test data
  train_test_data <- train_test(
    predicted = regression_data_train$predicted_val, 
    predictors = regression_data_train$predictors_val, 
    training_sample_size = training_sample_size,
    test_sample_size = test_sample_size)
  
    # Train random forest using full set of variables
    TrainData <- data.frame(train_test_data$trainingSet)
    names(TrainData) <- c('target', 
                          names(train_test_data$trainingSet$predictors))
    
    formula <- target ~ (s2lai 
    + mean 
    + cv 
    + rumple
    + vci 
    + lcv 
    + lskew
    + vdr
    + std
    + variance
    )
    rf_model <- randomForest::randomForest(formula, data = TrainData,
                                           ntree = ntree, mtry = mtry,
                                           importance = TRUE, do.trace = TRUE)
    
    # Metrics Path for testing (Full_Composition)
    metrics_path_test <- file.path(results_path[[site]], 
                                   'Metrics', 
                                   'Full_Composition')
    metrics_test <- list.files(metrics_path_test, pattern = ".tif$")
    
    # Create a DataFrame for predicted value's path: LiDAR LAI (for testing)
    lai_lidar_files_test <- grep("lidarlai", metrics_test, value = TRUE)
    predicted_path_test <- data.frame('lidar_LAI' = file.path(metrics_path_test,
                                                              lidar_lai_file))
    
    # Remove 'lai_lidar' files from the predictors' path
    metrics_test <- metrics_test[!metrics_test %in% lidar_lai_file]
    
    # Create a DataFrame for predictors' path (for testing)
    metrics_names_test <- unlist(
      lapply(stringr::str_split(string = metrics_test, pattern = '_'), '[[', 1))
    predictors_path_test <- data.frame('name' = metrics_names_test, 
                                       'file' = file.path(metrics_path_test, 
                                                          metrics_test))
    
    # Extract values to be used during regression (for testing)
    regression_data_test <- extract_raster_info(
      predicted_path = predicted_path_test, 
      predictors_path = predictors_path_test)
    
    # Create Test data
    test_data <- train_test(predicted = regression_data_test$predicted_val, 
                            predictors = regression_data_test$predictors_val, 
                            training_sample_size = training_sample_size,
                            test_sample_size = test_sample_size)
    
    TestData <- data.frame(test_data$testSet)
    names(TestData) <- c('target', names(test_data$testSet$predictors))
    predictions <- predict(rf_model, newdata = TestData)
    
    # Correlation test
    cor_result <- cor.test(predictions, TestData$target)
    print(cor_result)
    
    filename <- paste0('LiDAR_LAI_model_', site, '_', composition_mask)
    plot_density_scatterplot(c(predictions), TestData$target,
                             "Estimated LAI", "Measured LiDAR LAI",
                             title = paste("Estimated LAI vs LiDAR LAI:", site),
                             xlimits = c(0, 15), ylimits = c(0, 15), 
                             dirname = figure_path_site, filename = filename)
  
  # perform SFS to assess LAI using optimal feature selection
  # ---------------------------------- SFS RF ----------------------------------
  message('perform SFS RF')
  Vars_SFS <- names(train_test_data$training$predictors)
  outSFS <- SFS_Regression(algo = 'RF', 
                           TrainingData = TrainingData, 
                           TestData = TestData, 
                           Vars_SFS = Vars_SFS, 
                           ntree = ntree, 
                           mtry = mtry, 
                           site = site, 
                           figure_path_site = figure_path_site, 
                           nbWorkers = nbWorkers)
  SelectedVars_RF[[site]] <- outSFS$SelectedVars
  EvolCorr_RF[[site]] <- outSFS$EvolCorr
  }
}
