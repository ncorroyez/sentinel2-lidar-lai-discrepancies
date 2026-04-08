# ------------------------------ Code Info -------------------------------------
# title: "6a_explain_heterogeneity_unique_site.R"
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
library(doFuture)
# library(caret)
library(future)
library(parallel)
library(doParallel)
library(stringr)
library(randomForest)
# library(caret)

# ------ Define working dir as the directory where the script is located -------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# --------------------------- Import useful functions --------------------------
source("../libraries/functions_JBF.R")
source("../libraries/functions_plots.R")

# ----------------------------- RF hyperparameters ----------------------------
ntree <- 30

# ----------------------------- Plotting Parameters ----------------------------
default_par <- par(no.readonly = TRUE)
par(mar=c(5, 5, 4, 2) + 0.1) # Adjust margins as per your preference
par(oma=c(0, 0, 0, 0)) # Adjust outer margins as per your preference
par(cex=1.6, cex.axis=1.6, cex.names=1.6)

# --------------------------- Results Reproducibility --------------------------
training_sample_size <- 100000
test_sample_size <- 5000
set.seed(0)

# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
lidar_dir <- "LiDAR"
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
sites <- c('Aigoual')

# ---------------------------------- Process  ---------------------------------- 
# nb workers in parallel for SFS
nbWorkers <- 4

# create site directories
results_path <- figure_path_site <- list()

for (site in sites){
  results_path[[site]] <- file.path('../../03_RESULTS', site)
}

# --------------------------- Data processing ----------------------------------
regression_data <- train_test_data <- list()
SelectedVars <- EvolCorr <- EvolRMSE <- StabilityIndex <- list()
SelectedVars_RF <- EvolCorr_RF <- EvolRMSE_RF <- StabilityIndex_RF <- list()
force_s2lai_bools <- c(TRUE, FALSE)
# force_s2lai_bools <- c(FALSE)

# process each site independently
for (site in sites){
  
  if (site == "Aigoual"){
    # Define the composition masks
    composition_masks <- list(
      # Full_Composition = "Full_Composition",
      # Deciduous_Flex = "Deciduous_Flex",
      Deciduous_Only = "Deciduous_Only"
      # Coniferous_Flex = "Coniferous_Flex",
      # Coniferous_Only = "Coniferous_Only"
    )
  }
  if (site == "Blois"){
    # Define the composition masks
    composition_masks <- list(
      # Full_Composition = "Full_Composition",
      # Deciduous_Flex = "Deciduous_Flex",
      Deciduous_Only = "Deciduous_Only"
    )
  }
  if (site == "Mormal"){
    # Define the composition masks
    composition_masks <- list(
      # Full_Composition = "Full_Composition",
      # Deciduous_Flex = "Deciduous_Flex",
      Deciduous_Only = "Deciduous_Only"
    )
  }
  
  # Define the test data path (Full Composition)
  test_metrics_path <- file.path(results_path[[site]], 'Metrics', 
                                 'Deciduous_Only')
  test_metrics <- list.files(test_metrics_path, pattern = ".tif$")
  lidar_lai_file <- grep("lidarlai", test_metrics, value = TRUE)
  test_metrics <- test_metrics[!test_metrics %in% lidar_lai_file]
  test_predicted_path <- data.frame('lidar_LAI' = file.path(test_metrics_path, 
                                                            lidar_lai_file))
  test_predictors_path <- data.frame('name' = str_remove(test_metrics, ".tif"), 
                                     'file' = file.path(test_metrics_path, 
                                                        test_metrics))
  
  # Extract values for test data
  regression_test_data <- extract_raster_info(predicted_path = test_predicted_path, 
                                              predictors_path = test_predictors_path)
  
  # Keep the test data constant (Deciduous Only)
  testSet <- train_test(predicted = regression_test_data$predicted_val, 
                        predictors = regression_test_data$predictors_val, 
                        training_sample_size = training_sample_size,
                        test_sample_size = test_sample_size)$testSet
  TestData <- data.frame(testSet)
  names(TestData) <- c('target', names(testSet$predictors))
  
  for (mask_name in names(composition_masks)) {
    metrics_path <- file.path(results_path[[site]], 'Metrics', 
                              composition_masks[[mask_name]])
    metrics <- list.files(metrics_path, pattern = ".tif$")
    
    # Create a DataFrame for the training data
    lidar_lai_file <- grep("lidarlai", metrics, value = TRUE)
    metrics <- metrics[!metrics %in% lidar_lai_file]
    predicted_path <- data.frame('lidar_LAI' = file.path(metrics_path, 
                                                         lidar_lai_file))
    predictors_path <- data.frame('name' = str_remove(metrics, ".tif"), 
                                  'file' = file.path(metrics_path, metrics))
    
    # Extract values to be used during regression
    regression_data[[site]] <- extract_raster_info(predicted_path = predicted_path, 
                                                   predictors_path = predictors_path)
    
    # Create training and test data for current mask
    train_test_data[[site]] <- train_test(predicted = regression_data[[site]]$predicted_val, 
                                          predictors = regression_data[[site]]$predictors_val, 
                                          training_sample_size = training_sample_size,
                                          test_sample_size = test_sample_size)
    
    TrainingData <- data.frame(train_test_data[[site]]$trainingSet)
    names(TrainingData) <- c('target', 
                             names(train_test_data[[site]]$trainingSet$predictors))
    formula <- as.formula(paste("target ~", 
                                paste(names(TrainingData)[-1], collapse = " + ")))
    
    # Train random forest using the masked training data
    rf_model <- randomForest(formula, 
                             data = TrainingData,
                             ntree = ntree,
                             importance = TRUE,
                             do.trace = TRUE)
    
    # Apply random forest on test dataset
    predictions <- predict(rf_model, newdata = TestData)
    
    # Compute metrics
    correlation <- cor(predictions, TestData$target)
    rsquared <- caret::R2(predictions, TestData$target)
    rmse <- Metrics::rmse(predictions, TestData$target)
    bias <- Metrics::bias(predictions, TestData$target)
    mae <- Metrics::mae(predictions, TestData$target)
    
    # Print metrics
    cat("Metrics:",
        "\nCorrelation:", correlation, 
        "\nR-squared:", rsquared, 
        "\nRoot Mean Squared Error:", rmse,
        "\nBias:", bias, 
        "\nMAE:", mae)
    
    filename <- paste0('LiDAR_LAI_model_', site, '_', mask_name)
    
    for(force_s2lai in force_s2lai_bools){
      figure_path_site[[site]] <- file.path(results_path[[site]], 
                                            "ML_Figures", 
                                            mask_name,
                                            paste0("force_s2lai_",
                                                   force_s2lai))
      if (!dir.exists(figure_path_site[[site]])) {
        dir.create(figure_path_site[[site]],
                   showWarnings = FALSE, 
                   recursive = TRUE)
      }
      
      plot_density_scatterplot(c(predictions), 
                               TestData$target,
                               "Estimated LAI", 
                               "Measured LiDAR LAI",
                               title = paste("Estimated LAI vs LiDAR LAI:", 
                                             site, 
                                             mask_name),
                               xlimits = c(0,15), 
                               ylimits = c(0,15), 
                               dirname = figure_path_site[[site]],
                               filename = filename)
      
      # perform SFS to assess LAI using optimal feature selection
      # ---------------------------------- SFS RF --------------------------------
      message('**perform SFS RF for ', site, ' ', mask_name, ' forest mask**')
      Vars_SFS <- names(train_test_data[[site]]$training$predictors)
      outSFS <- SFS_Regression(algo = 'RF',
                               TrainingData = TrainingData,
                               TestData = TestData, 
                               Vars_SFS = Vars_SFS,
                               force_s2lai = force_s2lai,
                               ntree = ntree,
                               site = site, 
                               figure_path_site = figure_path_site[[site]], 
                               nbWorkers = nbWorkers)
      
      # Store the results
      selected_vars_key <- paste(site, mask_name, sep = "_")
      SelectedVars_RF[[selected_vars_key]] <- outSFS$SelectedVars
      EvolCorr_RF[[selected_vars_key]] <- outSFS$EvolCorr
      EvolRMSE_RF[[selected_vars_key]] <- outSFS$EvolRMSE
      StabilityIndex_RF[[selected_vars_key]] <- outSFS$StabilityIndex
      
      # Save the results to a text file
      save_results_to_file(SelectedVars_RF[[selected_vars_key]],
                           EvolCorr_RF[[selected_vars_key]],
                           EvolRMSE_RF[[selected_vars_key]],
                           StabilityIndex_RF[[selected_vars_key]],
                           figure_path_site[[site]],
                           site,
                           mask_name)
    }
  }
}