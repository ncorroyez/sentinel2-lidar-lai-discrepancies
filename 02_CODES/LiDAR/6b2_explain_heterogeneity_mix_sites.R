# ------------------------------ Code Info -------------------------------------
# title: "6b_explain_heterogeneity_mix_sites.R"
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
library(stringr)
library(randomForest)
library(pls)
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
test_sample_size <- 3000
set.seed(0)

# ---------------------------- Directories & Setup  ---------------------------- 
data <- "../../01_DATA"
lidar_dir <- "LiDAR"
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
# sites <- c('Blois')

# ---------------------------------- Process  ---------------------------------- 
# nb workers in parallel for SFS
nbWorkers <- 6

# create site directories
results_path <- figure_path_site <- list()

for (site in sites){
  results_path[[site]] <- file.path('../../03_RESULTS', site)
}

# --------------------------- Data processing ----------------------------------
regression_train_data <- regression_test_data <- train_data <- list()
SelectedVars <- EvolCorr <- EvolRMSE <- list()
SelectedVars_RF <- EvolCorr_RF <- EvolRMSE_RF <- list()
force_s2lai_bools <- c(TRUE, FALSE)
algo <- 'RF'

# Define the combinations
combinations <- list(
  list(training_sites = c("Mormal", "Aigoual", "Blois"), 
       test_site = c("Mormal", "Aigoual", "Blois")),
  list(training_sites = c("Mormal", "Aigoual", "Blois"), 
       test_site = c("Mormal", "Aigoual", "Blois")),
  list(training_sites = c("Mormal", "Aigoual", "Blois"), 
       test_site = c("Mormal", "Aigoual", "Blois")),
  list(training_sites = c("Aigoual"), test_site = "Aigoual"),
  list(training_sites = c("Blois"), test_site = "Blois"),
  list(training_sites = c("Mormal"), test_site = "Mormal")
)

combinations <- list(
  list(training_sites = c("Mormal", "Aigoual", "Blois"), 
       test_site = c("Mormal", "Aigoual", "Blois"))
)

# Iterate over combinations
for (combo in combinations) {
  training_sites <- combo$training_sites
  test_site <- combo$test_site
  
  # Always train on same number of samples
  if (length(training_sites) == 1){
    training_sample_size <- 30000
  }
  if (length(training_sites) == 2){
    training_sample_size <- 15000
  }
  if (length(training_sites) == 3){
    training_sample_size <- 10000
  }
  
  # Construct a name for the combination
  combo_name <- paste(training_sites, collapse = "_")
  
  # Aggregate training data from selected sites
  all_training_data <- all_test_data <- data.frame()
  
  for (site in training_sites){
    # Define composition masks based on the site
    if (site == "Aigoual"){
      composition_masks <- list(
        Deciduous_Only = "Deciduous_Only"
      )
    } else if (site == "Blois"){
      composition_masks <- list(
        Deciduous_Only = "Deciduous_Only"
      )
    } else if (site == "Mormal"){
      composition_masks <- list(
        Deciduous_Only = "Deciduous_Only"
      )
    }
    
    # Iterate over composition masks
    for (mask_name in names(composition_masks)) {
      # Define the metrics path for training data using the current mask
      metrics_path <- file.path(results_path[[site]], 'Metrics', 
                                composition_masks[[mask_name]])
      
      # Get the list of metrics files in the metrics path
      metrics <- list.files(metrics_path, pattern = ".tif$")
      lidar_lai_file <- grep("lidarlai", metrics, value = TRUE)
      metrics <- metrics[!metrics %in% lidar_lai_file]
      predicted_path <- data.frame('lidar_LAI' = file.path(metrics_path, 
                                                           lidar_lai_file))
      predictors_path <- data.frame('name' = str_remove(metrics, ".tif"), 
                                    'file' = file.path(metrics_path, metrics))
      
      # if(all(c("Mormal", "Aigoual", "Blois") %in% training_sites) && test_site == "Aigoual"){
      #   names_to_keep <- c(
      #     # "s2lai_res_10_m",
      #     "mean_res_10_m",
      #     "lcv_res_10_m"
      #     # "rumple_res_10_m"
      #     # "lskew_res_10_m", 
      #     # "max_res_10_m",
      #     # "cv_res_10_m", "vci_res_10_m", "vdr_res_10_m", "std_res_10_m", "variance_res_10_m",
      #     # "shade_res_10_m"
      #     )
      #   predictors_path <- subset(predictors_path, grepl(paste(names_to_keep, collapse="|"), name))
      # }
      if(all(c("Mormal", "Aigoual", "Blois") %in% training_sites) 
              && all(c("Mormal", "Aigoual", "Blois") %in% test_site)){
        names_to_keep <- c(
          "s2lai_res_10_m",
          "mean_res_10_m",
          # "vci_res_10_m",
          "lcv_res_10_m",
          "lskew_res_10_m"
          # "max_res_10_m"
          # "std_res_10_m"
          # "rumple_res_10_m"
          # "cv_res_10_m"
          )
        # Create a pattern to match exact names
        pattern <- paste0("^", names_to_keep, "$", collapse = "|")
        
        # Subset predictors_path using the exact name pattern
        predictors_path <- subset(predictors_path, grepl(pattern, name))
        # cat(1)
      }
      else if(all(c("Mormal", "Aigoual", "Blois") %in% training_sites) && test_site == "Mormal"){
        names_to_keep <- c(
          # "s2lai_res_10_m", 
          "mean_res_10_m", "rumple_res_10_m", "max_res_10_m", "lcv_res_10_m", "lskew_res_10_m")
        predictors_path <- subset(predictors_path, grepl(paste(names_to_keep, collapse="|"), name))
      }
      else if(all(c("Mormal") %in% training_sites) && test_site == "Mormal"){
        names_to_keep <- c(
          "s2lai_res_10_m",
          "lcv_res_10_m", 
          "mean_res_10_m" 
          # "vci_res_10_m", "lskew_res_10_m", "vdr_res_10_m"
          )
        predictors_path <- subset(predictors_path, grepl(paste(names_to_keep, collapse="|"), name))
      }
      else if(all(c("Blois") %in% training_sites) && test_site == "Blois"){
        names_to_keep <- c(
          # "s2lai_res_10_m",
          "mean_res_10_m", 
          "lcv_res_10_m"
          # "vci_res_10_m", "lskew_res_10_m" "vdr_res_10_m"
          )
        predictors_path <- subset(predictors_path, grepl(paste(names_to_keep, collapse="|"), name))
      }
      else if(all(c("Aigoual") %in% training_sites) && test_site == "Aigoual"){
        names_to_keep <- c(
          # "s2lai_res_10_m",
          "mean_res_10_m",
          "lcv_res_10_m"
          # "lskew_res_10_m", "rumple_res_10_m", "max_res_10_m"
          )
        predictors_path <- subset(predictors_path, grepl(paste(names_to_keep, collapse="|"), name))
      }
      else{
        stop("Error: Proposed Train-Test Setup is not valid.")
      }
      
      # Extract values for training data
      regression_train_data[[site]] <- extract_raster_info(predicted_path = predicted_path, 
                                                           predictors_path = predictors_path)
      
      # Create training data for current site and mask
      train_data[[site]] <- train_test(predicted = regression_train_data[[site]]$predicted_val, 
                                       predictors = regression_train_data[[site]]$predictors_val, 
                                       training_sample_size = training_sample_size,
                                       test_sample_size = test_sample_size)
      
      TrainingData <- data.frame(train_data[[site]]$trainingSet)
      names(TrainingData) <- c('target', 
                               names(train_data[[site]]$trainingSet$predictors))
      
      # Aggregate training data from selected sites and masks
      all_training_data <- rbind(all_training_data, TrainingData)
      
      # Shuffle
      all_training_data <- all_training_data[sample(nrow(all_training_data)), ]
      rownames(all_training_data) <- NULL
      
      TestData <- data.frame(train_data[[site]]$testSet)
      names(TestData) <- c('target',
                           names(train_data[[site]]$testSet$predictors))

      # Aggregate training data from selected sites and masks
      all_test_data <- rbind(all_test_data, TestData)

      # Shuffle
      all_test_data <- all_test_data[sample(nrow(all_test_data)), ]
      rownames(all_test_data) <- NULL
    }
  }
  
  # TestData <- data.frame(train_data$Mormal$testSet) #train_data[[site]
  # names(TestData) <- c('target', 
  #                      names(train_data$Mormal$testSet$predictors))
  
  # Train random forest using the aggregated training data
  formula <- as.formula(paste("target ~", 
                              paste(names(all_training_data)[-1], 
                                    collapse = " + ")))
  
  # rf_model <- randomForest(formula, 
  #                          data = all_training_data,
  #                          ntree = ntree,
  #                          importance = TRUE,
  #                          do.trace = TRUE)
  
  # Fit the multivariate regression model
  lm_model <- lm(formula, data = all_training_data)
  
  # Summary of the model
  print(summary(lm_model))
  
  # Perform stepwise regression
  stepwise_model <- step(lm_model, direction = "both")
  
  # Summary of the stepwise model (AIC)
  summary(stepwise_model)
  
  # Fit the PLS regression model
  # pls_model <- plsr(formula, data = all_training_data, ncomp = 4)
  
  # Summary of the PLS model
  # summary(pls_model)
  
  # Apply random forest on test dataset
  predictions <- predict(stepwise_model, newdata = TestData)
  # predictions <- as.data.frame(predictions)
  # predictions <- predictions$`target.3 comps`
  
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
  
  # Load the test data for the site
  # for (mask_name in names(composition_masks)) {
  #   # Define the test data path
  #   test_metrics_path <- file.path(results_path[[test_site]], 'Metrics', 
  #                                  composition_masks[[mask_name]])
  #   test_metrics <- list.files(test_metrics_path, pattern = ".tif$")
  #   lidar_lai_file <- grep("lidarlai", test_metrics, value = TRUE)
  #   test_metrics <- test_metrics[!test_metrics %in% lidar_lai_file]
  #   test_predicted_path <- data.frame('lidar_LAI' = file.path(test_metrics_path, 
  #                                                             lidar_lai_file))
  #   test_predictors_path <- data.frame('name' = str_remove(test_metrics, ".tif"), 
  #                                      'file' = file.path(test_metrics_path, 
  #                                                         test_metrics))
  #   
  #   # Extract values for test data
  #   regression_test_data <- extract_raster_info(predicted_path = test_predicted_path, 
  #                                               predictors_path = test_predictors_path)
  #   
  #   # Keep the test data constant (Deciduous Only)
  #   testSet <- train_test(predicted = regression_test_data$predicted_val, 
  #                         predictors = regression_test_data$predictors_val, 
  #                         training_sample_size = training_sample_size,
  #                         test_sample_size = test_sample_size)$testSet
  #   
  #   # Apply random forest on the site's test dataset
  #   TestData <- data.frame(testSet)
  #   names(TestData) <- c('target', names(testSet$predictors))
  #   
  #   # Apply random forest on test dataset
  #   predictions <- predict(rf_model, newdata = TestData)
  #   
  #   # Compute metrics
  #   correlation <- cor(predictions, TestData$target)
  #   rsquared <- caret::R2(predictions, TestData$target)
  #   rmse <- Metrics::rmse(predictions, TestData$target)
  #   bias <- Metrics::bias(predictions, TestData$target)
  #   mae <- Metrics::mae(predictions, TestData$target)
  #   
  #   # Print metrics
  #   cat("Metrics:",
  #       "\nCorrelation:", correlation, 
  #       "\nR-squared:", rsquared, 
  #       "\nRoot Mean Squared Error:", rmse,
  #       "\nBias:", bias, 
  #       "\nMAE:", mae)
  # }
}