# ---
# title: "functions_explain_heterogeneity.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-14"
# ---

standardize_train_data <- function(data) {
  standardized_data <- data
  means <- vector("numeric", length = ncol(data))
  sds <- vector("numeric", length = ncol(data))
  
  for (id_var in seq_along(data)) {
    train_mean <- mean(data[, id_var], na.rm = TRUE)
    train_sd <- sd(data[, id_var], na.rm = TRUE)
    
    means[id_var] <- train_mean
    sds[id_var] <- train_sd
    
    standardized_data[, id_var] <- (data[, id_var] - train_mean) / train_sd
  }
  
  return(list(standardized_data = standardized_data,
              means = means,
              sds = sds)
  )
}

standardize_test_data <- function(data, means, sds) {
  standardized_data <- data
  
  for (id_var in seq_along(data)) {
    standardized_data[, id_var] <- (data[, id_var] - means[id_var]) / sds[id_var]
  }
  
  return(standardized_data)
}

unstandardize_variable <- function(variable, mean_val, sd_val) {
  unstandardized_variable <- (variable * sd_val) + mean_val
  return(unstandardized_variable)
}

initialize_predicted_raster <- function(masks_dir){
  # Predicted raster will be LiDAR LAI
  lai_lidar_raster <- terra::rast(file.path(masks_dir,
                                            "lai_lidar_masked_res_10_m.envi")) 
  lai_lidar_matrix <- as.matrix(lai_lidar_raster)
  na_indices <- which(is.na(lai_lidar_matrix), arr.ind = TRUE)
  non_na_indices <- which(!is.na(lai_lidar_matrix), arr.ind = TRUE)
  na_indices <- na_indices[,1]
  non_na_indices <- non_na_indices[,1]
  predicted_raster <- rast(lai_lidar_raster)
  return(predicted_raster)
}

extract_predictors_values <- function(masks_dir, lidr_dir){
  
  # We want to predict LiDAR LAI
  lai_lidar_raster <- terra::rast(file.path(masks_dir,
                                            "lai_lidar_masked_res_10_m.envi"))
  # Predictors
  lai_s2_raster <- terra::rast(
    file.path(masks_dir, "lai_s2_masked_res_10_m.envi"))
  mean_heights_raster <- terra::rast(
    file.path(masks_dir, "mnc_mean_heights_res_10_m.envi"))
  cv_raster <- terra::rast(file.path(masks_dir,
                                     "mnc_coeff_variation_res_10_m.envi"))
  rumple_raster <- terra::rast(file.path(lidr_dir, 
                                         "rumple_res_10_m_non_norm.tif"))
  vci_raster <- terra::rast(file.path(lidr_dir, "vci_res_10_m_non_norm.tif"))
  lcv_raster <- terra::rast(file.path(lidr_dir, "lcv_res_10_m_non_norm.tif"))
  lskew_raster <- terra::rast(file.path(lidr_dir, 
                                        "lskew_res_10_m_non_norm.tif"))
  
  # Define indices: keep same number of values for each raster
  indices <- list(
    na_indices_lcv = countNA(lcv_raster),
    na_indices_std = countNA(cv_raster),
    na_indices_lai_lidar = countNA(lai_lidar_raster),
    na_indices_lai_s2 = countNA(lai_s2_raster)
  )
  
  # Loop through indices
  for (index_name in names(indices)) {
    na_index <- indices[[index_name]]
    lai_lidar_raster[na_index == 1] <- NA
    lai_s2_raster[na_index == 1] <- NA
    mean_heights_raster[na_index == 1] <- NA
    cv_raster[na_index == 1] <- NA
    rumple_raster[na_index == 1] <- NA
    vci_raster[na_index == 1] <- NA
    lcv_raster[na_index == 1] <- NA
    lskew_raster[na_index == 1] <- NA
  }
  
  # Values 
  lai_lidar <- values(lai_lidar_raster)
  lai_s2 <- values(lai_s2_raster)
  mean_h <- values(mean_heights_raster)
  cv <- values(cv_raster)
  rumple <- values(rumple_raster)
  vci <- values(vci_raster)
  lcv <- values(lcv_raster)
  lskew <- values(lskew_raster)
  
  # No NA
  lai_lidar <- lai_lidar[complete.cases(lai_lidar), ]
  lai_s2 <- lai_s2[complete.cases(lai_s2), ]
  mean_h <- mean_h[complete.cases(mean_h), ]
  cv <- cv[complete.cases(cv), ]
  rumple <- rumple[complete.cases(rumple), ]
  vci <- vci[complete.cases(vci), ]
  lcv <- lcv[complete.cases(lcv), ]
  lskew <- lskew[complete.cases(lskew), ]
  
  # print(length(lai_lidar))
  # print(length(lai_s2))
  
  return(list(
    lai_lidar = lai_lidar, 
    lai_s2 = lai_s2, 
    mean_h = mean_h, 
    cv = cv, 
    rumple = rumple, 
    vci = vci,
    lcv = lcv, 
    lskew = lskew)
  )
}

create_train_test_datasets <- function(predictors_values, 
                                       total_sample_size,
                                       training_sample_size,
                                       test_sample_size){
  
  # Sample size verification
  if(training_sample_size + test_sample_size > total_sample_size){
    stop(paste("Error: Total sample size is:", total_sample_size, ",",
               "it can't be smaller than the sum of the training sample size",
               "and the testing sample size (respectively", 
               training_sample_size, "and", test_sample_size, ".\n"))
  }
  # Visualize correlation matrix
  df_predictors <- as.data.frame(predictors_values)
  correlation_matrix <- cor(df_predictors[, 
                                          -which(
                                            names(df_predictors)
                                            == "lai_lidar")])
  corrplot::corrplot(correlation_matrix, method = "number", type = "upper")
  
  # Train-test split
  
  # Subset training_data using the row indices
  # train_prop <- 0.7
  # train_indices <- createDataPartition(sampled_data$lai_lidar, 
  #                                      p = train_prop, 
  #                                      list = FALSE)
  # training_data <- sampled_data[train_indices, ]
  # test_data <- sampled_data[-train_indices, ]
  # test_data <- standardized_data
  
  # Sampling
  sample_indices <- sample(nrow(df_predictors), total_sample_size)
  sampled_data <- df_predictors[sample_indices, ]
  
  # Split into training and testing sets
  indices <- sample(1:nrow(sampled_data), training_sample_size)
  training_data <- sampled_data[indices, ]
  test_data <- sampled_data[-indices, ]
  
  # Standardize training data
  train_data_info <- standardize_train_data(training_data)
  standardized_training_data <- train_data_info$standardized_data
  means <- train_data_info$means
  sds <- train_data_info$sds
  
  # Standardize test data using mean and sd from training data
  standardized_test_data <- standardize_test_data(test_data, means, sds)
  
  # Check dimensions
  cat("Dimensions of training_data:", dim(standardized_training_data),
      "\nDimensions of test_data:", dim(standardized_test_data), "\n")
  
  return(list(training_data = standardized_training_data,
              test_data = standardized_test_data,
              lai_lidar_train_mean = train_data_info$means[1],
              lai_lidar_train_sd = train_data_info$sds[1]
  )
  )
}

extract_raster_info <- function(predicted_path, predictors_path){
  
  # Read predicted value and predictors included in raster data
  predicted <- terra::rast(predicted_path[[1]]) 
  predictors <- terra::rast(predictors_path$file)
  names(predictors) <- predictors_path$name
  
  # Eliminate NA from all data when one layer has NA
  na_index_predictors <- terra::countNA(predictors)
  na_index_predicted <- terra::countNA(predicted)
  na_index <- na_index_predicted + na_index_predictors
  predictors[na_index > 0] <- NA
  predicted[na_index > 0] <- NA
  
  # Extract values and eliminate NA
  # Predicted
  predicted_val <- terra::values(predicted)
  predicted_val <- data.frame(predicted_val[complete.cases(predicted_val), ])
  names(predicted_val) <- names(predicted_path)
  
  # Predictors
  predictors_val <- terra::values(predictors)
  predictors_val <- data.frame(predictors_val[complete.cases(predictors_val), ])
  names(predictors_val) <- predictors_path$name
  
  return(list('predictors_val' = predictors_val, 
              'predicted_val' = predicted_val))
}

#' produce training & test sets in preparation for ML 
#' initial data includes predicted and predictors
#' - requires predicted and predictors values, 
#' - extract raster data and eliminates NA from extracted data
#' @param predicted data.frame predicted values
#' @param predictors data.frame predictors values
#' @param training_sample_size numeric number of training samples
#' @param test_sample_size numeric number of test samples
#' @return list 
#' @importFrom terra rast values countNA
#' @export

train_test <- function(predicted, predictors, 
                       training_sample_size = 10000, test_sample_size = NULL){
  
  nbsamples <- nrow(predicted)
  
  # check sample size
  if (nbsamples < training_sample_size)
    stop('Error: training sample size smaller than number of samples available')
  if (nbsamples == training_sample_size)
    stop('Error: training sample size equals number of samples available')
  if (is.null(test_sample_size)) 
    test_sample_size <- nbsamples - training_sample_size
  if (training_sample_size + test_sample_size > nbsamples){
    test_sample_size <- nbsamples - training_sample_size
    warning('Warning: test sample size too large for nb of samples available')
    warning(paste('test sample size reduced to', test_sample_size))
    warning('adjust training sample size tpo increase test sample size')
  }
  
  # Train-test split
  train_ind <- sample(nbsamples, training_sample_size)
  trainingSet <- testSet <- list()
  trainingSet[['predicted']] <- data.frame(predicted[train_ind,])
  trainingSet[['predictors']] <- predictors[train_ind,]
  names(trainingSet[['predicted']]) <- names(predicted)
  
  testSet[['predicted']] <- data.frame(predicted[-train_ind,])
  testSet[['predictors']] <- predictors[-train_ind,]
  
  # adjust to desired number of samples
  testSet[['predicted']] <- data.frame(testSet[['predicted']][seq_len(test_sample_size), ])
  names(testSet[['predicted']]) <- names(predicted)
  testSet[['predictors']] <- testSet[['predictors']][seq_len(test_sample_size), ]
  
  return(list('trainingSet' = trainingSet,
              'testSet' = testSet))
}


#' Forward sequential feature selection (SFS) using a machine learning 
#' algorithm for regression task: rank features optimizing a criterion 
#' (here correlation) for regression task, using training and test sets
#' expects training and test dataset, 
#' can also save density scatterplot for each combination of features tested   
#' @param algo character. currently RF or liquidSVM
#' @param TrainingData numeric. training dataset
#' @param TestData numeric. test dataset
#' @param Vars_SFS character. variable names to be tested from training and test sets
#' @param ntree numeric. nb of trees for RF algo
#' @param mtry numeric. nb of attempts for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @param nbWorkers numeric. number of workers for parallel computing
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression <- function(algo = 'RF', TrainingData, TestData, 
                           Vars_SFS, ntree = NULL, mtry = NULL, site = site, 
                           figure_path_site, nbWorkers = 1){
  
  NbPCs_To_Keep <- length(Vars_SFS)
  CorrSFS <- AssessSFS <- list()
  SelectedVars <- EvolCorr <- c()
  
  # multi-thread feature selection (SFS)
  doFuture::registerDoFuture()
  cl <- parallel::makeCluster(nbWorkers)
  plan("cluster", workers = cl)
  
  # progressbar
  pb <- progress::progress_bar$new(
    format = "Perform feature selection [:bar] :percent in :elapsedfull",
    total = NbPCs_To_Keep, clear = FALSE, width= 100)
  
  # explore each variable
  for (nbvars2select in seq_len(NbPCs_To_Keep)){
    NumVar_list <- as.list(seq_len(length(Vars_SFS)))
    subSFS <- SFS_Regression_step(algo = algo, 
                                  NumVar_list = NumVar_list, 
                                  SelectedVars = SelectedVars,
                                  Vars_SFS = Vars_SFS, 
                                  TrainingData = TrainingData,
                                  TestData = TestData, 
                                  ntree = ntree, 
                                  mtry = mtry,
                                  site = site, 
                                  figure_path_site = figure_path_site)
    pb$tick()
    
    CorrSFS[[nbvars2select]] <- data.frame(
      matrix(unlist(lapply(subSFS,'[[','coeffcorr')),ncol = 1))
    rownames(CorrSFS[[nbvars2select]]) <- Vars_SFS
    SelVar <- which(CorrSFS[[nbvars2select]] == max(
      CorrSFS[[nbvars2select]],na.rm = T))
    AssessSFS[[nbvars2select]] <- c(
      unlist((lapply(subSFS,'[[','AssessedVal'))[SelVar[1]]))
    
    # which criterion to maximize with SFS?
    WhichVar <- rownames(CorrSFS[[nbvars2select]])[SelVar[1]]
    
    # add selected component to selected vars
    SelectedVars <- c(SelectedVars,WhichVar)
    
    # delete selected component from Vars_SFS
    Vars_SFS <- Vars_SFS[-which(Vars_SFS==WhichVar)]
    EvolCorr <- c(EvolCorr,CorrSFS[[nbvars2select]][SelVar[1],])
  }
  parallel::stopCluster(cl)
  plan(sequential)
  return(list('SelectedVars' = SelectedVars, 
              'EvolCorr' = EvolCorr))
}


#' subprocess for forward sequential feature selection using a machine learning 
#' algorithm for regression task: tests all features for a given number of 
#' features already extracted
#' expects training and test dataset, 
#' can also save density scatterpolot for each combination of features tested   
#' @param algo character. currently RF or liquidSVM
#' @param NumVar_list list. list of variables
#' @param SelectedVars character. 
#' @param Vars_SFS character. 
#' @param TrainingData numeric. training dataset
#' @param TestData numeric. test dataset
#' @param TestData numeric. test dataset
#' @param ntree numeric. nb of trees for RF algo
#' @param mtry numeric. nb of attempts for RF algo
#' @param site character. name of site
#' @param figure_path_site character. path where to save scatterplot
#' @return list 
#' @importFrom terra rast values countNA
#' @export

SFS_Regression_step <- function(algo = 'RF', NumVar_list, SelectedVars, Vars_SFS, 
                                TrainingData, TestData, ntree = NULL, mtry = NULL, 
                                site = NULL, figure_path_site = NULL) {
  foreach(numvar = NumVar_list) %dopar% {
    SelFeat_tmp <- c(SelectedVars,Vars_SFS[[numvar]])
    Trainingsubset <- TrainingData[c('target',SelFeat_tmp)]
    Testsubset <- TestData[c('target',SelFeat_tmp)]
    # train random forest using full set of variables
    if (algo == 'RF'){
      formula <- target ~  .
      rf_model <- randomForest::randomForest(formula,
                                             data = Trainingsubset,
                                             ntree = ntree,
                                             mtry = mtry,
                                             importance = TRUE,
                                             do.trace = FALSE)
      
      # apply random forest on test dataset
      predictions <- predict(rf_model, newdata = Testsubset)
    }
    
    if (algo == 'liquidSVM'){
      svmModel <- liquidSVM::svmRegression(as.matrix(TrainingData[SelFeat_tmp]),
                                           as.matrix(TrainingData['target']))
      # apply SVR on test dataset
      predictions <- predict(svmModel, TestData[SelFeat_tmp])
    }
    # produce density scatterplot
    if (!is.na(figure_path_site)){
      filename <- paste0('LiDAR_LAI_model_', algo, '_vars_',paste(SelFeat_tmp, collapse = '_'))
      plot_density_scatterplot(c(predictions), Testsubset$target,
                               "Estimated LAI", "Measured LiDAR LAI",
                               title = paste("Estimated LAI vs LiDAR LAI:", site),
                               xlimits = c(0,15), ylimits = c(0,15), 
                               dirname = figure_path_site, filename = filename)
    }
    CorrVal <- cor.test(c(predictions),Testsubset$target)$estimate
    return(list('coeffcorr' = CorrVal, 
                'AssessedVal' = c(predictions)))
  }
}