# ---
# title: "1.calculate.pai_lidar.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-14"
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

source("libraries/functions_JBF.R")
source("libraries/functions_train_spatially.R")
source("libraries/functions_plots.R")

library("foreach")
library("doFuture")
library("doParallel")
library("future")
library("parallel")
library("terra")
library("randomForest")
library("caret")
library("MASS")
library("BlandAltmanLeh")
library("irr")
library("png")

# Setup
results_dir <- "../03_RESULTS"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual")
set.seed(123)

# Loop setup
sampling_sizes <- c(1000)
algos <- c('RF') #  'xGB', 'PolyRegr' , 'liquidSVM'
# algos <- c('RF')
study_name <- "mix_leaf_on_off" # mix_leaf_on_off # diff_LAI # leaf_on_only
i_test_zone <- c(1, 2, 3, 4)
if (study_name == "diff_LAI"){
  s2params <- c("s2lai")
} else {
  s2params <- c("s2lai", "evi", "ndvi", "savi")
}

# Define a pattern to exclude based on the current VI
excluded_metrics <- paste0(
  "/cv_res_10_m.tif$|",
  "/variance_res_10_m.tif$|", 
  "/shade_res_10_m.tif$|", 
  "/std_res_10_m.tif$|", 
  "/aspect_cos_res_10_m.tif$|",
  "/aspect_sin_res_10_m.tif$|",
  "/s2lai_res_10_m.tif$"
)
# excluded_metrics <- NULL

final_train_all_sites <- list()
final_test_sfs_all_sites <- list()
final_test_all_sites <- list()

for (sampling_size in sampling_sizes){
  for (ind in i_test_zone){
    for (s2param in s2params){
      for (site in sites){
        cat("Processing site:", site, "\n")
        
        metrics_files <- list.files(path = file.path(results_dir, site, metrics_dir), 
                                    pattern = "\\.tif$", full.names = TRUE)
        ranked_metrics_files <- rank_files(metrics_files, "lidarlai", "s2lai")
        
        # Apply the exclusion pattern to filter out unwanted files
        ranked_metrics_files <- ranked_metrics_files[!grepl(excluded_metrics, ranked_metrics_files)]
        
        ranked_metrics_files <- ranked_metrics_files[!grepl(
          paste(s2params[s2params != s2param], collapse = "|"), 
          ranked_metrics_files) | grepl(s2param, ranked_metrics_files)]
        
        # Spatial train-test
        zones <- create_spatial_train_test_values(ranked_metrics_files)
        n_zones <- length(zones)
        
        # Loop through each zone and rename the column
        for (i in 1:n_zones) {
          zones[[i]] <- zones[[i]][sample(1:nrow(zones[[i]]), sampling_size), ]
          names(zones[[i]])[names(zones[[i]]) == paste0(s2param, "_summer")] <- "s2param"
          names(zones[[i]])[names(zones[[i]]) == paste0(s2param, "_winter")] <- "s2param_winter"
          # s2lai_summer evi_summer ndvi_summer savi_summer
          # s2lai_winter evi_winter ndvi_winter savi_winter
        }
        
        # Train-test split
        # test_zone <- zones[[i]]
        # train_zones <- do.call(rbind, zones[-i])
        # ind <- 1
        test_zone <- zones[[ind]]
        train_zones <- do.call(rbind, zones[-ind])
        
        # Print zones
        cat("Train zones:", setdiff(1:n_zones, ind), "\n")
        cat("Test zone:", ind, "\n\n")
        
        # Discard winter S2 LAI or VI
        s2param_winter_train_zones <- train_zones$s2param_winter
        s2param_winter_test_zone <- test_zone$s2param_winter
        train_zones$s2param_winter <- NULL
        test_zone$s2param_winter <- NULL
        
        # Create duplicates and modify LAIs
        duplicate_train <- train_zones
        duplicate_test <- test_zone
        
        if ("lidarlai" %in% colnames(duplicate_train)) {
          duplicate_train$lidarlai <- rnorm(nrow(duplicate_train), mean=1, sd=0.25)
          duplicate_train$lidarlai[duplicate_train$lidarlai < 0.5] <- 0.5
          duplicate_train$lidarlai[duplicate_train$lidarlai > 1.5] <- 1.5
          
          duplicate_test$lidarlai <- rnorm(nrow(duplicate_test), mean=1, sd=0.25)
          duplicate_test$lidarlai[duplicate_test$lidarlai < 0.5] <- 0.5
          duplicate_test$lidarlai[duplicate_test$lidarlai > 1.5] <- 1.5
        }
        if ("s2param" %in% colnames(duplicate_train)) {
          # duplicate_train$s2lai <- rnorm(nrow(duplicate_train), mean=0.65, sd=0.125)
          # duplicate_train$s2lai[duplicate_train$s2lai < 0.4] <- 0.4
          # duplicate_train$s2lai[duplicate_train$s2lai > 0.9] <- 0.9
          # 
          # duplicate_test$s2lai <- rnorm(nrow(duplicate_test), mean=0.65, sd=0.13)
          # duplicate_test$s2lai[duplicate_test$s2lai < 0.4] <- 0.4
          # duplicate_test$s2lai[duplicate_test$s2lai > 0.9] <- 0.9
          
          duplicate_train$s2param <- s2param_winter_train_zones
          duplicate_test$s2param <- s2param_winter_test_zone
        }
        
        # Concatenate original and modified sets
        final_train_set <- rbind(train_zones, duplicate_train)
        final_train_set <- final_train_set[sample(1:nrow(final_train_set)), ]
        final_test_set_sfs <- rbind(test_zone, duplicate_test)
        final_test_set_sfs <- final_test_set_sfs[sample(1:nrow(final_test_set_sfs)), ]
        final_test_set <- test_zone
        final_test_set <- final_test_set[sample(1:nrow(final_test_set)), ]
        
        # Diff_LAI
        # final_train_set$diff_lai <- final_train_set$lidarlai - final_train_set$s2param
        # final_test_set_sfs$diff_lai <- final_test_set_sfs$lidarlai - final_test_set_sfs$s2param
        # final_test_set$diff_lai <- final_test_set$lidarlai - final_test_set$s2param
        
        # Remove vars if needed
        # final_train_set$s2param <- NULL
        # final_test_set_sfs$s2param <- NULL
        # final_test_set$s2param <- NULL
        # final_train_set$lidarlai <- NULL
        # final_test_set_sfs$lidarlai <- NULL
        # final_test_set$lidarlai <- NULL
        
        # Rename lidarlai or diff_lai to target
        colnames(final_train_set)[colnames(final_train_set) == "lidarlai"] <- "target"
        colnames(final_test_set)[colnames(final_test_set) == "lidarlai"] <- "target"
        colnames(final_test_set_sfs)[colnames(final_test_set_sfs) == "lidarlai"] <- "target"
        # colnames(final_train_set)[colnames(final_train_set) == "diff_lai"] <- "target"
        # colnames(final_test_set_sfs)[colnames(final_test_set_sfs) == "diff_lai"] <- "target"
        # colnames(final_test_set)[colnames(final_test_set) == "diff_lai"] <- "target"
        
        # Rename S2 LAI or S2 VI
        colnames(final_train_set)[colnames(final_train_set) == "s2param"] <- s2param
        colnames(final_test_set)[colnames(final_test_set) == "s2param"] <- s2param
        colnames(final_test_set_sfs)[colnames(final_test_set_sfs) == "s2param"] <- s2param
        
        # Move the target column to the first position
        final_train_set <- final_train_set[, c("target", setdiff(names(final_train_set), "target"))]
        final_test_set_sfs <- final_test_set_sfs[, c("target", setdiff(names(final_test_set_sfs), "target"))]
        final_test_set <- final_test_set[, c("target", setdiff(names(final_test_set), "target"))]
        
        # All sites
        final_train_all_sites[[site]] <- final_train_set
        final_test_all_sites[[site]] <- final_test_set
        final_test_sfs_all_sites[[site]] <- final_test_set_sfs
      }
      # a
      # Step 2: Train and evaluate models on each train and test set
      
      # Combine all sites' train and test sets for the current split
      combined_train_set <- do.call(rbind, 
                                    lapply(names(final_train_all_sites), 
                                           function(x)
                                             final_train_all_sites[[x]]))
      combined_test_set <- do.call(rbind, 
                                   lapply(names(final_test_all_sites), 
                                          function(x)
                                            final_test_all_sites[[x]]))
      combined_test_set_sfs <- do.call(rbind, 
                                       lapply(names(final_test_sfs_all_sites), 
                                              function(x)
                                                final_test_sfs_all_sites[[x]]))
      
      # Correlation matrix of train data
      train_correlation_matrix <- cor(combined_train_set)
      corrplot::corrplot(train_correlation_matrix,
                         method = "number",
                         type = "upper")
      
      # Step 1: Scale the response variable and predictors
      # scaled_combined_train_set <- scale(combined_train_set)
      scaled_combined_train_set <- combined_train_set %>%
        mutate(across(where(is.numeric), ~ (.-mean(.))/sd(.)))
      a
      # Step 2: Fit the multiple regression model on the scaled data
      multiple_regression_model <- lm(scaled_combined_train_set$target ~ ., 
                                      data = scaled_combined_train_set[, -which(names(scaled_combined_train_set) == "target")])
      
      # Step 3: Predict the scaled lidar LAI
      scaled_predicted_lidar_lai <- predict(multiple_regression_model, newdata = scaled_predictors)
      
      # Step 4: Unscale the predicted values
      # Extract the mean and standard deviation of the original lidar_lai_values
      original_mean <- attr(scaled_lidar_lai, "scaled:center")
      original_sd <- attr(scaled_lidar_lai, "scaled:scale")
      
      # Unscale the predicted values
      unscaled_predicted_lidar_lai <- scaled_predicted_lidar_lai * original_sd + original_mean
      

      # Convert variable importance to a data frame for plotting
      importance_df <- as.data.frame(importance(rf_model))
      importance_df$Variable <- rownames(importance_df)

      # Plot variable importance
      plot(ggplot(importance_df, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
             geom_bar(stat = "identity", fill = "steelblue") +
             coord_flip() +
             labs(title = paste("Variable Importance for Test Zone", ind),
                  x = "Variable",
                  y = "Importance (IncNodePurity)") +
             theme_minimal() +
             theme(
               axis.title.x = element_text(size = 14, face = "bold"),
               axis.title.y = element_text(size = 14, face = "bold"),
               axis.text = element_text(size = 12),
               plot.title = element_text(size = 16, face = "bold"),
               plot.subtitle = element_text(size = 14),
               panel.grid.major = element_line(color = "gray", size = 0.5),
               panel.grid.minor = element_line(color = "gray", size = 0.25)
             )
      )

      # Evaluate on the corresponding test set
      predictions <- predict(rf_model, newdata = combined_test_set)
      correlation <- cor(predictions, combined_test_set$target)
      rsquared <- caret::R2(predictions, combined_test_set$target)
      rmse <- Metrics::rmse(predictions, combined_test_set$target)
      bias <- Metrics::bias(predictions, combined_test_set$target)
      mae <- Metrics::mae(predictions, combined_test_set$target)

      cat("Metrics for combined test set (split", ind, "):",
          "\nCorrelation:", correlation,
          "\nR-squared:", rsquared,
          "\nRoot Mean Squared Error:", rmse,
          "\nBias:", bias,
          "\nMAE:", mae, "\n\n")

      # Evaluate on each individual site's test set
      for (site in sites){
        cat("Evaluating model on test set for site:", site, "\n")
        site_test_set <- do.call(rbind,
                                 lapply(names(final_test_all_sites),
                                        function(x) if (
                                          grepl(paste0(site,
                                                       "_", i, "$"), x))
                                          final_test_all_sites[[x]]
                                        else NULL))
        site_test_set <- final_test_all_sites[[site]]

        if (nrow(site_test_set) > 0) {
          predictions <- predict(rf_model, newdata = site_test_set)
          correlation <- cor(predictions, site_test_set$target)
          rsquared <- caret::R2(predictions, site_test_set$target)
          rmse <- Metrics::rmse(predictions, site_test_set$target)
          bias <- Metrics::bias(predictions, site_test_set$target)
          mae <- Metrics::mae(predictions, site_test_set$target)

          cat("Metrics for site", site, "test set (split", ind, "):",
              "\nCorrelation:", correlation,
              "\nR-squared:", rsquared,
              "\nRoot Mean Squared Error:", rmse,
              "\nBias:", bias,
              "\nMAE:", mae, "\n\n")

          # Plot density scatterplot
          plot_density_scatterplot(predictions,
                                   site_test_set$target,
                                   "Estimated LiDAR LAI",
                                   "Measured LiDAR LAI",
                                   title = paste0("LiDAR LAI: ",
                                                  site, " Test Zone ", ind)
          )
        }
      }
    }
  }
}

stop("End")





# Define features to retain
features_to_retain <- c("target", 
                        "lcv",
                        # "cv",
                        # "std",
                        # "max",
                        "mean",
                        # "vci",
                        "s2lai"
                        # "lskew"
                        # "rumple"
)

# Create an empty list to store filtered DataFrames
filtered_train_data_list <- filtered_test_data_list <- list()

# Loop through each element (DataFrame) in the list
for (i in 1:length(final_train_all_sites)) {
  # Select desired features from current DataFrame
  filtered_train_df <- final_train_all_sites[[i]][features_to_retain]
  filtered_test_df <- final_test_all_sites[[i]][features_to_retain]
  
  # Add the filtered DataFrame to the list
  filtered_train_data_list[[i]] <- filtered_train_df
  filtered_test_data_list[[i]] <- filtered_test_df
}

# Names
names(filtered_train_data_list) <- names(final_train_all_sites)
names(filtered_test_data_list) <- names(final_test_all_sites)

# Final: Train on full sites and evaluate models on each site, with SFS Features

# Combine all sites' train and test sets for the current split
combined_train_set <- do.call(rbind, 
                              lapply(names(filtered_train_data_list), 
                                     function(x) if (
                                       grepl(paste0("_", i, "$"), x)) 
                                       filtered_train_data_list[[x]] else NULL))

# Correlation matrix of train data
train_correlation_matrix <- cor(combined_train_set)
corrplot::corrplot(train_correlation_matrix, 
                   method = "number", 
                   type = "upper")

# Train Random Forest model
rf_model <- randomForest(target ~ ., 
                         data = combined_train_set, 
                         ntree = 30, 
                         importance = TRUE, 
                         do.trace = TRUE)

# Convert variable importance to a data frame for plotting
importance_df <- as.data.frame(importance(rf_model))
importance_df$Variable <- rownames(importance_df)

# Plot variable importance
plot(ggplot(importance_df, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
       geom_bar(stat = "identity", fill = "steelblue") +
       coord_flip() +
       labs(title = paste("Variable Importance for Test Zone", i),
            x = "Variable",
            y = "Importance (IncNodePurity)") +
       theme_minimal() +
       theme(
         axis.title.x = element_text(size = 14, face = "bold"),
         axis.title.y = element_text(size = 14, face = "bold"),
         axis.text = element_text(size = 12),
         plot.title = element_text(size = 16, face = "bold"),
         plot.subtitle = element_text(size = 14),
         panel.grid.major = element_line(color = "gray", size = 0.5),
         panel.grid.minor = element_line(color = "gray", size = 0.25)
       )
)

# Evaluate on each individual site's test set
for (site in sites){
  cat("Evaluating model on test set for site:", site, "\n")
  site_test_set <- do.call(rbind, 
                           lapply(names(filtered_test_data_list), 
                                  function(x) if (
                                    grepl(paste0(site,
                                                 "_", i, "$"), x)) 
                                    filtered_test_data_list[[x]] 
                                  else NULL))
  
  if (nrow(site_test_set) > 0) {
    predictions <- predict(rf_model, newdata = site_test_set)
    correlation <- cor(predictions, site_test_set$target)
    rsquared <- caret::R2(predictions, site_test_set$target)
    rmse <- Metrics::rmse(predictions, site_test_set$target)
    bias <- Metrics::bias(predictions, site_test_set$target)
    mae <- Metrics::mae(predictions, site_test_set$target)
    
    cat("Metrics for site", site, "test set (split", i, "):",
        "\nCorrelation:", correlation, 
        "\nR-squared:", rsquared, 
        "\nRoot Mean Squared Error:", rmse,
        "\nBias:", bias, 
        "\nMAE:", mae, "\n\n")
    
    # Plot density scatterplot
    plot_density_scatterplot(predictions,
                             site_test_set$target,
                             "Estimated LiDAR LAI",
                             "Measured LiDAR LAI",
                             title = paste0("LiDAR LAI: ", site, " Split ", i)
                             # dirname = "../04_FIGURES/winter",
                             # filename = paste0(site, "_winter")
                             # xlimits = c(0, 10),
                             # ylimits = c(-2.5, 10)
    )
  }
}

# Forcing an order to plot EvolCorr and EvolRMSE
study_name <- "mix_leaf_on_off" # mix_leaf_on_off # diff_LAI_on # leaf_on_only
figure_path_force <- file.path('../04_FIGURES', 
                               "ML_Figures",
                               "Stability",
                               "FS",
                               "RF",
                               study_name)

# If leaf-on
# ordered_vars <- c("lcv",
#                   "mean",
#                   "vci",
#                   "lskew",
#                   "vdr",
#                   "rumple",
#                   "max",
#                   "std",
#                   "shade",
#                   "s2lai",
#                   "variance",
#                   "cv")
# If leaf-on + leaf-off
ordered_vars <- c("s2lai",
                  "lcv",
                  "max",
                  "vci",
                  "mean",
                  "lskew",
                  "vdr",
                  "rumple",
                  "variance",
                  "std",
                  "shade",
                  "cv")
all_EvolCorr <- vector("list", length(ordered_vars))
all_EvolRMSE <- vector("list", length(ordered_vars))

for (num_vars in 1:length(ordered_vars)) {
  cat("Training with", num_vars, "variables\n")
  
  EvolCorr <- numeric(n_zones)
  EvolRMSE <- numeric(n_zones)
  
  # Combine all sites' train and test sets for the current split
  combined_train_set <- do.call(rbind, 
                                lapply(names(final_train_all_sites), 
                                       function(x) if (grepl(paste0("_", i, "$"), x)) 
                                         final_train_all_sites[[x]] else NULL))
  combined_test_set <- do.call(rbind, 
                               lapply(names(final_test_all_sites), 
                                      function(x) if (grepl(paste0("_", i, "$"), x)) 
                                        final_test_all_sites[[x]] else NULL))
  
  # Subset the data to include only the selected number of variables
  selected_vars <- ordered_vars[1:num_vars]
  combined_train_set_subset <- combined_train_set[, c(selected_vars, "target")]
  combined_test_set_subset <- combined_test_set[, c(selected_vars, "target")]
  
  # Train Random Forest model
  rf_model <- randomForest(target ~ ., 
                           data = combined_train_set_subset, 
                           ntree = 30, 
                           importance = TRUE, 
                           do.trace = TRUE)
  
  # Evaluate on the corresponding test set
  predictions <- predict(rf_model, newdata = combined_test_set_subset)
  correlation <- cor(predictions, combined_test_set_subset$target)
  rmse <- sqrt(mean((predictions - combined_test_set_subset$target)^2))
  
  
  # Store the metrics
  EvolCorr[i] <- correlation
  EvolRMSE[i] <- rmse
}

# Store the correlation and RMSE for the current number of variables
all_EvolCorr[[num_vars]] <- EvolCorr
all_EvolRMSE[[num_vars]] <- EvolRMSE

# Calculate mean correlation and RMSE across splits
mean_EvolCorr <- sapply(all_EvolCorr, mean)
mean_EvolRMSE <- sapply(all_EvolRMSE, mean)

# Prepare data for plotting
data_plot <- data.frame(
  NumVars = seq_along(ordered_vars),
  MeanEvolCorr = mean_EvolCorr,
  MeanEvolRMSE = mean_EvolRMSE
)

# Plot the mean evolution of correlation and RMSE
p <- ggplot(data_plot, aes(x = NumVars)) +
  geom_line(aes(y = MeanEvolCorr, color = 'Correlation')) +
  geom_line(aes(y = MeanEvolRMSE, color = 'RMSE')) +
  scale_y_continuous(
    name = "Correlation",
    sec.axis = sec_axis(~., name = "RMSE")
  ) +
  labs(title = "Mean Evolution of Correlation and RMSE in Mixed Leaf-On and Leaf-Off Period",
       x = "Number of Variables") +
  theme_bw() +
  scale_color_manual(name = "Metric", values = c("Correlation" = "blue", "RMSE" = "red"))

ggsave(filename = file.path(figure_path_force, "Mean_Evolution_Corr_RMSE.png"), plot = p, width = 8, height = 6)

cat("Plot saved to", file.path(figure_path_force, "Mean_Evolution_Corr_RMSE.png"), "\n")


















# Correlation matrix of train data
# train_correlation_matrix <- cor(train_predictors)
# corrplot::corrplot(train_correlation_matrix, method = "number", type = "upper")
# 
# # Correlation matrix of test data
# test_correlation_matrix <- cor(test_predicted)
# corrplot::corrplot(test_correlation_matrix, method = "number", type = "upper")

# Identify the dependent variable (assumed to be the first column for this example)
# names(all_train_data)[1] <- "target"
# dependent_variable <- names(all_train_data)[1]  # or the specific name of your dependent variable column
# 
# # Create the formula
# predictor_variables <- names(all_train_data)[-1]  # Exclude the dependent variable from predictors
# formula <- as.formula(paste(dependent_variable, "~", paste(predictor_variables, collapse = " + ")))
# 
# # Define initial linear model!
# initial_lm <- lm(formula, data = all_train_data)
# # summary(initial_lm)
# 
# rf_model <- randomForest(formula,
#                          data = all_train_data,
#                          ntree = 30,
#                          importance = TRUE,
#                          do.trace = TRUE)
# 
# # Convert variable importance to a data frame for plotting
# importance_df <- as.data.frame(importance(rf_model))
# importance_df$Variable <- rownames(importance_df)
# 
# # Plot variable importance
# ggplot(importance_df, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Variable Importance",
#        x = "Variable",
#        y = "Importance (IncNodePurity)") +
#   theme_minimal()
# 
# # Perform stepwise model selection
# stepwise_lm <- stepAIC(initial_lm, direction = "both")
# # summary(stepwise_lm)

# Set up train control with cross-validation
# train_control <- trainControl(method = "cv", number = 5)

# Train the random forest model

# rf_model <- train(formula, 
#                   data = train_predictors, 
#                   method = "rf", 
#                   trControl = train_control,
#                   do.trace = TRUE
#                   )

# Apply random forest on test dataset
# predictions <- predict(rf_model, newdata = test_predicted)


# for (site in sites){
#   cat("Processing site:", site, "\n")
#   # names(test_data[[site]])[1] <- "target"
#   # testing <- test_data[[site]]
#   
#   names(test_data_final[[site]])[1] <- "target"
#   testing <- test_data_final[[site]]
#   
#   predictions <- predict(rf_model, newdata = testing)
#   
#   # Compute metrics
#   correlation <- cor(predictions, testing$target)
#   rsquared <- caret::R2(predictions, testing$target)
#   rmse <- Metrics::rmse(predictions, testing$target)
#   bias <- Metrics::bias(predictions, testing$target)
#   mae <- Metrics::mae(predictions, testing$target)
#   
#   # Print metrics
#   cat("Metrics:",
#       "\nCorrelation:", correlation, 
#       "\nR-squared:", rsquared, 
#       "\nRoot Mean Squared Error:", rmse,
#       "\nBias:", bias, 
#       "\nMAE:", mae, "\n")
#   
#   # Plot density scatterplot
#   plot_density_scatterplot(predictions,
#                            testing$target,
#                            "Estimated LAI Difference",
#                            "Measured LAI Difference",
#                            title = paste0("LAI Difference: ", site),
#                            dirname = "../04_FIGURES/winter",
#                            filename = paste0(site, "_winter")
#                            # xlimits = c(0, 10),
#                            # ylimits = c(-2.5, 10)
#   )
# }

# xlimits = c(-2.5, 10),
# ylimits = c(-2.5, 10)