# ---
# title: "main_variogram_60_plots.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-22"
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
library("vegan")
library("spdep")
library("blockCV")
library("cluster")
library("mgcv")
library("geojsonio")
library("drf")
library("MASS")
library("ranger")
library("glmnet")
library("reshape2")
library("VSURF")
source("libraries/functions_general_tools.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Blois" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/variogram"
forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)

# Lists
all_predictors_filtered <- all_predictors_excluded <- all_observations <- list()
moran_filtered_results <- moran_variogram_results <- all_predictors_variogram <- list()
moran_all_obs <- list()
# Load data for all sites
for (site in sites) {
  site_data <- load_metrics_for_site(site, results_dir, forest_composition)
  all_observations[[site]] <- site_data
  
  # Load the GeoJSON file for the current site
  geojson_file <- file.path(results_dir, site, "data_utm31n.geojson")
  geo_data <- geojson_read(geojson_file, what = "sp") # Read GeoJSON
  
  # Extract coordinates
  coord_x <- geo_data@data$coord_x_utm31n
  coord_y <- geo_data@data$coord_y_utm31n
  
  # Create a data frame of the coordinates
  coord_df <- data.frame(coord_x = coord_x, coord_y = coord_y)
  
  # Initialize an empty list to store the closest matches
  closest_matches <- list()
  
  # Loop through each coordinate in coord_df
  for (j in 1:nrow(coord_df)) {
    # Extract the current target coordinate
    target_coord <- coord_df[j, ]
    
    # Calculate the distances to all observations in site_data
    distances <- sqrt((site_data$lat - target_coord$coord_y)^2 + 
                        (site_data$lon - target_coord$coord_x)^2)
    
    # Find the index of the closest match
    closest_index <- which.min(distances)
    
    # Store the closest match coordinates along with the original data
    closest_matches[[j]] <- cbind(site_data[closest_index, ], target_coord)
  }
  
  # Combine closest matches into a data frame
  closest_data <- do.call(rbind, closest_matches)
  closest_data <- closest_data[, !colnames(closest_data) %in% c("coord_x", "coord_y")]
  
  # Variogram 500m
  points_data <- read.csv(file.path(results_dir, site, "points_500m.csv"))
  coord_x_vario <- points_data$x
  coord_y_vario <- points_data$y
  filter_coords_vario <- data.frame(coord_x = coord_x_vario, coord_y = coord_y_vario)
  
  # matched_data <- site_data %>%
  #   semi_join(filter_coords_vario, by = c("lon" = "lon", "lat" = "lat"))  
  
  # Initialize an empty list to store the closest matches
  closest_matches_vario <- list()
  
  # Loop through each coordinate in coord_df
  for (j in 1:nrow(filter_coords_vario)) {
    # Extract the current target coordinate
    target_coord_vario <- filter_coords_vario[j, ]
    
    # Calculate the distances to all observations in site_data
    distances_vario <- sqrt((site_data$lat - target_coord_vario$coord_y)^2 + 
                        (site_data$lon - target_coord_vario$coord_x)^2)
    
    # Find the index of the closest match
    closest_index_vario <- which.min(distances_vario)
    
    # Store the closest match coordinates along with the original data
    closest_matches_vario[[j]] <- cbind(site_data[closest_index_vario, ], target_coord_vario)
  }
  # Combine closest matches into a data frame
  closest_data_vario <- do.call(rbind, closest_matches_vario)
  closest_data_vario <- closest_data_vario[, !colnames(closest_data_vario) %in% c("coord_x", "coord_y")]
  
  # Store the matched data in the list
  all_predictors_variogram[[site]] <- closest_data_vario
  
  # Moran
  moran_data <- closest_data
  coordinates(moran_data) <- ~ lon + lat  # Adjust based on your actual column names
  proj4string(moran_data) <- CRS("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")  # Set the correct CRS
  
  k <- 4
  neighbors <- knearneigh(coordinates(moran_data), k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  moran_variable <- moran_data$lidar_lai

  moran_test <- moran.test(moran_variable, listw = weights)
  moran_filtered_results[[site]] <- moran_test

  cat("Moran's I results (Field) for", site, ":\n")
  print(moran_filtered_results[[site]])
  cat("\n")
  
  # Moran
  moran_data <- closest_data_vario
  coordinates(moran_data) <- ~ lon + lat  # Adjust based on your actual column names
  proj4string(moran_data) <- CRS("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")  # Set the correct CRS
  
  k <- 4
  neighbors <- knearneigh(coordinates(moran_data), k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  moran_variable <- moran_data$lidar_lai
  
  moran_test <- moran.test(moran_variable, listw = weights)
  moran_variogram_results[[site]] <- moran_test
  
  cat("Moran's I results (Variogram 500m) for", site, ":\n")
  print(moran_variogram_results[[site]])
  cat("\n")
  
  # Moran
  # moran_data <- all_observations[[site]][sample(nrow(all_observations[[site]]),
  #                                               20000), ]
  # coordinates(moran_data) <- ~ lon + lat  # Adjust based on your actual column names
  # proj4string(moran_data) <- CRS("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")  # Set the correct CRS
  # 
  # k <- 4
  # neighbors <- knearneigh(coordinates(moran_data), k = k)
  # weights <- nb2listw(knn2nb(neighbors), style = "W")
  # moran_variable <- moran_data$lidar_lai
  # 
  # moran_test <- moran.test(moran_variable, listw = weights)
  # moran_all_obs[[site]] <- moran_test
  # 
  # cat("Moran's I results all for", site, ":\n")
  # print(moran_all_obs[[site]])
  # cat("\n")
  
  # Store the filtered data frame in the list
  all_predictors_filtered[[site]] <- closest_data
  
  # Exclude rows in site_data that are also in closest_data
  combined_data <- rbind(site_data, closest_data)  # Combine both data frames
  unique_data <- combined_data[!duplicated(combined_data), ]
  
  # Filter out rows in site_data that are found in closest_data
  combined_predictors_excluded <- anti_join(site_data, closest_data, by = colnames(site_data))
  all_predictors_excluded[[site]] <- as.data.frame(combined_predictors_excluded)
}

# Combine all filtered data frames into one
combined_predictors_filtered <- do.call(rbind, all_predictors_filtered)
combined_predictors_excluded <- do.call(rbind, all_predictors_excluded)
combined_predictors_variogram <- do.call(rbind, all_predictors_variogram)
combined_predictors_all <- do.call(rbind, all_observations)

# ---------------------------------- RF ----------------------------------------
# 1:102 103:244 245:682
# sample_indices <- sample(1:nrow(combined_predictors_all), size = 500, replace = FALSE)
# train_data <- combined_predictors_all[sample_indices, -c(16:53)]
# test_data <- combined_predictors_all[-sample_indices, -c(16:53)]

# train_data <- combined_predictors_variogram[, -c(16:53)]
# test_data <- combined_predictors_excluded[, -c(16:53)]

train_data <- combined_predictors_filtered[1:60, -c(18:55)]
# train_data <- combined_predictors_variogram[103:244, -c(18:55)]
# train_data <- train_data[sample(nrow(train_data), 60), ]

train_data$deltaLAI <- train_data$lidar_lai - train_data$s2_lai
# test_data$deltaLAI <- test_data$lidar_lai - test_data$s2_lai

# train_data$lon <- NULL
# train_data$lat <- NULL
train_data$lidar_lai <- NULL
train_data$s2_lai <- NULL
# test_data$lon <- NULL
# test_data$lat <- NULL
# test_data$lidar_lai <- NULL
# test_data$s2_lai <- NULL

corrplot::corrplot(cor(train_data),
                   method = "number",
                   type = "upper")

formula <- deltaLAI ~ . #s2_lai lidar_lai deltaLAI
rf_model <- ranger(
  formula = formula,
  data = train_data,
  num.trees = 500,
  # mtry = floor((ncol(train_data) - 1) / 3),
  mtry= 8,
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
  # splitrule = "extratrees",
  # num.random.splits = 10,
  # sample.fraction = 0.2
)
oob_mse <- round(rf_model$prediction.error, 2)
r_squared <- round(rf_model$r.squared, 2)
varImp <- sort(rf_model$variable.importance, decreasing = TRUE)
barplot(varImp)

# Output results
cat("OOB MSE:", round(oob_mse, 2), "\n")
cat("OOB R-squared:", round(r_squared, 2), "\n")
a

# predictions_excluded <- predict(rf_model, data = test_data)
# actual_excluded <- test_data$lidar_lai # lidar_lai s2_lai deltaLAI
# mse_excluded <- mean((predictions_excluded$predictions - actual_excluded) ^ 2)
# rss <- sum((predictions_excluded$predictions - actual_excluded) ^ 2)
# tss <- sum((actual_excluded - mean(actual_excluded)) ^ 2)
# r_squared_excluded <- 1 - rss/tss

# Output results
# cat("MSE on excluded data:", round(mse_excluded, 2), "\n")
# cat("R-squared on excluded data:", round(r_squared_excluded, 2), "\n")

# ---------------------------Sampling: all sites -------------------------------
show_excluded <- combined_predictors_all
show_filtered <- combined_predictors_filtered
show_variogram <- combined_predictors_variogram
show_all_sampled <- combined_predictors_all[sample(nrow(combined_predictors_all),
                                                   500), ]
show_excluded$Source <- 'All data'
show_filtered$Source <- 'Field'
show_variogram$Source <- 'Variogram 500m'
show_all_sampled$Source <- 'Sampled (500 obs)'

# Combine the two data frames
combined_data <- rbind(
                       #show_excluded, 
                       show_filtered,
                       show_variogram,
                       show_all_sampled)

# Melt the data frame for easier plotting
melted_data <- melt(combined_data, id.vars = 'Source', measure.vars = c('lcv', 'mean', 'lskew', 'gap_fraction'))

# Create boxplots
plot <- ggplot(melted_data, aes(x = value, fill = Source)) + 
  geom_histogram(bins = 30, alpha = 0.7, position = 'dodge') +  # Overlay histograms with transparency
  facet_wrap(~ variable, scales = 'free_x') +  # Facet for each variable, free x-axis scaling
  labs(title = "All Sites: Histogram Comparison of lcv, mean, lskew, max",
       x = "Value",
       y = "Frequency") +
  theme_bw() +
  theme(legend.position = "right")  # Keep the legend visible
ggsave(filename = paste0(output_dir, "/histogram_all.png"), plot = plot, width = 10, height = 6, dpi = 300)

# plot <- ggplot(melted_data, aes(x = Source, y = value, fill = Source)) +
#   geom_histogram() +
#   facet_wrap(~ variable, scales = 'free_y') +  # Create a facet for each variable
#   labs(title = "All sites Comparison of Variables: lcv, mean, lskew, max",
#        x = "Data Source",
#        y = "Value") +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave(filename = paste0(output_dir, "/boxplot_all.png"), plot = plot, width = 10, height = 6, dpi = 300)

# -------------------------Sampling: site per site -----------------------------
for (site in sites){
  show_excluded <- combined_predictors_all
  show_filtered <- combined_predictors_filtered
  show_variogram <- combined_predictors_variogram
  show_all_sampled <- combined_predictors_all[sample(nrow(combined_predictors_all),
                                                     500), ]
  show_excluded$Source <- 'All data'
  show_filtered$Source <- 'Field'
  show_variogram$Source <- 'Variogram 500m'
  show_all_sampled$Source <- 'Sampled (500 obs)'
  
  # Combine the two data frames
  combined_data <- rbind(show_excluded, 
                         show_filtered,
                         show_variogram,
                         show_all_sampled)
  
  # Melt the data frame for easier plotting
  melted_data <- melt(combined_data, id.vars = 'Source', measure.vars = c('lcv', 'mean', 'lskew', 'max'))
  
  # Create boxplots
  plot <- ggplot(melted_data, aes(x = Source, y = value, fill = Source)) +
    geom_boxplot() +
    facet_wrap(~ variable, scales = 'free_y') +  # Create a facet for each variable
    labs(title = paste(site, "Comparison of Variables: lcv, mean, lskew, max"),
         x = "Data Source",
         y = "Value") +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(filename = paste0(output_dir, "/boxplot_", site, ".png"), plot = plot, width = 10, height = 6, dpi = 300)
  # print(t.test(lcv ~ Source, data = rbind(show_excluded, show_filtered)))
  # print(t.test(lcv ~ Source, data = rbind(show_excluded, show_all_sampled)))
}

# ------------------------------------------------------------------------------
# formula <- lidar_lai ~ .
# # Fit the Random Forest model using ranger
# rf_model <- ranger(
#   formula = formula,
#   data = combined_predictors_excluded,
#   num.trees = 500,
#   mtry = floor((ncol(combined_predictors_filtered) - 1) / 3),
#   importance = "permutation",
#   min.node.size = 5,
#   num.threads = 12
# )
# # You can compute OOB MSE or R-squared here if needed
# oob_mse_excluded <- round(rf_model$prediction.error, 2)
# r_squared_excluded <- round(rf_model$r.squared, 2)
# 
# cat("OOB MSE for Random Forest on Excluded Data:", oob_mse_excluded, "\n")
# cat("R-squared for Random Forest on Excluded Data:", r_squared_excluded, "\n")
# varImp <- sort(rf_model$variable.importance, decreasing = TRUE)
# barplot(varImp)
# 
# # Predictions on the excluded data
# predictions <- predict(rf_model, data = combined_predictors_excluded)
# 
# # Evaluate predictions
# mse_excluded <- mean((combined_predictors_excluded$lidar_lai - predictions$predictions)^2)
# r_squared_excluded <- 1 - (sum((combined_predictors_excluded$lidar_lai - predictions$predictions)^2) / 
#                              sum((combined_predictors_excluded$lidar_lai - mean(predictions$predictions))^2))
# 
# cat("Mean Squared Error on Excluded Data:", round(mse_excluded, 2), "\n")
# cat("R-squared on Excluded Data:", round(r_squared_excluded, 2), "\n")

# --------------------------------- DRF ----------------------------------------

# site <- "Mormal"
# X = all_predictors_filtered[[site]][, -1]
# Y = as.data.frame(all_predictors_filtered[[site]][, 1])

# combined_predictors_variogram combined_predictors_filtered 
# drf_data <- combined_predictors_excluded[sample(nrow(combined_predictors_excluded),
#                                                 500), ]
drf_data <- combined_predictors_filtered[121:180,]
drf_data$deltaLAI <- drf_data$lidar_lai - drf_data$s2_lai
drf_data <- drf_data[, c("deltaLAI", setdiff(names(drf_data), "deltaLAI"))]

drf_data$lon <- NULL
drf_data$lat <- NULL
drf_data$lidar_lai <- NULL
# drf_data$aspect <- NULL
X = drf_data[, -c(1:2, 18:55)] # 16:53 17:54
Y = as.data.frame(drf_data[, 1]); 
names(Y) <- "deltaLAI"; # lidar_lai deltaLAI s2_lai
# W = as.data.frame(drf_data[, 19:56])

# combined_predictors_filtered$lon <- NULL
# combined_predictors_filtered$lat <- NULL
# combined_predictors_filtered$lidar_lai <- NULL
# combined_predictors_filtered$aspect <- NULL

# corrplot::corrplot(cor(W),
#                    method = "number",
#                    type = "upper")
# 
# heatmap(cor(W), 
#         main = "Correlation Matrix Heatmap", 
#         col = heat.colors(256), 
#         scale = "column", 
#         margins = c(10, 10))
# 
# # Reshape the correlation matrix into long format
# correlation_melted <- melt(cor(W))
# 
# # Create a heatmap
# ggplot(data = correlation_melted, aes(Var1, Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), name="Correlation") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = "Correlation Matrix Heatmap", x = "Variables", y = "Variables")

# floor((ncol(combined_predictors_filtered) - 1) / 3)
drf_modele <- drf(
  X = X,
  Y = Y,
  num.trees = 2000,
  # num.features = 5,
  # mtry = ncol(combined_predictors_filtered) - 1, 
  mtry = 5,
  honesty = FALSE,
  splitting.rule = "FourierMMD", # FourierMMD FastMMD
  # bandwidth = density_estimate$bw,
  num.threads = 12)
pred.oob <- predict(drf_modele, 
                    functional = "mean",
                    num.threads = 12)$mean
err.oob <- apply((pred.oob - Y)^2, 2, "mean")
varex <- 1 - err.oob/apply(Y, 2, "var")
cat("OOB MSE:", round(err.oob, 2), "\n")
cat("OOB R-squared:", round(varex, 2), "\n")

VI_MMD <- variableImportance(drf_modele)
names(VI_MMD) <- names(drf_data[, -c(1:2, 18:55)])
VI_MMD <- sort(VI_MMD, decreasing = T)
barplot(VI_MMD)

VI_count <- variable_importance(drf_modele)[,1]
names(VI_count) <- names(drf_data[, -c(1:2, 18:55)])
VI_count <- sort(VI_count, decreasing = T)

# par(mfrow=c(1,2))
# barplot(VI_MMD)
barplot(VI_count)
# par(mfrow=c(1,1))

# poids <- predict(drf_modele, newdata = combined_predictors_filtered[1, -c(1, 14:51)])$weights[1,]
poids <- predict(drf_modele)$weights[1,]

W_poly <- function(W,p){
  W_mat <- matrix(NA,length(W),p)
  for (k in 1:p){
    W_mat[,k] <- W^k
  }
  return(data.frame(W_mat))
}

# Modified W_poly function to handle data frames with multiple variables
Ws_poly <- function(W, p) {
  # Initialize an empty list to store polynomial features for each column
  W_poly_list <- list()
  
  # Loop over each column in W
  for (col in seq_along(W)) {
    W_col <- W[, col]  # Extract the current column
    W_mat <- matrix(NA, nrow = nrow(W), ncol = p)
    
    # Generate polynomial features for the current column
    for (k in 1:p) {
      W_mat[, k] <- W_col^k
    }
    
    # Store the result as a data frame and add to the list
    W_poly_list[[colnames(W)[col]]] <- data.frame(W_mat)
  }
  
  # Combine all polynomial features into one data frame
  W_poly_combined <- do.call(cbind, W_poly_list)
  return(W_poly_combined)
}

# Setup
# degree <- 2
# times <- 38
# depth <- seq(0,10,1)
# depth_poly <- W_poly(depth, degree)
# # depth_poly_repeated <- do.call(cbind, replicate(times, depth_poly, simplify = FALSE))
# newx <- as.matrix(depth_poly)
# 
# # s2_lai deltaLAI
# # 1 var
# s2_lai_poly5 <- W_poly(W$PAD_35.5, degree)
# s2_lai_fit_poly5 <- lm(Y$s2_lai~., data = s2_lai_poly5, weights = poids)
# 
# # Multiple vars
# # s2_lai_poly5s <- Ws_poly(W, degree)
# # s2_lai_fit_poly5s <- lm(Y$deltaLAI~., data = s2_lai_poly5s, weights = poids)
# 
# # Convert your response variable to a matrix
# Y_response <- as.vector(Y$s2_lai)  # Ensure Y is in the correct format
# X_response <- as.matrix(s2_lai_poly5)
# 
# # Perform cross-validation to find the best lambda
# set.seed(123)  # For reproducibility
# cv_lasso <- cv.glmnet(X_response, Y_response, alpha = 1, weights = poids)
# plot(cv_lasso)
# 
# # Get the best lambda value
# best_lambda <- cv_lasso$lambda.min
# cat("Best lambda value for Lasso:", best_lambda, "\n")
# 
# # Fit the final model using the best lambda
# final_lasso_model <- glmnet(X_response, Y_response, alpha = 1, lambda = best_lambda)
# print(final_lasso_model)
# 
# # Make predictions on new data if needed
# causal_effect_lasso <- predict(final_lasso_model, newx = newx)
# plot(depth, causal_effect_lasso, type = 'l', ylab = "delta LAI", xlab = "LiDAR PAD")





# degree <- 2
# depth <- seq(0, 10, 1)
# depth_poly <- W_poly(depth, degree)
# newx <- as.matrix(depth_poly)
# vars_to_model <- c("PAD_30.5", "PAD_25.5", "PAD_20.5",
#                    "PAD_15.5", "PAD_10.5", "PAD_5.5")
# colors <- rainbow(length(vars_to_model))
# plot(NULL, xlim=c(0, 10), ylim=c(2, 6), xlab="LiDAR PAD", ylab="Sentinel-2 LAI", type="n")
# for (i in seq_along(vars_to_model)) {
#   var_name <- vars_to_model[i]
#   s2_lai_poly <- W_poly(W[[var_name]], degree)
#   Y_response <- as.vector(Y$s2_lai)
#   X_response <- as.matrix(s2_lai_poly)
#   set.seed(123)
#   cv_lasso <- cv.glmnet(X_response, Y_response, alpha = 1, weights = poids)
#   best_lambda <- cv_lasso$lambda.min
#   final_lasso_model <- glmnet(X_response, Y_response, alpha = 1, lambda = best_lambda)
#   causal_effect_lasso <- predict(final_lasso_model, newx = newx)
#   
#   # s2_lai_fit_poly5 <- lm(Y$s2_lai~., data = s2_lai_poly, weights = poids)
#   # causal_effect_lasso <- predict(s2_lai_fit_poly5, depth_poly)
#   
#   lines(depth, causal_effect_lasso, col=colors[i], lty=i, lwd=2)
# }
# legend("bottomright", legend=vars_to_model, col=colors, lty=1:length(vars_to_model), lwd=2)


depth <- seq(0, 10, 1)
vars_to_model <- c("PAD_35.5", "PAD_30.5", "PAD_25.5", "PAD_20.5", "PAD_15.5", "PAD_10.5", "PAD_5.5")
colors <- rainbow(length(vars_to_model))
degrees <- 2:5
plot(NULL, xlim=c(0, 10), ylim=c(0, 6), xlab="LiDAR PAD", ylab="Sentinel-2 LAI", type="n")
for (i in seq_along(vars_to_model)) {
  var_name <- vars_to_model[i]
  best_degree <- NULL
  best_lambda <- NULL
  min_cv_error <- Inf
  final_lasso_model <- NULL
  
  for (degree in degrees) {
    s2_lai_poly <- W_poly(W[[var_name]], degree)
    
    Y_response <- as.vector(Y$s2_lai) # s2_lai deltaLAI
    X_response <- as.matrix(s2_lai_poly)
    
    set.seed(123)
    cv_lasso <- cv.glmnet(X_response, Y_response, alpha = 1, weights = poids)
    current_best_lambda <- cv_lasso$lambda.min
    current_min_cv_error <- min(cv_lasso$cvm)
    if (current_min_cv_error < min_cv_error) {
      best_degree <- degree
      best_lambda <- current_best_lambda
      min_cv_error <- current_min_cv_error
      final_lasso_model <- glmnet(X_response, Y_response, alpha = 1, lambda = best_lambda)
    }
  }
  
  depth_poly <- W_poly(depth, best_degree)
  newx <- as.matrix(depth_poly)
  causal_effect_lasso <- predict(final_lasso_model, newx = newx)
  
  # Plot the resulting curve using the color palette
  lines(depth, causal_effect_lasso, col=colors[i], lty=i, lwd=2)
  cat("Variable:", var_name, "- Best Degree:", best_degree, "- Best Lambda:", best_lambda, "- Min CV Error:", min_cv_error, "\n")
}
legend("bottomright", legend=vars_to_model, col=colors, lty=1:length(vars_to_model), lwd=2)








# Rename the columns to reflect duplication
names <- names(s2_lai_fit_poly5s$model)
names <- names[-c(1, length(names))]
colnames(depth_poly_repeated) <- names

# causal_effect <- predict(s2_lai_fit_poly5, depth_poly)
causal_effect <- predict(s2_lai_fit_poly5s, depth_poly_repeated)

plot(depth, causal_effect,
     type='l',
     ylab = "Sentinel-2 LAI")

# --
# List of variables for which to compute causal effects (replace with your variable names)
variables <- c(
  # "PAD_2.5", 
  # "PAD_3.5", 
  # "PAD_4.5", 
  "PAD_5.5",
  # "PAD_6.5", "PAD_7.5", "PAD_8.5", "PAD_9.5",
  "PAD_10.5",
  # "PAD_11.5",
  # "PAD_12.5", "PAD_13.5",
  # "PAD_14.5", 
  "PAD_15.5", 
  # "PAD_16.5", 
  # "PAD_17.5",
  # "PAD_18.5", "PAD_19.5", 
  "PAD_20.5", 
  # "PAD_21.5",
  # "PAD_22.5", "PAD_23.5", "PAD_24.5", 
  "PAD_25.5",
  # "PAD_26.5", "PAD_27.5", "PAD_28.5", "PAD_29.5",
  "PAD_30.5",
  # "PAD_31.5", "PAD_32.5", "PAD_33.5",
  # "PAD_34.5",
  "PAD_35.5"
  # "PAD_36.5", "PAD_37.5",
  # "PAD_38.5", "PAD_39.5"
)

# Create a list to store the causal effects for each variable
causal_effects_list <- list()
Depthh <- seq(0, 10, 1)
deg <- 2
Depth_polyy <- W_poly(Depthh, deg)

# Loop through each variable
for (var_name in variables) {
  var_poly5 <- W_poly(W[[var_name]], deg)
  fit_patientx <- lm(Y$deltaLAI ~ ., data = var_poly5, weights = poids)
  causal_effect_var <- predict(fit_patientx, Depth_polyy)
  causal_effects_list[[var_name]] <- causal_effect_var
}

# Plot the causal effects for all variables
plot(Depthh, causal_effects_list[[1]], type = 'l', ylab = "delta LAI", 
     xlab = "LiDAR PAD", col = "red", ylim = range(unlist(causal_effects_list)))

# Add lines for each variable's causal effect
for (i in seq_along(variables)) {
  lines(Depthh, causal_effects_list[[i]], col = i)
}

# Add a legend to distinguish the lines
legend("topright", legend = variables, col = seq_along(variables), lty = 1, cex = 0.7)


# modele.lin <- lm(Y$s2_lai~., data=drf_data[, 14:51])
# ## Pseudo R^2
# 1 - mean((modele.lin$fitted.values -
#             Y$s2_lai)^2)/var(Y$s2_lai)
# 
# modele.rf <- randomForest(X, Y$s2_lai)
# 1 - mean((predict(modele.rf) - Y$s2_lai)^2)/var(Y$s2_lai)
