# ---
# title: "functions_train_spatially.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-03-14"
# ---

# Function to create a raster from coordinates and values for each layer
create_raster_from_coords <- function(template, coords, values) {
  r <- terra::rast(template)
  cells <- terra::cellFromXY(r, coords)
  r[cells] <- values
  return(r)
}

# Function to subset values
subset_values <- function(values, non_na_cells, zone_start, zone_end) {
  if (ncol(values) == 1) {
    return(values[non_na_cells][zone_start:zone_end])
  } else {
    return(values[non_na_cells, ][zone_start:zone_end, ])
  }
}

create_spatial_train_test_values <- function(predictors_path) {
  # Read predicted value and predictors included in raster data
  predictors <- terra::rast(predictors_path)
  
  # Eliminate NA from all data when one layer has NA
  na_index_predictors <- terra::countNA(predictors)
  predictors[na_index_predictors > 0] <- NA
  
  # Extract values and eliminate NA
  predictors_val <- terra::values(predictors)
  
  # Identify non-NA cells
  non_na_cells <- rowSums(is.na(predictors_val)) == 0
  non_na_coords <- terra::xyFromCell(predictors, which(non_na_cells))
  
  # Calculate the total number of non-NA cells
  total_non_na <- sum(non_na_cells)
  quarter <- ceiling(total_non_na / 4)
  
  # Divide non-NA cells into four equal zones
  zone1_coords <- non_na_coords[1:quarter, ]
  zone2_coords <- non_na_coords[(quarter + 1):(2 * quarter), ]
  zone3_coords <- non_na_coords[(2 * quarter + 1):(3 * quarter), ]
  zone4_coords <- non_na_coords[(3 * quarter + 1):total_non_na, ]
  
  zone1_values <- subset_values(predictors_val, non_na_cells, 1, quarter)
  zone2_values <- subset_values(predictors_val, non_na_cells, quarter + 1, 2 * quarter)
  zone3_values <- subset_values(predictors_val, non_na_cells, 2 * quarter + 1, 3 * quarter)
  zone4_values <- subset_values(predictors_val, non_na_cells, 3 * quarter + 1, total_non_na)
  
  # Create rasters for each zone for predictors and predicted
  raster_zone1_predictors <- create_raster_from_coords(predictors, zone1_coords, zone1_values)
  raster_zone2_predictors <- create_raster_from_coords(predictors, zone2_coords, zone2_values)
  raster_zone3_predictors <- create_raster_from_coords(predictors, zone3_coords, zone3_values)
  raster_zone4_predictors <- create_raster_from_coords(predictors, zone4_coords, zone4_values)
  # writeRaster(x = raster_zone1_predictors[[1]], filename = "zone1.tif")
  # writeRaster(x = raster_zone2_predictors[[1]], filename = "zone2")
  # writeRaster(x = raster_zone3_predictors[[1]], filename = "zone3")
  # writeRaster(x = raster_zone4_predictors[[1]], filename = "zone4")
  
  # png(filename = "zone1.png", width = 800, height = 600)
  # plot(raster_zone1_predictors[[1]])
  # dev.off()
  # png(filename = "zone2.png", width = 800, height = 600)
  # plot(raster_zone2_predictors[[1]])
  # dev.off()
  # png(filename = "zone3.png", width = 800, height = 600)
  # plot(raster_zone3_predictors[[1]])
  # dev.off()
  # png(filename = "zone4.png", width = 800, height = 600)
  # plot(raster_zone4_predictors[[1]])
  # dev.off()
  # png(filename = "full.png", width = 800, height = 600)
  # plot(predictors[[1]])
  # dev.off()
  
  # Create data frames and set column names
  zone1_predictors_df <- as.data.frame(zone1_values)
  zone2_predictors_df <- as.data.frame(zone2_values)
  zone3_predictors_df <- as.data.frame(zone3_values)
  zone4_predictors_df <- as.data.frame(zone4_values)
  
  # Extract variable names from file paths and set as column names
  variable_names <- gsub("_res_10_m.tif", "", basename(predictors_path))
  colnames(zone1_predictors_df) <- variable_names
  colnames(zone2_predictors_df) <- variable_names
  colnames(zone3_predictors_df) <- variable_names
  colnames(zone4_predictors_df) <- variable_names
  
  return(list(
    'zone1_predictors' = zone1_predictors_df,
    'zone2_predictors' = zone2_predictors_df,
    'zone3_predictors' = zone3_predictors_df,
    'zone4_predictors' = zone4_predictors_df
  ))
}

compute_range_correlation <- function(raster1, raster2){
  a <- 1
}

index_of_agreement <- function(observed, predicted) {
  numerator <- sum((predicted - observed)^2)
  denominator <- sum((abs(predicted - mean(observed)) + abs(observed - mean(observed)))^2)
  ioa <- 1 - (numerator / denominator)
  return(ioa)
}

# Function to filter data based on specific features
filter_features <- function(data, features) {
  # Check if 'target' column exists and retain it
  if ("target" %in% colnames(data)) {
    target_col <- data$target
    # Filter data for only the specified features
    filtered_data <- data[, colnames(data) %in% features]
    # Add back the target column
    filtered_data <- cbind(filtered_data, target = target_col)
  } else {
    # Filter data for only the specified features
    filtered_data <- data[, colnames(data) %in% features]
  }
  return(filtered_data)
}


























# create_spatial_train_test_raster <- function(predicted_path, predictors_path) {
#   # Read predicted value and predictors included in raster data
#   predicted <- terra::rast(predicted_path) 
#   predictors <- terra::rast(predictors_path)
#   
#   # Eliminate NA from all data when one layer has NA
#   na_index_predictors <- terra::countNA(predictors)
#   na_index_predicted <- terra::countNA(predicted)
#   na_index <- na_index_predicted + na_index_predictors
#   predictors[na_index > 0] <- NA
#   predicted[na_index > 0] <- NA
#   
#   # Extract values and eliminate NA
#   predicted_val <- terra::values(predicted)
#   predictors_val <- terra::values(predictors)
#   
#   # Identify non-NA cells
#   non_na_cells <- rowSums(is.na(predicted_val)) == 0
#   non_na_coords <- terra::xyFromCell(predicted, which(non_na_cells))
#   
#   # Calculate the total number of non-NA cells
#   total_non_na <- sum(non_na_cells)
#   quarter <- ceiling(total_non_na / 4)
#   
#   # Divide non-NA cells into four equal zones
#   zone1_coords <- non_na_coords[1:quarter, ]
#   zone2_coords <- non_na_coords[(quarter + 1):(2 * quarter), ]
#   zone3_coords <- non_na_coords[(2 * quarter + 1):(3 * quarter), ]
#   zone4_coords <- non_na_coords[(3 * quarter + 1):total_non_na, ]
#   
#   zone1_values <- subset_values(predictors_val, non_na_cells, 1, quarter)
#   zone2_values <- subset_values(predictors_val, non_na_cells, quarter + 1, 2 * quarter)
#   zone3_values <- subset_values(predictors_val, non_na_cells, 2 * quarter + 1, 3 * quarter)
#   zone4_values <- subset_values(predictors_val, non_na_cells, 3 * quarter + 1, total_non_na)
#   
#   zone1_predicted_values <- subset_values(predicted_val, non_na_cells, 1, quarter)
#   zone2_predicted_values <- subset_values(predicted_val, non_na_cells, quarter + 1, 2 * quarter)
#   zone3_predicted_values <- subset_values(predicted_val, non_na_cells, 2 * quarter + 1, 3 * quarter)
#   zone4_predicted_values <- subset_values(predicted_val, non_na_cells, 3 * quarter + 1, total_non_na)
#   
#   # Create rasters for each zone for predictors and predicted
#   raster_zone1_predictors <- create_raster_from_coords(predictors, zone1_coords, zone1_values)
#   raster_zone2_predictors <- create_raster_from_coords(predictors, zone2_coords, zone2_values)
#   raster_zone3_predictors <- create_raster_from_coords(predictors, zone3_coords, zone3_values)
#   raster_zone4_predictors <- create_raster_from_coords(predictors, zone4_coords, zone4_values)
#   
#   raster_zone1_predicted <- create_raster_from_coords(predicted, zone1_coords, zone1_predicted_values)
#   raster_zone2_predicted <- create_raster_from_coords(predicted, zone2_coords, zone2_predicted_values)
#   raster_zone3_predicted <- create_raster_from_coords(predicted, zone3_coords, zone3_predicted_values)
#   raster_zone4_predicted <- create_raster_from_coords(predicted, zone4_coords, zone4_predicted_values)
#   
#   # Combine
#   combined_predictors_values <- rbind(zone1_values, zone2_values, zone3_values)
#   combined_predicted_values  <- as.vector(rbind(zone1_predicted_values, 
#                                                 zone2_predicted_values, 
#                                                 zone3_predicted_values))
#   
#   # Create data frames and set column names
#   combined_predictors_df <- as.data.frame(combined_predictors_values)
#   combined_predicted_df <- as.data.frame(combined_predicted_values)
#   zone4_predictors_df <- as.data.frame(zone4_values)
#   zone4_predicted_df <- as.data.frame(zone4_predicted_values)
#   
#   # Extract variable names from file paths and set as column names
#   variable_names <- gsub("_res_10_m.tif", "", basename(predictors_path))
#   colnames(combined_predictors_df) <- variable_names
#   names(zone4_predicted_df) <- "lidarlai"
#   
#   return(list(
#     'combined_predictors' = combined_predictors_df,
#     'combined_predicted' = combined_predicted_df,
#     'zone4_predictors' = zone4_predictors_df,
#     'zone4_predicted' = zone4_predicted_df
#   ))
# }
