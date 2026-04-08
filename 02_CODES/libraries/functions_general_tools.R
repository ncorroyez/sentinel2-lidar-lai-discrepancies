# ---
# title: "functions_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD,
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"
# ---
library(e1071)
# Function to format the computation time
format_time <- function(start, end) {
  elapsed <- as.numeric(difftime(end, start, units = "secs"))
  hours <- floor(elapsed / 3600)
  minutes <- floor((elapsed %% 3600) / 60)
  seconds <- elapsed %% 60
  
  sprintf("Elapsed time: %02d hours, %02d minutes, %.2f seconds", hours, minutes, seconds)
}

# Function to print raster size
print_raster_size <- function(raster) {
  size <- dim(raster)
  sprintf("%dx%d (rows x cols)", size[1], size[2])
}

# Function to extract all dates for a given site
get_site_dates <- function(site_name, data) {
  site_row <- data[data$Site == site_name, ]
  dates <- unlist(site_row[, -1], use.names = FALSE)
  dates <- dates[!is.na(dates)]
  return(as.Date(dates, origin = "1970-01-01"))
}

# get_site_full_dates <- function(site_name, data) {
#   site_row <- data[data$Site == site_name, -1, drop = FALSE]
#   raw_strings <- unlist(site_row, use.names = FALSE)
#   raw_strings <- raw_strings[!is.na(raw_strings) & raw_strings != ""]
#   print(raw_strings)
#   all_dates <- unlist(strsplit(as.character(raw_strings), "\\s*\\+\\s*"))
#   print(all_dates)
#   dates_clean <- trimws(all_dates)
#   print(dates_clean)
#   return(sort(as.Date(dates_clean)))
# }

get_site_full_dates <- function(site_name, data) {
  
  site_row <- data[data$Site == site_name, -1, drop = FALSE]
  raw_strings <- unlist(site_row, use.names = FALSE)
  raw_strings <- raw_strings[!is.na(raw_strings) & raw_strings != ""]
  
  # split on +
  all_dates <- unlist(strsplit(as.character(raw_strings), "\\s*\\+\\s*"))
  dates_clean <- trimws(all_dates)
  
  # ---- SEPARATE ISO vs Excel numeric dates ----
  
  is_numeric_date <- grepl("^\\d+$", dates_clean)
  
  dates <- as.Date(rep(NA, length(dates_clean)))
  
  # ISO dates
  dates[!is_numeric_date] <- as.Date(dates_clean[!is_numeric_date])
  
  # Excel serial dates
  dates[is_numeric_date] <- as.Date(
    as.numeric(dates_clean[is_numeric_date]),
    origin = "1970-01-01"
  )
  
  sort(unique(dates))
}

copy_and_move_folder <- function(source, destination) {
  # Check if the source folder exists
  if (!dir.exists(source)) {
    stop("The source folder does not exist.")
  }
  
  # Create the destination folder if it doesn't already exist
  if (!dir.exists(destination)) {
    dir.create(destination, recursive = TRUE)
  }
  
  # Copy the source folder to the destination
  copy_success <- file.copy(from = source, to = destination, recursive = TRUE)
  
  # If the copy is successful, delete the source folder
  if (copy_success) {
    unlink(source, recursive = TRUE)  # Delete the source folder
    message("The folder has been successfully copied and deleted.")
  } else {
    stop("The folder copy failed.")
  }
}

# Define a custom function to rank files
rank_files <- function(file_list, primary_metric, secondary_metric) {
  # Extract filenames
  filenames <- basename(file_list)
  
  # Rank files containing primary_metric and secondary_metric
  primary_files <- filenames[grep(primary_metric, filenames)]
  secondary_files <- filenames[grep(secondary_metric, filenames)]
  other_files <- filenames[!(filenames %in% c(primary_files, secondary_files))]
  
  # Combine files into the desired order
  ordered_files <- c(primary_files, secondary_files, other_files)
  
  # Return the full paths of ordered files
  file_list[match(ordered_files, filenames)]
}

load_metrics_for_site <- function(site, results_dir, forest_composition) {
  # Define the metrics directory for the given site
  metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
  
  # Define the names of the metrics
  metric_files <- c("lidarlai_res_10_m.tif",
                    "s2lai_summer_res_10_m.tif",
                    "lcv_res_10_m.tif",
                    "mean_res_10_m.tif",
                    "max_res_10_m.tif",
                    "lskew_res_10_m.tif",
                    # "shade_res_10_m.tif",
                    "hillshade_res_10_m.tif",
                    "slope_res_10_m.tif",
                    "aspect_res_10_m.tif",
                    # "northness_res_10_m.tif",
                    # "northness_slope_res_10_m.tif",
                    # "twi_res_10_m.tif",
                    "std_res_10_m.tif",
                    "cv_res_10_m.tif",
                    "rumple_res_10_m.tif",
                    "vci_res_10_m.tif",
                    # "vdr_res_10_m.tif",
                    # "gap_fraction_res_10_m.tif",
                    "fCover_res_10_m.tif",
                    "ndvi_summer_res_10_m.tif",
                    # "dtm.tif",
                    "cv_lad_res_10_m.tif",
                    "cv_lad_dtm_res_10_m.tif",
                    "stand_type_raster.tif"
  )
  
  # Function to load and extract raster values
  load_metric <- function(filename) {
    rast <- terra::rast(file.path(metrics_dir, filename))
    values(rast)
  }
  
  # Load all the metrics in a list
  metrics_values <- lapply(metric_files, load_metric)
  
  # Assign meaningful names to the metrics
  names(metrics_values) <- c("lidar_lai",
                             "s2_lai",
                             "lcv",
                             "mean",
                             "max",
                             "lskew",
                             # "shade",
                             "hillshade",
                             "slope",
                             "aspect",
                             # "northness",
                             # "northness_slope",
                             # "twi",
                             "std",
                             "cv",
                             "rumple",
                             "vci",
                             # "vdr",
                             # "gap_fraction",
                             "fCover",
                             "ndvi",
                             # "dtm",
                             "cv_lad_dsm",
                             "cv_lad_dtm",
                             "stands"
  )
  # print(str(metrics_values))
  
  # PAD profiles
  pad_files <- sprintf("PAD_Profiles_updated_modifminz/PAD_%.1f_40.tif", seq(2.5, 39.5, by = 1))
  pad_metrics_values <- lapply(pad_files, load_metric)
  names(pad_metrics_values) <- sprintf("PAD_%.1f", seq(2.5, 39.5, by = 1))
  metrics_values <- c(metrics_values, pad_metrics_values)
  
  # Extract coordinates (assuming all rasters have the same structure)
  lidar_lai <- terra::rast(file.path(metrics_dir, "lidarlai_res_10_m.tif"))
  coordinates <- xyFromCell(lidar_lai, 1:ncell(lidar_lai))
  
  # Add longitude and latitude to the list of metrics
  metrics_values$lon <- coordinates[, 1]
  metrics_values$lat <- coordinates[, 2]
  # metrics_values$stands <- rast(file.path(metrics_dir, "stand_type_raster.tif"))
  
  # Filter out NA values across all metrics
  valid_idx <- complete.cases(do.call(cbind, metrics_values))
  
  # Filter metrics to retain only valid (non-NA) rows
  valid_metrics <- lapply(metrics_values, function(x) x[valid_idx])
  valid_metrics$stands <- factor(valid_metrics$stands,
                                 levels = 0:3,
                                 labels = c("land", "deciduous", "coniferous", "mixtures"))
  # Combine all predictor variables into a data frame
  predictors <- as.data.frame(valid_metrics)
  # predictors$site <- site
  
  return(predictors)
}

# load_metrics_with_pad <- function(site, 
#                                   results_dir, 
#                                   forest_composition,
#                                   num_depths = 1) {
#   # Define the metrics directory for the given site
#   metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
#   
#   # Define the names of the static metrics
#   metric_files <- c("lidarlai_res_10_m.tif",
#                     "s2lai_summer_res_10_m.tif",
#                     "lcv_res_10_m.tif",
#                     "mean_res_10_m.tif",
#                     "max_res_10_m.tif",
#                     "lskew_res_10_m.tif",
#                     "hillshade_res_10_m.tif",
#                     "slope_res_10_m.tif",
#                     "aspect_res_10_m.tif",
#                     "std_res_10_m.tif",
#                     "cv_res_10_m.tif",
#                     "rumple_res_10_m.tif",
#                     "vci_res_10_m.tif",
#                     "fCover_res_10_m.tif",
#                     "ndvi_summer_res_10_m.tif",
#                     "cv_lad_res_10_m.tif",
#                     "cv_lad_dtm_res_10_m.tif",
#                     "stand_type_raster.tif"
#   )
#   
#   # Function to load and extract raster values
#   load_metric <- function(filename) {
#     rast <- terra::rast(file.path(metrics_dir, filename))
#     values(rast)
#   }
#   
#   # Load all the metrics in a list
#   metrics_values <- lapply(metric_files, load_metric)
#   
#   # Assign meaningful names to the metrics
#   names(metrics_values) <- c("lidar_lai",
#                              "s2_lai",
#                              "lcv",
#                              "mean",
#                              "max",
#                              "lskew",
#                              "hillshade",
#                              "slope",
#                              "aspect",
#                              "std",
#                              "cv",
#                              "rumple",
#                              "vci",
#                              "fCover",
#                              "ndvi",
#                              "cv_lad_dsm",
#                              "cv_lad_dtm",
#                              "stands"
#   )
#   
#   # Load PAD profiles and process them for both dtm and dsm
#   pad_types <- c("dtm", "dsm")
#   # pad_types <- c("dtm")
#   pad_results <- lapply(pad_types, function(pad_type, site) {
#     set.seed(42)
#     if (pad_type == "dtm") {
#       pad_files <- sprintf("PAD_Profiles_NA/PAD_%.1f_40.tif", 
#                            seq(2.5, 39.5, by = 1))
#       pad_stack_file <- "PAD_Profiles_NA/ladstack.tif"
#     } else if (pad_type == "dsm") {
#       pad_files <- sprintf("PAD_Profiles_updated_modifminz_NA/PAD_%.1f_40.tif", 
#                            seq(2.5, 39.5, by = 1))
#       pad_stack_file <- "PAD_Profiles_updated_modifminz_NA/ladstack.tif"
#     }
#     # Load PAD metrics
#     pad_metrics_values <- lapply(pad_files, load_metric)
#     # valid_cells <- complete.cases(do.call(cbind, pad_metrics_values))
#     # pad_metrics_values <- lapply(pad_metrics_values, function(x) x[valid_cells])
#     
#     # Extract coordinates from one PAD raster (assume same structure for all rasters)
#     pad_rast <- terra::rast(file.path(metrics_dir, pad_files[1]))
#     coordinates <- xyFromCell(pad_rast, 1:ncell(pad_rast))
#     
#     # Load PAD Stack
#     lad_stack_values <- lapply(pad_stack_file, load_metric)
#     # valid_stack <- complete.cases(do.call(cbind, lad_stack_values))
#     pad_stack <- terra::rast(file.path(metrics_dir, pad_stack_file))
#     valid_stack <- !is.na(pad_stack) 
#     
#     # Function to sample one valid PAD value per pixel
#     sample_valid_pad <- function(pad_metrics_values, n_cells, num_depths) {
#       sampled_values <- matrix(NA, nrow = n_cells, ncol = num_depths)
#       sampled_layers <- matrix(NA, nrow = n_cells, ncol = num_depths)
#       pb <- progress_bar$new(
#         format = "Processing cells [:bar] :percent in :elapsed ETA: :eta",
#         total = n_cells,  # Total number of cells to process
#         clear = FALSE,
#         width = 60
#       )
#       row_index <- 1
#       for (cell in seq_len(n_cells)) {
#         pixel_values <- sapply(pad_metrics_values, function(layer) layer[cell])
#         valid_layers <- which(!is.na(pixel_values))
#         # print(valid_layers)
#         
#         # if (length(valid_layers) > 0) {
#         #   sampled_layer <- sample(valid_layers, size = 1)
#         #   sampled_values[cell] <- pixel_values[sampled_layer]
#         #   sampled_layers[cell] <- sampled_layer
#         # } else {
#         #   sampled_values[cell] <- NA
#         #   sampled_layers[cell] <- NA
#         # }
#         if (length(valid_layers) >= num_depths) {
#           sampled_layers[cell, ] <- sample(valid_layers, size = num_depths, replace = FALSE)
#           sampled_values[cell, ] <- pixel_values[sampled_layers[cell, ]]
#         } else if (length(valid_layers) > 0) {
#           sampled_layers[cell, ] <- NA
#           sampled_values[cell, ] <- NA
#         }
#         pb$tick()
#       }
#       list(values = sampled_values, layers = sampled_layers)
#     }
#     n_cells <- ncell(pad_rast)
#     pad_sampling_result <- sample_valid_pad(pad_metrics_values, n_cells, num_depths)
#     random_pad_values <- pad_sampling_result$values
#     random_pad_layers <- pad_sampling_result$layers
#     depth <- apply(random_pad_layers, 
#                    2, function(layer) ifelse(!is.na(layer), 
#                                              40 - (2 + (layer - 1) * 1), NA))
#     print(random_pad_values)
#     print(depth)
#     # depth <- ifelse(!is.na(random_pad_layers), 
#     #                 40 - (2 + (random_pad_layers - 1) * 1), NA)
#     
#     # Create random PAD selection per pixel
#     # set.seed(42)
#     # random_pad_indices <- sample(seq_along(pad_files), size = ncell(pad_rast), replace = TRUE)
#     # print(random_pad_indices)
#     # random_pad_values <- sapply(seq_len(ncell(pad_rast)), function(cell) {
#     #   pad_metrics_values[[random_pad_indices[cell]]][cell]
#     # })
#     
#     # Calculate PAD LiDAR - S2 LAI and add depth
#     s2_lai <- metrics_values[["s2_lai"]]
#     pad_minus_lai <- random_pad_values - s2_lai
#     # depth <- 40 - (2 + (random_pad_indices - 1) * 1)
#     
#     list(random_pad = random_pad_values,
#          depth = depth,
#          pad_minus_lai = pad_minus_lai,
#          coordinates = coordinates)
#   })
#   
#   # Add both dtm and dsm results to metrics
#   metrics_values$random_pad_dtm <- pad_results[[1]]$random_pad
#   metrics_values$depth_dtm <- pad_results[[1]]$depth
#   metrics_values$pad_minus_lai_dtm <- pad_results[[1]]$pad_minus_lai
#   
#   metrics_values$random_pad_dsm <- pad_results[[2]]$random_pad
#   metrics_values$depth_dsm <- pad_results[[2]]$depth
#   metrics_values$pad_minus_lai_dsm <- pad_results[[2]]$pad_minus_lai
#   
#   # Add longitude and latitude to the list of metrics
#   metrics_values$lon <- pad_results[[1]]$coordinates[, 1]
#   metrics_values$lat <- pad_results[[1]]$coordinates[, 2]
#   
#   # Filter out NA values across all metrics
#   valid_idx <- complete.cases(do.call(cbind, metrics_values))
#   valid_metrics <- lapply(metrics_values, function(x) x[valid_idx])
#   valid_metrics$stands <- factor(valid_metrics$stands,
#                                  levels = 0:3,
#                                  labels = c("land", "deciduous", "coniferous", "mixtures"))
#   predictors <- as.data.frame(valid_metrics)
#   predictors$deltaLAI <- predictors$lidar_lai - predictors$s2_lai
#   return(predictors)
# }

load_metrics_with_pad <- function(site, 
                                  results_dir, 
                                  forest_composition,
                                  num_depths = 1) {
  # Define the metrics directory for the given site
  metrics_dir <- file.path(results_dir, site, "Metrics", forest_composition)
  
  # Define the names of the static metrics
  metric_files <- c(
                    "lidarlai_res_10_m.tif",
                    "lidarlai_optim_depth_res_10_m.tif",
                    "lidarlai_under_res_10_m.tif",
                    "s2lai_summer_atbd_res_10_m.tif",
                    # "s2lai_summer_depth_study_common_res_10_m.tif",
                    "s2lai_summer_best_indiv_res_10_m.tif",
                    "lcv_res_10_m.tif",
                    # "lcv_no_ground.tif",
                    "mean_res_10_m.tif",
                    "max_res_10_m.tif",
                    "lskew_res_10_m.tif",
                    # "lskew_no_ground.tif",
                    "hillshade_res_10_m.tif",
                    "slope_res_10_m.tif",
                    "aspect_res_10_m.tif",
                    "std_res_10_m.tif",
                    "cv_res_10_m.tif",
                    "rumple_res_10_m.tif",
                    "vci_res_10_m.tif",
                    "fCover_res_10_m.tif",
                    "ndvi_summer_res_10_m.tif",
                    "cv_lad_res_10_m.tif",
                    "cv_lad_dtm_res_10_m.tif",
                    "stand_trees_raster.tif"
  )
  
  # Function to load and extract raster values
  load_metric <- function(filename, is_factor = FALSE, mapping = NULL) {
    rast <- terra::rast(file.path(metrics_dir, filename))
    
    if (is_factor && !is.null(levels)) {
      # Map numeric values to factor labels
      values <- terra::values(rast)
      factor(values, levels = as.numeric(mapping), labels = names(mapping))
    } else {
      # Return numeric/continuous raster values
      terra::values(rast)
    }
  }
  
  stand_type_mapping <- c(
    "oak" = 0,
    "deciduous" = 1,
    "beech" = 2,
    "poplar" = 3,
    "coniferous" = 4,
    "mixed" = 5,
    "douglas" = 6,
    "larch" = 7,
    "fir/spruce" = 8,
    "scots pine" = 9,
    "laricio pine/black pine" = 10,
    "mixed pines" = 11,
    "nc" = 12,
    "nr" = 13,
    "land" = 14
  )
  
  # Load all the metrics in a list
  metrics_values <- lapply(metric_files, function(file) {
    if (file == "stand_trees_raster.tif") {
      # Load the "stands" raster as a factor
      load_metric(file, is_factor = TRUE, mapping = stand_type_mapping)
    } else {
      # Load other metrics as numeric rasters
      load_metric(file)
    }
  })
  
  # Assign meaningful names to the metrics
  names(metrics_values) <- c("lidar_lai",
                             "lidar_lai_optimD",
                             "lidar_lai_under",
                             "s2_lai_atbd",
                             "s2_lai_common",
                             "lcv",
                             "mean",
                             "max",
                             "lskew",
                             "hillshade",
                             "slope",
                             "aspect",
                             "std",
                             "cv",
                             "rumple",
                             "vci",
                             "fCover",
                             "ndvi",
                             "cv_lad_dsm",
                             "cv_lad_dtm",
                             "stands"
  )
  
  # Load PAD profiles and process them for both dtm and dsm
  pad_types <- c("dtm", "dsm")
  pad_results <- lapply(pad_types, function(pad_type) {
    if (pad_type == "dtm") {
      pad_files <- sprintf("PAD_Profiles_NA/PAD_%.1f_40.tif", seq(2.5, 39.5, by = 1))
    } else if (pad_type == "dsm") {
      pad_files <- sprintf("PAD_Profiles_updated_modifminz_NA/PAD_%.1f_40.tif", seq(2.5, 39.5, by = 1))
    }
    # Load PAD metrics into a list
    pad_metrics <- lapply(pad_files, load_metric)
    names(pad_metrics) <- sprintf("PAD_%s_%.1f", pad_type, seq(2.5, 39.5, by = 1))
    pad_metrics
  })
  
  # Combine PAD results into metrics_values
  metrics_values <- c(metrics_values, unlist(pad_results, recursive = FALSE))
  
  # Add longitude and latitude to the list of metrics
  pad_rast <- terra::rast(file.path(metrics_dir, sprintf("PAD_Profiles_NA/PAD_2.5_40.tif")))
  coordinates <- xyFromCell(pad_rast, 1:ncell(pad_rast))
  metrics_values$lon <- coordinates[, 1]
  metrics_values$lat <- coordinates[, 2]
  
  # Filter out rows where lidar_lai is NA
  valid_idx <- (!is.na(metrics_values$lidar_lai) & !is.na(metrics_values$s2_lai_atbd)
                & !is.na(metrics_values$lcv) & !is.na(metrics_values$std)
                & !is.na(metrics_values$hillshade) & !is.na(metrics_values$slope))
  valid_metrics <- lapply(metrics_values, function(x) x[valid_idx])
  # valid_metrics$stands <- factor(valid_metrics$stands,
  #                                levels = 0:3,
  #                                labels = c("land", "deciduous", "coniferous", "mixtures"))
  
  # Convert to data frame and calculate deltaLAI
  predictors <- as.data.frame(valid_metrics)
  # predictors$deltaLAI <- predictors$lidar_lai - predictors$s2_lai_dsm
  predictors$deltaLAI_atbd <- predictors$lidar_lai - predictors$s2_lai_atbd
  predictors$deltaLAI_common <- predictors$lidar_lai - predictors$s2_lai_common
  
  # Get PAD variable names
  pad_dtm_vars <- grep("^PAD_dtm_", names(predictors), value = TRUE)
  pad_dsm_vars <- grep("^PAD_dsm_", names(predictors), value = TRUE)
  
  depths <- seq(1, 38, by = 1)
  
  generate_depth_weights <- function(depth_seq, row) {
    non_na_depths <- depth_seq[!is.na(row) & !is.nan(row)]
    # print(row)
    
    # if (length(non_na_depths) < 2) {
    #   return(rep(0, length(depth_seq)))
    # }
    
    max_depth <- max(non_na_depths)
    min_depth <- min(non_na_depths)
    
    # Calculate raw weights: closer to 1 gives low weight, closer to 38 gives high weight
    raw_weights <- ((non_na_depths - min_depth) / (max_depth - min_depth))^3
    raw_weights <- 1 - raw_weights  # Inverse so that closer to 1 has higher weight
    # if (all(raw_weights == 0)) {
    #   raw_weights <- rep(1e-6, length(raw_weights))  # Assign a very small value if all weights are zero
    # }
    # scaled_depths <- (non_na_depths - min_depth) / (max_depth - min_depth)
    # raw_weights <- exp(-scaled_depths * 3)
    # print(non_na_depths)
    # print(raw_weights)
    
    # Determine the frequency of each depth in the full dataset
    # depth_frequencies <- table(non_na_depths)  # Frequencies of each depth
    # freq_weights <- 1 / as.numeric(depth_frequencies)  # Inverse of frequency (rarer depths get higher weight)
    # names(freq_weights) <- names(depth_frequencies)    # Assign names from depth_frequencies
    # print(depth_frequencies)
    # raw_weights <- freq_weights[as.character(non_na_depths)]
    # print(raw_weights)
    
    # if (sum(raw_weights) == 0) {
    #   return(rep(1, length(depth_seq)))  # Return zero weights if raw weights sum to 0
    # }
    # raw_weights[is.na(raw_weights)] <- 0
    # Normalize weights so that their sum is 1
    normalized_weights <- raw_weights / sum(raw_weights, na.rm = TRUE)
    # print(normalized_weights)
    # normalized_weights <- 1 - normalized_weights
    
    # Match the size of full depth_seq, fill in 0 for NAs
    full_weights <- rep(0, length(depth_seq))
    depth_match <- match(non_na_depths, depth_seq)
    full_weights[depth_match] <- normalized_weights
    
    # print(full_weights)
    # if (length(non_na_depths) < 2) {
    #   full_weights <- rep(0, length(depth_seq))
    #   full_weights[1] <- 1
    # }
    
    # print(full_weights)
    return(full_weights)
  }
  
  select_random_pad_and_depth <- function(row, depth_seq) {
    depth_weights <- generate_depth_weights(depth_seq, row)
    non_na_indices <- which(!is.na(row) & !is.null(row))
    if (length(non_na_indices) > 0) {
      # Ensuring that the probabilities correspond correctly to the non-NA indices
      depth_indices <- depth_seq[non_na_indices]  # Get corresponding depth indices
      prob_weights <- depth_weights[match(depth_indices, depth_seq)]  # Match with depth_weights
      if (length(non_na_indices) == 1) {
        prob_weights <- NULL
      }
      # Sample based on the depth weights
      selected_idx <- sample(non_na_indices, 1, prob = prob_weights)
      # selected_idx <- sample(non_na_indices, 1)
      random_pad <- row[selected_idx]
      depth <- 40 - (2 + (depth_seq[selected_idx]) - 1)
      return(c(random_pad, depth))
    } else {
      return(c(NA, NA))
    }
  }
  
  # Apply the function to dtm and dsm
  # set.seed(42)
  # dtm_results <- t(apply(predictors[pad_dtm_vars], 1, select_random_pad_and_depth, depth_seq = depths))
  # dsm_results <- t(apply(predictors[pad_dsm_vars], 1, select_random_pad_and_depth, depth_seq = depths))
  # 
  # # Extract values
  # predictors$random_pad_dtm <- dtm_results[, 1]
  # predictors$depth_dtm <- dtm_results[, 2]
  # predictors$random_pad_dsm <- dsm_results[, 1]
  # predictors$depth_dsm <- dsm_results[, 2]
  # 
  # predictors$deltaLAI_dtm <- predictors$random_pad_dtm - predictors$s2_lai
  # predictors$deltaLAI_dsm <- predictors$random_pad_dsm - predictors$s2_lai
  # predictors <- predictors[, !names(predictors) %in% c(pad_dtm_vars, pad_dsm_vars)]
  
  
  
  # select_all_valid_pad_and_depths <- function(row, depth_seq) {
  #   non_na_indices <- which(!is.na(row))  # Find indices with valid values
  #   if (length(non_na_indices) > 0) {
  #     random_pads <- row[non_na_indices]
  #     depths <- 40 - (2 + (depth_seq[non_na_indices] - 1))
  #     return(data.frame(random_pad = random_pads, depth = depths))
  #   } else {
  #     return(data.frame(random_pad = numeric(0), depth = numeric(0)))
  #   }
  # }
  # dtm_results_list <- apply(predictors[pad_dtm_vars], 1, select_all_valid_pad_and_depths, depth_seq = depths)
  # dsm_results_list <- apply(predictors[pad_dsm_vars], 1, select_all_valid_pad_and_depths, depth_seq = depths)
  # dtm_counts <- sapply(dtm_results_list, nrow)
  # dsm_counts <- sapply(dsm_results_list, nrow)
  # max_counts <- pmax(dtm_counts, dsm_counts)
  # predictors <- predictors[rep(1:nrow(predictors), times = max_counts), ]
  # dtm_results <- do.call(rbind, lapply(seq_along(dtm_results_list), function(i) {
  #   result <- dtm_results_list[[i]]
  #   if (nrow(result) < max_counts[i]) {
  #     result <- rbind(result, data.frame(random_pad = rep(NA, max_counts[i] - nrow(result)),
  #                                        depth = rep(NA, max_counts[i] - nrow(result))))
  #   }
  #   result
  # }))
  # dsm_results <- do.call(rbind, lapply(seq_along(dsm_results_list), function(i) {
  #   result <- dsm_results_list[[i]]
  #   if (nrow(result) < max_counts[i]) {
  #     result <- rbind(result, data.frame(random_pad = rep(NA, max_counts[i] - nrow(result)),
  #                                        depth = rep(NA, max_counts[i] - nrow(result))))
  #   }
  #   result
  # }))
  # predictors$random_pad_dtm <- dtm_results$random_pad
  # predictors$depth_dtm <- dtm_results$depth
  # predictors$random_pad_dsm <- dsm_results$random_pad
  # predictors$depth_dsm <- dsm_results$depth
  # predictors$deltaLAI_dtm <- predictors$random_pad_dtm - predictors$s2_lai
  # predictors$deltaLAI_dsm <- predictors$random_pad_dsm - predictors$s2_lai
  predictors <- predictors[, !names(predictors) %in% c(pad_dtm_vars, pad_dsm_vars)]
  predictors$site <- site
  return(predictors)
}

partialPlot <- function(obj, pred.data, xname, n.pt = 10, discrete.x = FALSE,
                        subsample = pmin(1, n.pt * 100 / nrow(pred.data)), which.class = NULL,
                        xlab = deparse(substitute(xname)), ylab = "", type = if (discrete.x) "p" else "b",
                        main = "", rug = TRUE, seed = NULL, show_plot = TRUE, target_var,
                        save_path = NULL, ...) {
  stopifnot(dim(pred.data) >= 1)
  
  # Subsample the data if needed
  if (subsample < 1) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    n <- nrow(pred.data)
    picked <- sample(n, trunc(subsample * n))
    pred.data <- pred.data[picked, , drop = FALSE]
  }
  
  # Extract and validate xname column
  xv <- pred.data[[xname]]
  if (is.list(xv)) {
    xv <- unlist(xv)  # Flatten list to vector
  }
  xv <- as.numeric(xv)  # Ensure it's numeric
  
  # Define quantiles from 0.1 to 0.9 (9 points)
  if (discrete.x) {
    x <- unique(xv)
  } else {
    x <- quantile(xv, probs = seq(0, 1, length.out = n.pt), names = FALSE)  # Quantiles from 0.1 to 0.9
  }
  
  y <- numeric(length(x))
  
  isRanger <- inherits(obj, "ranger")
  isLm <- inherits(obj, "lm") | inherits(obj, "lmrob") | inherits(obj, "lmerMod")
  
  # Calculate partial dependence by predicting for each x[i]
  for (i in seq_along(x)) {
    pred.data[, xname] <- x[i]
    
    if (isRanger) {
      if (!is.null(which.class)) {
        if (obj$treetype != "Probability estimation") {
          stop("Choose probability = TRUE when fitting ranger multiclass model")
        }
        preds <- predict(obj, pred.data)$predictions[, which.class]
      }
      else {
        preds <- predict(obj, pred.data)$predictions
      }
    } else if (isLm) {
      preds <- predict(obj, pred.data)
    } else {
      if (!is.null(which.class)) {
        preds <- predict(obj, pred.data, reshape = TRUE)[, which.class + 1]
      } else {
        preds <- predict(obj, pred.data)
      }
    }
    
    y[i] <- mean(preds)
  }
  
  if (show_plot) {
    # Save the plot if save_path is provided
    if (!is.null(save_path)) {
      png(filename = file.path(save_path, paste0("pp_", xname, "_target_", target_var, ".png")), width = 1920, height = 1080, res = 300)
    }
    
    # Create the plot
    plot(x, y, xlab = xlab, ylab = ylab,
         main = paste("Partial Plot for", xname, "target var", target_var),
         type = type, ...)
    
    # Cross a line at the variable value where y reaches 0 for specified xnames
    if (xname %in% c("depth_dtm", "depth_dsm", "mean", "max")) {
      zero_crossings <- which(diff(sign(y)) != 0)
      
      if (length(zero_crossings) > 0) {
        # Interpolate between the first crossing points
        first_crossing <- zero_crossings[1]
        x1 <- x[first_crossing]
        x2 <- x[first_crossing + 1]
        y1 <- y[first_crossing]
        y2 <- y[first_crossing + 1]
        zero_x <- x1 - y1 * (x2 - x1) / (y2 - y1)
        abline(v = zero_x, col = "red", lty = 2)
      }
    }
    if (!is.null(save_path)) {
      dev.off()  # Close the graphical device
    }
  }
  par(mfrow = c(1,1))
  # Return the partial dependence data for further use
  # data.frame(x = x, y = y)
}

biPartialPlot <- function(obj, pred.data, xname1, xname2, n.pt = 10, discrete.x1 = FALSE, discrete.x2 = FALSE,
                          subsample = pmin(1, n.pt * 100 / nrow(pred.data)), which.class = NULL,
                          xlab1 = deparse(substitute(xname1)), xlab2 = deparse(substitute(xname2)),
                          ylab = "", main = "", rug = TRUE, seed = NULL, show_plot = TRUE, target_var,
                          save_path = NULL, ...) {
  
  stopifnot(dim(pred.data) >= 1)
  
  # Subsample the data if needed
  if (subsample < 1) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    n <- nrow(pred.data)
    picked <- sample(n, trunc(subsample * n))
    pred.data <- pred.data[picked, , drop = FALSE]
  }
  
  xv1 <- pred.data[, xname1]
  xv2 <- pred.data[, xname2]
  
  if (is.list(xv1)) {
    xv1 <- unlist(xv1)
  }
  if (is.list(xv2)) {
    xv2 <- unlist(xv2)
  }
  xv2 <- as.numeric(xv2)
  
  # Define quantiles for each predictor
  if (discrete.x1) {
    x1 <- unique(xv1)
  } else {
    x1 <- quantile(xv1, probs = seq(0, 1, length.out = n.pt), names = FALSE)
  }
  
  if (discrete.x2) {
    x2 <- unique(xv2)
  } else {
    x2 <- quantile(xv2, probs = seq(0, 1, length.out = n.pt), names = FALSE)
  }
  
  # Create a grid for x1 and x2
  grid <- expand.grid(x1 = x1, x2 = x2)
  y <- numeric(nrow(grid))
  
  isRanger <- inherits(obj, "ranger")
  isLm <- inherits(obj, "lm") | inherits(obj, "lmrob") | inherits(obj, "lmerMod")
  
  # Calculate partial dependence by predicting for each combination of x1 and x2
  for (i in seq_len(nrow(grid))) {
    pred.data[, xname1] <- grid$x1[i]
    pred.data[, xname2] <- grid$x2[i]
    
    if (isRanger) {
      if (!is.null(which.class)) {
        if (obj$treetype != "Probability estimation") {
          stop("Choose probability = TRUE when fitting ranger multiclass model")
        }
        preds <- predict(obj, pred.data)$predictions[, which.class]
      }
      else {
        preds <- predict(obj, pred.data)$predictions
      }
    } else if (isLm) {
      preds <- predict(obj, pred.data)
    } else {
      if (!is.null(which.class)) {
        preds <- predict(obj, pred.data, reshape = TRUE)[, which.class + 1]
      } else {
        preds <- predict(obj, pred.data)
      }
    }
    
    y[i] <- mean(preds)
  }
  
  # Plotting the results
  if (show_plot) {
    if (!is.null(save_path)) {
      png(filename = file.path(save_path, paste0("bi_partial_", xname1, "_", xname2, "_target_", target_var, ".png")), width = 1920, height = 1080, res = 300)
    }
    
    # Create the plot
    filled.contour(x1, x2, matrix(y, nrow = length(x1)),
                   xlab = xlab1, ylab = xlab2,
                   main = paste("Bi-Partial Plot for", xname1, "and", xname2,
                                "target var", target_var), ...)
    
    if (!is.null(save_path)) {
      dev.off()  # Close the graphical device
    }
  }
  
  # Return the partial dependence data for further use
  data.frame(x1 = rep(x1, each = length(x2)), x2 = rep(x2, times = length(x1)), y = y)
}

# Define a function to handle each sampling case and save results to a text file
run_random_forest_analysis <- function(train_data_vario, data_indices, file_name) {
  # Initialize empty vectors to store results
  oob_mse_list <- numeric(100)
  r_squared_list <- numeric(100)
  oob_mse_opt_list <- numeric(100)
  r_squared_opt_list <- numeric(100)
  morans_i_list <- numeric(100)
  gearys_c_list <- numeric(100)
  
  # Prepare data for the specified indices
  train_data <- train_data_vario[data_indices, -c(16:55)]
  train_data$deltaLAI <- train_data$lidar_lai - train_data$s2_lai
  train_data <- train_data[, c("deltaLAI", setdiff(names(train_data), "deltaLAI"))]
  train_data$lidar_lai <- NULL
  
  # Initialize a vector to track feature selection counts
  feature_names <- names(train_data[, -1]) # Exclude deltaLAI
  feature_counts <- setNames(numeric(length(feature_names)), feature_names)
  
  # Create a spatial weights matrix
  coords <- cbind(train_data$lon, train_data$lat)  # Adjust to your coordinate columns
  nb <- 10 # dnearneigh(coords, 0, 1)  # Create neighbors based on distance (adjust as necessary)
  lw <- nb2listw(nb, style = "W")  # Create spatial weights list
  
  # Repeat the process 100 times
  for (i in 1:100) {
    print(paste("Iteration:", i))
    
    # Step 1: Sample 60 points from train_data
    sampled_data <- train_data[sample(nrow(train_data), 60), ]
    
    # Step 2: Train initial Random Forest model
    rf_model <- ranger(
      formula = deltaLAI ~ .,
      data = sampled_data,
      num.trees = 500,
      mtry = ceiling(ncol(sampled_data[, -1]) / 3),
      importance = "none",
      min.node.size = 5,
      num.threads = 12
    )
    
    # Store the OOB MSE and R-squared
    oob_mse_list[i] <- rf_model$prediction.error
    r_squared_list[i] <- rf_model$r.squared
    
    # Step 4: Feature selection with VSURF
    vsurf <- VSURF(deltaLAI ~ ., sampled_data, mtry = 6, parallel = TRUE)
    selected_vars <- names(sampled_data[, -1])[vsurf$varselect.interp]
    
    # Update feature counts
    for (var in selected_vars) {
      feature_counts[var] <- feature_counts[var] + 1
    }
    
    # Step 5: Train optimized Random Forest model on selected variables
    selected_vars <- c("deltaLAI", selected_vars) # Add response variable
    sampled_data_opt <- sampled_data[selected_vars]
    rf_model_opt <- ranger(
      formula = deltaLAI ~ .,
      data = sampled_data_opt,
      num.trees = 500,
      importance = "permutation",
      min.node.size = 5,
      num.threads = 12
    )
    
    # Store the optimized OOB MSE and R-squared
    oob_mse_opt_list[i] <- rf_model_opt$prediction.error
    r_squared_opt_list[i] <- rf_model_opt$r.squared
    
    # Calculate Moran's I and Geary's C for deltaLAI
    moran_result <- moran.test(sampled_data$deltaLAI, listw = lw)
    geary_result <- geary.test(sampled_data$deltaLAI, listw = lw)
    
    # Store Moran's I and Geary's C
    morans_i_list[i] <- moran_result$estimate[1]
    gearys_c_list[i] <- geary_result$estimate[1]
  }
  
  # Calculate averages across the 100 iterations
  mean_oob_mse <- mean(oob_mse_list)
  mean_r_squared <- mean(r_squared_list)
  mean_oob_mse_opt <- mean(oob_mse_opt_list)
  mean_r_squared_opt <- mean(r_squared_opt_list)
  mean_morans_i <- mean(morans_i_list)
  mean_gearys_c <- mean(gearys_c_list)
  
  # Calculate average feature selection frequency
  selection_frequency <- feature_counts / 100
  selection_frequency <- sort(selection_frequency, decreasing = TRUE)
  
  # Save results to a text file
  results <- c(
    paste("Average OOB MSE (original model):", round(mean_oob_mse, 2)),
    paste("Average R-squared (original model):", round(mean_r_squared, 2)),
    paste("Average OOB MSE (optimized model):", round(mean_oob_mse_opt, 2)),
    paste("Average R-squared (optimized model):", round(mean_r_squared_opt, 2)),
    paste("Average Moran's I:", round(mean_morans_i, 3)),
    paste("Average Geary's C:", round(mean_gearys_c, 3)),
    "Feature selection frequency (percentage):"
  )
  
  # Add selection frequencies to results
  results <- c(results, paste(names(selection_frequency), round(selection_frequency * 100, 2), sep = ": "))
  
  # Write results to disk
  writeLines(results, con = file_name)
  cat("Results saved to:", file_name, "\n")
}

# Function to plot the evolution of OOB MSE and R² with mtry
plot_rf_mtry_evolution <- function(train_data, 
                                   X,
                                   formula, 
                                   ntree = 500, 
                                   min_node_size = 5,
                                   num_threads = 12) {
  
  # Get the number of predictors
  mtry_values <- seq(1, X, by = 1) # - deltaLAI, lidar_lai
  
  # Initialize data storage for the results
  data_plot <- data.frame(mtry = integer(),
                          OOB_MSE = numeric(),
                          R_squared = numeric())
  
  # Loop over each mtry value to build the model
  for (m in mtry_values) {
    print(paste("Training with mtry =", m))
    
    # Train the Random Forest model
    rf_model <- ranger(
      formula = formula,
      data = train_data,
      num.trees = ntree,
      mtry = m,
      importance = "permutation",
      min.node.size = min_node_size,
      num.threads = num_threads
    )
    
    # Save the results
    data_plot <- rbind(data_plot, data.frame(
      mtry = m,
      OOB_MSE = rf_model$prediction.error,
      R_squared = rf_model$r.squared
    ))
  }
  
  # Find best mtry based on OOB MSE
  best_mtry <- data_plot[which.min(data_plot$OOB_MSE), "mtry"]
  best_oob_mse <- data_plot[which.min(data_plot$OOB_MSE), "OOB_MSE"]
  best_r_squared <- data_plot[data_plot$mtry == best_mtry, "R_squared"]
  
  # Plot the results using ggplot
  p <- ggplot(data_plot, aes(x = mtry)) +
    geom_line(aes(y = OOB_MSE, color = 'OOB MSE'), size = 1) +
    geom_line(aes(y = R_squared, color = 'OOB R²'), size = 1) +
    geom_point(aes(y = OOB_MSE, color = 'OOB MSE'), size = 2) +
    geom_point(aes(y = R_squared, color = 'OOB R²'), size = 2) +
    geom_vline(aes(xintercept = best_mtry, color = "Best mtry"), linetype = "dashed") +
    scale_y_continuous(
      name = "OOB MSE",
      sec.axis = sec_axis(~., name = "R²")
    ) +
    scale_x_continuous(breaks = seq(min(data_plot$mtry), max(data_plot$mtry), by = 1)) +
    labs(title = paste("Evolution of OOB MSE and R² with mtry"),
         x = "mtry") +
    theme_bw() +
    scale_color_manual(name = "", 
                       values = c("OOB MSE" = "blue", "OOB R²" = "red", "Best mtry" = "black")) +
    theme(legend.position = "top")
  print(p)
  
  return(best_mtry)
}

create_corrected_s2lai_scatterplot <- function(study_site,
                                               predictions,
                                               test_data, 
                                               lai_col, 
                                               s2_lai_col,
                                               depth = NULL) {
  # Predictions
  # predictions <- predict(rf_model, data = test_data)$predictions
  s2_plus_delta_lai <- predictions + test_data[[s2_lai_col]]
  test_data <- cbind(S2_plus_deltaLAI = s2_plus_delta_lai, test_data)
  
  # Create plot data
  plot_data <- data.frame(
    S2_plus_deltaLAI = s2_plus_delta_lai,
    LiDAR_LAI = test_data[[lai_col]]
  )
  if (!is.null(depth)) {
    plot_data$Depth <- test_data[[depth]]
  }
  # Linear model
  lm_fit <- lm(LiDAR_LAI ~ S2_plus_deltaLAI, data = plot_data)
  r_squared <- summary(lm_fit)$r.squared
  mse <- mean((plot_data$LiDAR_LAI - predict(lm_fit, plot_data))^2)
  range_LAI <- range(plot_data$LiDAR_LAI)
  nrmse <- sqrt(mse) / (range_LAI[2] - range_LAI[1])
  
  # Legend adaptation
  title <- switch(lai_col,
                  "random_pad_dtm" = paste(study_site, 
                                           "Scatterplot: S2 + deltaLAI vs Random LiDAR LAI (DTM)"),
                  "random_pad_dsm" = paste(study_site, 
                                           "Scatterplot: S2 + deltaLAI vs Random LiDAR LAI (DSM)"),
                  "lidar_lai" = paste(study_site, 
                                      "Scatterplot: S2 + deltaLAI vs True LiDAR LAI")
  )
  ylab = switch(lai_col,
                "random_pad_dtm" = "Random LiDAR LAI (DTM)",
                "random_pad_dsm" = "Random LiDAR LAI (DSM)",
                "lidar_lai" = "LiDAR LAI"
  )
  
  # Scatterplot with regression line
  scatter_plot <- ggplot(plot_data, aes(x = S2_plus_deltaLAI, y = LiDAR_LAI)) +
    geom_point(alpha = 0.7, aes(color = if (!is.null(Depth)) Depth else NULL)) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +
    scale_color_viridis_c() +
    labs(
      title = title,
      x = "Sentinel-2 LAI + deltaLAI",
      y = ylab,
      color = if (!is.null(depth)) "Depth" else NULL
    ) +
    annotate(
      "text", 
      x = max(plot_data$S2_plus_deltaLAI, na.rm = TRUE) * 0.8,
      y = max(plot_data$LiDAR_LAI, na.rm = TRUE) * 0.5,
      label = paste0("R² = ", round(r_squared, 2), 
                     "\nNRMSE = ", round(nrmse, 2)),
      color = "black",
      size = 5,
      hjust = 0
    ) +
    # xlim(-1, 11) +
    # ylim(0, 13) +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      legend.position = if (!is.null(depth)) "right" else "none"
    )
  print(scatter_plot)
  return(list(test_data = test_data
              # plot = scatter_plot, r_squared = r_squared, mse = mse
  ))
}

# Define a function to compute performance metrics
compute_rf_metrics <- function(model,
                               test_predictions, 
                               train_data, 
                               test_data, 
                               target_col) {
  # Predictions on the training data
  # train_predictions <- predict(rf_model, data = train_data)
  # train_actual <- train_data[[target_col]]
  # train_mse <- mean((train_predictions$predictions - train_actual) ^ 2)
  # train_rss <- sum((train_predictions$predictions - train_actual) ^ 2)
  # train_tss <- sum((train_actual - mean(train_actual)) ^ 2)
  # train_r_squared <- 1 - train_rss / train_tss
  
  # Predictions on the test data
  # test_predictions <- predict(rf_model, data = test_data)
  test_actual <- test_data[[target_col]]
  test_mse <- mean((test_prediction - test_actual) ^ 2)
  test_rss <- sum((test_prediction - test_actual) ^ 2)
  test_tss <- sum((test_actual - mean(test_actual)) ^ 2)
  test_r_squared <- 1 - test_rss / test_tss
  
  # Extract OOB metrics
  oob_mse <- model$prediction.error
  oob_r_squared <- model$r.squared
  
  # Output results
  cat("Training MSE:", round(train_mse, 2), "\n")
  cat("Training R-squared:", round(train_r_squared, 2), "\n")
  cat("Test MSE:", round(test_mse, 2), "\n")
  cat("Test R-squared:", round(test_r_squared, 2), "\n")
  cat("OOB MSE:", round(oob_mse, 2), "\n")
  cat("OOB R-squared:", round(oob_r_squared, 2), "\n")
  
  # Return results as a list
  return(list(
    train_mse = train_mse,
    train_r_squared = train_r_squared,
    test_mse = test_mse,
    test_r_squared = test_r_squared,
    oob_mse = oob_mse,
    oob_r_squared = oob_r_squared
  ))
}

# Function to perform Moran's I and Geary's C tests with =/= neighbor strategies
spatial_autocorrelation_tests <- function(moran_data, 
                                          variable_name, 
                                          site = NULL) {
  k_values = 1:20
  dist_values <- switch(site,
                        "Aigoual" = seq(700, 2000, by = 100),
                        "Blois" = seq(800, 2000, by = 100),
                        "Mormal" = seq(1200, 2000, by = 100)
  )
  
  coordinates(moran_data) <- ~ lon + lat
  moran_values_k <- numeric(length(k_values))
  geary_values_k <- numeric(length(k_values))
  moran_values_dist <- numeric(length(dist_values))
  geary_values_dist <- numeric(length(dist_values))
  
  for (k in k_values) {
    neighbors_k <- knearneigh(coordinates(moran_data), k = k)
    lw_weights_k <- nb2listw(knn2nb(neighbors_k), style = "W")
    moran_variable <- moran_data[[variable_name]]
    moran_test_k <- moran.test(moran_variable, listw = lw_weights_k)
    moran_values_k[k] <- moran_test_k$estimate[1]
    geary_test_k <- geary.test(moran_variable, listw = lw_weights_k)
    geary_values_k[k] <- geary_test_k$estimate[1]
  }
  
  for (dist in dist_values) {
    neighbors_dist <- dnearneigh(coordinates(moran_data), 0, dist)
    lw_weights_dist <- nb2listw(neighbors_dist, style = "W", zero.policy = T)
    moran_test_dist <- moran.test(moran_variable, listw = lw_weights_dist)
    moran_values_dist[which(dist_values == dist)] <- moran_test_dist$estimate[1]
    geary_test_dist <- geary.test(moran_variable, listw = lw_weights_dist)
    geary_values_dist[which(dist_values == dist)] <- geary_test_dist$estimate[1]
  }
  
  # Prepare results as a data frame
  results <- list(
    moran_values_k = moran_values_k,
    geary_values_k = geary_values_k,
    moran_values_dist = moran_values_dist,
    geary_values_dist = geary_values_dist
  )
  
  # Increase the size of the plotting window (if necessary)
  par(mfrow = c(2, 2)) # Reset to single plot layout for setting margins
  # dev.new(width = 20, height = 20)  # Increase the plot size if using RStudio or similar
  
  # Correct layout (2x2 grid for 4 plots)
  # layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
  
  # Set margins if needed to avoid "too large" error
  # par(mar = c(5, 4, 2, 1))  # Bottom, left, top, right margins
  plot(k_values, moran_values_k, type = "o", col = "blue", xlab = "Number of Neighbors (k)", 
       ylab = "Moran's I", main = "Moran's I with k-Nearest Neighbors")
  grid()
  plot(k_values, geary_values_k, type = "o", col = "red", xlab = "Number of Neighbors (k)", 
       ylab = "Geary's C", main = "Geary's C with k-Nearest Neighbors")
  grid()
  plot(dist_values, moran_values_dist, type = "o", col = "green", xlab = "Distance (m)", 
       ylab = "Moran's I", main = "Moran's I with Circular Distance")
  grid()
  plot(dist_values, geary_values_dist, type = "o", col = "purple", xlab = "Distance (m)", 
       ylab = "Geary's C", main = "Geary's C with Circular Distance")
  grid()
  mtext(paste(site, "Spatial Autocorrelation Tests"), 
        side = 3, line = -16, cex = 1.2, outer = TRUE)
  par(mfrow=c(1,1))
  # Return the results
  return(results)
}

spatial_moran_test <- function(moran_data, 
                               variable_name, 
                               site = NULL, 
                               k_values = 1:20, 
                               iterations = 100) {
  coordinates(moran_data) <- ~ lon + lat
  
  # Initialize vectors to store results
  moran_mean_k <- numeric(length(k_values))
  moran_std_k <- numeric(length(k_values))
  moran_p_value_mean_k <- numeric(length(k_values))
  moran_p_value_std_k <- numeric(length(k_values))
  
  # Running Moran's I test 'iterations' times for k-nearest neighbors
  for (k in k_values) {
    cat(paste("k =", k), "\n")
    moran_values_k <- numeric(iterations)
    p_values_k <- numeric(iterations)
    
    # Perform the specified number of iterations
    for (i in 1:iterations) {
      neighbors_k <- knearneigh(coordinates(moran_data), k = k)
      lw_weights_k <- nb2listw(knn2nb(neighbors_k), style = "W")
      moran_variable <- moran_data[[variable_name]]
      moran_test_k <- moran.test(moran_variable, listw = lw_weights_k)
      moran_values_k[i] <- moran_test_k$estimate[1]  # Moran's I value
      p_values_k[i] <- moran_test_k$p.value  # p-value
    }
    
    # Compute mean and standard deviation of Moran's I and p-values
    moran_mean_k[k] <- mean(moran_values_k)
    moran_std_k[k] <- sd(moran_values_k)
    moran_p_value_mean_k[k] <- mean(p_values_k)
    moran_p_value_std_k[k] <- sd(p_values_k)
  }
  
  # Prepare results as a data frame
  results <- list(
    moran_mean_k = moran_mean_k,
    moran_std_k = moran_std_k,
    moran_p_value_mean_k = moran_p_value_mean_k,
    moran_p_value_std_k = moran_p_value_std_k
  )
  
  # Plot results
  par(mfrow = c(1, 2))  # Setup 1x2 grid of plots
  
  # Plot Mean Moran's I for k-nearest neighbors with error bars (std)
  plot(k_values, moran_mean_k, type = "o", col = "blue", xlab = "Number of Neighbors (k)", 
       ylab = "Mean Moran's I", main = "Moran's I with k-Nearest Neighbors")
  arrows(k_values, moran_mean_k - moran_std_k, k_values, moran_mean_k + moran_std_k, 
         angle = 90, code = 3, length = 0.1, col = "blue")
  grid()
  
  # Plot Mean p-value for k-nearest neighbors with error bars (std)
  plot(k_values, moran_p_value_mean_k, type = "o", col = "green", xlab = "Number of Neighbors (k)", 
       ylab = "Mean p-value", main = "Mean p-value for k-Nearest Neighbors")
  arrows(k_values, moran_p_value_mean_k - moran_p_value_std_k, k_values, moran_p_value_mean_k + moran_p_value_std_k, 
         angle = 90, code = 3, length = 0.1, col = "green")
  grid()
  
  mtext(paste0(site), 
        side = 3, 
        # line = 12, 
        cex = 1.2, outer = TRUE)
  par(mfrow=c(1,1))  # Reset plot layout
  
  # Return the results
  return(results)
}

round_to_odd_multiple_of_5 <- function(x) {
  # Function to round a single number
  round_single <- function(num) {
    # Find the nearest multiple of 5
    rounded_to_5 <- round(num / 5) * 5
    
    # If the result is an even multiple of 5, adjust it to the nearest odd multiple
    if (rounded_to_5 %% 10 == 0) {
      if (rounded_to_5 > num) {
        # If the rounded value is larger than the input, decrease by 5
        return(rounded_to_5 - 5)
      } else {
        # Otherwise, increase by 5
        return(rounded_to_5 + 5)
      }
    }
    
    # If the result is already an odd multiple of 5, return it
    return(rounded_to_5)
  }
  
  # Apply the function to each element in the vector x
  sapply(x, round_single)
}

calculate_depth_vectors <- function(valid_weights_matrix, iterations = 100) {
  depth_matrix <- matrix(NA, nrow = nrow(valid_weights_matrix), ncol = iterations)
  selected_files_matrix <- matrix(NA, nrow = nrow(valid_weights_matrix), ncol = iterations)
  for (iter in 1:iterations) {
    selected_file_indices <- apply(valid_weights_matrix, 1, 
                                   function(x) sample(seq_along(x), 
                                                      size = 1, prob = x))
    depth_matrix[, iter] <- 40 - (2 + selected_file_indices) + 1
    selected_files_matrix[, iter] <- selected_file_indices
  }
  return(list(depth_matrix = depth_matrix, 
              selected_files_matrix = selected_files_matrix))
}

extract_pad_value <- function(pad_files, lon, lat, selected_file_index) {
  pad_raster <- rast(pad_files[selected_file_index])
  pad_value <- terra::extract(pad_raster, cbind(lon, lat))
  pad_value <- as.numeric(unlist(pad_value))
  return(pad_value)
}

# Function to process the test_data_dtm dataframe
process_depth_and_pad <- function(dataframe, pad_files, site,
                                  plot = FALSE, num_iterations = 100) {
  # Initialize necessary matrices
  coords <- data.frame(lon = dataframe$lon, lat = dataframe$lat)
  non_na_counts_matrix <- matrix(NA, nrow = nrow(coords), ncol = length(pad_files))
  
  # Extract values from pad_files
  for (i in seq_along(pad_files)) {
    raster <- rast(pad_files[i])
    # print(raster)
    extracted_values <- terra::extract(raster, coords)
    non_na_counts_matrix[, i] <- extracted_values[[2]]
  }
  
  # Calculate weights based on NA counts
  na_counts_per_column <- colSums(is.na(non_na_counts_matrix))
  na_counts_per_column <- ifelse(na_counts_per_column == 0, 1, na_counts_per_column)
  print(na_counts_per_column)
  epsilon <- 1e-20
  # weights <- (1/na_counts_per_column)^1
  # print(weights)
  # weights <- 1 - weights
  # weights <- 1 / (weights)
  # print(weights)
  # weights <- log(na_counts_per_column + 1)^1
  # weights <- ifelse(weights == 0, epsilon, weights)
  # Step 2: Normalize na_counts to range [0, 1]
  max_na <- max(na_counts_per_column)
  min_na <- min(na_counts_per_column)
  normalized_na <- (na_counts_per_column - min_na) / (max_na - min_na)
  normalized_na <- ifelse(normalized_na == 0, epsilon, normalized_na)
  print(normalized_na)
  coeff <- switch(site,
                  # "Aigoual" = 0.5,
                  # "Blois" = 0.1,
                  # "Mormal" = 0.8
                  "Aigoual" = 5,
                  "Blois" = 3,
                  "Mormal" = 3
  )
  sd_val <- sd(normalized_na, na.rm = TRUE)
  alpha <- 1.5 # Power to amplify the sensitivity
  scale_factor <- 1  # Multiplier to globally scale the coefficient
  coeff <- (1 / (sd_val + 1e-6))^alpha * scale_factor
  # coeff <- 1 / (sd_val^1.5 + epsilon)
  # weights <- normalized_na^coeff
  # weights <- ifelse(weights == 0, epsilon, weights)
  # weights <- weights / sum(weights)
  # weights <- normalized_na / sum(normalized_na)
  # weights <- 1 / (normalized_na * (1 - normalized_na) + epsilon)
  weights <- exp(normalized_na)^coeff
  weights <- weights / sum(weights, na.rm = TRUE)
  print(paste("Standard Deviation:", sd_val))
  print(paste("Enhanced Coefficient:", coeff))
  print(weights)
  
  # Create valid weights matrix
  valid_indices_matrix <- ifelse(!is.na(non_na_counts_matrix), 1, 0)
  print(str(valid_indices_matrix))
  valid_weights_matrix <- matrix(0, nrow = nrow(valid_indices_matrix), ncol = ncol(valid_indices_matrix))
  for (i in 1:nrow(valid_indices_matrix)) {
    valid_weights_matrix[i, ] <- valid_indices_matrix[i, ] * weights
    row_sum <- sum(valid_weights_matrix[i, ], na.rm = TRUE)
    if (row_sum > 0) {
      valid_weights_matrix[i, ] <- valid_weights_matrix[i, ] / row_sum
    } else {
      # Handle rows with all NA values or zero probabilities
      valid_weights_matrix[i, ] <- 0  # Set to 0 to avoid sampling errors
      valid_weights_matrix[i, ncol(valid_weights_matrix)] <- 1
    }
  }
  # print(valid_weights_matrix[, ])
  # print(sum(valid_weights_matrix[, ]))  # Ensure sum is 1
  # print(is.na(valid_weights_matrix[, ]))  # Check for NA values
  # print(all(valid_weights_matrix[, ] >= 0))  # Ensure no negative values
  # Select file indices based on the weights
  selected_file_indices <- apply(valid_weights_matrix, 1, function(x) {
    if (sum(x, na.rm = TRUE) > 0) {
      sample(seq_along(x), size = 1, prob = x)
    } else {
      print(x)
      NA  # Assign NA if all probabilities are 0
    }
  })
  
  # Calculate depth vectors
  depth_vectors <- matrix(NA, nrow = nrow(valid_weights_matrix), ncol = num_iterations)
  for (iter in 1:num_iterations) {
    depth_vectors[, iter] <- 40 - (2 + selected_file_indices) + 1
  }
  
  dataframe$depth <- floor(rowMeans(depth_vectors, na.rm = TRUE))
  selected_file_indices <- 39 - dataframe$depth
  
  dataframe$random_pad <- NA
  for (i in 1:nrow(dataframe)) {
    selected_file_index <- selected_file_indices[i]
    lon <- dataframe$lon[i]
    lat <- dataframe$lat[i]
    dataframe$random_pad[i] <- extract_pad_value(pad_files, lon, lat, selected_file_index)
  }
  dataframe$delta <- dataframe$random_pad - dataframe$s2_lai
  # Optionally plot histograms
  if (plot) {
    hist_bins <- seq(min(depth_vectors, na.rm = TRUE), max(depth_vectors, na.rm = TRUE), length.out = 11)
    hist_counts <- matrix(0, nrow = num_iterations, ncol = length(hist_bins) - 1)
    for (iter in 1:num_iterations) {
      hist_counts[iter, ] <- hist(depth_vectors[, iter], breaks = hist_bins, plot = FALSE)$counts
    }
    mean_hist_counts <- colMeans(hist_counts)
    
    # Set up plotting area
    # par(mfrow = c(1, 2))
    par(
      mfrow=c(1,2),
      mar=c(4,4,1,0)
    )
    # combined_range <- range(c(dataframe$depth, dataframe$mean + 1), na.rm = TRUE)
    combined_range <- c(0, 40)
    breaks <- seq(combined_range[1], combined_range[2], length.out = 9)
    hist(dataframe$depth, breaks = breaks, col = "lightblue", 
         main = "Histogram of Depth", xlab = "Depth", ylab = "Count",
         xlim = combined_range)
    hist(dataframe$mean, breaks = breaks, col = "lightgreen", 
         main = "Histogram of Mean", xlab = "Mean", ylab = "Count",
         xlim = combined_range)
  }
  par(mfrow = c(1, 1))
  print(cor(dataframe$depth, dataframe$mean))
  
  return(dataframe)
}

process_depth_and_pad_uniform <- function(dataframe, pad_files, target_per_depth = NULL) {
  # Number of points and depths
  n_points <- nrow(dataframe)
  n_depths <- length(pad_files)
  
  # If target_per_depth is not provided, aim for an even split
  if (is.null(target_per_depth)) {
    target_per_depth <- floor(n_points / n_depths)
  }

  # Build a validity matrix (rows = points, columns = pad files)
  coords <- data.frame(lon = dataframe$lon, lat = dataframe$lat)
  valid_matrix <- matrix(0, nrow = n_points, ncol = n_depths)
  for (i in seq_along(pad_files)) {
    raster <- rast(pad_files[i])
    extracted_values <- terra::extract(raster, coords)
    # Mark as 1 if value is non-NA, else 0
    valid_matrix[, i] <- ifelse(is.na(extracted_values[[2]]), 0, 1)
  }
  
  # Initialize vector to store selected pad file index (one per point)
  selected_depth_indices <- rep(NA, n_points)
  assigned <- rep(FALSE, n_points)
  
  # Process pad files in order from highest depth to lowest.
  # (Note: In your original code, depth is computed as 40 - (2 + index) + 1, so pad_files[1] gives the highest depth.)
  for (i in 1:n_depths) {
    eligible <- which(valid_matrix[, i] == 1 & !assigned)
    if (length(eligible) > 0) {
      n_assign <- min(length(eligible), target_per_depth)
      selected <- sample(eligible, n_assign)
      selected_depth_indices[selected] <- i
      assigned[selected] <- TRUE
    }
  }
  
  # For any unassigned rows, assign them a depth from any available pad file
  unassigned <- which(!assigned)
  if (length(unassigned) > 0) {
    depth_counts <- rep(0, n_depths)  # Track how many points are assigned to each depth
    
    for (j in unassigned) {
      available <- which(valid_matrix[j,] == 1)
      if (length(available) > 0) {
        # Choose the depth with the least assigned points among available ones
        selected_depth_indices[j] <- available[which.min(depth_counts[available])]
        depth_counts[selected_depth_indices[j]] <- depth_counts[selected_depth_indices[j]] + 1
      }
    }
  }
  # print(assigned)
  # print(selected_depth_indices)
  
  # Convert the pad file index into depth value as before.
  # (With your original formula: depth = 40 - (2 + selected_index) + 1)
  dataframe$depth <- 40 - (2 + selected_depth_indices) + 1
  
  # Extract the corresponding pad value for each point
  dataframe$random_pad <- NA
  for (i in 1:n_points) {
    selected_index <- selected_depth_indices[i]
    if (!is.na(selected_index)) {
      dataframe$random_pad[i] <- extract_pad_value(pad_files, dataframe$lon[i], dataframe$lat[i], selected_index)
    }
  }
  
  # Calculate delta (or any other derived quantity)
  dataframe$delta <- dataframe$random_pad - dataframe$s2_lai_dsm
  
  return(dataframe)
}

sample_sites <- function(all_sites_data, sites) {
  selected_data_list <- list()
  for (current_site in sites) {
    if (current_site %in% names(all_sites_data)) {
      current_site_data <- all_sites_data[[current_site]]
      if (nrow(current_site_data) > 0) {
        selected_data_list[[current_site]] <- current_site_data
      }
    } else {
      warning(paste("Site", current_site, "not found in sites data"))
    }
  }
  if (length(selected_data_list) > 0) {
    all_sampled_sites <- do.call(rbind, selected_data_list)
  } else {
    stop("No valid data found for the specified sites.")
  }
  return(all_sampled_sites)
}

# Function to process each site and combine results
process_sites_dtm <- function(all_sites_data, sites, results_dir, metrics_dir, plot = TRUE, num_iterations = 10) {
  # Initialize an empty list to store processed data for each site
  processed_data_list <- list()
  
  # Loop through each site
  for (site in sites) {
    cat("Processing site:", site, "\n")
    
    # Generate PAD file paths for the current site
    pad_files <- sprintf(file.path(results_dir, site, metrics_dir, 
                                   "PAD_Profiles_NA/PAD_%.1f_40.tif"), 
                         seq(2.5, 39.5, by = 1))
    
    # Extract data for the current site
    train_data_dtm <- all_sites_data[[site]]
    
    # Remove unnecessary variables
    remove_vars_dsm <- c("random_pad_dsm", "depth_dsm", 
                         "deltaLAI_dsm", "cv_lad_dsm", "closest_distance", 
                         "site", "strata", "type", "rule", "deltaLAI",
                         "slope", "aspect", "stands", "cv_lad_dsm", "ndvi")
    for (remove_var in remove_vars_dsm) {
      train_data_dtm[[remove_var]] <- NULL
    }
    
    # Process depth and PAD for the current site
    # train_data_dtm <- process_depth_and_pad(train_data_dtm, pad_files, plot = plot, num_iterations = num_iterations)
    train_data_dtm <- process_depth_and_pad_uniform(train_data_dtm, pad_files)
    
    # Rename columns
    colnames(train_data_dtm)[colnames(train_data_dtm) == "depth"] <- "depth_dtm"
    colnames(train_data_dtm)[colnames(train_data_dtm) == "random_pad"] <- "random_pad_dtm"
    colnames(train_data_dtm)[colnames(train_data_dtm) == "delta"] <- "deltaLAI_dtm"
    
    # Add the processed data to the list
    processed_data_list[[site]] <- train_data_dtm
    # print(str(train_data_dtm))
  }
  
  # Concatenate all processed site data
  combined_data <- do.call(rbind, processed_data_list)
  
  return(combined_data)
}


# Function to process each site and combine results
process_sites_dsm <- function(all_sites_data, sites, results_dir, metrics_dir, plot = TRUE, num_iterations = 10) {
  # Initialize an empty list to store processed data for each site
  processed_data_list <- list()
  
  # Loop through each site
  for (site in sites) {
    cat("Processing site:", site, "\n")
    
    # Generate PAD file paths for the current site
    pad_files <- sprintf(file.path(results_dir, site, metrics_dir,
                                   "PAD_Profiles_updated_modifminz_NA/PAD_%.1f_40.tif"), 
                         seq(2.5, 39.5, by = 1))
    
    # Extract data for the current site
    train_data_dsm <- all_sites_data[[site]]
    
    # Remove unnecessary variables
    remove_vars_dtm <- c("random_pad_dtm", "depth_dtm", 
                         "deltaLAI_dtm", "cv_lad_dtm", "closest_distance", 
                         "site", "strata", "type", "rule", "deltaLAI",
                         "slope", "aspect", "stands", "ndvi")
    for (remove_var in remove_vars_dtm) {
      train_data_dtm[[remove_var]] <- NULL
    }
    
    # Process depth and PAD for the current site
    # train_data_dsm <- process_depth_and_pad(train_data_dsm, pad_files, plot = plot, num_iterations = num_iterations)
    train_data_dsm <- process_depth_and_pad_uniform(train_data_dsm, pad_files)
    
    # Rename columns
    colnames(train_data_dsm)[colnames(train_data_dsm) == "depth"] <- "depth_dsm"
    colnames(train_data_dsm)[colnames(train_data_dsm) == "random_pad"] <- "random_pad_dsm"
    colnames(train_data_dsm)[colnames(train_data_dsm) == "delta"] <- "deltaLAI_dsm"
    
    # Add the processed data to the list
    processed_data_list[[site]] <- train_data_dsm
    # print(str(train_data_dsm))
  }
  
  # Concatenate all processed site data
  combined_data <- do.call(rbind, processed_data_list)
  
  return(combined_data)
}

process_sites_dsm_dtm <- function(all_sites_data, sites, results_dir, metrics_dir, plot = TRUE, num_iterations = 10) {
  # Initialize separate lists for DSM and DTM processed data
  processed_dtm_list <- list()
  processed_dsm_list <- list()
  
  # Loop through each site
  for (site in sites) {
    cat("Processing site:", site, "\n")
    
    # Generate PAD file paths for DTM and DSM
    pad_files_dtm <- sprintf(file.path(results_dir, site, metrics_dir, 
                                       "PAD_Profiles_NA/PAD_%.1f_40.tif"), 
                             seq(2.5, 39.5, by = 1))
    pad_files_dsm <- sprintf(file.path(results_dir, site, metrics_dir,
                                       "PAD_Profiles_updated_modifminz_NA/PAD_%.1f_40.tif"), 
                             seq(2.5, 39.5, by = 1))
    
    # Extract data for the current site
    train_data_dtm <- all_sites_data[[site]]
    train_data_dsm <- all_sites_data[[site]]
    
    # Remove unnecessary variables for DTM
    remove_vars_dtm <- c("random_pad_dtm", "depth_dtm", 
                         "deltaLAI_dtm", "cv_lad_dtm", "closest_distance", 
                         "site", "strata", "type", "rule", "deltaLAI",
                         "slope", "aspect", "stands", "ndvi")
    
    # Remove unnecessary variables for DSM
    remove_vars_dsm <- c("random_pad_dsm", "depth_dsm", 
                         "deltaLAI_dsm", "cv_lad_dsm", "closest_distance", 
                         "site", "strata", "type", "rule", "deltaLAI",
                         "slope", "aspect", "stands", "ndvi")
    
    for (remove_var in remove_vars_dsm) {
      train_data_dtm[[remove_var]] <- NULL
    }
    
    for (remove_var in remove_vars_dtm) {
      train_data_dsm[[remove_var]] <- NULL
    }
    
    # Process depth and PAD for DTM (use DTM PAD files for depth generation)
    # train_data_dtm <- process_depth_and_pad(train_data_dtm, pad_files_dtm,
    #                                         site, plot = plot,
    #                                         num_iterations = num_iterations)
    train_data_dtm <- process_depth_and_pad_uniform(train_data_dtm, pad_files_dtm)
    train_data_dsm <- process_depth_and_pad_uniform(train_data_dsm, pad_files_dsm)

    # Apply the same depths to DSM
    # train_data_dsm$depth <- train_data_dtm$depth
    # train_data_dsm$random_pad <- NA
    # for (i in 1:nrow(train_data_dsm)) {
    #   selected_file_index <- 39 - train_data_dsm$depth[i]
    #   lon <- train_data_dsm$lon[i]
    #   lat <- train_data_dsm$lat[i]
    #   train_data_dsm$random_pad[i] <- extract_pad_value(pad_files_dsm, lon, lat, selected_file_index)
    # }
    # train_data_dsm$delta <- train_data_dsm$random_pad - train_data_dsm$s2_lai_dsm
    
    # Rename columns for consistency
    colnames(train_data_dtm)[colnames(train_data_dtm) == "depth"] <- "depth_dtm"
    colnames(train_data_dtm)[colnames(train_data_dtm) == "random_pad"] <- "random_pad_dtm"
    colnames(train_data_dtm)[colnames(train_data_dtm) == "delta"] <- "deltaLAI_dtm"
    
    colnames(train_data_dsm)[colnames(train_data_dsm) == "depth"] <- "depth_dsm"
    colnames(train_data_dsm)[colnames(train_data_dsm) == "random_pad"] <- "random_pad_dsm"
    colnames(train_data_dsm)[colnames(train_data_dsm) == "delta"] <- "deltaLAI_dsm"
    
    # Check for NA rows in both datasets
    dtm_na_indices <- which(rowSums(is.na(train_data_dtm)) > 0)
    dsm_na_indices <- which(rowSums(is.na(train_data_dsm)) > 0)
    
    # Combine NA indices from both datasets
    na_indices <- unique(c(dtm_na_indices, dsm_na_indices))
    if (length(na_indices) > 0) {
      cat("Removing rows with NAs for site:", site, "at indices:", na_indices, "\n")
      train_data_dtm <- train_data_dtm[-na_indices, ]
      train_data_dsm <- train_data_dsm[-na_indices, ]
    }
    
    # Add the processed data to their respective lists
    processed_dtm_list[[site]] <- train_data_dtm
    processed_dsm_list[[site]] <- train_data_dsm
  }
  
  # Return a list with separate DTM and DSM datasets
  return(list(dtm = processed_dtm_list, dsm = processed_dsm_list))
}