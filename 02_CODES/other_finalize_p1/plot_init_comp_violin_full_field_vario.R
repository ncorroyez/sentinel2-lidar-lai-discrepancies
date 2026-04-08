# ---
# title: "plot_init_comp_violin_full_field_vario.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-10-01"
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
library("ggbeeswarm")
library("ggdist")
library("tidyquant")
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
library("caret")
library("VSURF")
library("corrplot")
library("spatialreg")
source("libraries/functions_general_tools.R")
source("libraries/functions_plots.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/variogram"
pp_dir <- "../04_FIGURES/pp"
forest_composition <- "Not_Masked" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
set.seed(42)

# Lists
all_predictors_filtered <- all_predictors_excluded <- all_observations <- list()
moran_filtered_results <- moran_variogram_results <- all_predictors_variogram <- list()
moran_all_obs <- list()
# Load data for all sites
for (site in sites) {
  site_data <- load_metrics_for_site(site, results_dir, forest_composition)
  all_observations[[site]] <- site_data

  geojson_file <- file.path(results_dir, site, "data_utm31n.geojson")
  geo_data <- geojson_read(geojson_file, what = "sp")
  coord_x <- geo_data@data$coord_x_utm31n
  coord_y <- geo_data@data$coord_y_utm31n
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

  # for (i in seq(1,10))
  k <- 8
  neighbors <- knearneigh(coordinates(moran_data), k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  moran_variable <- moran_data$lidar_lai - moran_data$s2_lai

  moran_test <- moran.test(moran_variable, listw = weights)
  moran_filtered_results[[site]] <- moran_test

  cat("Moran's I results (Field) for", site, ":\n")
  print(moran_filtered_results[[site]])
  cat("\n")

  # Moran
  moran_data <- closest_data_vario
  coordinates(moran_data) <- ~ lon + lat  # Adjust based on your actual column names
  proj4string(moran_data) <- CRS("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")  # Set the correct CRS

  # k <- 4
  neighbors <- knearneigh(coordinates(moran_data), k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  moran_variable <- moran_data$lidar_lai - moran_data$s2_lai

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
a
# ------------------------------------------------------------------------------
for (site in sites){
  plot_raincloud_one_site(site,
                          all_observations,
                          all_predictors_filtered,
                          all_predictors_variogram,
                          output_dir)
}
# ------------------------------------------------------------------------------
plot_violin(all_observations,
            all_predictors_filtered,
            all_predictors_variogram,
            sites,
            output_dir)
