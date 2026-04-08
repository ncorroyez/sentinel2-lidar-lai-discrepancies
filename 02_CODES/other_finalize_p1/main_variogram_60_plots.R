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
library("caret")
library("VSURF")
library("corrplot")
library("progress")
library("reshape2")
library("spatialreg")
library("Metrics")
library("pdp")
library("CCA")
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
  site_data <- load_metrics_with_pad(site, results_dir, forest_composition)
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
  k <- 1
  neighbors <- knearneigh(coordinates(moran_data), k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  moran_variable <- moran_data$pad_minus_lai
  
  moran_test <- moran.test(moran_variable, listw = weights)
  moran_filtered_results[[site]] <- moran_test
  
  cat("Moran's I results (Field) for", site, ":\n")
  print(moran_filtered_results[[site]])
  cat("\n")
  
  # Moran
  # moran_data <- closest_data_vario
  # coordinates(moran_data) <- ~ lon + lat  # Adjust based on your actual column names
  # proj4string(moran_data) <- CRS("+proj=utm +zone=31 +north +datum=WGS84 +units=m +no_defs")  # Set the correct CRS
  # 
  # neighbors <- knearneigh(coordinates(moran_data), k = k)
  # weights <- nb2listw(knn2nb(neighbors), style = "W")
  # moran_variable <- moran_data$lidar_lai
  # 
  # moran_test <- moran.test(moran_variable, listw = weights)
  # moran_variogram_results[[site]] <- moran_test
  # 
  # cat("Moran's I results (Variogram 500m) for", site, ":\n")
  # print(moran_variogram_results[[site]])
  # cat("\n")
  
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

# ------------------------------------------------------------------------------
all_balanced_samples <- list()
sites <- c("Aigoual")
for (site in sites){
  print(paste("Site:", site))
  site_data <- all_observations[[site]]
  site_data <- site_data %>% filter(stands == "deciduous")
  
  ####### Test
  # Sort by spatial or structural variables for systematic sampling
  site_data <- site_data %>% arrange(lon, lat)
  
  # Set desired sample size
  sample_size <- 10000
  
  # Calculate systematic sampling interval
  sampling_interval <- nrow(site_data) / sample_size
  print(sampling_interval)
  indices <- round(seq(1, nrow(site_data), by = sampling_interval))
  
  # Select rows based on indices
  site_data <- site_data[indices, ]
  
  
  
  
  
  
  # Normalize structural variables
  # site_data <- site_data[sample(nrow(site_data), 10000), ]
  site_data_scaled <- site_data %>%
    dplyr::select(-lon, -lat, -stands) %>%  # Exclude spatial columns
    scale()
  
  # Compute pairwise Euclidean distances for the structural variables
  structural_dist_matrix <- dist(site_data_scaled)
  
  # Compute spatial distance
  coords <- cbind(site_data$lon, site_data$lat)
  spatial_dist_matrix <- dist(coords)
  
  # Combine spatial and structural distances (adjust the weights)
  combined_dist_matrix_full <- 0.5 * as.matrix(structural_dist_matrix) + 0.5 * as.matrix(spatial_dist_matrix)
  
  # Initialize selected points
  selected_points <- numeric(0)
  remaining_points <- 1:nrow(site_data)
  
  # Select first point randomly
  first_point <- sample(remaining_points, 1)
  selected_points <- c(selected_points, first_point)
  remaining_points <- setdiff(remaining_points, first_point)
  
  # Set desired sample size
  sample_size <- 5000
  
  pb <- progress_bar$new(
    total = sample_size,
    format = "Progress [:bar] :percent :eta",
    clear = TRUE,
    width = 50
  )
  
  for (i in 2:sample_size) {
    # Compute distances from selected points to all remaining points
    dist_to_selected <- apply(combined_dist_matrix_full[remaining_points, ], 1, function(x) min(x[selected_points]))
    
    # Choose the point with the maximum combined distance
    next_point <- remaining_points[which.max(dist_to_selected)]
    selected_points <- c(selected_points, next_point)
    remaining_points <- setdiff(remaining_points, next_point)
    
    pb$tick()
  }
  
  # Extract the sampled data
  balanced_sample <- site_data[selected_points, ]
  all_balanced_samples[[site]] <- balanced_sample
  saveRDS(balanced_sample, paste0(output_dir, "/balanced_sample_", 
                                  sample_size, "_", site, ".rds"))
  coords_sampled <- cbind(balanced_sample$lon, balanced_sample$lat)
  
  # Create spatial weights matrix based on the sampled points
  neighbors_sampled <- dnearneigh(coords_sampled, 0, 1000)  # Adjust distance
  lw_sampled <- nb2listw(neighbors_sampled)
  
  # Spatial autocorrelation tests
  morans_I_sampled <- moran.test(balanced_sample$lidar_lai - balanced_sample$s2_lai, lw_sampled)  # Replace with actual variable
  print(paste("Moran's I for site", site, ":"))
  print(morans_I_sampled)
  geary_sampled <- geary.test(balanced_sample$lidar_lai - balanced_sample$s2_lai, lw_sampled)  # Replace with actual variable
  print(paste("Geary C for site", site, ":"))
  print(geary_sampled)
  
  # Plot
  balanced_sample_plot <- ggplot() +
    geom_point(data = site_data, aes(x = lon, y = lat), color = "gray", alpha = 0.5) +
    # Highlight selected points
    geom_point(data = balanced_sample, aes(x = lon, y = lat), color = "red", size = 3) +
    labs(title = "Spatially Balanced Sample") +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(
    filename = paste0(output_dir, "/balanced_sample_", sample_size, "_", site, ".png"),
    plot = balanced_sample_plot,
    width = 8,  # Width in inches
    height = 6,  # Height in inches
    dpi = 300  # High resolution
  )
}
combined_balanced_samples <- do.call(rbind, all_balanced_samples)

# ------------------------------------ Test ------------------------------------
# Initialize an empty list to store results
results <- list()

# Loop through sites and sample sizes
sample_sizes <- c(300, 500, 2000, 5000)  # Define your desired sample sizes

sites <- c("Blois")
for (site in sites) {
  print(paste("Site:", site))
  site_data <- all_observations[[site]]
  site_data <- site_data %>% filter(stands == "deciduous")
  site_data <- site_data[sample(nrow(site_data), 15000), ]
  
  site_data_scaled <- site_data %>%
    dplyr::select(-lon, -lat, -stands) %>%
    scale()
  
  structural_dist_matrix <- dist(site_data_scaled)
  coords <- cbind(site_data$lon, site_data$lat)
  spatial_dist_matrix <- dist(coords)
  
  combined_dist_matrix_full <- 0.5 * as.matrix(structural_dist_matrix) +
    0.5 * as.matrix(spatial_dist_matrix)
  
  for (sample_size in sample_sizes) {
    print(paste("Processing sample size:", sample_size))
    
    selected_points <- numeric(0)
    remaining_points <- 1:nrow(site_data)
    
    first_point <- sample(remaining_points, 1)
    selected_points <- c(selected_points, first_point)
    remaining_points <- setdiff(remaining_points, first_point)
    
    pb <- progress_bar$new(
      total = sample_size,
      format = "Progress [:bar] :percent :eta",
      clear = TRUE,
      width = 50
    )
    
    for (i in 2:sample_size) {
      dist_to_selected <- apply(combined_dist_matrix_full[remaining_points, ], 1, function(x) min(x[selected_points]))
      next_point <- remaining_points[which.max(dist_to_selected)]
      selected_points <- c(selected_points, next_point)
      remaining_points <- setdiff(remaining_points, next_point)
      pb$tick()
    }
    
    balanced_sample <- site_data[selected_points, ]
    coords_sampled <- cbind(balanced_sample$lon, balanced_sample$lat)
    
    neighbors_sampled <- dnearneigh(coords_sampled, 0, 1000)
    lw_sampled <- nb2listw(neighbors_sampled)
    
    morans_I_sampled <- moran.test(balanced_sample$lidar_lai - balanced_sample$s2_lai, lw_sampled)
    geary_sampled <- geary.test(balanced_sample$lidar_lai - balanced_sample$s2_lai, lw_sampled)
    
    # Save results for this site and sample size
    results[[paste(site, sample_size, sep = "_")]] <- list(
      sample_size = sample_size,
      site = site,
      morans_I = morans_I_sampled,
      geary = geary_sampled
    )
    
    # Save plot
    balanced_sample_plot <- ggplot() +
      geom_point(data = site_data, aes(x = lon, y = lat), color = "gray", alpha = 0.5) +
      geom_point(data = balanced_sample, aes(x = lon, y = lat), color = "red", size = 3) +
      labs(title = paste("Sample Size:", sample_size, "at Site:", site)) +
      theme_minimal() +
      theme(legend.position = "none")
    ggsave(
      filename = paste0(output_dir, "/balanced_sample_", site, "_", sample_size, ".png"),
      plot = balanced_sample_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

# Write results to a text file
output_file <- paste0(output_dir, "/spatial_autocorrelation_results.txt")
sink(output_file)
for (result_name in names(results)) {
  cat("\n--- Results for:", result_name, "---\n")
  cat("Sample Size:", results[[result_name]]$sample_size, "\n")
  cat("Site:", results[[result_name]]$site, "\n")
  cat("Moran's I:\n")
  print(results[[result_name]]$morans_I)
  cat("\nGeary's C:\n")
  print(results[[result_name]]$geary)
  cat("\n")
}
sink()

print("All results saved to text file.")


# ----------------------------------- SYS --------------------------------------
all_balanced_sys_samples <- list()
for (site in sites) {
  # site <- "Blois"
  print(paste("Site:", site))
  site_data <- all_observations[[site]]
  site_data <- site_data %>% filter(stands == "deciduous")
  
  # Sort by spatial or structural variables for systematic sampling
  site_data <- site_data %>% arrange(lon, lat)
  
  # Set desired sample size
  sample_size <- 50 # 500 1000 2000 5000 10000
  
  # Calculate systematic sampling interval
  sampling_interval <- nrow(site_data) / sample_size
  print(sampling_interval)
  indices <- round(seq(1, nrow(site_data), by = sampling_interval))
  
  # Select rows based on indices
  balanced_sys_sample <- site_data[indices, ]
  all_balanced_sys_samples[[site]] <- balanced_sys_sample
  
  # detrended_data <- balanced_sys_sample
  # for (var in names(detrended_data)) {
  #   if (is.numeric(detrended_data[[var]])) {
  #     detrended_data[[var]] <- residuals(lm(detrended_data[[var]] ~ detrended_data$lon + detrended_data$lat))
  #   }
  # }
  
  # Save the sampled data
  # saveRDS(balanced_sample, paste0(output_dir, "/balanced_sample_", site, ".rds"))
  
  # Spatial autocorrelation tests
  coords_sampled <- cbind(balanced_sys_sample$lon, balanced_sys_sample$lat)
  # 500m
  # neighbors_sampled <- dnearneigh(coords_sampled, 0, 500) 
  # lw_sampled <- nb2listw(neighbors_sampled)
  # morans_I_sampled <- moran.test(balanced_sys_sample$lidar_lai - balanced_sys_sample$s2_lai, lw_sampled)
  # print(paste("Moran's I (500m) for site", site, ":"))
  # print(morans_I_sampled)
  # 1000m
  # neighbors_sampled <- dnearneigh(coords_sampled, 0, 1000)
  # lw_sampled <- nb2listw(neighbors_sampled)
  # morans_I_sampled <- moran.test(balanced_sys_sample$lidar_lai - balanced_sys_sample$s2_lai, lw_sampled)
  # print(paste("Moran's I (1000m) for site", site, ":"))
  # print(morans_I_sampled)
  
  k <- 10
  neighbors <- knearneigh(coords_sampled, k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  # morans_I_sampled <- localmoran(sampled_points$lidar_lai - sampled_points$s2_lai, lw_sampled)
  morans_I_sampled <- moran.test(balanced_sys_sample$random_pad - balanced_sys_sample$s2_lai, weights)
  # print(paste("Moran's I (500m) for site", site, ":"))
  print(morans_I_sampled)
  
  # Plot
  # balanced_sample_plot <- ggplot() +
  #   geom_point(data = site_data, aes(x = lon, y = lat), color = "gray", alpha = 0.5) +
  #   # Highlight selected points
  #   geom_point(data = balanced_sys_sample, aes(x = lon, y = lat), color = "red", size = 3) +
  #   labs(title = "Systematic Sampling") +
  #   theme_bw() +
  #   theme(legend.position = "none")
  # print(balanced_sample_plot)
  # ggsave(
  #   filename = paste0(output_dir, "/balanced_sys_sample_", site, ".png"),
  #   plot = balanced_sample_plot,
  #   width = 8,  # Width in inches
  #   height = 6,  # Height in inches
  #   dpi = 300  # High resolution
  # )
}
combined_balanced_sys_samples <- do.call(rbind, all_balanced_sys_samples)

# for (site in sites){
#   plot_raincloud_one_site(site,
#                           all_observations[[site]],
#                           all_predictors_filtered[[site]],
#                           all_balanced_sys_samples[[site]],
#                           output_dir)  
# }

# ---------------------------------- SBS ---------------------------------------
all_sbs_samples <- list()
sites <- c("Aigoual", "Blois")
for (site in sites) {
  sampled_points <- NA
  cat("Site:", site, "\n")
  # Filter data for the current site and prepare spatial data
  site_data <- all_observations[[site]] %>% filter(stands == "deciduous")
  df_sf <- st_as_sf(site_data, coords = c("lon", "lat"), crs = 32631)
  df_sf <- df_sf %>%
    mutate(lai_diff = lidar_lai - s2_lai)
  lai_diff_shuffled <- sample(df_sf$lai_diff, nrow(df_sf))
  df_sf <- df_sf %>% mutate(lai_diff = lai_diff_shuffled)
  df_sf <- df_sf %>%
    mutate(lai_diff_bin = cut(lai_diff, breaks = seq(min(lai_diff), max(lai_diff), length.out = 10), include.lowest = TRUE)) %>%
    group_by(lai_diff_bin) %>%
    mutate(bin_weight = 1 / n()) %>%
    ungroup()
  
  # target_lai_diff <- all_predictors_filtered[[site]]$lidar_lai - all_predictors_filtered[[site]]$s2_lai
  # target_density <- density(target_lai_diff)
  # df_sf <- df_sf %>%
  #   mutate(lai_diff_density = approx(x = target_density$x,
  #                                    y = target_density$y,
  #                                    xout = lai_diff,
  #                                    rule = 2)$y)
  # df_sf <- df_sf %>%
  #   mutate(lai_diff_weight = lai_diff_density / sum(lai_diff_density, na.rm = TRUE))
  # print(summary(df_sf$lai_diff_density))
  
  # Compute binning for lai_diff to approximate a uniform distribution
  n_bins <- 10  # Number of bins for the uniform distribution
  bin_edges <- seq(min(df_sf$lai_diff, na.rm = TRUE), max(df_sf$lai_diff, na.rm = TRUE), length.out = n_bins + 1)

  df_sf <- df_sf %>%
    mutate(lai_diff_bin = cut(lai_diff, breaks = bin_edges, include.lowest = TRUE)) %>%
    group_by(lai_diff_bin) %>%
    mutate(bin_weight = 1 / n()) %>%  # Inverse of the count in each bin
    ungroup() %>%
    mutate(lai_diff_weight = bin_weight / sum(bin_weight, na.rm = TRUE))  # Normalize weights
  
  # Create 500m grid
  grid <- st_make_grid(df_sf, cellsize = 500)
  intersections <- st_intersects(grid, df_sf, sparse = FALSE)
  valid_grid_cells <- grid[apply(intersections, 1, any), ]
  grid_centers <- st_centroid(valid_grid_cells)
  
  # Compute global weights based on the difference
  # Discretize
  # df_sf <- df_sf %>%
  #   mutate(lai_diff_weight = 1 / (table(lai_diff)[as.character(lai_diff)] + 1))
  # unique(table(df_sf$lai_diff))
  
  # Sampling with global weights across valid grid cells
  sampled_points <- lapply(1:length(valid_grid_cells), function(i) {
    cell <- valid_grid_cells[i]
    points_in_cell <- df_sf[st_within(df_sf, cell, sparse = FALSE), ]

    if (nrow(points_in_cell) > 0) {
      # Sample one point per grid cell using global weights
      sampled_point <- points_in_cell[sample(1:nrow(points_in_cell), 1, prob = points_in_cell$lai_diff_weight), ]
      # sampled_point <- points_in_cell[sample(1:nrow(points_in_cell), 1), ]
      # points_in_cell <- points_in_cell %>%
      #   mutate(lai_diff_transformed = scales::rescale(lai_diff, to = c(0, 1)))
      # sampled_point <- points_in_cell[sample(1:nrow(points_in_cell), 1), ]
      return(sampled_point)
    }
    return(NULL)
  })
  
  # sampled_points <- lapply(1:length(grid_centers), function(i) {
  #   center <- grid_centers[i]
  #   distances <- st_distance(df_sf, center)
  #   closest_point <- df_sf[which.min(distances), ]
  #   return(closest_point)
  # })
  
  # Combine sampled points into a single sf object for the current site
  sampled_points <- do.call(rbind, sampled_points)
  sampled_points_sf <- st_as_sf(sampled_points)
  
  # Save the sampled points for the current site as an RDS file
  saveRDS(sampled_points_sf, file = file.path(output_dir, paste0(site, "_uni3_sampled_points.rds")))
  
  # Create and save the plot
  plot <- ggplot() +
    geom_sf(data = df_sf, aes(color = lai_diff), size = 0.5, alpha = 0.6) +
    geom_sf(data = sampled_points_sf, color = "red", size = 1, alpha = 0.9) +
    scale_color_viridis_c() +
    theme_minimal() +
    ggtitle(paste("Sampled Points for", site)) +
    theme(legend.position = "bottom")
  
  ggsave(filename = file.path(output_dir, paste0(site, "_uni3_sampled_points_plot.png")),
         plot = plot, width = 8, height = 6, dpi = 300)
  
  all_sbs_samples[[site]] <- sampled_points_sf
}
combined_sbs_samples <- do.call(rbind, all_sbs_samples)

site_test <- "Blois" # Aigoual Blois Mormal
coords <- all_sbs_samples[[site_test]]$geometry
neighbors_sampled <- dnearneigh(coords, 0, 700)
lw_sampled <- nb2listw(neighbors_sampled)
k <- 10
neighbors <- knearneigh(coords, k = k)
lw_sampled <- nb2listw(knn2nb(neighbors), style = "W")
# morans_I_sampled <- localmoran(sampled_points$lidar_lai - sampled_points$s2_lai, lw_sampled)
morans_I_sampled <- moran.test(all_sbs_samples[[site_test]]$lidar_lai - all_sbs_samples[[site_test]]$s2_lai, lw_sampled)
# print(paste("Moran's I (500m) for site", site, ":"))
print(morans_I_sampled)



# sampled_points_sf <- st_sample(sampled_points, size = length(sampled_points))

Metrics::bias(all_observations$Blois$lidar_lai - all_observations$Blois$s2_lai,
              sampled_points$lidar_lai - sampled_points$s2_lai)
Metrics::bias(all_observations$Blois$vci, sampled_points$vci)


valid_points <- st_intersection(grid, st_convex_hull(df_sf))
# sampled_points <- st_sample(valid_points, size = 3)
sampled_points <- st_centroid(valid_points)
st_crs(sampled_points) <- st_crs(df_sf)
ggplot() +
  geom_sf(data = df_sf, color = "blue", size = 3) +  # Original points in blue
  geom_sf(data = sampled_points_sf, color = "red", size = 4) +  # Sampled points in red
  theme_minimal() +
  ggtitle("Sampled Points on Original Data") +
  theme(axis.text = element_blank(), axis.title = element_blank())



set1 <- train_data[, c("lcv", "mean", "vci", "s2_lai")]
set2 <- train_data[, c("deltaLAI")]
cca_result <- cancor(set1, set2)
cca_result$cor
canonical_x <- as.data.frame(as.matrix(set1) %*% cca_result$xcoef)
canonical_y <- as.data.frame(as.matrix(set2) %*% cca_result$ycoef)
ggplot(data = data.frame(Can1_X = canonical_x$V, Can1_Y = canonical_y$V1),
       aes(x = Can1_X, y = Can1_Y)) +
  geom_point() +
  labs(title = "First Canonical Variables",
       x = "Canonical Variable 1 (Set 1)",
       y = "Canonical Variable 1 (Set 2)") +
  theme_minimal()

ggplot(train_data, aes(x = mean, y = deltaLAI)) + geom_point() + geom_smooth(method = "loess") +
  labs(title = "Relationship between LCV and deltaLAI", x = "LCV", y = "deltaLAI")

all_sbs_samples <- list()
for (site in sites){
  all_sbs_samples[[site]] <- readRDS(file.path(output_dir, 
                                               paste0(site, "_sampled_points.rds")))
}

# ---------------------------------- RF ----------------------------------------
site <- "Mormal" # Aigoual Blois Mormal
# 1:102 103:244 245:682 Aigoual Blois Mormal all_predictors_variogram
# train_data <- all_sbs_samples[[site]]
# train_data <- all_balanced_sys_samples[[site]]
train_data <- combined_balanced_sys_samples
# train_data <- all_balanced_sys_samples[[site]]
# train_data <- all_predictors_variogram[[site]]
# mat <- dist(train_data[, c("lon", "lat")]) 
# mat <- as.matrix(mat)
# train_data <- train_data[train_data$stands == "deciduous", ]
# train_data$lat = NULL
# train_data$lon = NULL
train_data$stands <- NULL
train_data$geometry <- NULL
train_data$lai_diff <- NULL
train_data$lai_diff_transformed <- NULL
train_data$lai_diff_weight <- NULL
train_data$lidar_strata <- NULL
train_data$lai_diff_bin <- NULL
train_data$lai_diff_density <- NULL
train_data$bin_weight <- NULL
train_data$pad_minus_lad <- NULL
# train_data <- train_data[, -c(18:55)] # 1:60 61:120 121:180

# train_data <- all_observations$Aigoual[, -c(16:55)]
# train_data <- combined_predictors_variogram[245:682, -c(16:55)]
# train_data <- train_data[sample(nrow(train_data), 60), ]
train_data$deltaLAI <- train_data$random_pad - train_data$s2_lai
train_data <- train_data[,
                         c("deltaLAI", setdiff(names(train_data), "deltaLAI"))]
# train_data$lidar_lai <- NULL
# train_data$s2_lai <- NULL

corrplot::corrplot(cor(train_data),
                   method = "number",
                   type = "upper")

cor_matrix <- abs(cor(train_data))
cor_hclust <- hclust(as.dist(1 - cor_matrix))
cor_matrix <- cor_matrix[cor_hclust$order, cor_hclust$order]
corrplot(cor_matrix,
         method = "number",
         type = "full",
         order = "hclust",
         hclust.method = "ward.D",
         addrect = 3) # "complete" "ward.D"

# test_data <- all_predictors_filtered[[site]] # -c(18:20)
# test_data <- all_predictors_filtered$Aigoual # -c(18:20)
test_data <- combined_predictors_filtered
test_data$stands = NULL
# test_data <- combined_predictors_filtered # all_predictors_excluded
test_data$deltaLAI <- test_data$random_pad - test_data$s2_lai
test_data <- test_data[,
                       c("deltaLAI", setdiff(names(test_data), "deltaLAI"))]
# test_data$lidar_lai <- NULL
# test_data$s2_lai <- NULL

X <- train_data[, c(-1, -2)]
Y <- train_data[, 1]

# best_mtry <- plot_rf_mtry_evolution(train_data, X, formula)
mtry = ceiling((ncol(X)) / 2)

# formula <- as.formula(paste("deltaLAI ~ . - lidar_lai"))
formula <- as.formula(paste("deltaLAI ~ . - pad_minus_lai - random_pad"))
formula <- deltaLAI ~ lcv + mean + depth
rf_model <- ranger(
  formula = formula,
  data = train_data,
  num.trees = 500,
  # mtry = best_mtry,
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
)
# model <- rf_spatial(data = train_data, 
#                     dependent.variable.name = "deltaLAI", 
#                     predictor.variable.names = setdiff(names(train_data), "deltaLAI"),
#                     distance.matrix = mat)

oob_mse <- round(rf_model$prediction.error, 2)
oob_r_squared <- round(rf_model$r.squared, 2)
print(oob_mse)
print(oob_r_squared)
varImp <- sort(rf_model$variable.importance, decreasing = TRUE)
barplot(varImp)

pred <- predict(rf_model, data = test_data)
new_s2_plus_diff <- pred$predictions + test_data$s2_lai #[test_data$max > 2]

plot_data <- data.frame(
  S2_plus_deltaLAI = new_s2_plus_diff,
  LiDAR_LAI = test_data$random_pad
)
lm_fit <- lm(LiDAR_LAI ~ S2_plus_deltaLAI, data = plot_data)
r_squared <- summary(lm_fit)$r.squared
mse <- mean((plot_data$LiDAR_LAI - predict(lm_fit, plot_data))^2)
scatter_plot <- ggplot(plot_data, aes(x = S2_plus_deltaLAI, y = LiDAR_LAI)) +
  geom_point(alpha = 0.6, color = "blue") +  # Points with slight transparency
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +  # Optional regression line
  labs(
    title = "Scatterplot: S2 + deltaLAI vs LiDAR LAI",
    x = "S2 + deltaLAI",
    y = "LiDAR LAI"
  ) +
  annotate(
    "text", 
    x = max(plot_data$S2_plus_deltaLAI) * 0.8,  # Adjust the x-position
    y = max(plot_data$LiDAR_LAI) * 0.5,         # Adjust the y-position
    label = paste0("R² = ", round(r_squared, 2), "\nMSE = ", round(mse, 2)),
    color = "black",
    size = 5,
    hjust = 0
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "none"
  )
print(scatter_plot)

# plot_density_scatterplot(pred$predictions + test_data$s2_lai,
#                          test_data$lidar_lai[test_data$max > 2],
#                          "S2 LAI",
#                          "LiDAR LAI")

# Predictions on the training data
train_predictions <- predict(rf_model, data = train_data)
train_actual <- train_data$deltaLAI
train_mse <- mean((train_predictions$predictions - train_actual) ^ 2)
train_rss <- sum((train_predictions$predictions - train_actual) ^ 2)
train_tss <- sum((train_actual - mean(train_actual)) ^ 2)
train_r_squared <- 1 - train_rss / train_tss

# Predictions on the test data
predictions_excluded <- predict(rf_model, data = test_data)
actual_excluded <- test_data$deltaLAI # lidar_lai s2_lai deltaLAI
mse_excluded <- mean((predictions_excluded$predictions - actual_excluded) ^ 2)
rss <- sum((predictions_excluded$predictions - actual_excluded) ^ 2)
tss <- sum((actual_excluded - mean(actual_excluded)) ^ 2)
r_squared_excluded <- 1 - rss/tss

# Output results
cat("Training MSE:", round(train_mse, 2), "\n")
cat("Training R-squared:", round(train_r_squared, 2), "\n")
cat("MSE on excluded data:", round(mse_excluded, 2), "\n")
cat("R-squared on excluded data:", round(r_squared_excluded, 2), "\n")
cat("OOB MSE:", round(oob_mse, 2), "\n")
cat("OOB R-squared:", round(oob_r_squared, 2), "\n")

# PDP
for (var in names(rf_model$variable.importance)) {
  partialPlot(obj = rf_model,
              pred.data = train_data[, setdiff(names(train_data), "deltaLAI")],
              xname = var,
              target_var = "deltaLAI",
              show_plot = TRUE,
              save_path = file.path(output_dir, site)) # NULL pp_dir
}

# ------------------------------------------------------------------------------
# Curves
# oob_error_curve <- rf_model$prediction.error
# train_error_curve <- numeric(rf_model$num.trees)
# test_error_curve <- numeric(rf_model$num.trees)
# 
# # Loop through tree counts and calculate errors
# for (i in seq_len(rf_model$num.trees)) {
#   train_predictions <- predict(rf_model, data = train_data, num.trees = i)$predictions
#   test_predictions <- predict(rf_model, data = test_data, num.trees = i)$predictions
#   
#   # Calculate MSE for train and test sets
#   train_error_curve[i] <- mean((train_predictions - train_data$deltaLAI)^2)
#   test_error_curve[i] <- mean((test_predictions - test_data$deltaLAI)^2)
# }
# 
# # Combine errors into a single data frame for plotting
# error_data <- data.frame(
#   Trees = seq_len(rf_model$num.trees),
#   TrainError = train_error_curve,
#   TestError = test_error_curve,
#   OOBError = rep(oob_error_curve, length.out = rf_model$num.trees)
# )
# 
# error_data_melted <- melt(error_data, id.vars = "Trees", 
#                           variable.name = "ErrorType", 
#                           value.name = "MSE")
# ggplot(error_data_melted, aes(x = Trees, y = MSE, color = ErrorType)) +
#   geom_line(size = 1) +
#   labs(title = "Random Forest Error Curves",
#        x = "Number of Trees",
#        y = "Mean Squared Error (MSE)",
#        color = "Error Type") +
#   theme_minimal() +
#   scale_color_manual(values = c("blue", "red", "green"),
#                      labels = c("Train Error", "Test Error", "OOB Error")) +
#   theme(legend.position = "top")

# --------------------------------- LOOCV --------------------------------------
# Prepare test data
test_data <- all_predictors_filtered$Aigoual
test_data$stands = NULL
test_data$deltaLAI <- test_data$lidar_lai - test_data$s2_lai
test_data <- test_data[,
                       c("deltaLAI", setdiff(names(test_data), "deltaLAI"))]
# test_data$lidar_lai <- NULL
X_test <- test_data[, c(-1, -2)]  # Exclude deltaLAI and any non-predictor columns
Y_test <- test_data[, 1]          # Target variable: deltaLAI

predictions <- numeric(nrow(test_data))
formula <- deltaLAI ~ .
formula <- deltaLAI ~ lcv + mean
for (i in 1:nrow(test_data)) {
  train_subset <- test_data[-i, ]
  
  rf_model_loocv <- ranger(
    formula = formula,
    data = train_subset,
    num.trees = 500,
    importance = "permutation",
    min.node.size = 5,
    num.threads = 12
  )
  predictions[i] <- predict(rf_model_loocv, data = test_data[i, ])$predictions
}
mse <- mean((Y_test - predictions)^2)
r_squared <- 1 - sum((Y_test - predictions)^2) / sum((Y_test - mean(Y_test))^2)
print(paste("LOOCV MSE:", round(mse, 2)))
print(paste("LOOCV R-squared:", round(r_squared, 2)))

test_pred <- all_predictors_filtered$Aigoual
pred <- predict(rf_model_loocv, data = test_pred)
new_s2_plus_diff <- predictions + test_pred$s2_lai #[test_data$max > 2]

plot_data <- data.frame(
  S2_plus_deltaLAI = new_s2_plus_diff,
  LiDAR_LAI = test_pred$lidar_lai
)
lm_fit <- lm(LiDAR_LAI ~ S2_plus_deltaLAI, data = plot_data)
r_squared <- summary(lm_fit)$r.squared
mse <- mean((plot_data$LiDAR_LAI - predict(lm_fit, plot_data))^2)
scatter_plot <- ggplot(plot_data, aes(x = S2_plus_deltaLAI, y = LiDAR_LAI)) +
  geom_point(alpha = 0.6, color = "blue") +  # Points with slight transparency
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = TRUE) +  # Optional regression line
  labs(
    title = "Scatterplot: S2 + deltaLAI vs LiDAR LAI",
    x = "S2 + deltaLAI",
    y = "LiDAR LAI"
  ) +
  annotate(
    "text", 
    x = max(plot_data$S2_plus_deltaLAI) * 0.8,  # Adjust the x-position
    y = max(plot_data$LiDAR_LAI) * 0.5,         # Adjust the y-position
    label = paste0("R² = ", round(r_squared, 2), "\nMSE = ", round(mse, 2)),
    color = "black",
    size = 5,
    hjust = 0
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),  # Center the title
    legend.position = "none"
  )
print(scatter_plot)

plot(Y_test, predictions, 
     xlab = "Observed deltaLAI", 
     ylab = "Predicted deltaLAI", 
     main = "LOOCV Results",
     pch = 19, col = "blue")
abline(a = 0, b = 1, col = "red")

varImp <- sort(rf_model_loocv$variable.importance, decreasing = TRUE)
barplot(varImp)

vsurf_result <- VSURF(X_test, Y_test, mtry = 8, parallel = TRUE) 
interp_vars <- names(X_test)[vsurf_result$varselect.interp]
print(interp_vars)
# ------------------------------------------------------------------------------
# Prepare test data
test_data <- all_predictors_filtered$Blois
test_data$stands = NULL
test_data$deltaLAI <- test_data$lidar_lai - test_data$s2_lai
test_data <- test_data[,
                       c("deltaLAI", setdiff(names(test_data), "deltaLAI"))]
test_data$lidar_lai <- NULL

X_test <- test_data[, c(-1, -2)]  # Exclude deltaLAI and any non-predictor columns
Y_test <- test_data[, 1]          # Target variable: deltaLAI

# Initialize variables
mtry_values <- 1:ncol(X_test)  # Test mtry from 1 to the number of predictors
mtry_loocv_mse <- numeric(length(mtry_values))  # Store LOOCV MSE for each mtry

# Loop over each mtry value
for (mtry in mtry_values) {
  predictions <- numeric(nrow(test_data))
  
  # Perform LOOCV for the current mtry
  for (i in 1:nrow(test_data)) {
    train_subset <- test_data[-i, ]
    
    # Train Random Forest with the current mtry
    rf_model_loocv <- ranger(
      formula = deltaLAI ~ .,
      data = train_subset,
      num.trees = 500,
      mtry = mtry,
      importance = "none",
      min.node.size = 5,
      num.threads = 12
    )
    
    # Predict for the left-out observation
    predictions[i] <- predict(rf_model_loocv, data = test_data[i, ])$predictions
  }
  
  # Compute LOOCV MSE for the current mtry
  mtry_loocv_mse[mtry] <- mean((Y_test - predictions)^2)
}

# Find the best mtry
optimal_mtry <- which.min(mtry_loocv_mse)
print(paste("Optimal mtry (LOOCV):", optimal_mtry))

# Plot LOOCV MSE vs mtry
plot(mtry_values, 
     mtry_loocv_mse, 
     type = "b", 
     xlab = "mtry", 
     ylab = "LOOCV MSE", 
     main = "Optimization of mtry with LOOCV",
     pch = 19)
abline(v = optimal_mtry, col = "red", lty = 2)

# Train final model with optimal mtry
predictions_final <- numeric(nrow(test_data))
for (i in 1:nrow(test_data)) {
  train_subset <- test_data[-i, ]
  
  # Train Random Forest with the optimal mtry
  rf_model_loocv <- ranger(
    formula = deltaLAI ~ .,
    data = train_subset,
    num.trees = 500,
    mtry = optimal_mtry,
    importance = "permutation",
    min.node.size = 5,
    num.threads = 12
  )
  
  # Predict for the left-out observation
  predictions_final[i] <- predict(rf_model_loocv, data = test_data[i, ])$predictions
}

# Evaluate performance of the final model
final_mse <- mean((Y_test - predictions_final)^2)
final_r_squared <- 1 - sum((Y_test - predictions_final)^2) / sum((Y_test - mean(Y_test))^2)
print(paste("Final LOOCV MSE with optimal mtry:", round(final_mse, 2)))
print(paste("Final LOOCV R-squared with optimal mtry:", round(final_r_squared, 2)))

# Plot observed vs. predicted (final model)
plot(Y_test, predictions_final, 
     xlab = "Observed deltaLAI", 
     ylab = "Predicted deltaLAI", 
     main = "Final LOOCV Results with Optimal mtry",
     pch = 19, col = "blue")
abline(a = 0, b = 1, col = "red")

vsurf_result <- VSURF(X_test, Y_test, mtry = 8, parallel = TRUE) 
interp_vars <- names(X_test)[vsurf_result$varselect.interp]
print(interp_vars)
# ------------------------------------------------------------------------------
# VSURF
vsurf_formula <- deltaLAI ~ .
formula <- as.formula(paste("deltaLAI ~ . - pad_minus_lai - random_pad"))
train_data_vsurf <- train_data[, -2]
# train_data_vsurf <- test_data[, -2]
vsurf <- VSURF(vsurf_formula, train_data, ntree.interp = 500)
# plot(vsurf)
# length(vsurf$varselect.pred)
# names(X)[vsurf$varselect.pred]
# names(X)[-vsurf$varselect.pred]
vsurf_vars <- names(X)[vsurf$varselect.interp]
print(vsurf_vars)

# Optimal
# selected_vars <- names(train_data[, -1])[vsurf$varselect.interp]
selected_vars <- c("deltaLAI", vsurf_vars)
train_data_opt <- train_data[selected_vars]

# best_mtry_vsurf <- plot_rf_mtry_evolution(train_data_opt, formula)
formula <- deltaLAI ~ .
rf_model_opt <- ranger(
  formula = formula,
  data = train_data_opt,
  num.trees = 500,
  # mtry = ncol(train_data_opt) - 1,
  importance = "permutation",
  num.threads = 12
)

for (vsurf_var in vsurf_vars) {
  partialPlot(obj = rf_model_opt,
              pred.data = train_data_opt[, setdiff(names(train_data_opt), "deltaLAI")],
              xname = vsurf_var,
              target_var = "deltaLAI",
              show_plot = TRUE,
              save_path = file.path(output_dir, site)) # NULL pp_dir
}

oob_mse_opt <- round(rf_model_opt$prediction.error, 2)
r_squared_opt <- round(rf_model_opt$r.squared, 2)
varImp_opt <- sort(rf_model_opt$variable.importance, decreasing = TRUE)
barplot(varImp_opt)

# Predictions on the training data
opt_train_predictions <- predict(rf_model_opt, data = train_data_opt)
opt_train_actual <- train_data_opt$deltaLAI
opt_train_mse <- mean((opt_train_predictions$predictions - opt_train_actual) ^ 2)
opt_train_rss <- sum((opt_train_predictions$predictions - opt_train_actual) ^ 2)
opt_train_tss <- sum((opt_train_actual - mean(opt_train_actual)) ^ 2)
opt_train_r_squared <- 1 - opt_train_rss / opt_train_tss

# Predictions on the test data
opt_predictions_excluded <- predict(rf_model_opt, data = test_data)
opt_actual_excluded <- test_data$deltaLAI # lidar_lai s2_lai deltaLAI
opt_mse_excluded <- mean((opt_predictions_excluded$predictions - opt_actual_excluded) ^ 2)
opt_rss <- sum((opt_predictions_excluded$predictions - opt_actual_excluded) ^ 2)
opt_tss <- sum((opt_actual_excluded - mean(opt_actual_excluded)) ^ 2)
opt_r_squared_excluded <- 1 - opt_rss/opt_tss

# Output results
cat("Opt Training MSE:", round(opt_train_mse, 2), "\n")
cat("Opt Training R-squared:", round(opt_train_r_squared, 2), "\n")
cat("Opt MSE on excluded data:", round(opt_mse_excluded, 2), "\n")
cat("Opt R-squared on excluded data:", round(opt_r_squared_excluded, 2), "\n")
cat("Opt OOB MSE:", round(oob_mse_opt, 2), "\n")
cat("Opt OOB R-squared:", round(r_squared_opt, 2), "\n")

# Scatterplot
ggplot(data = train_data_opt, aes(x = deltaLAI, y = rf_model_opt$predictions)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +  # Scatter plot points
  geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Actual vs. Predicted Values of deltaLAI",
    x = "Actual deltaLAI",
    y = "Predicted deltaLAI"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Mean
# Determine the common x and y limits
pred <- predict(rf_model_opt, data = test_data)
new_s2_plus_diff <- pred$predictions + test_data$s2_lai
x_limits <- c(min(c(test_data$s2_lai, new_s2_plus_diff)), 
              max(c(test_data$s2_lai, new_s2_plus_diff)))

y_limits <- c(min(c(test_data$lidar_lai, test_data$lidar_lai)), 
              max(c(test_data$lidar_lai, test_data$lidar_lai)))
# c(0, 1/3, 2/3, 1) c(0, 1/5, 2/5, 3/5, 4/5, 1) [sample(nrow(all_predictors_filtered$Blois), 200), ]
train_data_qt <- test_data %>%
  mutate(class = cut(vci, breaks = quantile(vci, probs = c(0, 1/3, 2/3, 1)),
                     labels = c("Low", "Medium", "High"), include.lowest = TRUE))
before <- ggplot(data = train_data_qt, aes(x = s2_lai, y = test_data$lidar_lai, color = class)) +
  geom_point(size = 3, alpha = 0.7) +
  # geom_bin2d(bins = 200, aes(fill = class)) +
  geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  labs(
    title = "LiDAR LAI vs Sentinel-2 LAI with Quantile-based Classes",
    x = "Sentinel-2 LAI",
    y = "LiDAR LAI",
    color = "Class (Quantiles of vci)"
  ) +
  scale_color_manual(values = c("Low" = "steelblue", "Medium" = "goldenrod", "High" = "tomato")) +
  xlim(x_limits) +
  ylim(y_limits) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

# Lcv
train_data_qt <- test_data %>%
  mutate(class = cut(vci, breaks = quantile(vci, probs = c(0, 1/3, 2/3, 1)),
                     labels = c("Low", "Medium", "High"), include.lowest = TRUE))
after <- ggplot(data = train_data_qt, aes(x = new_s2_plus_diff, y = test_data$lidar_lai, color = class)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter plot points with color by class
  geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  labs(
    title = "LiDAR LAI vs Corrected Sentinel-2 LAI with Quantile-based Classes",
    x = "Corrected Sentinel-2 LAI",
    y = "LiDAR LAI",
    color = "Class (Quantiles of vci)"
  ) +
  scale_color_manual(values = c("Low" = "steelblue", "Medium" = "goldenrod", "High" = "tomato")) +
  xlim(x_limits) +
  ylim(y_limits) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

grid.arrange(before, after, ncol = 2)
# Add a site identifier column to differentiate between sites
# train_data_opt <- train_data_opt %>%
#   mutate(
#     site = case_when(
#       row_number() >= 1 & row_number() <= 60 ~ "Aigoual",
#       row_number() >= 61 & row_number() <= 120 ~ "Blois",
#       row_number() >= 121 & row_number() <= 180 ~ "Mormal"
#     ),
#     class = cut(mean, breaks = quantile(mean, probs = c(0, 1/3, 2/3, 1)),
#                 labels = c("Low", "Medium", "High"), include.lowest = TRUE)
#   )
#
# # Create scatter plot with color by class and shape by site
# ggplot(data = train_data_opt, aes(x = deltaLAI, y = rf_model_opt$predictions, color = class, shape = site)) +
#   geom_point(size = 2, alpha = 0.7) +  # Scatter plot points with color by class and shape by site
#   geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
#   labs(
#     title = "Actual vs. Predicted Values of deltaLAI with Site and Quantile-based Classes",
#     x = "Actual deltaLAI",
#     y = "Predicted deltaLAI",
#     color = "Class (Quantiles of mean)",
#     shape = "Site"
#   ) +
#   scale_color_manual(values = c("Low" = "steelblue", "Medium" = "goldenrod", "High" = "tomato")) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     legend.position = "right"
#   )

# ------------------------------------------------------------------------------
# 1:102 103:244 245:682
# Run the function for each data subset and save results to text files
run_random_forest_analysis(combined_predictors_variogram,
                           1:102,
                           file.path(output_dir, "results_vsurf_sampled_vario_Aigoual.txt"))
run_random_forest_analysis(combined_predictors_variogram,
                           103:244,
                           file.path(output_dir, "results_vsurf_sampled_vario_Blois.txt"))
run_random_forest_analysis(combined_predictors_variogram,
                           245:682,
                           file.path(output_dir, "results_vsurf_sampled_vario_Mormal.txt"))

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
X = drf_data[, -c(1, 16:55)] # 16:53 17:54
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
  # mtry = 5,
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
names(VI_MMD) <- names(drf_data[, -c(1, 16:55)])
VI_MMD <- sort(VI_MMD, decreasing = T)
barplot(VI_MMD)

VI_count <- variable_importance(drf_modele)[,1]
names(VI_count) <- names(drf_data[, -c(1, 16:55)])
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
