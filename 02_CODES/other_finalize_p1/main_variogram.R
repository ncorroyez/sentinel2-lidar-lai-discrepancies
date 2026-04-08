# ---
# title: "main_variogram.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-27"
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
library("sf")
library("vegan")
library("spdep")
library("blockCV")
library("cluster")
library("mgcv")
source("libraries/functions_general_tools.R")

# Pre-processing Parameters
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/variogram"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)

# Function to check distance
check_distance <- function(points, min_distance) {
  n <- nrow(points)
  to_keep <- logical(n)  # Logical vector to keep track of points to retain
  to_keep[1] <- TRUE      # Always keep the first point

  # Iterate through points starting from the second one
  for (i in 2:n) {
    # Calculate distances from the current point to all previously kept points
    dists <- spDists(points[which(to_keep), , drop = FALSE], points[i, , drop = FALSE], longlat = FALSE)

    # Keep the point if all distances are greater than or equal to min_distance
    to_keep[i] <- all(dists >= min_distance)
  }

  return(points[to_keep, , drop = FALSE])  # Ensure to return a matrix
}

select_max_dist_points <- function(points, n_select = 60) {
  n <- nrow(points)

  # Check if there are enough points
  if (n <= n_select) {
    warning("The number of points is less than or equal to the selection count; returning all points.")
    return(points)
  }

  # Start by randomly selecting the first point
  selected_indices <- integer(n_select)
  selected_indices[1] <- sample(1:n, 1)
  selected_points <- points[selected_indices[1], , drop = FALSE]

  for (i in 2:n_select) {
    # Calculate distances from each unselected point to the selected points
    dists <- spDists(points, selected_points, longlat = FALSE)
    min_dists <- apply(dists, 1, min)  # Minimum distance to any selected point

    # Exclude already selected points
    min_dists[selected_indices[1:(i-1)]] <- -Inf

    # Select the point with the maximum of the minimum distances
    selected_indices[i] <- which.max(min_dists)
    selected_points <- rbind(selected_points, points[selected_indices[i], , drop = FALSE])
  }

  return(points[selected_indices, , drop = FALSE])  # Return the selected points
}

select_fixed_distance_points <- function(raster, spacing = 500, n_select = 60) {
  # Ensure raster has a valid CRS
  raster_crs <- crs(raster, proj = TRUE)

  # Get raster extent and convert it to sf geometry
  raster_extent <- terra::ext(raster)  # Extract extent
  raster_bbox <- st_bbox(c(
    xmin = raster_extent$xmin,
    ymin = raster_extent$ymin,
    xmax = raster_extent$xmax,
    ymax = raster_extent$ymax
  ), crs = raster_crs)
  raster_sf_extent <- st_as_sfc(raster_bbox)  # Create sf object for extent

  # Generate a grid of points spaced `spacing` meters apart
  grid_points <- st_make_grid(raster_sf_extent, cellsize = spacing, what = "centers")

  # Convert grid to sf object
  grid_sf <- st_sf(geometry = grid_points, crs = raster_crs)

  # Convert raster to sf polygons and mask grid
  raster_polygons <- as.polygons(raster, na.rm = TRUE)  # Remove NA areas
  raster_sf <- st_as_sf(raster_polygons)  # Convert to sf object

  # Clip grid points to raster area
  valid_points <- grid_sf[st_intersects(grid_sf, raster_sf, sparse = FALSE), ]

  # Check if there are enough points
  if (nrow(valid_points) < n_select) {
    stop("Not enough valid points available to select the required number.")
  }

  # Randomly select exactly `n_select` points
  selected_points <- valid_points[sample(nrow(valid_points), n_select), ]

  return(selected_points)
}

select_equidistant_points <- function(points, n_select = 60, min_dist = 500) {
  # Ensure points is a matrix or data.frame with proper coordinates
  if (!is.matrix(points) && !is.data.frame(points)) {
    stop("Input 'points' must be a matrix or data.frame with coordinates.")
  }

  n <- nrow(points)

  # Check if there are enough points to select from
  if (n < n_select) {
    stop("Not enough points to select the desired number while maintaining the distance constraint.")
  }

  # Randomly pick the first point
  selected_indices <- integer(n_select)
  selected_indices[1] <- sample(1:n, 1)

  selected_points <- points[selected_indices[1], , drop = FALSE]

  for (i in 2:n_select) {
    # Compute distances of unselected points to the current set of selected points
    dists <- spDists(points, selected_points, longlat = FALSE)
    min_dists <- apply(dists, 1, min)  # Minimum distance to any selected point

    # Exclude already selected points and those too close
    valid_indices <- which(min_dists >= min_dist & !(1:n %in% selected_indices[1:(i - 1)]))

    if (length(valid_indices) == 0) {
      stop("Cannot find enough points meeting the minimum distance constraint.")
    }

    # Randomly select the next point from the valid candidates
    selected_indices[i] <- sample(valid_indices, 1)
    selected_points <- rbind(selected_points, points[selected_indices[i], , drop = FALSE])
  }

  return(points[selected_indices, , drop = FALSE])
}

for (site in sites){
  lai_lidar_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                            "lidarlai_res_10_m.tif"))
  lai_s2_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                         "s2lai_summer_res_10_m.tif"))
  lcv_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                      "lcv_res_10_m.tif"))
  lskew_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                        "lskew_res_10_m.tif"))
  mean_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                       "mean_res_10_m.tif"))
  rumple_raster <- terra::rast(file.path(results_dir, site, metrics_dir,
                                         "rumple_res_10_m.tif"))

  # Convert raster to points (remove NA values)
  r_points <- as.data.frame(c(
    lai_lidar_raster - lai_s2_raster
    # lai_s2_raster,
    # lcv_raster,
    # mean_raster,
    # rumple_raster,
    # lskew_raster
  ), xy = TRUE, na.rm = TRUE)

  coordinates(r_points) <- ~ x + y

  g <- gstat(id = "lai_model", formula = V1 ~ 1, data = r_points)
  v <- variogram(g, width = 10) # width 10m

  # Create the variogram plot using ggplot2
  variogram_plot <- ggplot(as.data.frame(v), aes(x = dist, y = gamma)) +
    geom_point() +
    geom_line() +
    labs(title = paste(site, "Empirical Variogram"), x = "Distance", y = "Semivariance") +
    theme_bw()
  print(variogram_plot)
  # # Save the plot to a file
  ggsave(file.path(output_dir, paste0(site, "_variogram_deciduous_plot.png")),
         plot = variogram_plot,
         width = 1920, height = 1080, units = "px", dpi = 300)
  # a
  # coords <- as.matrix(coordinates(r_points))
  # filtered_points <- check_distance(coords, 500)
  # filtered_points <- select_max_dist_points(coords, 60)
  # filtered_points <- select_equidistant_points(coords, n_select = 60, min_dist = 500)
  # geojson_path <- file.path(results_dir, site, "points_60_maxdist.geojson")
  # st_write(st_as_sf(as.data.frame(filtered_points),
  #                   coords = c("x", "y"),
  #                   crs = st_crs(lai_lidar_raster)),
  #          geojson_path, driver = "GeoJSON", delete_layer = TRUE)
  # write.csv(filtered_points, file.path(results_dir, site, "points_60_maxdist.csv"))
}
a
# random_indices <- sample(1:nrow(r_points), size = 2000, replace = FALSE)

# Subset the data frame to get the randomly selected points
# r_points <- r_points[random_indices, ]

# Create spatial data frame
coordinates(r_points) <- ~x + y

# Ensure V1 is the correct column for LAI
colnames(r_points@data)[1] <- "lai_lidar"  # Rename the first column to 'lai' for clarity
# colnames(r_points@data)[2:6] <- c("lai_s2", "lcv", "mean", "rumple", "lskew")  # Rename remaining columns

# Set the number of folds and the size of blocks (e.g., 5-fold CV)
# Assuming 'r_points' is your SpatialPointsDataFrame from the previous code
# num_folds <- 5
# folds <- kmeans(as.matrix(coordinates(r_points)), centers = num_folds)
#
# # Add the fold information to the data frame
# r_points$fold <- as.factor(folds$cluster)
#
# # List to store GAM results
# gam_results <- list()
#
# # Loop through each fold for GAM modeling
# for (i in unique(r_points$fold)) {
#   # Subset the data for the current fold
#   fold_data <- r_points[r_points$fold == i, ]
#
#   # Check if there are enough points
#   if (nrow(fold_data) < 3) {
#     next
#   }
#
#   # Fit the GAM with spatial smooths (s() indicates smooth terms)
#   gam_model <- gam(lai_lidar ~ s(x, y) + lai_s2 + lcv + mean + rumple + lskew, data = fold_data)
#
#   # Store the model summary
#   gam_results[[as.character(i)]] <- summary(gam_model)
#
#   # Check residuals for spatial autocorrelation
#   residuals_gam <- residuals(gam_model)
#
#   # Create a neighbors list for the fold
#   neighbors <- knn2nb(knearneigh(coordinates(fold_data), k = 2))
#
#   # Create a listw object for spatial weights
#   listw <- nb2listw(neighbors, style = "W")
#
#   # Calculate Moran's I for the residuals
#   moran_result <- moran.test(residuals_gam, listw)
#
#   # Store the Moran's I result
#   gam_results[[paste("Fold", i, "Moran's I")]] <- moran_result
#
#   # Calculate Geary's C on the residuals
#   geary_test <- geary.test(residuals_gam, listw)
#
#   # Print Geary's C test results
#   print(geary_test)
# }
#
# # Print GAM results and Moran's I for each fold
# for (i in seq_along(gam_results)) {
#   cat("GAM Model Summary for Fold:", names(gam_results)[i], "\n")
#   print(gam_results[[i]])
# }
#
# # Print Moran's I results
# for (i in unique(r_points$fold)) {
#   cat("Moran's I for Fold:", i, "\n")
#   print(gam_results[[paste("Fold", i, "Moran's I")]])
# }
#
# a





# Use gstat to model LAI as a function of LCV, mean, etc.
g <- gstat(id = "lai_model", formula = lai_lidar ~ 1, data = r_points)
v <- variogram(g, width = 10) # width 10m

# Create the variogram plot using ggplot2
variogram_plot <- ggplot(as.data.frame(v), aes(x = dist, y = gamma)) +
  geom_point() +
  geom_line() +
  labs(title = paste(site, "Empirical Variogram"), x = "Distance", y = "Semivariance") +
  theme_bw()
print(variogram_plot)
# Save the plot to a file
# ggsave(file.path(output_dir, paste0(site, "_variogram_plot.png")),
#        plot = variogram_plot,
#        width = 1920, height = 1080, units = "px", dpi = 300)
# }

# Get the coordinates
coords <- as.matrix(coordinates(r_points))
# nb <- knn2nb(knearneigh(coords, k = 4))  # Nearest neighbor (4 nearest)
# listw <- nb2listw(nb)

# Moran's I for LAI (or any variable of interest)
# moran_i <- moran.test(r_points$V1, listw)
# print(moran_i)

# Filter points that are at least x meters apart
coords <- as.matrix(coordinates(r_points))
filtered_points <- check_distance(coords, 500)
a
# Create a new SpatialPointsDataFrame
filtered_r_points <- SpatialPointsDataFrame(filtered_points,
                                            data = as.data.frame(r_points@data[rownames(filtered_points),]))

# Create a distance matrix for spatial coordinates and another for LAI values
coords_filtered <- as.matrix(filtered_points)

# Extract the filtered row indices based on coordinates
filtered_indices <- match(rownames(coords_filtered), rownames(r_points@data))

# Ensure these indices are valid
filtered_indices <- na.omit(filtered_indices)  # Remove any NA values in case of mismatch

# Extract LAI values for these filtered indices
lai_values_filtered <- r_points@data$V1[filtered_indices]

# Create distance matrices for the filtered points
geo_dist_filtered <- dist(coords_filtered)
lai_dist_filtered <- dist(lai_values_filtered)

# Perform Mantel test
mantel_test <- mantel(geo_dist_filtered, lai_dist_filtered)
print(mantel_test)

nb <- knn2nb(knearneigh(coords_filtered, k = 4))  # Nearest neighbor (4 nearest)
listw <- nb2listw(nb)

# Moran's I for LAI (or any variable of interest)
moran_i <- moran.test(filtered_r_points@data$`r_points@data[rownames(filtered_points), ]`, listw)
print(moran_i)
# }
