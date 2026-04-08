# ---
# title: "main_RF_sampling_depth.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-11-29"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

library(terra)
library(sgsR)
library(ggplot2)
library(sf)
library(sp)
library(spdep)
library(dplyr)
library(geojsonio)
library(ranger)
library(VSURF)
library(progress)
library(corrplot)
library(patchwork)
library(ncf)
library(miscTools)

source("libraries/functions_general_tools.R")
source("libraries/functions_plots.R")

# Pre-processing Parameters
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual", "Blois")
# sites <- "Aigoual" # Mormal Blois Aigoual
results_dir <- "../03_RESULTS"
output_dir <- "../04_FIGURES/variogram"
pp_dir <- "../04_FIGURES/pp"
forest_composition <- "Deciduous_Only" # Full_Composition Deciduous_Only Not_Masked
metrics_dir <- file.path("Metrics", forest_composition)
set.seed(42)

# Lists
all_observations_one_random <- list()
all_observations_multiple_randoms <- list()
all_observations <- list()
all_predictors_filtered <- list()

# --------------------------------Load data ------------------------------------
for (site in sites) {
  cat("Site:", site, "\n")
  site_data <- load_metrics_with_pad(site, results_dir, "Deciduous_Only")
  all_observations[[site]] <- site_data 
  
  geojson_file <- file.path(results_dir, site, "data_utm31n.geojson")
  geo_data <- geojson_read(geojson_file, what = "sp")
  # coord_x <- geo_data@data$coord_x_utm31n
  # coord_y <- geo_data@data$coord_y_utm31n
  coord_x <- round_to_odd_multiple_of_5(geo_data@data$coord_x_utm31n)
  coord_y <- round_to_odd_multiple_of_5(geo_data@data$coord_y_utm31n)
  coord_df <- data.frame(coord_x = coord_x, coord_y = coord_y)
  
  closest_matches <- list()
  lidarlai_raster <- terra::rast(file.path(results_dir, site, metrics_dir, "lidarlai_res_10_m.tif"))
  for (j in 1:nrow(coord_df)) {
    target_coord <- coord_df[j, ]
    distances <- sqrt((site_data$lat - target_coord$coord_y)^2 +
                        (site_data$lon - target_coord$coord_x)^2)
    closest_index <- which.min(distances)
    closest_match <- cbind(site_data[closest_index, ], target_coord)
    lidarlai_value <- terra::extract(lidarlai_raster, data.frame(lon = closest_match$lon, lat = closest_match$lat))
    if (is.na(lidarlai_value[[2]])) {
      closest_match[] <- NA
    }
    closest_matches[[j]] <- closest_match
  }
  closest_data <- do.call(rbind, closest_matches)
  # closest_data <- closest_data[, !colnames(closest_data) %in% c("coord_x", "coord_y")]
  
  # closest_data <- merge(coord_df, site_data, 
  #                       by.x = c("coord_x", "coord_y"), 
  #                       by.y = c("lon", "lat"), 
  #                       all.x = TRUE)
  
  # closest_data <- site_data[site_data$lon %in% coord_df$coord_x & site_data$lat %in% coord_df$coord_y, ]
  closest_data <- closest_data[
    # closest_data$stands %in% c("oak", "beech", "deciduous", "poplar") &
    !is.na(closest_data$lidar_lai) &
      closest_data$mean >= 2,
  ]
  all_predictors_filtered[[site]] <- closest_data
}
combined_predictors_filtered <- do.call(rbind, all_predictors_filtered)
combined_predictors_all <- do.call(rbind, all_observations)

# -------------------------------- Sampling ------------------------------------
all_sites_closest_data <- list()
# sites <- c("Aigoual", "Blois", "Mormal")
for (site in sites) {
  cat("Site:", site, "\n")
  nSamp <- switch(site,
                  "Aigoual" = 120,
                  "Blois" = 120,
                  "Mormal" = 120)
  mindist <- switch(site,
                    "Aigoual" = 320,
                    "Blois" = 410,
                    "Mormal" = 500)
  # nSamp <- 500
  # mindist <- NULL
  tif_files <- list.files(file.path(results_dir, site, metrics_dir), 
                          pattern = "\\.tif$", full.names = TRUE)
  files <- c(file.path(results_dir, site, metrics_dir, "mean_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "lcv_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "lskew_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "hillshade_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "vci_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "slope_res_10_m.tif"),
             file.path(results_dir, site, metrics_dir, "stand_trees_raster.tif"),
             file.path(results_dir, site, metrics_dir, "lidarlai_res_10_m.tif"))
  
  mraster <- terra::rast(files)
  names(mraster[[1]]) <- 'mean'
  names(mraster[[2]]) <- 'lcv'
  names(mraster[[3]]) <- 'lskew'
  names(mraster[[4]]) <- 'hillshade'
  names(mraster[[5]]) <- 'vci'
  names(mraster[[6]]) <- 'slope'
  names(mraster[[7]]) <- 'stands'
  names(mraster[[8]]) <- 'lidar_lai'
  
  br.mean <- c(15, 25)
  br.lskew <- c(-0.3, -0.1)
  br.mean <- switch(site,
                    "Aigoual" = c(12, 16),
                    "Blois" = c(18, 28),
                    "Mormal" = c(18, 27)
                    # "Aigoual" = c(11, 14, 17),
                    # "Blois" = c(12, 24, 31),
                    # "Mormal" = c(18, 25, 27),
                    # "Aigoual" = 14,
                    # "Blois" = 24,
                    # "Mormal" = 25
  )
  br.lcv <- switch(site,
                   # "Aigoual" = c(0.0015, 0.0025),
                   # "Blois" = c(0.02, 0.04),
                   # "Mormal" = c(0.015, 0.03)
                   "Aigoual" = c(0.0009, 0.0013, 0.0018),
                   "Blois" = c(0.02, 0.028, 0.034),
                   "Mormal" = c(0.014, 0.019, 0.024),
                   # "Aigoual" = 0.0013,
                   # "Blois" = 0.028,
                   # "Mormal" = 0.019
  )
  br.lskew <- switch(site,
                     "Aigoual" = c(-0.35, -0.2),
                     "Blois" = c(-0.35, -0.2),
                     "Mormal" = c(-0.3, -0.1)
                     # "Aigoual" = c(-0.35, -0.15, -0.1),
                     # "Blois" = c(-0.33, -0.3, -0.27),
                     # "Mormal" = c(-0.34, -0.28, -0.21)
                     # "Aigoual" = -0.17,
                     # "Blois" = -0.25,
                     # "Mormal" = -0.26
  )
  br.lidarlai <- switch(site,
                        # "Aigoual" = c(4, 7),
                        # "Blois" = c(3.5, 4.5),
                        # "Mormal" = c(4, 5.5)
                        "Aigoual" = 5.5,
                        "Blois" = 4,
                        "Mormal" = 5
  )
  sraster <- strat_breaks(mraster[[c(1,3)]],
                          breaks = list(br.mean, br.lskew),
                          plot = F, details = F)
  composite_strata <- sraster[[1]] * 10 + sraster[[2]]
  # composite_strata <- sraster
  # nStrata = 3
  totalStrata <- (length(br.mean) + 1) * (length(br.lskew) + 1) # nStrata * nStrata
  # totalStrata <- length(br.mean) + 1
  # totalStrata <- nStrata # * nStrata # * nStrata
  # totalStrata <- nStrata * nStrata # * nStrata
  # totalStrata <- 8
  # weights <- rep(1 / totalStrata, totalStrata)
  # sraster <- strat_quantiles(mraster[[c(1,2,3)]],
  #                            nStrata = list(nStrata, nStrata, nStrata),
  #                            plot = F)
  # sraster <- strat_quantiles(mraster[[8]], nStrata = 2,
  #                            plot = F, details = F)
  # composite_strata <- sraster[[1]]
  # composite_strata <- sraster[[1]] * 10 + sraster[[2]]
  # composite_strata <- sraster[[1]] * 100 + sraster[[2]] * 10 + sraster[[3]]
  # composite_strata <- sraster[[1]] * 1000 + sraster[[2]] * 100 + sraster[[3]] * 10 + sraster[[4]]
  # totalStrata <- 2
  weights <- rep(1 / totalStrata, totalStrata)
  # print(global(composite_strata == 211, "sum", na.rm = TRUE))
  # print(global(sraster == 24, "sum", na.rm = TRUE))
  
  # moran_i_values <- replicate(100, {
  #   set.seed(sample(1:1000, 1))
  #   sampled <- sample_strat(sraster = composite_strata,
  #                           nSamp = nSamp,
  #                           allocation = "manual",
  #                           weights = weights,
  #                           mindist = mindist,
  #                           force = TRUE, plot = FALSE)
  # 
  #   # Find closest matches
  #   site_data <- all_observations[[site]]
  #   site_data_sf <- st_as_sf(site_data, coords = c("lon", "lat"), crs = st_crs(sampled))
  #   sampled_coords <- st_coordinates(sampled)
  #   site_coords <- st_coordinates(site_data_sf)
  #   closest_matches <- vector("list", nrow(sampled_coords))
  #   for (i in 1:nrow(sampled_coords)) {
  #     target_coord <- sampled_coords[i, ]
  #     distances <- sqrt((site_coords[, 1] - target_coord[1])^2 +
  #                         (site_coords[, 2] - target_coord[2])^2)
  #     closest_index <- which.min(distances)
  #     closest_lon_lat <- as.numeric(site_coords[closest_index, ])
  #     closest_matches[[i]] <- cbind(
  #       sampled[i, ] %>% st_drop_geometry(),
  #       site_data_sf[closest_index, ] %>% st_drop_geometry(),
  #       lon = closest_lon_lat[1],
  #       lat = closest_lon_lat[2],
  #       closest_distance = distances[closest_index],
  #       site = site
  #     )
  #   }
  #   closest_data <- do.call(rbind, closest_matches)
  #   coordinates(closest_data) <- ~ lon + lat
  #   # Define k-nearest neighbors
  #   k <- 1  # Adjust k as needed
  #   neighbors <- knearneigh(coordinates(closest_data), k = k)
  #   lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
  # 
  #   # Extract Moran's I variable and calculate Moran's I
  #   moran_variable <- closest_data$lidar_lai
  #   moran_test <- moran.test(moran_variable, listw = lw_weights)
  #   print(moran_test$estimate[1])
  #   return(c(moran_test$estimate[1], moran_test$p.value))
  # })
  
  # Calculate mean and standard deviation of Moran's I and p-values
  # moran_i_values_mean <- mean(moran_i_values[1, ], na.rm = TRUE)
  # moran_i_values_sd <- sd(moran_i_values[1, ], na.rm = TRUE)
  # p_values_mean <- mean(moran_i_values[2, ], na.rm = TRUE)
  # p_values_sd <- sd(moran_i_values[2, ], na.rm = TRUE)
  # print(paste("Site:", site, "Mean Moran's I:", moran_i_values_mean, "SD:", moran_i_values_sd))
  # print(paste("Site:", site, "Mean p-value:", p_values_mean, "SD:", p_values_sd))
  
  sampled <- sample_strat(sraster = composite_strata, 
                          nSamp = nSamp, 
                          allocation = "manual",
                          weights = weights,
                          mindist = mindist, 
                          force = T,
                          plot = T)
  # sampled <- sample_systematic(mraster,
  # cellsize = 400,
  # location = "random",
  # force = TRUE)
  # sampled_ind <- sample(1:nrow(all_observations[[site]]), size = nSamp)
  # sampled <- all_observations[[site]][sampled_ind, ]
  # all_sites_closest_data[[site]] <- sampled
  
  # Find closest matches
  site_data <- all_observations[[site]]
  site_data <- site_data[site_data$stands %in% c("oak", "beech", "deciduous", "poplar"), ]
  site_data_sf <- st_as_sf(site_data, coords = c("lon", "lat"), 
                           crs = st_crs(sampled))
  sampled_coords <- st_coordinates(sampled)
  site_coords <- st_coordinates(site_data_sf)
  closest_matches <- vector("list", nrow(sampled_coords))
  for (i in 1:nrow(sampled_coords)) {
    target_coord <- sampled_coords[i, ]
    distances <- sqrt((site_coords[, 1] - target_coord[1])^2 +
                        (site_coords[, 2] - target_coord[2])^2)
    closest_index <- which.min(distances)
    # print(closest_index)
    closest_lon_lat <- as.numeric(site_coords[closest_index, ])
    closest_matches[[i]] <- cbind(
      sampled[i, ] %>% st_drop_geometry(),
      site_data_sf[closest_index, ] %>% st_drop_geometry(),
      lon = closest_lon_lat[1],
      lat = closest_lon_lat[2],
      closest_distance = distances[closest_index],
      site = site
    )
  }
  closest_data <- do.call(rbind, closest_matches)
  all_sites_closest_data[[site]] <- closest_data
  
  all_sites_closest_data[[site]]$strata = NULL
  all_sites_closest_data[[site]]$type = NULL
  all_sites_closest_data[[site]]$rule = NULL
  all_sites_closest_data[[site]]$closest_distance = NULL
  all_sites_closest_data[[site]]$site = NULL
  matching_rows <- all_observations[[site]][all_observations[[site]]$lat %in% all_sites_closest_data[[site]]$lat &
                                              all_observations[[site]]$lon %in% all_sites_closest_data[[site]]$lon &
                                              all_observations[[site]]$lidar_lai %in% all_sites_closest_data[[site]]$lidar_lai, ]
  all_sites_closest_data[[site]] <- rbind(all_sites_closest_data[[site]], matching_rows)
  all_sites_closest_data[[site]] <- all_sites_closest_data[[site]][!duplicated(all_sites_closest_data[[site]]), ]
}
final_closest_data <- do.call(rbind, all_sites_closest_data)

# saveRDS(all_sites_closest_data$Aigoual, "Aigoual_sgsR_0312.rds")
# saveRDS(all_sites_closest_data$Blois, "Blois_sgsR_0312.rds")
# saveRDS(all_sites_closest_data$Mormal, "Mormal_sgsR_0312.rds")
# all_sites_closest_data$Mormal <- readRDS("Mormal_sgsR_0312.rds")

# Moran
moran_data <- final_closest_data # all_predictors_filtered all_sites_closest_data
moran_data <- all_sites_closest_data$Blois # Aigoual Blois Mormal
coordinates(moran_data) <- ~ lon + lat
k <- 8
neighbors <- knearneigh(coordinates(moran_data), k = k)
lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
# lw_weights$weights[2:60] <- 0

# neighbors <- dnearneigh(coordinates(moran_data), 0, 1000)
# lw_weights <- nb2listw(neighbors, style = "W", zero.policy = T)

moran_variable <- moran_data$lidar_lai
moran_test <- moran.test(moran_variable, listw = lw_weights
                         # alternative = "two.sided"
)
print(moran_test)
cat("\n")

# Correlog

# Define the coordinates and variable of interest
lon <- moran_data$lon
lat <- moran_data$lat
lidar_lai <- moran_data$lidar_lai

# Run the correlogram
correlog_result <- correlog(
  x = lon,
  y = lat,
  z = lidar_lai,
  increment = 500,  # Distance lag increments (adjust as needed)
  resamp = 1000     # Number of resampling iterations for significance testing
)

# Print results
print(correlog_result)
plot(correlog_result)

# ---------------------------------- Test --------------------------------------

center_pixel_index <- 1
k = 8
neighbors <- knearneigh(coordinates(moran_data), k = k)
neighbor_indices <- neighbors$nn[center_pixel_index,]
weights <- rep(0, length(moran_data$lidar_lai))
weights[neighbor_indices] <- 1
nb_custom <- list(neighbor_indices)
class(nb_custom) <- "nb"
lw_custom <- nb2listw(nb_custom, style = "W")
moran_variable <- moran_data$lidar_lai
moran_variable <- moran_variable[neighbor_indices]
moran_test <- moran_t(moran_variable, listw = lw_custom)




# Coordinates for the pixels (simple 5x5 grid)
coordinates <- expand.grid(x = 1:5, y = 1:5)
coordinates <- as.matrix(coordinates)

# Create a simple Moran variable (just for testing)
moran_variable <- rnorm(nrow(coordinates))  # Random values for the variable

# Define the center pixel index and its neighbors
center_pixel_index <- 13  # Center of the grid (example: 13th pixel in a 5x5 grid)
pixel_neighbors <- c(8, 9, 10, 12, 14, 16, 17, 18)  # Example neighbors for the center pixel

# Create a neighbor list using knearneigh (for simplicity)
k <- 8
neighbors <- knearneigh(coordinates, k = k)
lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
lw_weights$neighbours <- lw_weights$neighbours[center_pixel_index]
lw_weights$weights <- lw_weights$weights[center_pixel_index]
new_region_ids <- as.character(center_pixel_index)
attr(lw_weights, "region.id") <- new_region_ids

moran_test <- localmoran(moran_variable, lw_weights)

moran_test <- moran_t_custom(moran_variable, listw = lw_weights, pixel_neighbors = pixel_neighbors)




moran_t_custom <- function(x, listw, randomisation = TRUE, 
                           zero.policy = attr(listw, "zero.policy"), 
                           alternative = "greater", rank = FALSE, na.action = na.fail, 
                           spChk = NULL, adjust.n = TRUE, drop.EI2 = FALSE, 
                           pixel_neighbors = NULL) {
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  wname <- deparse(substitute(listw))
  if (!inherits(listw, "listw")) 
    stop(wname, "is not a listw object")
  
  xname <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop(xname, " is not a numeric vector")
  
  if (is.null(zero.policy)) 
    zero.policy <- get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  
  n <- length(listw$neighbours)
  
  # Ensure listw has only the selected pixel and its neighbors
  if (!is.null(pixel_neighbors)) {
    # pixel_neighbors should be a vector of indices for the selected pixel and its neighbors
    listw$neighbours <- listw$neighbours[pixel_neighbors]
  }
  
  print(n)
  print(x)
  print(length(x))
  # Ensure the number of regions matches the length of the data vector
  if (n != length(x)) 
    stop("Length of data and number of regions in weights do not match")
  
  wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
  S02 <- wc$S0 * wc$S0
  res <- moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy, NAOK = NAOK)
  
  I <- res$I
  K <- res$K
  if (rank) 
    K <- (3 * (3 * wc$n^2 - 7))/(5 * (wc$n^2 - 1))
  
  EI <- (-1)/wc$n1
  if (randomisation) {
    VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n * wc$S2 + 3 * S02)
    tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 + 6 * S02)
    if (tmp > VI) 
      warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
    VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
    if (!drop.EI2) 
      VI <- (VI - EI^2)
    if (VI < 0) 
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  else {
    VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 * (wc$nn - 1))
    if (!drop.EI2) 
      VI <- (VI - EI^2)
    if (VI < 0) 
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Moran I statistic standard deviate"
  
  if (alternative == "two.sided") 
    PrI <- 2 * pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater") 
    PrI <- pnorm(ZI, lower.tail = FALSE)
  else PrI <- pnorm(ZI)
  
  if (!is.finite(PrI) || PrI < 0 || PrI > 1) 
    warning("Out-of-range p-value: reconsider test arguments")
  
  vec <- c(I, EI, VI)
  names(vec) <- c("Moran I statistic", "Expectation", "Variance")
  
  method <- paste("Moran I test under", ifelse(randomisation, 
                                               "randomisation", "normality"))
  data.name <- paste(xname, ifelse(rank, "using rank correction", ""), "\nweights:", wname, 
                     ifelse(is.null(na.act), "", paste("\nomitted:", paste(na.act, collapse = ", "))), 
                     ifelse(adjust.n && isTRUE(any(sum(card(listw$neighbours) == 0L))), "\n reduced by no-neighbour observations", ""), 
                     ifelse(drop.EI2, "\nEI^2 term dropped in VI", ""), "\n")
  
  res <- list(statistic = statistic, p.value = PrI, estimate = vec, 
              alternative = alternative, method = method, data.name = data.name)
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  class(res) <- "htest"
  return(res)
}







# ---------------------------------- SYS ---------------------------------------
all_balanced_sys_samples <- list()
for (site in sites) {
  # site <- "Blois"
  print(paste("Site:", site))
  site_data <- load_metrics_with_pad(site, results_dir, "Not_Masked")
  site_data <- site_data[
    # site_data$stands %in% c("oak", "beech", "deciduous", "poplar") &
    !is.na(site_data$lidar_lai) &
      site_data$mean >= 2,
  ]
  
  # Sort by spatial or structural variables for systematic sampling
  site_data <- site_data %>% arrange(lon, lat)
  
  # Set desired sample size
  sample_size <- 120 # 500 1000 2000 5000 10000
  
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
  
  k <- 8
  neighbors <- knearneigh(coords_sampled, k = k)
  weights <- nb2listw(knn2nb(neighbors), style = "W")
  # morans_I_sampled <- localmoran(sampled_points$lidar_lai - sampled_points$s2_lai, lw_sampled)
  morans_I_sampled <- moran.test(balanced_sys_sample$lidar_lai, weights)
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

# ---------------------------------- Plot --------------------------------------
vars <- list(c("lidar_lai", "s2_lai"), 
             "lcv", "lskew",
             "vci",
             "mean"
             )
# all_observations[[site]]
for (site in sites){
  for (var in vars)
    plot_raincloud_one_site(site,
                            all_balanced_sys_samples[[site]],
                            all_predictors_filtered[[site]],
                            all_sites_closest_data[[site]],
                            output_dir,
                            variables = var
                            )
}
par(mfrow = c(1, 1))

for (site in sites){
  # all_predictors_filtered[[site]] all_sites_closest_data[[site]] all_balanced_sys_samples[[site]]
  spatial_autocorrelation_tests(moran_data = all_sites_closest_data[[site]],
                                site = site,
                                variable_name = "lidar_lai")
}

for (site in sites) {
  sampled_lidar_lai <- all_sites_closest_data[[site]]$deltaLAI_dsm
  field_lidar_lai <- all_predictors_filtered[[site]]$deltaLAI_dsm
  combined_lai <- c(sampled_lidar_lai, field_lidar_lai)
  x_range <- range(combined_lai, na.rm = TRUE)
  nb_bins <- 10 
  breaks <- seq(x_range[1], x_range[2], length.out = nb_bins + 1)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 1, 0), oma = c(0, 0, 2, 0))
  hist_sampled <- hist(sampled_lidar_lai, breaks = breaks, plot = FALSE)
  hist_sampled$counts <- hist_sampled$counts / sum(hist_sampled$counts)
  plot(hist_sampled, 
       xlab = "Sampled", ylab = "Normalized Frequency", 
       main = "", col = "lightblue", xlim = x_range)
  hist_field <- hist(field_lidar_lai, breaks = breaks, plot = FALSE)
  hist_field$counts <- hist_field$counts / sum(hist_field$counts)
  plot(hist_field, 
       xlab = "Field", ylab = "Normalized Frequency", 
       main = "", col = "lightgreen", xlim = x_range)
  mtext(paste(site, "Comparison of Normalized Histograms"),
        outer = TRUE, cex = 1.5, line = -0.5)
}

# ------------------------------------------------------------------------------
par(mfrow = c(1,1))
# site <- "Blois" # Aigoual Blois Mormal All Sites Combined
sites <- c("Aigoual", "Blois", "Mormal") # Aigoual Blois Mormal
# sites <- c("Blois", "Mormal")
site <- "Mormal" # Aigoual Blois Mormal
# for (site in sites){
# -------------------------------- RF DTM --------------------------------------
train_loadings <- process_sites_dsm_dtm(all_sites_closest_data, site, 
                                        results_dir, metrics_dir, 
                                        plot = TRUE, num_iterations = 100)
# train_data_dtm <- train_loadings$dtm[[site]]
# train_data_dsm <- train_loadings$dsm[[site]]
train_data_dtm <- do.call(rbind, train_loadings$dtm)
train_data_dsm <- do.call(rbind, train_loadings$dsm)

test_loadings <- process_sites_dsm_dtm(all_predictors_filtered, site, 
                                       results_dir, metrics_dir, 
                                       plot = TRUE, num_iterations = 100)
# test_data_dtm <- test_loadings$dtm[[site]]
# test_data_dsm <- test_loadings$dsm[[site]]
test_data_dtm <- do.call(rbind, test_loadings$dtm)
test_data_dsm <- do.call(rbind, test_loadings$dsm)

# all_sites_closest_data all_balanced_sys_samples
# train_data_dtm <- process_sites_dtm(
#   all_sites_data = all_sites_closest_data, # all_sites_closest_data all_balanced_sys_samples
#   sites = site,
#   results_dir = results_dir,
#   metrics_dir = metrics_dir,
#   plot = T,
#   num_iterations = 100
# )
# 
# test_data_dtm <- process_sites_dtm(
#   all_sites_data = all_predictors_filtered,
#   sites = site,
#   results_dir = results_dir,
#   metrics_dir = metrics_dir,
#   plot = T,
#   num_iterations = 100
# )

# pad_files <- sprintf(file.path(results_dir, site, metrics_dir,
#                                "PAD_Profiles_NA/PAD_%.1f_40.tif"), 
#                      seq(2.5, 39.5, by = 1))
# train_data_dtm <- all_sites_closest_data[[site]]
# train_data_dtm <- final_closest_data # all_balanced_sys_samples all_sites_closest_data
# train_data_dtm <- sample_sites(all_balanced_sys_samples, sites)

# remove_vars_dsm <- c("random_pad_dsm", "depth_dsm", 
#                      "deltaLAI_dsm", "cv_lad_dsm", "closest_distance", 
#                      "site", "strata", "type", "rule", "deltaLAI",
#                      # "lon", "lat", 
#                      "slope", "aspect",
#                      "stands",
#                      "cv_lad_dsm", "ndvi"
# )
# for (remove_var in remove_vars_dsm){
#   train_data_dtm[[remove_var]] <- NULL
# }
# 
# train_data_dtm <- process_depth_and_pad(train_data_dtm, pad_files, plot = TRUE, num_iterations = 10)
# colnames(train_data_dtm)[colnames(train_data_dtm) == "depth"] <- "depth_dtm"
# colnames(train_data_dtm)[colnames(train_data_dtm) == "random_pad"] <- "random_pad_dtm"
# colnames(train_data_dtm)[colnames(train_data_dtm) == "delta"] <- "deltaLAI_dtm"

# train_data_dtm <- na.omit(train_data_dtm)
# train_data_dtm <- train_data_dtm %>%
#   group_by(lat, lon) %>%
#   slice_sample(n = 10) %>%
#   ungroup()
# 
# k <- 1
# coordinates(train_data_dtm) <- ~ lon + lat
# neighbors <- knearneigh(coordinates(train_data_dtm), k = k)
# lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
# moran_variable <- train_data_dtm$deltaLAI_dtm
# moran_test <- moran.test(moran_variable, listw = lw_weights)
# print(moran_test)

# train_control <- trainControl(method = "cv", number = 10)
# rf_model_cv <- train(as.formula(paste("deltaLAI_dtm ~ . - lidar_lai - random_pad_dtm")), 
#                      data = train_data_dtm, 
#                      method = "rf", trControl = train_control)
# print(rf_model_cv)

# test_data_dtm <- all_predictors_filtered[[site]]
# test_data_dtm <- combined_predictors_filtered
# 
# test_data_dtm <- process_depth_and_pad(test_data_dtm, pad_files, plot = TRUE, num_iterations = 10)
# colnames(test_data_dtm)[colnames(test_data_dtm) == "depth"] <- "depth_dtm"
# colnames(test_data_dtm)[colnames(test_data_dtm) == "random_pad"] <- "random_pad_dtm"
# colnames(test_data_dtm)[colnames(test_data_dtm) == "delta"] <- "deltaLAI_dtm"
# 
# subset_data <- test_data_dtm[, c("lidar_lai", "s2_lai", "deltaLAI", "mean", "depth_dtm", "deltaLAI_dtm")]

# RF
formula <- as.formula(paste("deltaLAI_dtm ~ . - lidar_lai - random_pad_dtm"))
formula <- deltaLAI_dtm ~ depth_dtm + lcv + vci + lskew + fCover + rumple # + cv_lad_dtm
# formula <- deltaLAI_dtm ~ depth_dtm + lcv
rf_model_dtm <- ranger(
  formula = formula,
  data = train_data_dtm,
  num.trees = 500,
  # mtry = 13,
  mtry = 5,
  # mtry = 2,
  # case.weights =  1 / table(train_data_dtm$deltaLAI_dtm)[train_data_dtm$deltaLAI_dtm],
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
)

# varImp
varImp_dtm <- sort(rf_model_dtm$variable.importance, decreasing = TRUE)
# png("/home/corroyez/Pictures/paper/Mormal/dtmImp.png", width = 1920, height = 1080)
barplot(varImp_dtm)
# dev.off()
oob_mse <- round(rf_model_dtm$prediction.error, 2)
range_LAI <- range(train_data_dtm$deltaLAI_dtm)
oob_nrmse <- round(sqrt(oob_mse) / IQR(train_data_dtm$deltaLAI_dtm, na.rm = TRUE), 2)
oob_r_squared <- round(rf_model_dtm$r.squared, 2)
print(oob_nrmse)
print(oob_r_squared)

# Pred
# dtm_train_pred <- predict(rf_model_dtm, newdata = train_data_dtm)$predictions
dtm_test_pred <- predict(rf_model_dtm, data = test_data_dtm)$predictions
test_data_dtm$s2_plus_deltaLAI <- dtm_test_pred + test_data_dtm$s2_lai
subset_data_dtm <- test_data_dtm[, c("random_pad_dtm", "s2_lai", "s2_plus_deltaLAI")]

# Test
residuals_test <- test_data_dtm$deltaLAI_dtm - dtm_test_pred
sst_test <- sum((test_data_dtm$deltaLAI_dtm - mean(test_data_dtm$deltaLAI_dtm))^2)
sse_test <- sum(residuals_test^2)
r2_test <- round(1 - (sse_test / sst_test), 2)
mse_test <- mean(residuals_test^2)
range_test_LAI <- range(test_data_dtm$deltaLAI_dtm)
nrmse_test <- round(sqrt(mse_test) / IQR(test_data_dtm$deltaLAI_dtm, na.rm = TRUE), 2)
print(paste("Test R²:", r2_test))
print(paste("Test NRMSE:", nrmse_test))

# Viz
dtm_result <- create_corrected_s2lai_scatterplot(
  study_site = site,
  predictions = dtm_test_pred,
  test_data = test_data_dtm,
  lai_col = "random_pad_dtm", # deltaLAI_dtm random_pad_dtm
  s2_lai_col = "s2_lai",
  depth = "depth_dtm"
)

# Metrics
# dtm_metrics <- compute_rf_metrics(
#   model = rf_model_dtm,
#   test_predictions = dtm_test_pred,
#   train_data = train_data_dtm,
#   test_data = test_data_dtm,
#   target_col = "deltaLAI_dtm"
# )

# Correlation Matrix
# cor_matrix_dtm <- abs(cor(train_data_dtm))
# cor_hclust_dtm <- hclust(as.dist(1 - cor_matrix_dtm))
# cor_matrix_dtm <- cor_matrix_dtm[cor_hclust_dtm$order, cor_hclust_dtm$order]
# corrplot(cor_matrix_dtm,
#          method = "number",
#          type = "full",
#          order = "hclust",
#          hclust.method = "ward.D",
#          addrect = 3) # "complete" "ward.D"

# PP
# for (var in names(rf_model_dtm$variable.importance)) {
for (var in 1) {
  partialPlot(obj = rf_model_dtm,
              pred.data = train_data_dtm,
              xname = "depth_dtm", # var,
              target_var = "deltaLAI_dtm",
              show_plot = TRUE,
              save_path = NULL
              # save_path = file.path(output_dir, site, "final/dtm")
              # save_path = "/home/corroyez/Pictures/paper/Mormal/"
  ) # NULL pp_dir
}

# Bi-PDP
list <- list(c("depth_dtm", "lcv")
)
for (pair in list) {
  xname1 <- pair[1]
  xname2 <- pair[2]
  biPartialPlot(
    obj = rf_model_dtm,
    pred.data = train_data_dtm,
    xname1 = xname1,
    xname2 = xname2,
    xlab1 = xname1,
    xlab2 = xname2,
    target_var = "deltaLAI_dtm",
    show_plot = TRUE,
    save_path = NULL
  )
}

# Corrplot
depth_values <- unique(train_data_dtm$depth_dtm)
correlations <- numeric(length(depth_values))

for (i in seq_along(depth_values)) {
  subset_data <- train_data_dtm[train_data_dtm$depth_dtm <= depth_values[i], ]
  correlations[i] <- cor(subset_data$s2_lai, subset_data$random_pad_dtm, use = "complete.obs")
}
correlation_data <- data.frame(depth_dtm = depth_values, correlation = correlations)
max_corr_depth <- correlation_data$depth_dtm[which.max(correlation_data$correlation)]
ggplot(correlation_data, aes(x = depth_dtm, y = correlation)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  # geom_vline(xintercept = max_corr_depth, color = "darkgreen", linetype = "dashed", size = 1) +
  labs(
    title = paste(site, "Evolution of Correlation Between s2_lai and pad_dtm for Train Set"),
    x = "Depth DTM",
    y = "Correlation"
  ) +
  theme_minimal()


# VSURF
deltaLAI_dtm <- train_data_dtm$deltaLAI_dtm
train_data_dtm_vsurf <- train_data_dtm[, c(-1, -14)]
train_data_dtm_vsurf <- train_data_dtm_vsurf[, 
                                             !(colnames(train_data_dtm_vsurf) %in% "deltaLAI_dtm")]
vsurf_formula <- deltaLAI_dtm ~ .
vsurf <- VSURF(vsurf_formula, train_data_dtm_vsurf, mtry = 13, ntree.interp = 500)
vsurf_vars <- colnames(train_data_dtm_vsurf)[vsurf$varselect.interp]
print(vsurf_vars)
# -------------------------------- RF DSM --------------------------------------
# sites <- c("Aigoual", "Blois", "Mormal") # Aigoual Blois Mormal
# sites <- "Blois" # Aigoual Blois Mormal
# # all_sites_closest_data all_balanced_sys_samples
# 
# train_data_dsm <- process_sites_dsm(
#   all_sites_data = all_balanced_sys_samples,
#   sites = sites,
#   results_dir = results_dir,
#   metrics_dir = metrics_dir,
#   plot = T,
#   num_iterations = 1000
# )
# 
# test_data_dsm <- process_sites_dsm(
#   all_sites_data = all_predictors_filtered,
#   sites = sites,
#   results_dir = results_dir,
#   metrics_dir = metrics_dir,
#   plot = T,
#   num_iterations = 1000
# )


# pad_files <- sprintf(file.path(results_dir, site, metrics_dir,
#                                "PAD_Profiles_updated_modifminz_NA/PAD_%.1f_40.tif"), 
#                      seq(2.5, 39.5, by = 1))
# train_data_dsm <- all_sites_closest_data[[site]]
# train_data_dsm <- final_closest_data # all_balanced_sys_samples all_sites_closest_data

# remove_vars_dtm <- c("random_pad_dtm", "depth_dtm", 
#                      "deltaLAI_dtm", "cv_lad_dtm", "closest_distance", 
#                      "site", "strata", "type", "rule", "deltaLAI",
#                      # "lon", "lat", 
#                      "slope", "aspect",
#                      "stands",
#                      "cv_lad_dtm", "ndvi"
# )
# for (remove_var in remove_vars_dtm){
#   train_data_dsm[[remove_var]] <- NULL
# }
# 
# train_data_dsm <- process_depth_and_pad(train_data_dsm, pad_files, plot = TRUE, num_iterations = 100)
# colnames(train_data_dsm)[colnames(train_data_dsm) == "depth"] <- "depth_dsm"
# colnames(train_data_dsm)[colnames(train_data_dsm) == "random_pad"] <- "random_pad_dsm"
# colnames(train_data_dsm)[colnames(train_data_dsm) == "delta"] <- "deltaLAI_dsm"


# train_data_dsm <- na.omit(train_data_dsm)
# train_data_dsm <- train_data_dsm %>%
#   group_by(lat, lon) %>%
#   slice_sample(n = 2) %>%
#   ungroup()

# test_data_dsm <- all_predictors_filtered[[site]]
# test_data_dsm <- process_depth_and_pad(test_data_dsm, pad_files, plot = TRUE, num_iterations = 100)
# colnames(test_data_dsm)[colnames(test_data_dsm) == "depth"] <- "depth_dsm"
# colnames(test_data_dsm)[colnames(test_data_dsm) == "random_pad"] <- "random_pad_dsm"
# colnames(test_data_dsm)[colnames(test_data_dsm) == "delta"] <- "deltaLAI_dsm"
# 
# subset_data <- test_data_dsm[, c("lidar_lai", "s2_lai", "deltaLAI", "mean", "depth_dsm", "deltaLAI_dsm")]
# test_data_dsm <- combined_predictors_filtered
# test_data_dsm <- all_observations$Aigoual

# train_control <- trainControl(method = "cv", number = 10)
# rf_model_cv <- train(as.formula(paste("deltaLAI_dsm ~ . - lidar_lai - random_pad_dsm")), 
#                      data = train_data_dsm, 
#                      method = "rf", trControl = train_control)
# print(rf_model_cv)

# RF
formula <- as.formula(paste("deltaLAI_dsm ~ . - lidar_lai - random_pad_dsm"))
formula <- deltaLAI_dsm ~ depth_dsm + lcv + lskew + vci + fCover + rumple
# formula <- deltaLAI_dsm ~ depth_dsm + lskew + vci
# formula <- deltaLAI_dsm ~ depth_dsm + lskew + lcv
rf_model_dsm <- ranger(
  formula = formula,
  data = train_data_dsm,
  num.trees = 500,
  # mtry = 13,
  mtry = 5,
  # mtry = 2,
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
)

# varImp
varImp_dsm <- sort(rf_model_dsm$variable.importance, decreasing = TRUE)
barplot(varImp_dsm)
# png("/home/corroyez/Pictures/paper/Mormal/dsmImp.png", width = 1920, height = 1080)
# barplot(varImp_dsm)
# dev.off()
oob_mse <- round(rf_model_dsm$prediction.error, 2)
range_LAI <- range(train_data_dsm$deltaLAI_dsm)
oob_nrmse <- oob_nrmse <- round(sqrt(oob_mse) / IQR(train_data_dsm$deltaLAI_dsm, na.rm = TRUE), 2)
oob_r_squared <- round(rf_model_dsm$r.squared, 2)
print(oob_nrmse)
print(oob_r_squared)

# Pred
# dsm_train_pred <- predict(rf_model_dsm, newdata = train_data_dsm)$predictions
dsm_test_pred <- predict(rf_model_dsm, data = test_data_dsm)$predictions
test_data_dsm$s2_plus_deltaLAI <- dsm_test_pred + test_data_dsm$s2_lai
subset_data_dsm <- test_data_dsm[, c("random_pad_dsm", "s2_lai", "s2_plus_deltaLAI")]

# Test
residuals_test <- test_data_dsm$deltaLAI_dsm - dsm_test_pred
sst_test <- sum((test_data_dsm$deltaLAI_dsm - mean(test_data_dsm$deltaLAI_dsm))^2)
sse_test <- sum(residuals_test^2)
r2_test <- round(1 - (sse_test / sst_test), 2)
mse_test <- mean(residuals_test^2)
range_test_LAI <- range(test_data_dsm$deltaLAI_dsm)
nrmse_test <- round(sqrt(mse_test) / IQR(test_data_dsm$deltaLAI_dsm, na.rm = TRUE), 2)
print(paste("Test R²:", r2_test))
print(paste("Test NRMSE:", nrmse_test))


# Viz
dsm_result <- create_corrected_s2lai_scatterplot(
  study_site = site,
  predictions = dsm_test_pred,
  test_data = test_data_dsm,
  lai_col = "random_pad_dsm", # deltaLAI_dsm random_pad_dsm
  s2_lai_col = "s2_lai",
  depth = "depth_dsm"
)

# Metrics
# dsm_metrics <- compute_rf_metrics(
#   model = rf_model_dsm,
#   test_predictions = dsm_test_pred,
#   train_data = train_data_dsm,
#   test_data = test_data_dsm,
#   target_col = "deltaLAI_dsm"
# )

# Correlation Matrix
# cor_matrix_dsm <- abs(cor(train_data_dsm))
# cor_hclust_dsm <- hclust(as.dist(1 - cor_matrix_dsm))
# cor_matrix_dsm <- cor_matrix_dsm[cor_hclust_dsm$order, cor_hclust_dsm$order]
# corrplot(cor_matrix_dsm,
#          method = "circle",
#          type = "upper",
#          order = "hclust",
#          hclust.method = "ward.D",
#          addrect = 3) # "complete" "ward.D"
# print(cor(train_data_dsm, use = "complete.obs"))

# PP
# for (var in names(rf_model_dsm$variable.importance)) {
for (var in 1) {
  partialPlot(obj = rf_model_dsm,
              pred.data = train_data_dsm,
              xname = "depth_dsm", # var,
              target_var = "deltaLAI_dsm",
              show_plot = TRUE,
              # save_path = file.path(output_dir, site, "final/dsm")
              save_path = NULL
              # save_path = "/home/corroyez/Pictures/paper/Mormal/"
  ) # NULL pp_dir
}

# Bi-PDP
list <- list(c("depth_dsm", "lcv")
)
for (pair in list) {
  xname1 <- pair[1]
  xname2 <- pair[2]
  biPartialPlot(
    obj = rf_model_dsm,
    pred.data = train_data_dsm,
    xname1 = xname1,
    xname2 = xname2,
    xlab1 = xname1,
    xlab2 = xname2,
    target_var = "deltaLAI_dsm",
    show_plot = TRUE,
    save_path = NULL
  )
}

# Corrplot
depth_values <- unique(train_data_dsm$depth_dsm)
correlations <- numeric(length(depth_values))

for (i in seq_along(depth_values)) {
  subset_data <- train_data_dsm[train_data_dsm$depth_dsm <= depth_values[i], ]
  correlations[i] <- cor(subset_data$s2_lai, subset_data$random_pad_dsm, use = "complete.obs")
}
correlation_data <- data.frame(depth_dsm = depth_values, correlation = correlations)
max_corr_depth <- correlation_data$depth_dsm[which.max(correlation_data$correlation)]
ggplot(correlation_data, aes(x = depth_dsm, y = correlation)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  # geom_vline(xintercept = max_corr_depth, color = "darkgreen", linetype = "dashed", size = 1) +
  labs(
    title = paste(site, "Evolution of Correlation Between s2_lai and pad_dsm for Train Set"),
    x = "Depth DSM",
    y = "Correlation"
  ) +
  theme_minimal()

# VSURF
deltaLAI_dsm <- train_data_dsm$deltaLAI_dsm
train_data_dsm_vsurf <- train_data_dsm[, c(-1, -14)]
train_data_dsm_vsurf <- train_data_dsm_vsurf[, 
                                             !(colnames(train_data_dsm_vsurf) %in% "deltaLAI_dsm")]
vsurf_formula <- deltaLAI_dsm ~ .
vsurf <- VSURF(vsurf_formula, train_data_dsm_vsurf, mtry = 13, ntree.interp = 500)
vsurf_vars <- colnames(train_data_dsm_vsurf)[vsurf$varselect.interp]
print(vsurf_vars)

# -------------------------------- RF PRED -------------------------------------
par(mfrow = c(1,1))
site <- "Aigoual" # Aigoual Blois Mormal
train_data_pred <- all_sites_closest_data[[site]] # all_sites_closest_data all_balanced_sys_samples
# train_data_pred <- all_sites_closest_data
# train_data_pred <- rbind(all_sites_closest_data$Aigoual,
#                          all_sites_closest_data$Blois,
#                          all_sites_closest_data$Mormal)

remove_vars_pred <- c("random_pad_dtm", "depth_dtm", 
                      "deltaLAI_dtm", "random_pad_dsm", "depth_dsm", 
                      "deltaLAI_dsm", "cv_lad_dtm", "cv_lad_dsm",
                      "closest_distance", "site", "strata", "type", "rule",
                      "stands", "lon", "lat", "ndvi", "hillshade",
                      "slope", "aspect")
for (remove_var in remove_vars_pred){
  train_data_pred[[remove_var]] <- NULL
}
test_data_pred <- all_predictors_filtered[[site]]
# test_data_pred <- combined_predictors_filtered
# test_data_pred <- all_observations$Aigoual

# train_size <- 0.8
# train_index <- sample(seq_len(nrow(train_data_pred)), size = floor(train_size * nrow(train_data_pred)))
# test_data_pred <- train_data_pred[-train_index, ]
# train_data_pred <- train_data_pred[train_index, ]

# train_control <- trainControl(method = "cv", number = 10)
# rf_model_cv <- train(deltaLAI ~ lcv + lskew + mean + vci ,
#                      data = train_data_pred,
#                      method = "rf", trControl = train_control)
# print(rf_model_cv)

# RF
# formula <- as.formula(paste("deltaLAI_atbd ~ . - lidar_lai"))
formula <- deltaLAI_atbd ~ lcv + mean + vci + lskew + rumple + fCover
# formula <- lidar_lai ~ s2_lai_common + lcv + mean + vci + lskew + rumple + fCover
# formula <- deltaLAI_atbd ~ mean + lcv + vci + lskew
rf_model_pred <- ranger(
  formula = formula,
  data = train_data_pred,
  num.trees = 500,
  # mtry = 13,
  # mtry = 6,
  # mtry = 4,
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
)

# varImp
varImp_pred <- sort(rf_model_pred$variable.importance, decreasing = TRUE)
barplot(varImp_pred)
# png("/home/corroyez/Pictures/paper/Mormal/dsmImp.png", width = 1920, height = 1080)
# barplot(varImp_pred)
# dev.off()
oob_mse <- round(rf_model_pred$prediction.error, 2)
# range_LAI <- range(train_data_pred$deltaLAI)
oob_nrmse <- oob_nrmse <- round(sqrt(oob_mse) / IQR(train_data_pred$deltaLAI_atbd, na.rm = TRUE), 2)
oob_r_squared <- round(rf_model_pred$r.squared, 2)
print(oob_nrmse)
print(oob_r_squared)

# Pred
# train_pred <- predict(rf_model_pred, newdata = train_data_pred)$predictions
test_pred <- predict(rf_model_pred, data = test_data_pred)$predictions
test_data_pred$s2_plus_deltaLAI <- test_pred + test_data_pred$s2_lai_atbd
subset_data_pred <- test_data_pred[, c("lidar_lai", "s2_lai_atbd", "s2_plus_deltaLAI")]

# Test
residuals_test <- test_data_pred$deltaLAI_atbd - test_pred
sst_test <- sum((test_data_pred$deltaLAI_atbd - mean(test_data_pred$deltaLAI_atbd))^2)
sse_test <- sum(residuals_test^2)
r2_test <- round(1 - (sse_test / sst_test), 2)
mse_test <- mean(residuals_test^2)
range_test_LAI <- range(test_data_pred$deltaLAI_atbd)
nrmse_test <- nrmse_test <- round(sqrt(mse_test) / IQR(test_data_pred$deltaLAI_atbd, na.rm = TRUE), 2)
print(paste("Test R²:", r2_test))
print(paste("Test NRMSE:", nrmse_test))

r2 <- rSquared(test_data_pred$deltaLAI_atbd, resid = residuals_test)

# Viz
rf_result <- create_corrected_s2lai_scatterplot(
  study_site = site,
  predictions = test_pred,
  test_data = test_data_pred,
  lai_col = "lidar_lai",
  s2_lai_col = "s2_lai_atbd",
  depth = "mean"
)
# res <- rf_result$test_data
# df_above_6 <- subset(res, S2_plus_deltaLAI > 6)
# df_below_6 <- subset(res, S2_plus_deltaLAI <= 6)
# 
# # Moran
# moran_data <- test_data_pred
# coordinates(moran_data) <- ~ lon + lat
# k <- 8
# neighbors <- knearneigh(coordinates(moran_data), k = k)
# lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
# moran_test <- moran.test(residuals_test, listw = lw_weights)
# print(moran_test)


# Metrics
# pred_metrics <- compute_rf_metrics(
#   model = rf_model_pred,
#   test_predictions = test_pred,
#   train_data = train_data_pred,
#   test_data = test_data_pred,
#   target_col = "deltaLAI"
# )

# # Correlation Matrix
# cor_matrix_pred <- abs(cor(train_data_pred))
# cor_hclust_pred <- hclust(as.dist(1 - cor_matrix_pred))
# cor_matrix_pred <- cor_matrix_pred[cor_hclust_pred$order, cor_hclust_pred$order]
# corrplot(cor_matrix_pred,
#          method = "circle",
#          type = "upper",
#          order = "hclust",
#          hclust.method = "ward.D",
#          addrect = 3) # "complete" "ward.D"
# 

# PP
# for (var in names(rf_model_pred$variable.importance)) {
for (var in 1) {
  partialPlot(obj = rf_model_pred,
              pred.data = train_data_pred,
              xname = "mean", # var,
              # xlab = var,
              target_var = "deltaLAI_atbd",
              show_plot = TRUE,
              # save_path = file.path(output_dir, site, "final/pred")
              save_path = NULL
              # save_path = "/home/corroyez/Pictures/paper/Mormal/"
  ) # NULL pp_dir
}

# Bi-PDP
list <- list(
             # c("mean", "lskew"),
             c("mean", "lcv")
             # c("lcv", "lskew"),
             # c("mean", "vci"),
             # c("lcv", "vci")
)
for (pair in list) {
  xname1 <- pair[1]
  xname2 <- pair[2]
  biPartialPlot(
    obj = rf_model_pred,
    pred.data = train_data_pred,
    xname1 = xname1,
    xname2 = xname2,
    xlab1 = xname1,
    xlab2 = xname2,
    target_var = "deltaLAI_atbd",
    show_plot = TRUE,
    save_path = NULL
  )
}

# best_mtry <- plot_rf_mtry_evolution(train_data_pred, 
#                                     X = 6,
#                                     # formula = as.formula(paste("deltaLAI ~ . - lidar_lai")), 
#                                     formula = deltaLAI ~ lcv + mean + vci + lskew + rumple + fCover,
#                                     ntree = 500, 
#                                     min_node_size = 5,
#                                     num_threads = 12)

# VSURF
# deltaLAI <- train_data_pred$deltaLAI
# train_data_pred_vsurf <- train_data_pred[, -1]
# train_data_pred_vsurf <- train_data_pred_vsurf[, !(colnames(train_data_pred_vsurf) %in% "deltaLAI")]
# vsurf_formula <- deltaLAI ~ lcv + mean + vci + lskew + rumple + fCover
# vars_in_formula <- all.vars(vsurf_formula)
# train_data_subset <- train_data_pred[, vars_in_formula]
# train_data_subset$deltaLAI <- NULL
# vsurf <- VSURF(vsurf_formula, train_data_subset, mtry = best_mtry, ntree.interp = 500)
# vsurf_vars <- colnames(train_data_subset)[vsurf$varselect.interp]
# print("Selected variables:")
# print(vsurf_vars)

# ---------------------------- RF PRED (optim) ---------------------------------
par(mfrow = c(1, 1))
site <- "Blois" # Aigoual Blois Mormal
train_data_pred <- all_sites_closest_data[[site]] # all_sites_closest_data all_balanced_sys_samples
# train_data_pred <- rbind(all_sites_closest_data$Aigoual,
#                          all_sites_closest_data$Blois,
#                          all_sites_closest_data$Mormal)

remove_vars_pred <- c("random_pad_dtm", "depth_dtm", 
                      "deltaLAI_dtm", "random_pad_dsm", "depth_dsm", 
                      "deltaLAI_dsm", "cv_lad_dtm", "cv_lad_dsm",
                      "closest_distance", "site", "strata", "type", "rule",
                      "stands", "lon", "lat", "ndvi", "hillshade",
                      "slope", "aspect")
for (remove_var in remove_vars_pred){
  train_data_pred[[remove_var]] <- NULL
}
test_data_pred <- all_predictors_filtered[[site]]
# test_data_pred <- combined_predictors_filtered
# test_data_pred <- all_observations$Aigoual

# train_size <- 0.8
# train_index <- sample(seq_len(nrow(train_data_pred)), size = floor(train_size * nrow(train_data_pred)))
# test_data_pred <- train_data_pred[-train_index, ]
# train_data_pred <- train_data_pred[train_index, ]

# train_control <- trainControl(method = "cv", number = 10)
# rf_model_cv <- train(deltaLAI ~ lcv + lskew + mean + vci ,
#                      data = train_data_pred,
#                      method = "rf", trControl = train_control)
# print(rf_model_cv)

# RF
# formula <- as.formula(paste("deltaLAI_common ~ . - lidar_lai"))
# formula <- deltaLAI_common ~ lcv + mean + vci + lskew + rumple + fCover
formula <- s2_lai_atbd ~ mean + lcv + vci + lskew + rumple + fCover # s2_lai_atbd s2_lai_common
# formula <- deltaLAI_atbd ~ mean + lcv + vci + lskew + rumple + fCover # deltaLAI_atbd
# formula <- lidar_lai_under ~ mean + lcv + vci + lskew # + rumple + fCover # s2_lai_atbd
rf_model_pred <- ranger(
  formula = formula,
  data = train_data_pred,
  num.trees = 500,
  # mtry = 13,
  # mtry = 5,
  # mtry = 7,
  importance = "permutation",
  min.node.size = 5,
  num.threads = 12
)

# varImp
varImp_pred <- sort(rf_model_pred$variable.importance, decreasing = TRUE)
barplot(varImp_pred)
# png("/home/corroyez/Pictures/paper/Mormal/dsmImp.png", width = 1920, height = 1080)
# barplot(varImp_pred)
# dev.off()
oob_mse <- round(rf_model_pred$prediction.error, 2)
range_LAI <- range(train_data_pred$deltaLAI)
oob_nrmse <- oob_nrmse <- round(sqrt(oob_mse) / IQR(train_data_pred$lidar_lai, na.rm = TRUE), 2)
oob_r_squared <- round(rf_model_pred$r.squared, 2)
print(oob_nrmse)
print(oob_r_squared)

# Pred
# train_pred <- predict(rf_model_pred, newdata = train_data_pred)$predictions
test_pred <- predict(rf_model_pred, data = test_data_pred)$predictions
test_data_pred$s2_plus_deltaLAI <- test_pred + test_data_pred$s2_lai_common
subset_data_pred <- test_data_pred[, c("lidar_lai", "s2_lai_common", "s2_plus_deltaLAI")]

# x <- test_data_pred$lidar_lai
# y <- test_pred + test_data_pred$s2_lai_atbd # s2_lai_atbd
x <- train_data_pred$lidar_lai
y <- train_data_pred$lidar_lai_under + train_data_pred$s2_lai_common # s2_lai_atbd s2_lai_common
# x <- train_data_pred$lidar_lai_optimD
# y <- train_data_pred$s2_lai_atbd # s2_lai_atbd s2_lai_common
# x <- test_data_pred$lidar_lai_under
# y <- test_pred
# x <- test_data_pred$s2_lai_atbd
# y <- test_pred

plot(x, y,
     xlab = "Optimized S2 LAI (target)",
     ylab = "Predicted LiDAR total LAI",
     xlim = c(0, 20), ylim = c(0, 20))
abline(0, 1, col="red", lty=2)
cor(x, y)
# rSquared(x, resid = test_data_pred$lidar_lai_under - (test_pred))
Metrics::rmse(x, y)
fit <- lm(y ~ x); coeffs <- coef(fit)
coeffs[2] # Slope
coeffs[1] # Intercept
summary(fit)$r.squared # R2

# Test
# residuals_test <- test_data_pred$deltaLAI_common - test_pred
# sst_test <- sum((test_data_pred$deltaLAI_common - mean(test_data_pred$deltaLAI_common))^2)
# sse_test <- sum(residuals_test^2)
# r2_test <- round(1 - (sse_test / sst_test), 2)
# mse_test <- mean(residuals_test^2)
# range_test_LAI <- range(test_data_pred$deltaLAI_common)
# nrmse_test <- nrmse_test <- round(sqrt(mse_test) / IQR(test_data_pred$deltaLAI_common, na.rm = TRUE), 2)
# print(paste("Test R²:", r2_test))
# print(paste("Test NRMSE:", nrmse_test))

residuals_test <- test_data_pred$lidar_lai - test_pred
sst_test <- sum((test_data_pred$lidar_lai - mean(test_data_pred$lidar_lai))^2)
sse_test <- sum(residuals_test^2)
r2_test <- round(1 - (sse_test / sst_test), 2)
mse_test <- mean(residuals_test^2)
range_test_LAI <- range(test_data_pred$lidar_lai)
nrmse_test <- nrmse_test <- round(sqrt(mse_test) / IQR(test_data_pred$lidar_lai, na.rm = TRUE), 2)
print(paste("Test R²:", r2_test))
print(paste("Test NRMSE:", nrmse_test))

r2 <- rSquared(test_data_pred$lidar_lai, resid = residuals_test)

# Viz
rf_result <- create_corrected_s2lai_scatterplot(
  study_site = site,
  predictions = test_pred,
  test_data = test_data_pred,
  lai_col = "lidar_lai",
  s2_lai_col = "s2_lai_common",
  depth = "mean"
)
# res <- rf_result$test_data
# df_above_6 <- subset(res, S2_plus_deltaLAI > 6)
# df_below_6 <- subset(res, S2_plus_deltaLAI <= 6)
# 
# # Moran
# moran_data <- test_data_pred
# coordinates(moran_data) <- ~ lon + lat
# k <- 8
# neighbors <- knearneigh(coordinates(moran_data), k = k)
# lw_weights <- nb2listw(knn2nb(neighbors), style = "W")
# moran_test <- moran.test(residuals_test, listw = lw_weights)
# print(moran_test)


# Metrics
# pred_metrics <- compute_rf_metrics(
#   model = rf_model_pred,
#   test_predictions = test_pred,
#   train_data = train_data_pred,
#   test_data = test_data_pred,
#   target_col = "deltaLAI"
# )

# # Correlation Matrix
# cor_matrix_pred <- abs(cor(train_data_pred))
# cor_hclust_pred <- hclust(as.dist(1 - cor_matrix_pred))
# cor_matrix_pred <- cor_matrix_pred[cor_hclust_pred$order, cor_hclust_pred$order]
# corrplot(cor_matrix_pred,
#          method = "circle",
#          type = "upper",
#          order = "hclust",
#          hclust.method = "ward.D",
#          addrect = 3) # "complete" "ward.D"
# 

# PP
# for (var in names(rf_model_pred$variable.importance)) {
for (var in 1) {
  partialPlot(obj = rf_model_pred,
              pred.data = train_data_pred,
              xname = "mean", # var,
              # xlab = var,
              target_var = "deltaLAI_common",
              show_plot = TRUE,
              # save_path = file.path(output_dir, site, "final/pred")
              save_path = NULL
              # save_path = "/home/corroyez/Pictures/paper/Mormal/"
  ) # NULL pp_dir
}

# Bi-PDP
list <- list(
  # c("mean", "lskew"),
  c("mean", "lcv")
  # c("lcv", "lskew"),
  # c("mean", "vci"),
  # c("lcv", "vci")
)
for (pair in list) {
  xname1 <- pair[1]
  xname2 <- pair[2]
  biPartialPlot(
    obj = rf_model_pred,
    pred.data = train_data_pred,
    xname1 = xname1,
    xname2 = xname2,
    xlab1 = xname1,
    xlab2 = xname2,
    target_var = "deltaLAI_common",
    show_plot = TRUE,
    save_path = NULL
  )
}

# --------------------------------- LOOCV --------------------------------------
site <- "Blois" # Aigoual Blois Mormal
test_data_loocv <- all_predictors_filtered[[site]]
observed <- numeric(nrow(test_data_loocv))
predicted <- numeric(nrow(test_data_loocv))

for (i in 1:nrow(test_data_loocv)) {
  train_loocv <- test_data_loocv[-i, ]
  test_loocv <- test_data_loocv[i, , drop = FALSE]
  
  # RF
  formula <- as.formula(paste("deltaLAI ~ . - lidar_lai"))
  formula <- deltaLAI ~ mean + lcv
  rf_model_loocv <- ranger(
    formula = formula,
    data = train_loocv,
    num.trees = 500,
    importance = "permutation",
    min.node.size = 5,
    num.threads = 12
  )
  
  # Predict for the excluded observation
  prediction <- predict(rf_model_loocv, test_loocv)$predictions
  observed[i] <- test_loocv$lidar_lai
  predicted[i] <- prediction + test_data_loocv[-i, ]$s2_lai
}
loocv_results <- data.frame(Observed = observed, Predicted = predicted)
correlation <- cor(loocv_results$Observed, loocv_results$Predicted)
rmse <- sqrt(mean((loocv_results$Observed - loocv_results$Predicted)^2))

cat("LOOCV Results:\n")
cat("Correlation:", correlation, "\n")
cat("RMSE:", rmse, "\n")

# Visualization: Observed vs Predicted
plot(loocv_results$Observed, loocv_results$Predicted,
     xlab = "Observed deltaLAI", ylab = "Predicted deltaLAI",
     main = "LOOCV Observed vs Predicted",
     pch = 16, col = "blue")
abline(0, 1, col = "red", lwd = 2)  # 1:1 line

# loocv_result <- create_corrected_s2lai_scatterplot(
#   study_site = study_site,
#   rf_model = rf_model_loocv,
#   test_data = test_data_pred,
#   lai_col = "lidar_lai",
#   s2_lai_col = "s2_lai"
# )

# ------------------------------------------------------------------------------
pad_files <- sprintf("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/Metrics/Deciduous_Only/PAD_Profiles_NA/PAD_%.1f_40.tif", seq(2.5, 39.5, by = 1))
non_na_counts_matrix <- matrix(NA, nrow = nrow(coords), ncol = length(pad_files))
coords <- data.frame(lon = test_data_dtm$lon, lat = test_data_dtm$lat)
for (i in seq_along(pad_files)) {
  raster <- rast(pad_files[i])
  extracted_values <- extract(raster, coords)
  non_na_counts_matrix[, i] <- extracted_values[[2]]
}

na_counts_per_column <- colSums(is.na(non_na_counts_matrix)) + 1
epsilon <- 1e-6
weights <- log(na_counts_per_column)
weights <- ifelse(weights == 0, epsilon, weights)
weights <- weights / sum(weights)
valid_indices_matrix <- ifelse(!is.na(non_na_counts_matrix), 1, 0) 
for (i in 1:nrow(valid_indices_matrix)) {
  valid_weights_matrix[i, ] <- valid_indices_matrix[i, ] * weights
}
selected_file_indices <- apply(valid_weights_matrix, 1, function(x) sample(seq_along(x), size = 1, prob = x))
depth_vector <- 40 - (2 + selected_file_indices) + 1

num_iterations <- 100
depth_vectors <- calculate_depth_vectors(valid_weights_matrix, num_iterations)
depth_matrix <- depth_vectors$depth_matrix
# depth_selection <- depth_matrix[, 1]
# selected_file_indices <- depth_vectors$selected_files_matrix[, 1]

# hist_bins <- seq(min(depth_matrix, na.rm = TRUE), max(depth_matrix, na.rm = TRUE), length.out = 11)
# hist_counts <- matrix(0, nrow = num_iterations, ncol = length(hist_bins) - 1)
# for (iter in 1:num_iterations) {
#   hist_counts[iter, ] <- hist(depth_matrix[, iter], breaks = hist_bins, plot = FALSE)$counts
# }
# mean_hist_counts <- colMeans(hist_counts)
# par(mfrow = c(1, 2))
# barplot(mean_hist_counts, names.arg = round(hist_bins[-1], 2), col = "lightblue",
#         main = "Mean Histogram of Depths", xlab = "Depth", ylab = "Mean Count")
# barplot(test_data_hist$counts, names.arg = round(hist_bins[-1], 2), col = "lightgreen",
#         main = "Histogram of test_data_dsm$mean", xlab = "Mean Value", ylab = "Count")
# par(mfrow = c(1, 1))
# cor(mean_hist_counts, test_data_hist$counts)

test_data_dtm$depth_dtm <- floor(rowMeans(depth_matrix, na.rm = TRUE))
selected_file_indices <- 39 - test_data_dtm$depth_dtm
test_data_dtm$random_pad_dtm <- NA

for (i in 1:nrow(test_data_dtm)) {
  selected_file_index <- selected_file_indices[i]
  lon <- test_data_dtm$lon[i]
  lat <- test_data_dtm$lat[i]
  test_data_dtm$random_pad_dtm[i] <- extract_pad_value(pad_files, lon, lat, selected_file_index)
}

non_na_counts_matrix <- matrix(NA, nrow = nrow(coords), ncol = length(pad_files))
coords <- data.frame(lon = test_data_dsm$lon, lat = test_data_dsm$lat)
for (i in seq_along(pad_files)) {
  raster <- rast(pad_files[i])
  extracted_values <- extract(raster, coords)
  non_na_counts_matrix[, i] <- extracted_values[[2]]
}
na_counts_per_column <- colSums(is.na(non_na_counts_matrix)) + 1
epsilon <- 1e-6
weights <- (1 / (na_counts_per_column + 1))^1
# weights <- (na_counts_per_column + epsilon) / sum(na_counts_per_column + epsilon)
weights <- log(na_counts_per_column)^1
weights <- ifelse(weights == 0, epsilon, weights)
weights <- 1 - weights
weights <- weights / sum(weights)
weights
valid_indices_matrix <- ifelse(!is.na(non_na_counts_matrix), 1, 0) 

for (i in 1:nrow(valid_indices_matrix)) {
  valid_weights_matrix[i, ] <- valid_indices_matrix[i, ] * weights
}

num_iterations <- 100
depth_matrix <- calculate_depth_vectors(num_iterations)
hist_bins <- seq(min(depth_matrix, na.rm = TRUE), max(depth_matrix, na.rm = TRUE), length.out = 11)
hist_counts <- matrix(0, nrow = num_iterations, ncol = length(hist_bins) - 1)
for (iter in 1:num_iterations) {
  hist_counts[iter, ] <- hist(depth_matrix[, iter], breaks = hist_bins, plot = FALSE)$counts
}
mean_hist_counts <- colMeans(hist_counts)
barplot(mean_hist_counts, names.arg = round(hist_bins[-1], 2), col = "lightblue",
        main = "Mean Histogram of Depths", xlab = "Depth", ylab = "Mean Count")