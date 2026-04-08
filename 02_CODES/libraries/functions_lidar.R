# ---
# title: "functions_lidar"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

library("lidR")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("sf")
library("dplyr")
library("rayshader")

#' Vertical Complexity Index (VCI)
#'
#' Calculates the Vertical Complexity Index (VCI) based on input data.
#'
#' @param z Numeric vector representing vertical-related data.
#' @return VCI value.
#' @export
VCI_local <- function(z){
  zmax <- max(z)
  vci <- VCI(z, zmax=zmax)
}

VCI_dsm <- function(z, zmax = 40, z_opt = 5){
  z_filtered <- z[z >= zmax - z_opt]
  zmax <- max(z_filtered)
  vci <- VCI(z_filtered, zmax = zmax)
  # return(vci)
}

#' Vertical Density Ratio (VDR)
#'
#' Calculates the Vertical Density Ratio (VDR) based on input data.
#'
#' @param z Numeric vector representing vertical-related data.
#' @return VDR value.
#' @export
VDR <- function(z){
  zmax <- max(z)
  vdr <- (zmax - median(z)) / zmax
}

#' Coefficient of Variation (CV)
#'
#' Calculates the Coefficient of Variation (CV) based on input data.
#'
#' @param z Numeric vector representing vertical-related data.
#' @return CV value.
#' @export
CV <- function(z){
  cv <- (sd(z) / mean(z)) * 100
}

#' Rumple Index Surface
#'
#' Calculates the Rumple Index for a given LiDAR point cloud.
#'
#' @param las LiDAR point cloud data.
#' @param res Resolution of the output raster in meters.
#' @return Rumple Index raster.
#' @export
rumple_index_surface = function(las, res)
{
  las <- filter_surfacepoints(las, 1)
  rumple <- pixel_metrics(las, ~rumple_index(X,Y,Z), res)
  return(rumple)
}

#' Plant Area Index (PAI)
#'
#' Calculates the Plant Area Index (PAI) based on input data.
#'
#' @param z Numeric vector representing vegetation-related data.
#' @param zmin Minimum threshold value.
#' @param k Coefficient for PAI calculation.
#' @return PAI value.
#' @export
myPAI <- function(z, zmin, k=0.5){
  Nout = sum(z<zmin)
  Nin = length(z)
  PAI = -log(Nout/Nin)/k
}

wrapper_gap_fraction <- function(z){
  max_z <- max(z)
  z <- z + 40 - max_z
  min_z <- min(z) + 2
  gf_i <- gap_fraction_profile(z, z0 = min_z)
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, gf_i$z)
  missing_data <- data.frame(z = missing_z, gf = rep(0, length(missing_z)))
  gf <- rbind(gf_i, missing_data)
  gf <- gf[order(gf$z), ]
  gf$gf[is.na(gf$gf)] <- 0
  row.names(gf) <- NULL
  GFs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    # Store GF profiles
    GFs[k] <- gf$gf[gf$z == i]
    
    # CV_LAD
    # CV_LADs[k] <- 100 * sd(lad_i$lad[lad_i$z >= i], na.rm = T) / mean(lad_i$lad[lad_i$z >= i], na.rm = T)
    
    # Inc
    k <- k + 1
  }
  
  overall_GF <- sum(GFs)
}

#' Plant Area Density (PAD) Layers
#'
#' Calculates the Plant Area Density (PAD) for different layers
#' based on input data.
#'
#' @param z Numeric vector representing vegetation-related data.
#' @param zmin Minimum threshold value.
#' @param zmax Maximum threshold value.
#' @return List of PAD values for each layer.
#' @export
myPAD <- function(z) {
  
  lad_i <- LAD(z)
  lad_i$z <- floor(lad_i$z) + 0.5
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  # print(lad)
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    # LADs[k] <- lad$lad[lad$z == i]
    # PAIs[k] <- sum(lad$lad[lad$z >= i])
    # CV_LAD
    # CV_LADs[k] <- 100 * sd(lad_i$lad[lad_i$z >= i], na.rm = T) / mean(lad_i$lad[lad_i$z >= i], na.rm = T)
    
    # Inc
    k <- k + 1
  }
  
  # Calculate overall CV_LAD (Coefficient of Variation for all LADs)
  overall_CV_LAD <- 100 * sd(lad$lad, na.rm = TRUE)/ mean(lad$lad, na.rm = TRUE)
  
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  # CV_LADs_list <- setNames(as.list(CV_LADs), paste0("CV_LAD_", new_z))
  
  # Return PAIs and CV_LAD
  # combined_list <- c(PAIs_list, LADs_list, CV_LADs_list)
  combined_list <- c(PAIs_list, LADs_list, cv_lad = overall_CV_LAD)
}

#' Plant Area Density (PAD) Layers with DTM-shifted normalization
#'
#' Calculates the Plant Area Density (PAD) for different layers
#' based on input data.
#'
#' @param z Numeric vector representing vegetation-related data.
#' @param zmin Minimum threshold value.
#' @param zmax Maximum threshold value.
#' @return List of PAD values for each layer.
#' @export
myPAD_dtm <- function(z) {
  
  max_z <- max(z)
  z <- z + 40 - max_z
  # min_z <- min(z) + 2
  # print(z)
  lad_i <- LAD(z)
  lad_i$z <- floor(lad_i$z) + 0.5
  # print(lad_i)
  # print(sum(lad_i$lad), na.rm = T)
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  # print(lad)
  # print(sum(as.numeric(lad$lad), na.rm = TRUE))
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    # LADs[k] <- lad$lad[lad$z == i]
    # PAIs[k] <- sum(lad$lad[lad$z >= i])
    # CV_LAD
    # CV_LADs[k] <- 100 * sd(lad_i$lad[lad_i$z >= i], na.rm = T) / mean(lad_i$lad[lad_i$z >= i], na.rm = T)
    
    # Inc
    k <- k + 1
  }
  
  # Calculate overall CV_LAD (Coefficient of Variation for all LADs)
  # overall_CV_LAD <- 100 * sd(lad$lad, na.rm = TRUE)/ mean(lad$lad, na.rm = TRUE)
  
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  # CV_LADs_list <- setNames(as.list(CV_LADs), paste0("CV_LAD_", new_z))
  
  # Return PAIs and CV_LAD
  combined_list <- c(PAIs_list, LADs_list)
  # combined_list <- c(PAIs_list, LADs_list, cv_lad = overall_CV_LAD)
  # combined_list <- c(LADs_list)
}






myPAD_dtm2 <- function(z) {
  
  max_z <- max(z)
  # z <- z + 40 - max_z
  min_z <- min(z) + 2
  lad_i <- LAD(z, z0 = 2)
  lad_i$z <- floor(lad_i$z) + 0.5
  print(z)
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    # LADs[k] <- lad$lad[lad$z == i]
    # PAIs[k] <- sum(lad$lad[lad$z >= i])
    # CV_LAD
    # CV_LADs[k] <- 100 * sd(lad_i$lad[lad_i$z >= i], na.rm = T) / mean(lad_i$lad[lad_i$z >= i], na.rm = T)
    
    # Inc
    k <- k + 1
  }
  
  # Calculate overall CV_LAD (Coefficient of Variation for all LADs)
  # overall_CV_LAD <- 100 * sd(lad$lad, na.rm = TRUE)/ mean(lad$lad, na.rm = TRUE)
  
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  # CV_LADs_list <- setNames(as.list(CV_LADs), paste0("CV_LAD_", new_z))
  
  # Return PAIs and CV_LAD
  # combined_list <- c(PAIs_list, LADs_list, CV_LADs_list)
  # combined_list <- c(PAIs_list, LADs_list, cv_lad = overall_CV_LAD)
  combined_list <- c(LADs_list)
}












#' Plant Area Density (PAD) Layers with DSM-shifted normalization
#'
#' Calculates the Plant Area Density (PAD) for different layers
#' based on input data.
#'
#' @param class Class of z value.
#' @param z Numeric vector representing vegetation-related data.
#' @return List of PAD values for each layer.
#' @export
myPAD_dsm <- function(z) {
  
  # Associate ground z values to under vegatation value
  # min_z_veg <- NA  # Initialize min_z
  # if (sum(class == 3) > 0) {
  #   min_z_veg <- min(z[class == 3])
  # } else if (sum(class == 4) > 0) {
  #   min_z_veg <- min(z[class == 4])
  # } else if (sum(class == 5) > 0) {
  #   min_z_veg <- min(z[class == 5])
  # } else {
  #   min_z_veg <- min(z)  # Use all values if none of the classes match
  # }
  # print(min_z_veg)
  # z[class == 2 & (z > min_z_veg | z == Inf | z == -Inf)] <- min_z_veg - 3
  # if (sum(class == 1) > 0){
  #   z[class == 1 & z > min_z_veg | z == Inf | z == -Inf] <- min_z_veg - 3
  # }
  # if (sum(class == 0) > 0){
  #   z[class == 0 & z > min_z_veg | z == Inf | z == -Inf] <- min_z_veg - 3
  # }
  
  # if (any(class %in% c(3, 4, 5))) {
  #   min_veg_z <- min(z[class %in% c(3, 4, 5)])
  # }
  # if (any(class %in% 2)) {
  #   max_ground_z <- max(z[class == 2])
  # }
  # z[class == 2] <- max_ground_z
  # if (sum(class == 1) > 0){
  #   z[class == 1] <- max_ground_z
  # }
  # if (sum(class == 0) > 0){
  #   z[class == 0] <- max_ground_z
  # }
  
  # a voir
  # min_z <- min_z_veg + 2
  # min_z <- min(z) + 2 # comme dtm ?
  # min_z_veg <- max(z[class == 2]) + 2
  # print(min_z_veg)
  # print(min(z[class == 2]))
  # print(max(z[class == 2]))
  lad_i <- LAD(z)
  lad_i$z <- floor(lad_i$z) + 0.5
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  
  print(lad)
  print(sum(as.numeric(lad$lad), na.rm = TRUE))
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    
    # CV_LAD
    # CV_LADs[k] <- 100 * sd(lad_i$lad[lad_i$z >= i], na.rm = T) / mean(lad_i$lad[lad_i$z >= i], na.rm = T)
    
    # Inc
    k <- k + 1
  }
  
  # Calculate overall CV_LAD (Coefficient of Variation for all LADs)
  overall_CV_LAD <- 100 * sd(lad$lad, na.rm = TRUE)/ mean(lad$lad, na.rm = TRUE)
  
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  # CV_LADs_list <- setNames(as.list(CV_LADs), paste0("CV_LAD_", new_z))
  
  # Return PAIs and CV_LAD
  combined_list <- c(PAIs_list, LADs_list)
  # combined_list <- c(LADs_list)
  # combined_list <- c(PAIs_list, LADs_list, cv_lad = overall_CV_LAD)
  # LADs
}

myPAD_dsm_new <- function(z) {
  
  # z <- z + 40
  # las_i <- classify_poi(las_i, class = int(2), poi = ~las_i@data$Z < 2)
  # las_i$Z[las_i$Classification == 2] <- 0
  
  lad_i <- LAD(z)
  lad_i$z <- floor(lad_i$z) + 0.5
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    k <- k + 1
  }
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  
  # Return PAIs and CV_LAD
  # combined_list <- c(PAIs_list, LADs_list, CV_LADs_list)
  combined_list <- c(PAIs_list, LADs_list)
  combined_list <- c(LADs_list)
  return(combined_list)
}

myPAD_dtm_new <- function(z, class) {
  
  # max_z <- max(z)
  z[class == 4] <- z[class == 4] + 40 - max(z[class == 4])
  
  lad_i <- LAD(z)
  lad_i$z <- floor(lad_i$z) + 0.5
  
  # Fill with 0s
  new_z <- seq(2.5, 39.5)
  missing_z <- setdiff(new_z, lad_i$z)
  missing_data <- data.frame(z = missing_z, lad = rep(NA, length(missing_z)))
  lad <- rbind(lad_i, missing_data)
  lad <- lad[order(lad$z), ]
  # lad$lad[is.na(lad$lad)] <- 0
  row.names(lad) <- NULL
  
  # Create numeric vectors
  LADs <- numeric(length(new_z))
  PAIs <- numeric(length(new_z))
  # CV_LADs <- numeric(length(new_z))
  
  k <- 1
  for (i in new_z) {
    if (!is.na(lad$lad[lad$z == i])){
      LADs[k] <- lad$lad[lad$z == i]
      PAIs[k] <- sum(lad$lad[lad$z >= i], na.rm = T)
    }
    else {
      LADs[k] <- NA
      PAIs[k] <- NA
    }
    k <- k + 1
  }
  # Convert numeric vectors in lists
  LADs_list <- setNames(as.list(LADs), paste0("LAD_Layer_", new_z))
  PAIs_list <- setNames(as.list(PAIs), paste0("PAI_", new_z))
  
  # Return PAIs and CV_LAD
  # combined_list <- c(PAIs_list, LADs_list, CV_LADs_list)
  combined_list <- c(PAIs_list, LADs_list)
  combined_list <- c(LADs_list)
  return(combined_list)
}

end_lad_dtm <- function(z, class){
  z[class == 4] <- z[class == 4] + 40 - max(z[class == 4])
  df <- LAD(z)
  full_bins <- seq(2.5, 39.5, by = 1)
  lad_values <- setNames(rep(NA, length(full_bins)), 
                         paste0("LAD_Layer_", full_bins))
  matched <- match(df$z, full_bins)
  valid <- !is.na(matched)
  lad_values[matched[valid]] <- df$lad[valid]
  return(as.list(lad_values))
}

end_lad_dsm <- function(z){
  df <- LAD(z)
  full_bins <- seq(2.5, 39.5, by = 1)
  lad_values <- setNames(rep(NA, length(full_bins)), 
                         paste0("LAD_Layer_", full_bins))
  matched <- match(df$z, full_bins)
  valid <- !is.na(matched)
  lad_values[matched[valid]] <- df$lad[valid]
  return(as.list(lad_values))
}


end_dtm_norm <- function(las){
  if (is.null(las) || npoints(las) == 0) return(NULL)
  dtm <- rasterize_terrain(las, res = 1, pkg = "terra", algorithm = tin()) #kriging(k = 40)
  # dtm <- rasterize_terrain(las, res = 1, pkg = "terra", knnidw(k = 10L, p = 2))
  # nlas <- las - dtm
  # nlas <- normalize_height(las, algo = dtm)
  # nlas <- normalize_height(las, algo = rasterize_terrain(las, res = 1, pkg = "terra", algorithm = tin()))
  nlas <- normalize_height(las, dtm)
  rm(dtm)
  # print(summary(nlas$Z))
  nlas <- filter_poi(nlas, Z <= 40)
  # print(summary(nlas$Z))
  nlas <- classify_poi(nlas, class = 2L, poi = ~nlas@data$Z < 2) # 0.5 2
  nlas$Z[nlas$Classification == 2L] <- 0
  # las$Z <- las$Z + 40 - max(las$Z)
  # print(summary(las@data$Classification))
  # nlas <- filter_poi(nlas, buffer == 0)
  gc()
  return(nlas)
}

end_dsm_norm <- function(las){
  if (is.null(las) || npoints(las) == 0) return(NULL)
  dsm <- rasterize_canopy(las, res = 1, pkg = "terra", algorithm = pitfree())
  nlas <- normalize_height(las, algo = dsm)
  rm(dsm)
  nlas <- filter_poi(nlas, !is.na(Z))
  # print(summary(nlas$Z))
  nlas$Z[nlas$Classification == 4 & nlas$Z > 0] <- 0
  nlas$Z[nlas$Classification == 4] <- nlas$Z[nlas$Classification == 4] + 40 - max(nlas$Z[nlas$Classification == 4])
  # Reclassify ground points with Z < 2
  nlas <- classify_poi(nlas, class = int(2), poi = ~Z < 2)
  nlas$Z[nlas$Classification == 2L] <- 0
  # print(summary(nlas$Z))
  gc()
  return(nlas)
}

# las_dsm$Z[las_dsm$Classification == 4 & las_dsm$Z > 0] <- 0
# las_dsm$Z[las_dsm$Classification == 4] <- las_dsm$Z[las_dsm$Classification == 4] + 40 - max(las_dsm$Z[las_dsm$Classification == 4])
# las_dsm <- classify_poi(las_dsm, class = int(2), poi = ~las_dsm@data$Z < thr)
# las_dsm$Z[las_dsm$Classification == 2] <- 0
# summary(las_dsm$Z)


is.empty <- function(x) {
  inherits(x, "LAS") && (is.null(x) || npoints(x) == 0)
}


fill.na <- function(x, i=5) { 
  if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }
  }







# Function to compute shadows analysis

#' Perform Shadows Analysis
#'
#' This function performs shadows analysis on a DSM (Digital Surface Model)
#' raster to determine shaded areas based on sun angle information.
#'
#' @param datapath Path to the data directory.
#' @param site Site identifier.
#'
#' @details This function opens a DSM raster, calculates sun angle based on
#' site location and time, and computes shadows using ray tracing techniques.
#' The resulting shaded areas are saved as raster files.
#'
#' @return 10-meters shadows raster + save files.
#'
#' @examples
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Mormal")
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Blois")
#' perform_shadows_analysis(datapath = "/path/to/data", site = "Aigoual")
#'
#' @export
perform_shadows_analysis <- function(dsm,
                                     metrics_dir){
  
  # Results storage
  metrics_dir <- file.path("../../03_RESULTS",
                           site,
                           "Metrics")
  results_dir <- file.path(metrics_dir,
                           "Raw")
  
  if (site == "Aigoual"){
    t <- as.POSIXct("2021-07-11 10:40:31", tz = "UTC")
    lat <- 44.12
    lon <- 3.52
  } else if (site == "Blois"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 47.57
    lon <- 1.29
  } else if (site == "Mormal"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 50.20
    lon <- 3.74
  } else{
    stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  }
  
  # Get sun angle and altitude
  df_sun <- oce::sunAngle(t, longitude = lon, latitude = lat)
  
  # Convert terra::rast in matrix
  dsm_mat <- terra::as.matrix(dsm, wide=TRUE)
  
  # Calculate shadows
  shade <- ray_shade(dsm_mat,
                     sunaltitude = df_sun$altitude,
                     sunangle = df_sun$azimuth,
                     lambert=FALSE)
  
  # Flip-ud to revert the one done in the function
  shade <- shade[nrow(shade):1, , drop = FALSE]
  
  # Final
  shade_rast <- terra::rast(shade,
                            ext = terra::ext(dsm),
                            crs = terra::crs(dsm))
  return(shade_rast)
}

my_hillshade <- function(dsm){
  
  # Results storage
  # metrics_dir <- file.path("../../03_RESULTS",
  #                          site,
  #                          "Metrics")
  # results_dir <- file.path(metrics_dir,
  #                          "Raw")
  
  if (site == "Aigoual"){
    t <- as.POSIXct("2021-07-11 10:40:31", tz = "UTC")
    lat <- 44.12
    lon <- 3.52
  } else if (site == "Blois"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 47.57
    lon <- 1.29
  } else if (site == "Mormal"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 50.20
    lon <- 3.74
  } else{
    stop("Error: Site must be Mormal, Blois, or Aigoual.\n")
  }
  
  df_sun <- oce::sunAngle(t, longitude = lon, latitude = lat)
  slope_dsm <- raster::terrain(dsm, v = "slope", unit = "radians")
  aspect_dsm <- raster::terrain(dsm, v = "aspect", unit = "radians")
  shade_dsm <- hillShade(raster::raster(slope_dsm),
                         raster::raster(aspect_dsm),
                         angle = df_sun$altitude,
                         direction = df_sun$azimuth)
  # Final
  shade_rast <- terra::rast(shade_dsm,
                            ext = terra::ext(dsm),
                            crs = terra::crs(dsm))
  return(1 - shade_rast)
}

# Function to calculate solar angles
calculate_solar_angles <- function(lat, long, date_time, slope, azimuth) {
  
  # Extract solar info
  df_sun <- oce::sunAngle(t, longitude = lon, latitude = lat)
  
  # Solar elevation angle
  solar_elevation <- df_sun$altitude
  
  # Solar zenith angle
  solar_zenith <- 90 - solar_elevation
  
  # Solar azimuth angle
  solar_azimuth <- df_sun$azimuth
  
  # Calculate solar incidence angle
  theta_i <- acos(sin(solar_zenith * pi/180) * sin(slope * pi/180) +
                    cos(solar_zenith * pi/180) * cos(slope * pi/180) *
                    cos((solar_azimuth - azimuth) * pi/180)) * 180/pi
  
  return(list(zenith_angle = solar_zenith, incidence_angle = theta_i))
}
