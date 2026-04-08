# ---
# title: "functions_shadows_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-27"
# ---

# Rotate the matrix 90 degrees clockwise
rotate90_clockwise <- function(mat) {
  t(apply(mat, 2, rev))
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
perform_shadows_analysis <- function(datapath,
                                     site){
  # Open DSM
  dsm <- terra::rast(file.path(datapath, 
                               site, 
                               "LiDAR/2-las_utm/dsm/rasterize_canopy.vrt"))
  
  # Results storage
  results <- "../../03_RESULTS"
  metrics_dir <- file.path(results,
                           site,
                           "Metrics")
  results_dir <- file.path(metrics_dir,
                           "Raw")
  
  # Open another raster to get the right dimensions
  tmp_raster <- terra::rast(file.path(results_dir,
                                      "lidarlai_res_10_m.tif"))
  
  if (site == "Mormal"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 50.20
    lon <- 3.74
  } else if (site == "Blois"){
    t <- as.POSIXct("2021-06-14 10:50:31", tz = "UTC")
    lat <- 47.57
    lon <- 1.29
  } else if (site == "Aigoual"){
    t <- as.POSIXct("2021-07-11 10:40:31", tz = "UTC")
    lat <- 44.12
    lon <- 3.52
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
  
  # 1m
  shade_rast <- terra::rast(shade, 
                            ext = terra::ext(dsm),
                            crs = terra::crs(dsm))
  shade_rast <- terra::project(shade_rast, 
                               dsm, 
                               method = 'average')
  writeRaster(shade_rast, 
              file.path(metrics_dir, "shade_res_1_m.tif"),
              overwrite=TRUE)
  
  # 10m
  # shade_rast <- terra::rast(shade, 
  #                           ext = terra::ext(tmp_raster),
  #                           crs = terra::crs(tmp_raster))
  shade_rast <- terra::project(shade_rast,
                               tmp_raster,
                               method = 'average')
  writeRaster(shade_rast, 
              file.path(results_dir, "shade_res_10_m.tif"),
              overwrite=TRUE)
  
  return(shade_rast)
}
