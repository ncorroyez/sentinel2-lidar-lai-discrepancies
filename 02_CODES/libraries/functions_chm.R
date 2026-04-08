# ---
# title: "functions_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-06-03"
# ---

#' Format the Computation Time
#'
#' This function formats the computation time from start to end into a string.
#'
#' @param start A POSIXct object representing the start time.
#' @param end A POSIXct object representing the end time.
#' @return A string representing the elapsed time (hours, minutes, and seconds).
#' @export
format_time <- function(start, end) {
  elapsed <- as.numeric(difftime(end, start, units = "secs"))
  hours <- floor(elapsed / 3600)
  minutes <- floor((elapsed %% 3600) / 60)
  seconds <- elapsed %% 60
  
  sprintf("Elapsed time: %02d hours, %02d minutes, %.2f seconds", 
          hours, minutes, seconds)
}

#' Print Raster Size
#'
#' This function returns the dimensions of a raster object as a string.
#'
#' @param raster A raster object.
#' @return A string representing the dimensions of the raster (rows x cols).
#' @export
print_raster_size <- function(raster) {
  size <- dim(raster)
  sprintf("%dx%d (rows x cols)", size[1], size[2])
}

dsm_tile <- function(las, resolution) {
  dsm <- rasterize_canopy(las,
                          res = 1,
                          algorithm = p2r(subcircle = 0.15),
                          pkg = "terra")
  # dsm <- rasterize_canopy(las,
  #                         res = resolution,
  #                         pkg = "terra",
  #                         algorithm = pitfree(subcircle = 0.15)
  
  # ) # p2r()
  # print(dsm)
  if (is.null(dsm)) {
    message("rasterize_canopy returned NULL — skipping.")
    return(NULL)
  }
  # dsm <- terra::focal(rasterize_canopy(las,
  #                                      res = 1,
  #                                      algorithm = p2r(subcircle = 0.15),
  #                                      pkg = "terra"),
  #                     w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
  # return(dsm)
  w <- matrix(1, 3, 3)
  # filled <- terra::focal(dsm, w = w, fun = fill.na)
  filled <- terra::focal(dsm, w = w, fun = mean, na.rm = T)
  return(filled)
}

#' Create a Canopy Height Model (CHM)
#'
#' This function creates a Canopy Height Model (CHM) for a given site.
#'
#' @param data_path A string representing the path to the data directory.
#' @param site A string representing the site name.
#' @param results_path A string representing the path to the results directory.
#' @return None. The function writes the CHM raster to a file.
#' @export
create_chm <- function(data_path, 
                       results_dir, 
                       site, 
                       s2_rast, 
                       leaf_state,
                       resolution) {
  
  results_path <- file.path(results_dir, site, "LiDAR", leaf_state)
  
  # Create DTM, DSM, and CHM directories if they do not exist
  dirs <- c("dtm", "dsm", "chm")
  for (dir in dirs) {
    dir_path <- file.path(results_path, dir, paste0("res_", resolution, "_m"))
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
    }
  }
  
  # Manage LAS data
  las_dir <- file.path(data_path, site, "LiDAR", leaf_state, "2-las_utm")
  # print(las_dir)
  las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Mask specific filenames for Mormal and Aigoual site
  # if (site == "Mormal") { # Acquisition problems
  #   filenames_to_mask <- c("LAS_754000_7016000.las", "LAS_754000_7016500.las", 
  #                          "LAS_754500_7016000.las", "LAS_754500_7016500.las")
  #   if (any(basename(las_files) %in% filenames_to_mask)) {
  #     las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
  #   } 
  # }
  # if (site == "Aigoual") { # DSM (pitfall) computation problem
  #   filenames_to_mask <- c("LAS_745500_6337000.las")
  #   if (any(basename(las_files) %in% filenames_to_mask)) {
  #     las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
  #   } 
  # }
  # if (site == "Blois") { # DSM (pitfall) computation problem
  #   filenames_to_mask <- c("LAS_567000_6724500.las")
  #   if (any(basename(las_files) %in% filenames_to_mask)) {
  #     las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
  #   } 
  # }
  
  # Mask specific filenames for problematic sites
  mask_filenames <- list(
    "Aigoual" = c("LAS_745500_6337000.las"),
    "Blois" = c("LAS_567000_6724500.las"),
    "Mormal" = c("LAS_754000_7016000.las", "LAS_754000_7016500.las", 
                 "LAS_754500_7016000.las", "LAS_754500_7016500.las"
                 # "LAS_747500_7013000.las", "LAS_748500_7014500.las"
    )
  )
  
  if (site %in% names(mask_filenames)) {
    las_files <- las_files[!basename(las_files) %in% mask_filenames[[site]]]
  }
  # print(las_files)
  # LAS Catalog
  ctg <- readLAScatalog(las_files)
  
  # LiDR optimization
  opt_chunk_size(ctg) <- 0 # Processing by files
  opt_chunk_buffer(ctg) <- 10
  opt_select(ctg) <- "xyzc"
  # opt_filter(ctg) <- "-keep_first"
  set_lidr_threads(1)
  # options(lidR.catalog_check = FALSE)
  
  st_crs(ctg) <- st_crs(s2_rast)
  # lidR.check.nested.parallelism = FALSE
  # plot(ctg)
  # opt_stop_early(ctg) <- FALSE
  
  # DTM
  # dtm <- rast(file.path(results_dir, site,
  #                       "LiDAR/dtm/res_1_m/rasterize_terrain.vrt"))
  
  opt_output_files(ctg) <- file.path(results_path,
                                     # "pitfree",
                                     "dtm",
                                     paste0("res_", resolution, "_m"),
                                     "{XLEFT}_{YBOTTOM}")
  start_time <- Sys.time()
  dtm <- rasterize_terrain(ctg,
                           res = resolution,
                           pkg = "terra",
                           algorithm = tin()
  ) # knnidw() # kriging(k = 40) tin()
  end_time <- Sys.time()
  cat("DTM computation time:", format_time(start_time, end_time), "\n")
  cat("DTM raster size:", print_raster_size(dtm), "\n")
  cat("DTM algorithm: TIN \n\n")
  
  # dtm <- terra::rast(paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/",
  #                           site, "/LiDAR/", leaf_state,
  #                           "/dtm/res_1_m/rasterize_terrain.vrt"))
  
  # Define a common CRS for consistency
  common_crs <- st_crs(ctg)$wkt
  
  # Ensure the CRS of DTM matches the common CRS
  crs(dtm) <- common_crs
  
  # DSM
  cat("Start DSM algorithm \n\n")
  opt_output_files(ctg) <- file.path(results_path,
                                     # "pitfree",
                                     "dsm",
                                     paste0("res_", resolution, "_m"),
                                     "{XLEFT}_{YBOTTOM}")
  start_time <- Sys.time()
  dsm <- rasterize_canopy(ctg,
                          # res = dtm, # 1,
                          res <- resolution,
                          pkg = "terra",
                          algorithm = p2r(subcircle = 0.15)
                          # algorithm = pitfree()
                          # 
  ) # p2r()
  # dsm <- terra::focal(dsm, w = matrix(1, 3, 3), fun = mean)
  # dsm <- terra::rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/LiDAR/dsm/res_10_m/rasterize_canopy.vrt")
  # dsm <- catalog_apply(ctg, function(las) dsm_tile(las, resolution))
  # raster_files <- list.files(
  #   path = file.path(results_path, "dsm", paste0("res_", resolution, "_m")),
  #   pattern = "\\.tif$", full.names = TRUE)
  # 
  # dsm <- terra::vrt(raster_files,
  #                   filename = file.path(results_path, "dsm", paste0("res_", resolution, "_m"), "rasterize_canopy.vrt"),
  #                   overwrite = TRUE)
  terra::crs(dsm) <- common_crs
  dsm <- terra::project(dsm, dtm)
  end_time <- Sys.time()
  cat("DSM computation time:", format_time(start_time, end_time), "\n")
  cat("DSM raster size:", print_raster_size(dsm), "\n")
  cat("DSM algorithm: p2r\n\n")
  
  # Ensure the CRS of DSM matches the common CRS
  # terra::crs(dsm) <- common_crs
  print(dtm)
  print(dsm)
  # CHM
  chm <- dsm - dtm
  # chm <- rast(file.path(results_dir, site,
  #                       "LiDAR/dsm/res_1_m/rasterize_canopy.vrt")) -
  #   rast(file.path(results_dir, site,
  #                  "LiDAR/dtm/res_1_m/rasterize_terrain.vrt"))
  chm_dir <- file.path(results_path, 
                       # "pitfree",
                       "chm", paste0("res_", resolution, "_m"))
  writeRaster(chm, filename = file.path(chm_dir, "chm.tif"), overwrite = TRUE)
}
