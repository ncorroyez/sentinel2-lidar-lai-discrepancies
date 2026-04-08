# ---
# title: "functions_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"
# ---

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

# Function to create a CHM for a given site
create_chm <- function(data_path, site, results_path) {
  
  # Create DTM, DSM, and CHM directories
  if (!dir.exists(file.path(results_path, "dtm"))) {
    dir.create(file.path(results_path, "dtm"), showWarnings = FALSE)
  }
  if (!dir.exists(file.path(results_path, "dsm"))) {
    dir.create(file.path(results_path, "dsm"), showWarnings = FALSE)
  }
  if (!dir.exists(file.path(results_path, "chm"))) {
    dir.create(file.path(results_path, "chm"), showWarnings = FALSE)
  }
  
  # Manage LAS data
  las_dir <- file.path(data_path, site, "LiDAR/2-las_utm")
  las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)
  
  # Mask specific filenames for Mormal and Aigoual site
  if (site == "Mormal") { # Acquisition problems
    filenames_to_mask <- c("LAS_754000_7016000.las", "LAS_754000_7016500.las", 
                           "LAS_754500_7016000.las", "LAS_754500_7016500.las")
    if (any(basename(las_files) %in% filenames_to_mask)) {
      las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
    } 
  }
  if (site == "Aigoual") { # DSM (pitfall) computation problem
    filenames_to_mask <- c("LAS_745500_6337000.las")
    if (any(basename(las_files) %in% filenames_to_mask)) {
      las_files <- las_files[!basename(las_files) %in% filenames_to_mask]
    } 
  }
  
  # LAS Catalog
  ctg <- readLAScatalog(las_files)
  plot(ctg)
  
  # LiDR optimization
  opt_chunk_size(ctg) <- 0 # Processing by files
  opt_chunk_buffer(ctg) <- 10
  set_lidr_threads(1)
  
  # DTM
  opt_output_files(ctg) <- paste0(file.path(results_path,
                                            "dtm/{XLEFT}_{YBOTTOM}"))
  start_time <- Sys.time()
  dtm <- rasterize_terrain(ctg, 
                           res = 1, 
                           pkg = "terra",
                           algorithm = tin(extrapolate = kriging())
                           ) # knnidw() # kriging(k = 40) tin()
  end_time <- Sys.time()
  cat("DTM computation time:", format_time(start_time, end_time), "\n")
  cat("DTM raster size:", print_raster_size(dtm), "\n")
  cat("DTM algorithm: TIN (extrapolation: kriging)\n\n")
  
  # Define a common CRS for consistency
  common_crs <- st_crs(ctg)$wkt
  
  # Ensure the CRS of DTM matches the common CRS
  terra::crs(dtm) <- common_crs
  
  # DSM
  opt_output_files(ctg) <- paste0(file.path(results_path,
                                            "dsm/{XLEFT}_{YBOTTOM}"))
  start_time <- Sys.time()
  dsm <- rasterize_canopy(ctg,
                          res = dtm, # 1,
                          pkg = "terra",
                          algorithm = pitfree(thresholds = c(0, 2, 5, 10, 15),
                                              max_edge = c(0, 1),
                                              subcircle = 0.5)
                          ) # p2r()
  end_time <- Sys.time()
  cat("DSM computation time:", format_time(start_time, end_time), "\n")
  cat("DSM raster size:", print_raster_size(dsm), "\n")
  cat("DSM algorithm: pitfree\n\n")
  
  # Ensure the CRS of DSM matches the common CRS
  terra::crs(dsm) <- common_crs
  
  # CHM
  chm <- dsm - dtm
  chm_dir <- file.path(results_path, "chm")
  writeRaster(chm, filename = file.path(chm_dir, "chm.tif"), overwrite=TRUE)
}