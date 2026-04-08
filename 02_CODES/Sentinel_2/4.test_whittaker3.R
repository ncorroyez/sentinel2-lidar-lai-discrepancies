# ---
# title: "4.test_whittaker3.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2026-01-14"
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
library(terra)
library(zoo)
library(readODS)
library(stringr)
library(signal)
library(ptw)
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_general_tools.R")

# PARAMETERS ---------------------------------------------------------------
set.seed(42)
data_dir <- '/media/corroyez/MyPassport/01_DATA'
results_dir <- '../../03_RESULTS'

sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Blois")
distribs <- c("atbd", "atbd_optim_common")
distribs <- c("atbd_optim_common_brownmodif")
distribs <- c("optim_Blois_brownmodif")

distribs <- c('atbd', 
              'atbd_optim_common',
              'atbd_optim_common_brownmodif',
              'optim_Blois'
              # 'optim_Blois_brownmodif'
              )

# TARGET GRID (Output for Machine Learning)
# Strictly June to September
common_dates <- seq(from = as.Date("2021-06-05"), 
                    to   = as.Date("2021-09-30"), 
                    by   = "10 days")

# LOADING WINDOW (Input)
# Extended to May and Dec/Jan to anchor the interpolation at boundaries
load_start <- as.Date("2021-05-01") 
load_end   <- as.Date("2022-01-31") 

# LiDAR Dates
lidar_dates <- c(
  "Aigoual" = "2021-07-18",
  "Blois"   = "2021-06-17",
  "Mormal"  = "2021-06-16"
)

# Load S2 Dates file (if needed for other tasks, otherwise unused in new loading logic)
dates_file <- file.path(data_dir, "S2_Dates_Summer2021.ods")
data_dates <- read_ods(dates_file)

# MAIN LOOP ----------------------------------------------------------------

for (site in sites) {
  cat("\n========================================\n")
  cat("PROCESSING SITE:", site, "\n")
  
  res_path <- file.path(results_dir, site, "Metrics/Not_Masked")
  lidar_date_site <- as.Date(lidar_dates[site])
  
  # Load Reference LiDAR
  # path_lidar <- file.path(res_path, "lidarlai_res_10_m.tif")
  path_lidar <- file.path(res_path, "lidarlai_optim_depth_res_10_m.tif")
  
  if(!file.exists(path_lidar)) {
    cat("[ERROR] LiDAR not found:", path_lidar, "\n")
    next
  }
  r_lidar <- rast(path_lidar)
  
  for (distrib in distribs) {
    cat("\n  -> Distribution:", distrib, "\n")
    
    # Load Files (Extended Window Strategy)
    # ----------------------------------------------------------------
    # 1. Match ANY date format to capture May/Dec files
    pattern <- paste0("^s2lai_\\d{4}-\\d{2}-\\d{2}_", distrib, "_res_10_m\\.tif$")
    all_files <- list.files(res_path, pattern = pattern, full.names = TRUE)
    
    if(length(all_files) == 0) { cat("     [SKIP] No files found.\n"); next }
    
    # 2. Filter dates based on load_start / load_end
    all_dates <- as.Date(str_extract(basename(all_files), "\\d{4}-\\d{2}-\\d{2}"))
    keep_idx  <- which(all_dates >= load_start & all_dates <= load_end)
    
    files <- all_files[keep_idx]
    current_dates <- all_dates[keep_idx] # Irregular dates (May -> Dec)
    
    if(length(files) < 2) { cat("     [SKIP] Not enough S2 files in window.\n"); next }
    
    r_stack <- rast(files)
    terra::time(r_stack) <- current_dates
    
    # Harmonization (Irregular -> Regular 10-day Grid)
    # ----------------------------------------------------------------
    cat("     1. Harmonization (Gap-Filling + Resampling to 10-day)...\n")
    
    dates_irregular_num <- as.numeric(current_dates)
    dates_regular_num   <- as.numeric(common_dates)
    
    harmonize_fun <- function(y) {
      # Need at least 2 points to draw a line
      if (sum(!is.na(y)) < 2) return(rep(NA, length(dates_regular_num)))
      
      # Interpolate from Irregular (x) to Regular Target (xout)
      # Because 'x' includes May/Dec, 'xout' (June/Sept) is interpolated, not extrapolated.
      return(approx(x = dates_irregular_num, 
                    y = y, 
                    xout = dates_regular_num, 
                    rule = 2)$y)
    }
    
    r_regular <- app(r_stack, fun = harmonize_fun, cores = 1)
    terra::time(r_regular) <- common_dates
    names(r_regular) <- paste0("LAI_Reg_", common_dates)
    
    # Save Intermediate (Optional)
    out_name_regular <- file.path(res_path,
                                  "Smooth_TS",
                                  paste0("LAI_Gapfilled_10d_",
                                         distrib, 
                                         ".tif"))
    if (!dir.exists(dirname(out_name_regular))) {
      dir.create(dirname(out_name_regular), showWarnings = FALSE, recursive = TRUE)
    }
    writeRaster(r_regular, out_name_regular, overwrite = TRUE)
    
    # Whittaker
    # ---------------------------------------------------------------
    cat("     2. Application of Whittaker Smoothing...\n")
    
    best_lambda <- 1
    whittaker_fun <- function(y) {
      if (all(is.na(y))) return(y)
      return(ptw::whit2(y, lambda = best_lambda)) 
    }
    
    r_smoothed <- app(r_regular, fun = whittaker_fun, cores = 1)
    terra::time(r_smoothed) <- common_dates
    names(r_smoothed) <- paste0("LAI_Smooth_", common_dates)
    
    # Save Intermediate (Optional)
    out_name_smooth <- file.path(res_path,
                                 "Smooth_TS",
                                 "smooth", 
                                 "gap",
                                 paste0("LAI_Gapfilled_10d_",
                                        distrib, 
                                        ".tif"))
    writeRaster(r_smoothed, out_name_smooth, overwrite = TRUE)
    
    # LiDAR Correction (Ratio)
    # ----------------------------------------------------------------
    cat("     3. LiDAR Correction (Ratio)...\n")
    
    # 1. Get S2 value at exact LiDAR date
    lidar_date_num <- as.numeric(lidar_date_site)
    
    get_s2_at_lidar <- function(y) {
      if (all(is.na(y))) return(NA)
      return(approx(x = dates_regular_num, y = y, xout = lidar_date_num, rule = 2)$y)
    }
    
    r_s2_ref_date <- app(r_smoothed, fun = get_s2_at_lidar, cores = 1)
    
    # 2. Geometry check
    if (!compareGeom(r_lidar, r_s2_ref_date, stopOnError = FALSE)) {
      r_lidar_proj <- project(r_lidar, r_s2_ref_date, method = "bilinear")
    } else {
      r_lidar_proj <- r_lidar
    }
    
    # 3. Calculate Ratio
    r_ratio <- r_lidar_proj / r_s2_ref_date
    # r_ratio <- r_s2_ref_date / r_lidar_proj
    
    # FIX: Clamp ratio to prevent outliers (shadows etc)
    # r_ratio <- clamp(r_ratio, lower = 0.1, upper = 10, values = FALSE)
    # r_ratio <- ifel(is.finite(r_ratio), r_ratio, 1)
    
    # 4. Apply Ratio
    r_corrected <- r_smoothed * r_ratio
    names(r_corrected) <- paste0("LAI_Corr_", common_dates)
    terra::time(r_corrected) <- common_dates
    
    # Save
    out_name_corr <- file.path(res_path,
                               "Smooth_TS",
                               "smooth", 
                               "corr",
                               paste0("LAI_Corrected_",
                                      distrib, 
                                      ".tif"))
    
    # Create dir if not exists
    if (!dir.exists(dirname(out_name_corr))) {
      dir.create(dirname(out_name_corr), showWarnings = FALSE, recursive = TRUE)
    }
    
    writeRaster(r_corrected, out_name_corr, overwrite = TRUE)
    cat("     [OK] Saved:", basename(out_name_corr), "\n")
    
  } # End Loop Distrib
} # End Loop Site