# ---
# title: "4.test_whittaker2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2026-01-13"
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

# 2. PARAMETRES ---------------------------------------------------------------
set.seed(42)
data_dir <- '/media/corroyez/MyPassport/01_DATA'
results_dir <- '../../03_RESULTS'

# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Blois")
distribs <- c("atbd", "atbd_optim_common")
distribs <- c("atbd")
# distribs <- c("atbd_optim_common")

# Common temporal grid for Machine Learning (Regular 10-day step)
common_dates <- seq(from = as.Date("2021-06-05"), 
                    to   = as.Date("2021-09-30"), 
                    by   = "10 days")
load_start <- as.Date("2021-05-01") 
load_end   <- as.Date("2021-12-31")

# LiDAR acquisition dates per site
lidar_dates <- c(
  "Aigoual" = "2021-07-18",
  "Blois"   = "2021-06-17",
  "Mormal"  = "2021-06-16"
)

# Load S2 dates file
dates_file <- file.path(data_dir, "S2_Dates_Summer2021.ods")
data_dates <- read_ods(dates_file)

# 3. MAIN LOOP ----------------------------------------------------------------

for (site in sites) {
  cat("\n========================================\n")
  cat("PROCESSING SITE:", site, "\n")
  
  # Paths and Dates
  res_path <- file.path(results_dir, site, "Metrics/Not_Masked")
  dateAcqs <- get_site_full_dates(site, data_dates) 
  lidar_date_site <- as.Date(lidar_dates[site])
  
  # Load Reference LiDAR
  path_lidar <- file.path(res_path, "lidarlai_res_10_m.tif")
  # path_lidar <- file.path(res_path, "lidarlai_optim_depth_res_10_m.tif")
  
  if(!file.exists(path_lidar)) {
    cat("[ERROR] LiDAR not found for", site, ":", path_lidar, "\n")
    next
  }
  r_lidar <- rast(path_lidar)
  
  for (distrib in distribs) {
    cat("\n  -> Distribution:", distrib, "\n")
    
    # A. Load & Sort S2 Images (Irregular Grid)
    # ----------------------------------------------------------------
    pattern <- paste0("^s2lai_(", paste(dateAcqs, collapse = "|"), ")_", distrib, "_res_10_m\\.tif$")
    files <- list.files(res_path, pattern = pattern, full.names = TRUE)
    
    if(length(files) < 2) { cat("     [SKIP] Not enough S2 files.\n"); next }
    
    # Sort chronologically
    dates_found <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
    ord <- order(dates_found)
    files <- files[ord]
    current_dates <- dates_found[ord] # Irregular dates
    
    r_stack <- rast(files)
    terra::time(r_stack) <- current_dates
    
    # B. Harmonization (Irregular -> Regular 10-day Grid)
    # ----------------------------------------------------------------
    # We interpolate FIRST so that Whittaker (Step D) gets equidistant data.
    # poids du temps pour l'interpolation --> OK, 
    # validation microclimat/GEDI,
    # on irregular dates is ok,
    # interrogation nuages/brume et ecart atbd/opt_common
    cat("     1. Harmonization (Gap-Filling + Resampling to 10-day)...\n")
    
    # Numeric dates for approx
    dates_irregular_num <- as.numeric(current_dates)
    dates_regular_num   <- as.numeric(common_dates)
    
    # Function to interpolate pixel-wise
    harmonize_fun <- function(y) {
      if (sum(!is.na(y)) < 2) {
        return(rep(NA, length(dates_regular_num)))
      }
      # Linear interpolation (Rule 2 = constant extrapolation if needed)
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
    # writeRaster(r_regular, out_name_regular, overwrite = TRUE)
    
    # C. Lambda Calibration (Optimizing Whittaker on Regular Grid vs LiDAR)
    # ----------------------------------------------------------------
    cat("     2. Calibration of Lambda (Whittaker vs LiDAR)...\n")
    
    # Stratified Sampling setup
    val_p98 <- global(r_lidar, fun = function(x) quantile(x, 0.98, na.rm=TRUE))[[1]]
    breaks <- seq(from = 2, to = floor(val_p98), by = 1)
    breaks <- c(breaks, floor(val_p98) + 1)
    breaks <- unique(breaks)
    r_strata <- classify(r_lidar, rcl = breaks, include.lowest = TRUE)
    r_strata <- mask(r_strata, r_regular[[1]])
    
    n_bins <- length(breaks) - 1
    target_total <- 5000
    pts_per_bin <- ceiling(target_total / n_bins)
    pts_sample <- spatSample(r_strata, 
                             size = pts_per_bin, 
                             method = "stratified", 
                             xy = TRUE, na.rm = TRUE, values = FALSE)
    
    if(nrow(pts_sample) > 0) {
      # Extract values (ID=FALSE is crucial)
      vals_lidar_sample <- extract(r_lidar, pts_sample)
      vals_s2_sample <- extract(r_regular, pts_sample)
      
      # Filter valid rows
      valid_rows <- complete.cases(vals_lidar_sample) & complete.cases(vals_s2_sample)
      vals_lidar_vec <- vals_lidar_sample[valid_rows, ]
      vals_s2_mat    <- as.matrix(vals_s2_sample[valid_rows, ])
      
      lambdas_test <- c(0.01, 0.1, 1, 5, 15, 50, 100, 500)
      # lambdas_test <- c(5, 15, 50, 100, 500)
      rmse_results <- numeric(length(lambdas_test))
      lidar_date_num <- as.numeric(lidar_date_site)
      
      # Initialize storage dataframe
      metrics_res <- data.frame(
        Lambda = lambdas_test,
        RMSE = numeric(length(lambdas_test)),
        Correlation = numeric(length(lambdas_test)),
        Bias = numeric(length(lambdas_test)),
        Slope = numeric(length(lambdas_test))
      )
      
      for(i in seq_along(lambdas_test)) {
        lam <- lambdas_test[i]
        
        # 1. Smooth the regular series
        s2_smooth_mat <- t(apply(vals_s2_mat, 1, function(y) {
          tryCatch(ptw::whit2(y, lambda = lam), error = function(e) y)
        }))
        
        # 2. Get value at exact LiDAR date
        preds_at_lidar <- apply(s2_smooth_mat, 1, function(y) {
          approx(dates_regular_num, y, xout = lidar_date_num)$y
        })
        
        # 3. Compute Metrics
        metrics_res$RMSE[i] <- sqrt(mean((preds_at_lidar - vals_lidar_vec)^2, na.rm = TRUE))
        metrics_res$Correlation[i] <- cor(preds_at_lidar, vals_lidar_vec, use = "complete.obs")
        metrics_res$Bias[i] <- mean(preds_at_lidar - vals_lidar_vec, na.rm = TRUE)
        metrics_res$Slope[i] <- coef(lm(preds_at_lidar ~ vals_lidar_vec))[2]
      }
      
      # --- DIAGNOSTIC PRINT ---
      cat("\n        --- Optimization Results ---\n")
      # Round for cleaner display
      print(format(metrics_res, digits = 2))
      cat("        ----------------------------\n")
      
      # Select best Lambda based on RMSE (standard approach)
      best_lambda <- metrics_res$Lambda[which.min(metrics_res$RMSE)]
      best_rmse   <- min(metrics_res$RMSE)
      
      cat("        => Optimal Lambda selected:", best_lambda, "(RMSE:", round(best_rmse, 3), ")\n")
      
    } else {
      cat("        [WARN] No common points. Default Lambda = 15.\n")
      best_lambda <- 15
    }
    best_lambda <- 5
    # stop()
    # D. Apply Smoothing (Whit2 on Regular Grid)
    # ----------------------------------------------------------------
    cat("     3. Application of Whittaker Smoothing...\n")
    
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
                                 paste0("LAI_Smoothed_10d_",
                                        "Lambda_", best_lambda, "_",
                                        distrib, 
                                        ".tif"))
    writeRaster(r_smoothed, out_name_smooth, overwrite = TRUE)
    
    # E. LiDAR Correction (Ratio / Prorata)
    # ----------------------------------------------------------------
    cat("     4. LiDAR Correction (Ratio)...\n")
    
    # 1. Get the S2 Smoothed value at the exact LiDAR date
    lidar_date_num <- as.numeric(lidar_date_site)
    
    get_s2_at_lidar <- function(y) {
      if (all(is.na(y))) return(NA)
      return(approx(x = dates_regular_num, 
                    y = y, 
                    xout = lidar_date_num, 
                    rule = 2)$y)
    }
    
    r_s2_ref_date <- app(r_smoothed, fun = get_s2_at_lidar, cores = 1)
    
    # 2. Align geometry if needed
    if (!compareGeom(r_lidar, r_s2_ref_date, stopOnError = FALSE)) {
      print("cc")
      r_lidar_proj <- project(r_lidar, r_s2_ref_date, method = "bilinear")
    } else {
      r_lidar_proj <- r_lidar
    }
    
    # 3. Calculate Ratio
    r_ratio <- r_lidar_proj / r_s2_ref_date
    # r_ratio <- clamp(r_ratio, lower = 0.1, upper = 10, values = FALSE)
    r_ratio <- ifel(is.finite(r_ratio), r_ratio, 1) # Handle NA/Inf
    
    # 4. Apply Ratio to the whole time series
    r_corrected <- r_smoothed * r_ratio
    names(r_corrected) <- paste0("LAI_Corr_", common_dates)
    terra::time(r_corrected) <- common_dates
    
    # Save Final
    out_name_corr <- file.path(res_path,
                               "Smooth_TS",
                               paste0("LAI_Corrected_",
                                      "Lambda_", best_lambda, "_",
                                      distrib, 
                                      ".tif"))
    writeRaster(r_corrected, out_name_corr, overwrite = TRUE)
    
    cat("     [OK] Saved:", basename(out_name_corr), "\n")
    
  } # End Loop Distrib
} # End Loop Site