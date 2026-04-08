# ---
# title: "4.test_RF.R"
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
library(randomForest)
library(gstat)
library(sf)
library(ape)
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

distribs <- c(
  'atbd',
  'atbd_optim_common',
  'atbd_optim_common_brownmodif',
  'optim_Blois'
  # 'optim_Blois_brownmodif'
)
# distribs <- c('optim_Blois')
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
  path_lidar <- file.path(res_path, "lidarlai_res_10_m.tif")
  # path_lidar <- file.path(res_path, "lidarlai_optim_depth_res_10_m.tif")
  
  if(!file.exists(path_lidar)) {
    cat("[ERROR] LiDAR not found:", path_lidar, "\n")
    next
  }
  r_lidar <- rast(path_lidar)
  
  # Load Static LiDAR Metrics (Predictors)
  # ----------------------------------------------------------------
  # These metrics are static for the season and used as features for the RF
  # Adjust the pattern or list to match your exact filenames
  metrics_names <- c("mean", "std", "vci") # Example keywords
  lidar_metrics_files <- list.files(res_path, pattern = "_res_10_m\\.tif$", full.names = TRUE)
  
  # Filter to keep only the structural metrics you want (exclude the target 'lidarlai')
  # This logic depends on your naming convention. Here implies files like 'lcv_res_10_m.tif'
  keep_metrics <- grepl(paste(metrics_names, collapse = "|"), basename(lidar_metrics_files)) &
    !grepl("lidarlai", basename(lidar_metrics_files))
  
  metrics_files <- lidar_metrics_files[keep_metrics]
  
  if(length(metrics_files) == 0) {
    cat("[WARNING] No LiDAR metrics found for RF. Skipping site.\n")
    next
  }
  
  r_lidar_metrics <- rast(metrics_files)
  names(r_lidar_metrics) <- metrics_names
  
  # Ensure geometry matches the target LiDAR LAI
  if (!compareGeom(r_lidar, r_lidar_metrics, stopOnError = FALSE)) {
    r_lidar_metrics <- project(r_lidar_metrics, r_lidar, method = "bilinear")
  }
  cat("   Loaded LiDAR Metrics:", names(r_lidar_metrics), "\n")
  
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
                                  "smooth",
                                  paste0("LAI_Gap_10d_",
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
    
    # LiDAR Correction (RF)
    # ----------------------------------------------------------------
    cat("     3. LiDAR Correction (Random Forest)...\n")
    
    # A. Data Preparation for Training (at LiDAR Date)
    # ------------------------------------------------
    lidar_date_num <- as.numeric(lidar_date_site)
    
    # Function to get S2 value at exact LiDAR date
    get_s2_at_lidar <- function(y) {
      if (all(is.na(y))) return(NA)
      return(approx(x = dates_regular_num, y = y, xout = lidar_date_num, rule = 2)$y)
    }
    
    r_s2_ref_date <- app(r_smoothed, fun = get_s2_at_lidar, cores = 1)
    names(r_s2_ref_date) <- "S2_LAI" # Important: Fix name for the model
    
    # Align geometries if needed (Target vs S2)
    if (!compareGeom(r_lidar, r_s2_ref_date, stopOnError = FALSE)) {
      r_lidar_proj <- project(r_lidar, r_s2_ref_date, method = "bilinear")
      r_metrics_proj <- project(r_lidar_metrics, r_s2_ref_date, method = "bilinear")
    } else {
      r_lidar_proj <- r_lidar
      r_metrics_proj <- r_lidar_metrics
    }
    
    # Create the Training Stack: Target (LiDAR) + Features (S2_Ref + Metrics)
    # We rename r_lidar_proj to "Target" for clarity in the formula
    names(r_lidar_proj) <- "Target"
    
    train_stack <- c(r_lidar_proj, r_s2_ref_date, r_metrics_proj)
    
    # Convert to dataframe for Random Forest
    # Using na.rm=TRUE removes pixels where ANY layer is NA (edge effects, clouds)
    # df_train <- as.data.frame(train_stack, na.rm = TRUE)
    n_samples <- 10000
    valid_indices <- cells(train_stack)
    num_valid <- length(valid_indices)
    idx_seq <- as.integer(seq(1, num_valid, length.out = n_samples))
    sample_indices <- valid_indices[idx_seq]
    
    df_train <- train_stack[sample_indices]
    df_train <- na.omit(df_train)
    df_train$Bias <- df_train$Target - df_train$S2_LAI
    df_train_model <- df_train[, !names(df_train) %in% c("Target", "S2_LAI")]
    # stop()
    # df_train <- spatSample(train_stack, 
    #                        size = n_samples, 
    #                        method = "regular", # "random"
    #                        na.rm = TRUE,
    #                        as.df = TRUE, 
    #                        values = TRUE)
    
    # B. Train Random Forest Model
    # ------------------------------------------------
    cat("      -> Training RF model (n =", nrow(df_train), "pixels)...\n")
    # stop()
    
    # Model: Target ~ S2_LAI + Metric1 + Metric2 ...
    rf_model <- randomForest(Bias ~ ., data = df_train_model, 
                             ntree = 100,
                             importance = TRUE, 
                             do.trace = FALSE)
    imp_perm <- importance(rf_model, type = 1)
    imp <- sort(rf_model$importance[,1], decreasing=TRUE)
    barplot(imp)
    
    # C. Prediction on Time Series
    # ------------------------------------------------
    cat("      -> Predicting on full time series...\n")
    
    # We loop through each time step of the harmonized S2 data.
    # For each step, we use the specific S2 LAI + the STATIC LiDAR metrics.
    corrected_layers <- list()
    for (i in 1:nlyr(r_smoothed)) {
      r_s2_t <- r_smoothed[[i]]
      pred_bias <- predict(r_metrics_proj, rf_model, na.rm = TRUE)
      r_pred_t <- r_s2_t + pred_bias
      r_pred_t_cl <- clamp(r_pred_t, lower = 0.001, upper = 15, values = TRUE)
      corrected_layers[[i]] <- r_pred_t_cl
    }
    
    r_corrected <- rast(corrected_layers)
    terra::time(r_corrected) <- common_dates
    names(r_corrected) <- paste0("LAI_Corr_RF_", common_dates)
    
    # Save
    out_name_corr <- file.path(res_path,
                               "Smooth_TS",
                               "smooth", 
                               "RF/LAI_ALS",
                               paste0("LAI_Corrected_RF_",
                                      distrib, 
                                      ".tif"))
    
    if (!dir.exists(dirname(out_name_corr))) {
      dir.create(dirname(out_name_corr), showWarnings = FALSE, recursive = TRUE)
    }
    
    writeRaster(r_corrected, out_name_corr, overwrite = TRUE)
    cat("     [OK] Saved:", basename(out_name_corr), "\n")
    
    # ----------------------------------------------------------------
    # D. Check Spatial Autocorrelation (Residuals)
    # ----------------------------------------------------------------
    
    # 1. Retrieve Coordinates for the training samples
    # We need to know WHERE the points are to calculate distances
    coords_train <- xyFromCell(r_lidar_proj, sample_indices)
    
    # 2. Calculate Model Residuals (OOB)
    # rf_model$predicted contains the Out-Of-Bag predictions for the training set
    residuals <- df_train$Bias - rf_model$predicted
    
    # Create a temporary dataframe for gstat
    df_spatial <- data.frame(
      x = coords_train[,1],
      y = coords_train[,2],
      resid = residuals
    )
    
    # 3. Compute Empirical Variogram
    # formula: resid ~ 1 means we assume constant mean (no trend)
    # cutoff: max distance to look at (e.g., 2000 meters). Adjust based on your site size.
    # width: bin size (e.g., calculate correlation every 50 meters)
    cat("      -> Computing variogram on residuals...\n")
    
    var_geo <- variogram(resid ~ 1, 
                         locations = ~x+y, 
                         data = df_spatial, 
                         cutoff = 2000,  # Look up to 2km
                         width = 50)     # Lag size of 50m
    
    print(plot(var_geo, 
               main = paste("Variogram of RF Residuals -", site),
               xlab = "Distance (m)", 
               ylab = "Semivariance"))
    
    # Moran
    idx_sub <- sample(1:nrow(df_spatial), size = nrow(df_spatial))
    dists_inv <- 1 / as.matrix(dist(df_spatial[idx_sub, c("x", "y")]))
    diag(dists_inv) <- 0
    
    # Calculate Moran's I
    moran_res <- Moran.I(df_spatial$resid[idx_sub], dists_inv)
    
    cat("      -> Moran's I (p-value):", moran_res$p.value, "\n")
    cat("      -> Moran's I (observed):", moran_res$observed, "\n")
    
    # ----------------------------------------------------------------
    # E. Visual Check (Prediction vs Target at LiDAR Date)
    # ----------------------------------------------------------------
    
    # 1. Find the layer closest to the LiDAR date in the result
    # We look for the index in common_dates that is closest to lidar_date_site
    # idx_closest <- which.min(abs(as.numeric(common_dates) - as.numeric(lidar_date_site)))
    # 
    # # Extract the predicted map for that specific date
    # r_pred_at_lidar <- r_corrected[[idx_closest]]
    # 
    # # 2. Calculate the Difference map (Residuals)
    # # Difference = Predicted (RF) - Observed (LiDAR)
    # # Ideally, this should be close to 0 (green) everywhere.
    # r_diff <- r_pred_at_lidar - r_lidar_proj
    # 
    # # Colors
    # col_lai <- colorRampPalette(c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))(100)
    # col_resid <- colorRampPalette(c("blue", "white", "red"))(100)
    # 
    # # 3. Plot Comparison
    # par(mfrow = c(1, 3)) # 3 plots side by side
    # 
    # # Plot 1: Target (LiDAR)
    # plot(r_lidar_proj, main = "Target (LiDAR)", range = c(0, 7), col = col_lai, mar=c(3,3,3,4))
    # 
    # # Plot 2: Prediction (RF)
    # plot(r_pred_at_lidar, main = "Prediction (RF)", range = c(0, 7), col = col_lai, mar=c(3,3,3,4))
    # 
    # # Plot 3: Residuals (Diff)
    # plot(r_diff, main = "Residuals (Pred - Target)", range = c(-2, 2), col = col_resid, mar=c(3,3,3,4))
    
    # ----------------------------------------------------------------
    # F. Model Evaluation (Quantitative Metrics)
    # ----------------------------------------------------------------
    cat("      -> Computing evaluation metrics (OOB)...\n")
    
    # TEST ON SUBSET
    
    # 1. Get Observations and OOB Predictions
    # rf_model$predicted contains predictions for samples NOT used in the tree construction
    # This acts as an independent validation set.
    obs_lai       <- df_train$Target
    pred_bias_oob <- rf_model$predicted
    pred_lai_oob  <- df_train$S2_LAI + pred_bias_oob
    
    # 2. Calculate Metrics
    # RMSE (Root Mean Square Error) - Magnitude of error
    rmse_val <- sqrt(mean((pred_lai_oob - obs_lai)^2))
    
    # R-squared (Coefficient of Determination)
    # 1 - (Sum of Squared Residuals / Total Sum of Squares)
    sse <- sum((pred_lai_oob - obs_lai)^2)
    sst <- sum((obs_lai - mean(obs_lai))^2)
    r2_val <- 1 - (sse / sst)
    
    # Bias (Mean Error) - Systematic over/under-estimation
    # Positive = Model overestimates; Negative = Model underestimates
    bias_val <- mean(pred_lai_oob - obs_lai)
    
    # MAE (Mean Absolute Error)
    mae_val <- mean(abs(pred_lai_oob - obs_lai))
    
    # 3. Save Metrics to CSV
    metrics_df <- data.frame(
      Site = site,
      Model = "RandomForest",
      RMSE = round(rmse_val, 4),
      R2 = round(r2_val, 4),
      Bias = round(bias_val, 4),
      MAE = round(mae_val, 4),
      N_Samples = nrow(df_train)
    )
    cat("      [RESULT] RMSE:", round(rmse_val, 3),
        "| R2:", round(r2_val, 3), 
        "| Bias:", round(bias_val, 3), "\n")
    
    par(mfrow = c(1, 1))
    # Density scatter plot using standard 'graphics' (lighter than ggplot2)
    # Use a simple color density representation
    smoothScatter(x = obs_lai, y = pred_lai_oob,
                  xlab = "Observed LAI (LiDAR)",
                  ylab = "Predicted LAI (S2 + RF Bias)",
                  main = paste("Model Performance -", site),
                  colramp = colorRampPalette(c("white", "lightblue", "blue", "darkblue", "red")))
    
    # Add 1:1 line (Perfect fit)
    abline(0, 1, col = "black", lwd = 2, lty = 2)
    
    # Add Fit line (Regression)
    abline(lm(pred_lai_oob ~ obs_lai), col = "red", lwd = 1)
    
    # Add metrics text on plot
    legend("topleft", 
           legend = c(paste("RMSE =", round(rmse_val, 2)),
                      paste("R2 =", round(r2_val, 2)),
                      paste("Bias =", round(bias_val, 3))),
           bty = "n", cex = 1.2)
    
  } # End Loop Distrib
} # End Loop Site