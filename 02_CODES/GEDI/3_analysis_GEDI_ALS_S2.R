# ---
# title: "3_analysis_GEDI_ALS_S2.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-12-09"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)
library(randomForest)
library(caret)
library(broom)
library(stringr)
library(ggplot2)
library(tidyr)
source("../libraries/GEDI/functions_GEDI.R")

set.seed(42)
data_dir <- "../../01_DATA"
metrics_dir <- "Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois")  # uncomment to test only one site

all_sites_metrics_list <- list()
bins <- seq(0, 35, 5)
gedi_old_names <- paste0("pai_", bins, "_", bins + 5, "m")
als_org_old_names <- paste0("LAD_", bins, "_", bins + 5, "_org")
als_cor_old_names <- paste0("LAD_", bins, "_", bins + 5, "_cor")
gedi_new_names <- paste0("PAI_", bins, "_", bins + 5, "m")

# --------------------------------------------------------------
# MAIN PROCESSING LOOP
# --------------------------------------------------------------
for (site in sites) {
  message("Processing site: ", site)
  results_dir <- file.path("../../03_RESULTS", site, "Metrics", metrics_dir)
  
  # A. Load GEDI
  gedi_als_data <- readRDS(file.path(data_dir, site, "GEDI", 
                                     paste0(site, "_GEDI_LiDAR_Metrics.rds")))
  gedi_als_data$dist <- sqrt(gedi_als_data$x_offset^2 + gedi_als_data$y_offset^2)
  
  # stop()
  
  gedi_als_data <- gedi_als_data %>%
    filter(max_accum > 1000) %>%
    filter(dist < 50)
  
  s2_path <- file.path(results_dir, "stack", paste0(site, "_lai_stacked_raster.tif"))
  message("  > Extracting Sentinel-2 LAI...")
  
  # Corrected Locations (Apply Offsets)
  gedi_sf_org <- st_as_sf(sf::st_drop_geometry(gedi_als_data), 
                          coords = c("x_utm", "y_utm"), 
                          crs = NA) 
  gedi_buf_org <- st_buffer(gedi_sf_org, dist = 12.5)
  
  # B. Corrected Locations
  # Directly using the columns present in your dataframe
  gedi_sf_cor <- st_as_sf(sf::st_drop_geometry(gedi_als_data), 
                          coords = c("x_corrected", "y_corrected"), 
                          crs = NA)
  gedi_buf_cor <- st_buffer(gedi_sf_cor, dist = 12.5)
  
  # 2. Run Extraction
  # LAI S2 for Original Footprints
  gedi_als_data$lai_s2_org <- extract_s2_temporal2(
    gedi_geometry = st_geometry(gedi_buf_org), 
    gedi_dates    = gedi_als_data$date, 
    s2_stack_path = s2_path
  )
  
  # LAI S2 for Corrected Footprints
  gedi_als_data$lai_s2_cor <- extract_s2_temporal2(
    gedi_geometry = st_geometry(gedi_buf_cor), 
    gedi_dates    = gedi_als_data$date, 
    s2_stack_path = s2_path
  )
  gedi_als_data <- gedi_als_data %>%
    drop_na(lai_s2_cor)
  
  # --- Original (GEDI vs ALS Original LAD/PAI) ---
  df_org <- gedi_als_data %>% 
    dplyr::select(
      # Copy GEDI PAI columns
      all_of(gedi_old_names),
      # Select original ALS LAD/PAI columns and rename them for the function
      all_of(als_org_old_names),
      lai_gedi,
      LAI_lidar_org,
      lai_s2_org
    ) %>%
    mutate(across(all_of(gedi_old_names), ~ .x / 5))
  df_org$geometry = NULL
  
  names(df_org) <- c(
    gedi_new_names,
    paste0("mean_als_pai_", bins, "_", bins + 5, "m"),
    "lai_gedi",
    "lai_lidar",
    "lai_s2"
  )
  
  # Add other required columns for plotting/stratification (e.g., site, a height metric like rh95)
  df_org$site <- site
  df_org$mean_h_max_p95 <- gedi_als_data$rh95 # Using rh95 as a proxy for stratification
  
  # --- Corrected (GEDI vs ALS Corrected LAD/PAI) ---
  df_cor <- gedi_als_data %>% 
    dplyr::select(
      # Copy GEDI PAI columns
      all_of(gedi_old_names),
      # Select corrected ALS LAD/PAI columns and rename them for the function
      all_of(als_cor_old_names),
      lai_gedi,
      LAI_lidar_cor,
      lai_s2_cor
    ) %>%
    mutate(across(all_of(gedi_old_names), ~ .x / 5))
  df_cor$geometry = NULL
  
  names(df_cor) <- c(
    gedi_new_names,
    paste0("mean_als_pai_", bins, "_", bins + 5, "m"),
    "lai_gedi",
    "lai_lidar",
    "lai_s2"
  )
  
  df_cor$site <- site
  df_cor$mean_h_max_p95 <- gedi_als_data$rh95
  
  # Calculate Shape Metrics (Vectorized)
  metrics_org <- calculate_shape_metrics_vectorized(df_org)
  metrics_cor <- calculate_shape_metrics_vectorized(df_cor)
  
  ## Calculate LAI/PAI Scalar Statistics (R2, RMSE, Bias)
  stats_lai_org <- get_scalar_stats(gedi_als_data$LAI_lidar_org, gedi_als_data$lai_gedi)
  stats_lai_cor <- get_scalar_stats(gedi_als_data$LAI_lidar_cor, gedi_als_data$lai_gedi)
  
  # Calculate Mean Profile Correlation
  mean_prof_cor_org <- mean(metrics_org$prof_cor, na.rm = TRUE)
  mean_prof_cor_cor <- mean(metrics_cor$prof_cor, na.rm = TRUE)
  
  # Store in a temporary data frame for this site
  site_results <- data.frame(
    Site = site,
    LAI_Org_r    = stats_lai_org["r"],
    LAI_Org_R2   = stats_lai_org["R2"],
    LAI_Org_RMSE = stats_lai_org["RMSE"],
    LAI_Org_Bias = stats_lai_org["Bias"],
    Prof_Cor_Org_Mean = mean_prof_cor_org,
    LAI_Cor_r    = stats_lai_cor["r"],
    LAI_Cor_R2   = stats_lai_cor["R2"],
    LAI_Cor_RMSE = stats_lai_cor["RMSE"],
    LAI_Cor_Bias = stats_lai_cor["Bias"],
    Prof_Cor_Cor_Mean = mean_prof_cor_cor,
    stringsAsFactors = FALSE
  )
  
  # Append to list
  all_sites_metrics_list[[site]] <- site_results
  
  # --- Profiles ---
  # p_avg_org <- plot_average_profiles(metrics_org, 
  #                                    site_name = site, 
  #                                    stratify_col = "mean_h_max_p95")
  # p_avg_cor <- plot_average_profiles(metrics_cor, 
  #                                    site_name = site, 
  #                                    stratify_col = "mean_h_max_p95")
  # p_strat_org <- plot_stratified_profiles(metrics_org, 
  #                                         site_name = site, 
  #                                         n_per_class = 3, 
  #                                         stratify_col = "mean_h_max_p95")
  # p_strat_cor <- plot_stratified_profiles(metrics_cor, 
  #                                         site_name = site,
  #                                         n_per_class = 3,
  #                                         stratify_col = "mean_h_max_p95")
  # plot(p_avg_org)
  # plot(p_avg_cor)
  # plot(p_strat_org)
  # plot(p_strat_cor)
  df_plot <- data.frame(site = site, mean_h_max_p95 = df_org$mean_h_max_p95)
  
  # 1. GEDI Columns (keep as is or rename to ensure clarity)
  df_plot <- bind_cols(df_plot, df_org %>% dplyr::select(starts_with("PAI_")))
  
  # 2. ALS Original Columns (Rename to 'ALS_Org_')
  als_org_cols <- df_org %>% dplyr::select(starts_with("mean_als_pai_")) %>% 
    rename_with(~gsub("mean_als_pai", "ALS_Org_PAI", .))
  
  # 3. ALS Corrected Columns (Rename to 'ALS_Cor_')
  als_cor_cols <- df_cor %>% dplyr::select(starts_with("mean_als_pai_")) %>% 
    rename_with(~gsub("mean_als_pai", "ALS_Cor_PAI", .))
  
  # Bind all together
  df_plot <- bind_cols(df_plot, als_org_cols, als_cor_cols)
  
  # --- [PLOT PROFILES] ---
  p_avg <- plot_average_profiles2(df_plot, 
                                  site_name = site, 
                                  stratify_col = "mean_h_max_p95")
  print(p_avg)
}

final_metrics_df <- dplyr::bind_rows(all_sites_metrics_list)
row.names(final_metrics_df) <- NULL

print(final_metrics_df)













summer_data <- gedi_als_data %>%
  sf::st_drop_geometry() %>%            # <--- THIS FIXES YOUR ERROR
  mutate(month = month(date)) %>%
  filter(month >= 5 & month <= 9) %>%   # Filter Summer
  mutate(shot_id = row_number())        # Create ID

# Reshape GEDI (Target)
gedi_long <- summer_data %>%
  select(shot_id, date, lai_s2_org, rh95, starts_with("pai_")) %>%
  pivot_longer(
    cols = starts_with("pai_"), 
    names_to = "layer_name", 
    values_to = "gedi_lad_raw"
  ) %>%
  mutate(height_bin = as.numeric(gsub("pai_([0-9]+)_.*", "\\1", layer_name)))

# Reshape ALS (Static Predictor)
als_long <- summer_data %>%
  select(shot_id, starts_with("LAD_")) %>%
  pivot_longer(
    cols = starts_with("LAD_"), 
    names_to = "als_layer_name", 
    values_to = "als_lad_static"
  ) %>%
  mutate(height_bin = as.numeric(gsub("LAD_([0-9]+)_.*", "\\1", als_layer_name)))

# Combine (Now they are just regular dataframes, so left_join works)
model_data <- left_join(gedi_long, als_long, by = c("shot_id", "height_bin"))


# 2. CALIBRATE (Simple Ratio)
# ---------------------------
# Since we are in summer, we assume the bias is consistent.
# We calculate one correction factor per height bin.

calibration_factors <- model_data %>%
  group_by(height_bin) %>%
  summarise(
    # K = Mean ALS / Mean GEDI
    # Adding small epsilon (0.0001) to avoid division by zero if GEDI is 0
    k_factor = mean(als_lad_static, na.rm=TRUE) / (mean(gedi_lad_raw, na.rm=TRUE) + 0.0001)
  )

print("Calibration Factors per Height Bin:")
print(calibration_factors)

# Apply the correction
model_data <- model_data %>%
  left_join(calibration_factors, by = "height_bin") %>%
  mutate(
    gedi_lad_corrected = gedi_lad_raw * k_factor,
    # Safety: If the correction is extreme (e.g. > 5x), cap it or trust ALS
    gedi_lad_corrected = ifelse(k_factor > 5, als_lad_static, gedi_lad_corrected)
  )

# 3. TRAIN RANDOM FOREST (Summer Model)
# -------------------------------------
# Even in summer, S2 LAI varies (e.g., 3.5 in May vs 5.0 in July vs 4.0 in Sept).
# The RF will learn how this S2 fluctuation affects the profile density.

rf_summer <- ranger(
  formula = gedi_lad_corrected ~ lai_s2_org + als_lad_static + height_bin + rh95,
  data = model_data,
  num.trees = 500,
  importance = "impurity",
  seed = 42
)

print(rf_summer)
print(rf_summer$variable.importance)

# 4. PREDICTION TEST (e.g., A dry August day)
# -------------------------------------------
new_pixel <- data.frame(
  lai_s2_org = 4.2,      # S2 LAI in August
  als_lad_static = 0.6,  # The "potential" density at this height
  height_bin = 20,       # 20-25m layer
  rh95 = 28              # Tree height
)

pred <- predict(rf_summer, data = new_pixel)
print(paste("Predicted Summer LAD:", pred$predictions))

