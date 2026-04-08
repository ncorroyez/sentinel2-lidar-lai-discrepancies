# ---
# title: "2b.GEDI_zonal_temporal_LAI_analysis.R"
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
# sites <- c("Mormal")  # uncomment to test only one site

all_sites_data <- list()
bin_breaks <- seq(0, 40, by = 5)
gedi_radius = 12.5

# --------------------------------------------------------------
# MAIN PROCESSING LOOP
# --------------------------------------------------------------
for (site in sites) {
  message("Processing site: ", site)
  results_dir <- file.path("../../03_RESULTS", site, "Metrics", metrics_dir)
  
  # A. Load GEDI
  gedi_als_data <- readRDS(file.path(data_dir, site, "GEDI", 
                                     paste0(site, "_GEDI_LiDAR_Metrics.rds")))
  df <- gedi_als_data %>%
    dplyr::select(x_utm, y_utm, x_corrected, y_corrected, date,
                  rh95, rh98, rh100,
                  lai_gedi, pai_0_5m, pai_5_10m, pai_10_15m, pai_15_20m,
                  pai_20_25m, pai_25_30m, pai_30_35m, pai_35_40m,
                  LAI_lidar_org, LAD_0_5_org, LAD_5_10_org, LAD_10_15_org,
                  LAD_15_20_org, LAD_20_25_org, LAD_25_30_org, LAD_30_35_org,
                  LAD_35_40_org, LAI_lidar_cor, LAD_0_5_cor, LAD_5_10_cor,
                  LAD_10_15_cor, LAD_15_20_cor, LAD_20_25_cor, LAD_25_30_cor,
                  LAD_30_35_cor, LAD_35_40_cor) %>%
    mutate(doy = yday(date),
           sin_doy = sin(2 * pi * doy / 365),
           cos_doy = cos(2 * pi * doy / 365)
    )
  stop()
  # B. Prepare ALS Rasters
  als_stack <- prepare_als_rasters(site, results_dir, bin_breaks)
  
  # C. Extract Sentinel-2 Time Series
  s2_path <- file.path(results_dir, "stack", paste0(site, "_lai_stacked_raster.tif"))
  gedi_sf$mean_lai_s2 <- extract_s2_temporal(gedi_sf, gedi_buf, s2_path)
  
  # D. Extract ALS Metrics
  zonal_df <- extract_als_metrics(gedi_buf, als_stack)
  zonal_df$check <- zonal_df$mean_lai_als - rowSums(
    zonal_df[, grep("mean_als_pai_", 
                    names(zonal_df), value = TRUE)], na.rm = T)
  
  # E. Combine
  site_data <- cbind(
    st_coordinates(gedi_sf),
    st_drop_geometry(gedi_sf),
    zonal_df
  )
  site_data$site <- site
  all_sites_data[[site]] <- site_data[complete.cases(site_data), ]
}

zonal_als <- exactextractr::exact_extract(als_stack, gedi_buf)
first_fp_data <- zonal_als[[1]]
manual_mean <- sum(first_fp_data$lai_als * first_fp_data$coverage_fraction) / 
  sum(first_fp_data$coverage_fraction)
# -----------------------------
# Combine all sites
combined_data <- bind_rows(all_sites_data)
stop()
months <- c("05", "06", "07", "08", "09")
# months <- c("01", "02", "03", "04", "05", "06",
#             "07", "08", "09", "10", "11", "12")
# months <- c("06", "07", "08", "09")
# months <- c("11", "12", "01", "02", "03", "04")
# stop()
# 1️⃣ Filter for May–September
df_filtered <- combined_data %>%
  filter(format(date, "%m") %in% months) %>%
  mutate(
    als_fraction_sub = mean_als_pai_0_5m / mean_lai_als,
    als_fraction_can = (mean_als_pai_5_10m + mean_als_pai_10_15m + mean_als_pai_15_20m +
                          mean_als_pai_20_25m + mean_als_pai_25_30m +
                          mean_als_pai_30_35m + mean_als_pai_35_40m) / mean_lai_als,
    als_fraction_sub = ifelse(is.na(als_fraction_sub), 0, als_fraction_sub),
    als_fraction_can = ifelse(is.na(als_fraction_can), 0, als_fraction_can),
    doy = yday(date),
    sin_doy = sin(2 * pi * doy / 365),
    cos_doy = cos(2 * pi * doy / 365)
  )
# df_filtered <- combined_data %>%
#   filter(format(date, "%m") %in% months) %>%
#   filter(site == "Mormal")

# Add Shape Metrics (Correlation, EMD)
df_advanced <- calculate_shape_metrics_vectorized(df_filtered)

# 2️⃣ For each site, compute correlations and simple linear models
site_stats <- df_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    cor_lai_gedi_lai_als = cor(lai_gedi, mean_lai_als, use = "complete.obs"),
    cor_lai_gedi_lai_s2  = cor(lai_gedi, mean_lai_s2, use = "complete.obs"),
    slope_lai_als = coef(lm(lai_gedi ~ mean_lai_als, na.action = na.omit))[2],
    slope_lai_s2  = coef(lm(lai_gedi ~ mean_lai_s2, na.action = na.omit))[2],
    intercept_lai_als = coef(lm(lai_gedi ~ mean_lai_als, na.action = na.omit))[1],
    intercept_lai_s2  = coef(lm(lai_gedi ~ mean_lai_s2, na.action = na.omit))[1],
    cor_rh100_hmax = cor(RH100, mean_h_max_p100, use = "complete.obs"),
    rmse_lai_als          = Metrics::rmse(lai_gedi, mean_lai_als),
    rmse_lai_s2           = Metrics::rmse(lai_gedi, mean_lai_s2),
    rmse_hmax             = Metrics::rmse(RH100, mean_h_max_p100)
  )

print(site_stats)

rh_stats <- df_filtered %>%
  group_by(site) %>%
  summarise(
    cor_rh100_hmax        = cor(RH100, mean_h_max_p100, use = "complete.obs"),
    rmse_rh100_hmax       = Metrics::rmse(RH100, mean_h_max_p100),
    cor_rh95_hmax         = cor(RH95, mean_h_max_p95, use = "complete.obs"),
    rmse_rh95_hmax        = Metrics::rmse(RH95, mean_h_max_p95),
    cor_rh98_hmax         = cor(RH98, mean_h_max_p98, use = "complete.obs"),
    rmse_rh98_hmax        = Metrics::rmse(RH98, mean_h_max_p98)
  )

# Print Correlations
site_stats <- df_advanced %>%
  group_by(site) %>%
  summarise(
    n = n(),
    r_GEDI_ALS = cor(lai_gedi, mean_lai_als, use = "complete.obs"),
    r_GEDI_S2  = cor(lai_gedi, mean_lai_s2, use = "complete.obs"),
    mean_prof_cor = mean(prof_cor, na.rm=T),
    med_prof_cor = median(prof_cor, na.rm=T),
    mean_emd = mean(emd_dist, na.rm=T)
  )
print(site_stats)

# --------------------------------------------------------------
# 4) MODELING (Random Forest)
# --------------------------------------------------------------
rf_results <- run_random_forest(df_advanced)

cat("--- Random Forest Results ---\n")
print(rf_results$metrics)
varImpPlot(rf_results$model)

# Prediction Plot
ggplot(rf_results$test_data, aes(x = lai_gedi, y = Predicted, color = site)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_bw() +
  labs(title = "RF Predictions vs Observed GEDI LAI")

for (site in unique(df_filtered$site)) {
  p <- plot_stratified_profiles(
    data = df_filtered, 
    site_name = site, 
    n_per_class = 3, 
    stratify_col = "mean_h_max_p95" 
  )
  print(p)
  
  p_avg <- plot_average_profiles(
    data = df_filtered, 
    site_name = site, 
    stratify_col = "mean_h_max_p95"
  )
  print(p_avg)
}

stop()



rf_canopee <- randomForest(
  y = df_advanced$lai_gedi,
  x = df_advanced %>% 
    dplyr::select(
      mean_lai_s2,
      mean_lai_als,
      als_fraction_can,
      sin_doy, cos_doy
    )
)






# ix <- 500005
# stack <- als_stack[ix]; stack[is.na(stack)] <- 0
# lai_als <- stack$lai_als
# lai_sum <- round(stack$als_pai_0_5m, 2) +
#   round(stack$als_pai_5_10m, 2) + 
#   round(stack$als_pai_10_15m, 2) + 
#   round(stack$als_pai_15_20m, 2) + 
#   round(stack$als_pai_20_25m, 2) + 
#   round(stack$als_pai_25_30m, 2) + 
#   round(stack$als_pai_30_35m, 2) + 
#   round(stack$als_pai_35_40m, 2)
# 
# ladstack <- rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/Metrics/Deciduous_Only/ladstack_classic.tif")
# lidarlai <- rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/Metrics/Deciduous_Only/lidarlai_res_10_m.tif")
# 
# ladstack_pix <- ladstack[ix]
# lidarlai_pix <- lidarlai[ix]
# 
# lidarlai_pix_tot <- round(lidarlai_pix$V1, 2)
# ladstack_pix_tot <- sum(round(ladstack_pix, 2), na.rm = T)




# Standard Stats
site_stats <- df_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    cor_lai_gedi_lai_als = cor(lai_gedi, mean_lai_als, use = "complete.obs"),
    rmse_lai_als         = Metrics::rmse(lai_gedi, mean_lai_als)
  )

print("--- Standard Metrics ---")
print(site_stats)

# Profile Stats
message("\n--- Vertical Profile Comparison ---")
profile_stats_list <- list()

for (i in 1:(length(bin_breaks)-1)) {
  lower <- bin_breaks[i]
  upper <- bin_breaks[i+1]
  
  gedi_col <- paste0("PAI_", lower, "_", upper, "m")
  als_col  <- paste0("mean_als_pai_", lower, "_", upper, "m")
  
  if(gedi_col %in% names(df_filtered) && als_col %in% names(df_filtered)){
    bin_stat <- df_filtered %>%
      group_by(site) %>%
      summarise(
        Height_Bin = paste0(lower, "-", upper, "m"),
        Cor = cor(.data[[gedi_col]], .data[[als_col]], use = "complete.obs"),
        RMSE = Metrics::rmse(.data[[gedi_col]], .data[[als_col]])
      )
    profile_stats_list[[paste0(lower)]] <- bin_stat
  }
}

profile_stats_df <- bind_rows(profile_stats_list) %>%
  arrange(site, Height_Bin)

print(profile_stats_df)


set.seed(123)
samples <- df_filtered %>%
  ungroup() %>%
  sample_n(min(9, n())) %>%  # Pick 9 samples (or fewer if data is small)
  mutate(sample_id = paste0(site, "_ID", row_number())) 

# 2. Reshape to long format
# We extract columns starting with "PAI_" (GEDI) and "als_pai_" (ALS)
long_df <- samples %>%
  dplyr::select(sample_id, site, matches("^PAI_\\d+_\\d+m$"), 
                matches("mean_als_pai_\\d+_\\d+m$")) %>%
  pivot_longer(
    cols = -c(sample_id, site),
    names_to = "original_col",
    values_to = "PAI_Value"
  ) %>%
  mutate(
    # Identify Sensor
    Sensor = ifelse(grepl("^mean_als_", original_col), "ALS", "GEDI"),
    
    # Extract Height Bin string (e.g., "0_5")
    # Remove prefixes "PAI_" or "als_pai_" and suffix "m"
    bin_str = original_col,
    bin_str = str_remove(bin_str, "mean_als_pai_"), # Remove ALS prefix
    bin_str = str_remove(bin_str, "^PAI_"),         # Remove GEDI prefix
    bin_str = str_remove(bin_str, "m$")             # Remove suffix "m"
  ) %>%
  # Split "0_5" into numeric min/max
  separate(bin_str, into = c("h_min", "h_max"), sep = "_", convert = TRUE) %>%
  mutate(h_mid = (h_min + h_max) / 2) # Calculate midpoint for plotting

# 3. Plot
p_profiles <- ggplot(long_df, aes(x = PAI_Value, y = h_mid, color = Sensor)) +
  # Draw lines connecting the points
  geom_path(aes(group = Sensor), linewidth = 1) +
  # Add points at the bin centers
  geom_point(size = 2) +
  # Facet by sample to see individual comparisons
  facet_wrap(~ sample_id, scales = "free_x") +
  labs(
    title = "Vertical PAI Profiles: GEDI vs ALS",
    subtitle = "Comparison of 5m vertical bins for random footprints",
    x = "PAI (m²/m²)",
    y = "Height (m)",
    color = "Sensor"
  ) +
  scale_y_continuous(breaks = seq(0, 40, 5)) + # Ticks every 5m
  theme_bw() +
  theme(legend.position = "bottom")

print(p_profiles)


# --------------------------------------------------------------
# 10) Advanced Profile Shape Comparison
# --------------------------------------------------------------

# Define the list of profile columns for both sensors
# (Adjust names if your dataframe uses different indices)
gedi_cols <- paste0("PAI_", seq(0, 35, 5), "_", seq(5, 40, 5), "m")
als_cols  <- paste0("mean_als_pai_", seq(0, 35, 5), "_", seq(5, 40, 5), "m")

# Check if columns exist
if (all(c(gedi_cols, als_cols) %in% names(df_filtered))) {
  
  # Function to calculate shape metrics for a single row
  calc_shape_metrics <- function(idx, data) {
    
    # Extract vectors
    y_gedi <- as.numeric(data[idx, gedi_cols])
    y_als  <- as.numeric(data[idx, als_cols])
    
    # Handle empty profiles (sum = 0) to avoid division by zero
    if (sum(y_gedi) == 0 || sum(y_als) == 0) {
      return(data.frame(
        prof_cor = NA, 
        emd_dist = NA, 
        delta_cog = NA,
        total_pai_bias = sum(y_gedi) - sum(y_als)
      ))
    }
    
    # 1. Profile Correlation (Shape similarity)
    # ----------------------------------------
    prof_cor <- cor(y_gedi, y_als, method = "pearson")
    
    # 2. Earth Mover's Distance (EMD) via CDFs
    # ----------------------------------------
    # Normalize to create Probability Density Function (PDF)
    pdf_gedi <- y_gedi / sum(y_gedi)
    pdf_als  <- y_als  / sum(y_als)
    
    # Calculate Cumulative Distribution Function (CDF)
    cdf_gedi <- cumsum(pdf_gedi)
    cdf_als  <- cumsum(pdf_als)
    
    # EMD in 1D is the Sum of absolute differences between CDFs 
    # (multiplied by bin width, here 5m)
    emd_dist <- sum(abs(cdf_gedi - cdf_als)) * 5
    
    # 3. Center of Gravity (Weighted Mean Height)
    # ----------------------------------------
    # Midpoints of bins: 2.5, 7.5, ..., 37.5
    mid_heights <- seq(2.5, 37.5, by = 5)
    
    cog_gedi <- sum(y_gedi * mid_heights) / sum(y_gedi)
    cog_als  <- sum(y_als * mid_heights) / sum(y_als)
    
    return(data.frame(
      prof_cor = prof_cor,
      emd_dist = emd_dist,      # Lower is better match
      delta_cog = cog_gedi - cog_als, # Positive = GEDI looks "higher"
      total_pai_bias = sum(y_gedi) - sum(y_als)
    ))
  }
  
  message("Calculating shape metrics (Correlation, EMD, CoG)...")
  
  # Apply function to all rows (using lapply for speed, then bind)
  shape_metrics_list <- lapply(1:nrow(df_filtered), calc_shape_metrics, data = df_filtered)
  shape_metrics_df   <- do.call(rbind, shape_metrics_list)
  
  # Combine back to main dataframe
  df_advanced <- bind_cols(df_filtered, shape_metrics_df)
  
  # --------------------------------------------------------------
  # 11) Analyze Shape Results
  # --------------------------------------------------------------
  
  shape_summary <- df_advanced %>%
    group_by(site) %>%
    summarise(
      n = n(),
      # Median Profile Correlation (how well do shapes match?)
      median_profile_cor = median(prof_cor, na.rm = TRUE),
      
      # Mean EMD (Average "distance" between shapes in meters)
      mean_emd = mean(emd_dist, na.rm = TRUE),
      
      # Bias in Center of Gravity
      mean_delta_cog = mean(delta_cog, na.rm = TRUE),
      
      # Total PAI Bias
      mean_pai_bias = mean(total_pai_bias, na.rm = TRUE)
    )
  
  print("--- Advanced Profile Shape Statistics ---")
  print(shape_summary)
  
  # --------------------------------------------------------------
  # 12) Visualization: Correlation vs EMD
  # --------------------------------------------------------------
  # This plot helps you see: 
  # Top Left: Good shape match, low shift.
  # Bottom Right: Terrible match.
  
  p_shape <- ggplot(df_advanced, aes(x = prof_cor, y = emd_dist)) +
    geom_point(alpha = 0.3, aes(color = site)) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed") +
    labs(
      title = "Profile Similarity Analysis",
      subtitle = "Profile Correlation (Shape) vs EMD (Vertical Shift magnitude)",
      x = "Profile Correlation (Pearson r)",
      y = "Earth Mover's Distance (EMD) [meters]"
    ) +
    theme_bw()
  
  print(p_shape)
}










df_ml <- df_filtered %>%
  mutate(
    doy = yday(date),
    sin_doy = sin(2 * pi * doy / 365),
    cos_doy = cos(2 * pi * doy / 365)
  ) %>%
  dplyr::select(
    lai_gedi, 
    mean_lai_s2,
    mean_lai_als, 
    mean_h_max_p98, 
    mean_dsm_sd,
    sin_doy, 
    cos_doy,
    site
  ) %>%
  na.omit()

# n_min <- df_ml %>% 
#   count(site) %>% 
#   pull(n) %>% 
#   min()
# 
# df_ml <- df_ml %>%
#   group_by(site) %>%
#   slice_sample(n = n_min) %>%   # equal samples per site
#   ungroup()

train_index <- createDataPartition(df_ml$lai_gedi, p = 0.7, list = FALSE)
train_data <- df_ml[train_index, ]
test_data  <- df_ml[-train_index, ]

rf_model <- randomForest(
  lai_gedi ~ . - site, 
  data = train_data, 
  ntree = 500, 
  mtry = 2,
  importance = TRUE
)

print(rf_model)

# Predictions
predictions <- predict(rf_model, test_data)
rf_rmse <- Metrics::rmse(test_data$lai_gedi, predictions)
rf_r    <- cor(test_data$lai_gedi, predictions)
rf_r2   <- broom::glance(lm(predictions ~ test_data$lai_gedi, test_data))$r.squared
rf_bias <- mean(predictions - test_data$lai_gedi)

cat("\n--- RÉSULTATS DU MODÈLE SUR LE JEU DE TEST ---\n")
cat("RMSE :", round(rf_rmse, 3), "\n")
cat("R²   :", round(rf_r2, 3), "\n")
cat("Biais:", round(rf_bias, 3), "\n")


plot_data <- data.frame(
  Observed = test_data$lai_gedi,
  Predicted = predictions,
  mean_lai_s2 = test_data$mean_lai_s2,
  Site = test_data$site
)

p_rf <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = Site), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = paste0("Random Forest: GEDI PAI Prediction (Test Set)"),
    subtitle = paste0("RMSE = ", round(rf_rmse, 2), " | R2 = ", round(rf_r2, 2)),
    x = "GEDI PAI (Observé)",
    y = "RF PAI (Prédit)"
  ) +
  theme_bw() +
  xlim(0, 10) + ylim(0, 10)

print(p_rf)

varImpPlot(rf_model, main = "Importance des variables (IncNodePurity)")
importance_df <- as.data.frame(importance(rf_model))
importance_df$Variable <- rownames(importance_df)

ggplot(importance_df, aes(x = reorder(Variable, IncNodePurity), y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Importance des Variables", x = "", y = "Gain de Pureté (Importance)") +
  theme_minimal()


ggplot(plot_data, aes(x = mean_lai_s2, y = Predicted)) +
  geom_point(aes(color = Site), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = paste0("GEDI PAI Prediction vs S2 LAI (Test Set)"),
    x = "S2 LAI (Observed)",
    y = "RF PAI (Predicted)"
  ) +
  theme_bw() +
  xlim(0, 10) + ylim(0, 10)



stop()













# Plot GEDI vs ALS and GEDI vs S2 per site
ggplot(df_filtered, aes(x = mean_lai_als, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~site) +
  labs(
    x = "ALS LAI",
    y = "GEDI LAI",
    title = "GEDI vs ALS LAI (May–September)"
  ) +
  xlim(c(0, 15)) + ylim(c(0, 15)) +
  theme_bw()

ggplot(df_filtered, aes(x = mean_lai_s2, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "Sentinel-2 LAI",
    y = "GEDI LAI",
    title = "GEDI vs Sentinel-2 LAI (May–September)"
  ) +
  xlim(c(0, 15)) + ylim(c(0, 15)) +
  theme_bw()

ggplot(df_filtered, aes(x = mean_h_max_p98, y = RH98)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "ALS Hmax98",
    y = "GEDI RH98",
    title = "GEDI RH98 vs ALS Hmax98 (May–September)"
  ) +
  theme_bw()

# Create classes of mean_h_max
df_filtered_classes <- df_filtered %>%
  mutate(h_class = cut(
    mean_h_max_p95,
    breaks = c(0, 10, 20, 30, 40, Inf),
    labels = c("0–10 m", "10–20 m", "20–30 m", "30–40 m", ">40 m"),
    right = FALSE
  ))

cor_stats_hclass <- df_filtered_classes %>%
  group_by(site, h_class) %>%
  summarise(
    n = n(),
    cor_lai_s2_gedi = cor(mean_lai_s2, lai_gedi, use = "complete.obs")
  ) %>%
  ungroup()

print(cor_stats_hclass)

ggplot(df_filtered_classes, aes(x = mean_lai_s2, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_grid(h_class ~ site) +
  labs(
    x = "Sentinel-2 LAI",
    y = "GEDI LAI",
    title = "GEDI vs Sentinel-2 LAI (May–September)"
  ) +
  theme_bw()
