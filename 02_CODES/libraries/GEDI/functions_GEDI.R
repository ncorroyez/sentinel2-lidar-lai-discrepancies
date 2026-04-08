# ---
# title: "functions_GEDI"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-11-13"
# ---

load_gedi_data <- function(site, 
                           data_dir,
                           results_dir,
                           bin_breaks = seq(0, 40, by = 5),
                           gedi_radius = 12.5) {
  # Construct path
  csv_path <- file.path(data_dir, site, "Geo_Files", paste0("gedi_shots_", site, ".csv"))
  als_path <- file.path(results_dir, "lidarlai_res_10_m.tif")
  
  if (!file.exists(csv_path)) stop("GEDI file not found: ", csv_path)
  if (!file.exists(als_path)) stop("ALS LAI file not found: ", als_path)
  
  gedi <- read.csv(csv_path)
  lai_als <- terra::rast(als_path)
  
  # formatting
  gedi$month <- sprintf("%02d", gedi$month)
  gedi$day   <- sprintf("%02d", gedi$day)
  
  # Spatial conversion
  gedi_sf <- st_as_sf(gedi, coords = c("X_utm", "Y_utm"), crs = 32631) %>%
    mutate(date = make_date(year, as.integer(month), as.integer(day)),
           lai_gedi = pai,
           month_year = format(date, "%Y-%m"))
  
  # Remove non-deciduous pixels
  als_values <- terra::extract(lai_als, gedi_sf, ID = FALSE)
  valid_indices <- !is.na(als_values[, 1])
  n_removed <- sum(!valid_indices)
  if (n_removed > 0) {
    message(paste("Removing", n_removed, "/", nrow(gedi_sf),
                  "GEDI shots falling on NA values in ALS LAI raster."))
  }
  gedi_sf <- gedi_sf[valid_indices, ]
  
  # Select columns
  pai_cols <- grep("^PAI_\\d+_\\d+m$", names(gedi), value = TRUE)
  gedi_sf <- gedi_sf %>%
    dplyr::select(lai_gedi, date, month_year, RH100, RH95, RH98, all_of(pai_cols))
  
  # Create Buffer
  gedi_buf <- st_buffer(gedi_sf, dist = gedi_radius)
  
  return(list(sf = gedi_sf, buffer = gedi_buf))
}

prepare_als_rasters <- function(site, results_dir, bin_breaks) {
  # 1. Load Standard Metrics
  r_list <- list(
    lai_als    = rast(file.path(results_dir, "lidarlai_res_10_m.tif")),
    h_max_p100 = rast(file.path(results_dir, "max_res_10_m.tif")),
    h_max_p95  = rast(file.path(results_dir, "max_p95_res_10_m.tif")),
    h_max_p98  = rast(file.path(results_dir, "max_p98_res_10_m.tif")),
    dsm_sd     = rast(file.path(results_dir, "dsm_sd_res_10_m.tif"))
  )
  
  # Standardize names
  names(r_list) <- names(r_list) # ensures raster names match list keys
  
  # 2. Process LAD Stack (The complex binning logic)
  ladstack <- rast(file.path(results_dir, "ladstack_classic.tif"))
  als_pai_list <- list()
  
  layer_names <- names(ladstack)
  # Extract height from "LAD_Layer_X.5" -> returns X
  layer_heights <- as.numeric(stringr::str_extract(layer_names, "(?<=Layer_)\\d+"))
  # layer_heights <- as.numeric(stringr::str_extract(layer_names, 
  #                                                  "(?<=Layer_)[0-9.]+"))
  # print(layer_heights)
  
  for(i in 1:(length(bin_breaks)-1)) {
    lower <- bin_breaks[i]
    upper <- bin_breaks[i+1]
    bin_name <- paste0("als_pai_", lower, "_", upper, "m")
    
    # Find layers belonging to this bin
    indices <- which(layer_heights >= lower & layer_heights < upper)
    # print(which(layer_heights >= lower & layer_heights < upper))
    # print(layer_heights >= lower & layer_heights < upper)
    
    if(length(indices) > 0) {
      # als_pai_list[[bin_name]] <- sum(ladstack[[indices]], na.rm = T)
      als_pai_list[[bin_name]] <- mean(ladstack[[indices]], na.rm = T)
    } else {
      # Handle empty bins
      if(terra::nlyr(ladstack) > 0) {
        als_pai_list[[bin_name]] <- ladstack[[1]] * 0
      }
    }
  }
  
  # Combine everything into one stack
  pai_stack <- rast(als_pai_list)
  # pai_stack <- terra::classify(pai_stack, cbind(NA, 0))
  std_stack <- rast(r_list)
  
  valid_mask <- !is.na(std_stack$lai_als) # & !is.na(ladstack[[1]])
  pai_stack_filled <- terra::classify(pai_stack, cbind(NA, 0))
  pai_stack_aligned <- terra::mask(pai_stack_filled, valid_mask)
  final_stack <- c(std_stack, pai_stack_aligned)
  # final_stack <- c(std_stack, pai_stack)
  # final_stack <- terra::classify(final_stack, cbind(NA, 0))
  
  return(final_stack)
}

load_gedi_and_geogedi_data <- function(site, data_dir, results_dir, bin_breaks, gedi_radius) {
  
  csv_path <- file.path(data_dir, site, "Geo_Files", paste0("gedi_joined_", site, ".csv"))
  als_path <- file.path(results_dir, "lidarlai_res_10_m.tif")
  
  if (!file.exists(csv_path)) stop("GEDI file not found: ", csv_path)
  # if (!file.exists(als_path)) stop("ALS LAI file not found: ", als_path) # Optional check
  
  gedi <- read.csv(csv_path)
  
  # Formatting dates
  gedi$month <- sprintf("%02d", gedi$month)
  gedi$day   <- sprintf("%02d", gedi$day)
  
  # ---------------------------------------------------------
  # KEY UPDATE: remove = FALSE to keep X_utm/Y_utm columns
  # ---------------------------------------------------------
  gedi_sf <- st_as_sf(gedi, coords = c("x_utm", "y_utm"), crs = 32631, remove = FALSE) %>%
    mutate(date = make_date(year, as.integer(month), as.integer(day)),
           lai_gedi = pai,
           month_year = format(date, "%Y-%m"))
  
  # Optional: Filter based on ALS overlap (using ORIGINAL coords)
  if (file.exists(als_path)) {
    lai_als <- terra::rast(als_path)
    als_values <- terra::extract(lai_als, gedi_sf, ID = FALSE)
    valid_indices <- !is.na(als_values[, 1])
    gedi_sf <- gedi_sf[valid_indices, ]
  }
  
  # ---------------------------------------------------------
  # KEY UPDATE: Select x_corrected and y_corrected
  # ---------------------------------------------------------
  # pai_cols <- grep("^PAI_\\d+_\\d+m$", names(gedi), value = TRUE)
  
  # Ensure columns exist before selecting to avoid crash
  # cols_to_keep <- c("lai_gedi", "date", "month_year", "rh100", "rh95", "rh98", 
  #                   "X_utm", "Y_utm", "x_corrected", "y_corrected")
  
  # Add PAI cols
  # cols_to_keep <- c(cols_to_keep, pai_cols)
  
  # Only select columns that actually exist in the dataframe
  # cols_to_keep <- intersect(cols_to_keep, names(gedi_sf))
  
  # gedi_sf <- gedi_sf %>% dplyr::select(all_of(cols_to_keep))
  
  # Create Buffer (Geometry is based on Original Coords)
  gedi_buf <- st_buffer(gedi_sf, dist = gedi_radius)
  
  return(list(sf = gedi_sf, buffer = gedi_buf))
}

extract_s2_temporal <- function(gedi_sf, gedi_buf, s2_stack_path) {
  # Load S2 Stack
  s2_stack <- rast(s2_stack_path)
  layer_dates <- as.Date(substr(names(s2_stack), 1, 8), "%Y%m%d")
  names(s2_stack) <- as.character(layer_dates)
  s2_month_year <- format(layer_dates, "%Y-%m")
  
  # Container for results
  gedi_sf$mean_lai_s2 <- NA
  unique_months <- unique(gedi_sf$month_year)
  
  # Temporal loop
  for (m in unique_months) {
    idx_layers <- which(s2_month_year == m)
    if (length(idx_layers) == 0) next
    
    s2_stack_month <- s2_stack[[idx_layers]]
    gedi_idx <- which(gedi_sf$month_year == m)
    
    # Extract
    if(length(gedi_idx) > 0){
      gedi_buf_month <- gedi_buf[gedi_idx, ]
      s2_values <- exactextractr::exact_extract(s2_stack_month, gedi_buf_month, fun = "mean")
      
      # Handle multi-layer average if multiple S2 images exist for that month
      if(is.matrix(s2_values) || is.data.frame(s2_values)) {
        s2_mean <- rowMeans(s2_values, na.rm = TRUE)
      } else {
        s2_mean <- s2_values
      }
      gedi_sf$mean_lai_s2[gedi_idx] <- s2_mean
    }
  }
  return(gedi_sf$mean_lai_s2)
}

extract_s2_temporal2 <- function(gedi_geometry, gedi_dates, s2_stack_path) {
  
  # 1. Load S2 Stack and Parse Dates
  s2_stack <- rast(s2_stack_path)
  
  # Assuming layer names are like "20210801" or similar. Adjust substr index if needed.
  layer_dates <- as.Date(substr(names(s2_stack), 1, 8), "%Y%m%d")
  s2_month_year <- format(layer_dates, "%Y-%m") # e.g., "2021-08"
  
  # Prepare GEDI dates
  gedi_month_year <- format(as.Date(gedi_dates), "%Y-%m")
  unique_months <- unique(gedi_month_year)
  
  # Result container (vector same length as input)
  results <- rep(NA, length(gedi_geometry))
  
  message(paste("Extracting S2 Data for", length(unique_months), "unique month-year combinations..."))
  
  # 2. Temporal Loop
  for (m in unique_months) {
    
    # A. Find S2 layers for this specific Month-Year
    idx_s2 <- which(s2_month_year == m)
    
    if (length(idx_s2) == 0) {
      # warning(paste("No S2 data found for month:", m))
      next
    }
    
    # Subset Raster Stack for this month (could be 1 or multiple images)
    s2_subset <- s2_stack[[idx_s2]]
    
    # B. Find GEDI footprints for this specific Month-Year
    idx_gedi <- which(gedi_month_year == m)
    
    # C. Extract
    # subset geometry
    subset_geom <- gedi_geometry[idx_gedi]
    
    # exact_extract returns a matrix/df if multiple layers, or vector if single layer
    vals <- exactextractr::exact_extract(s2_subset, subset_geom, fun = "mean")
    
    # If multiple S2 images exist for this month (e.g., 2021-08-01 and 2021-08-15),
    # vals will be a data.frame. We average them to get one monthly value.
    if (is.data.frame(vals) || is.matrix(vals)) {
      if (ncol(vals) > 1) {
        monthly_mean <- rowMeans(vals, na.rm = TRUE)
      } else {
        monthly_mean <- vals[,1]
      }
    } else {
      monthly_mean <- vals
    }
    
    # Assign back to results vector
    results[idx_gedi] <- monthly_mean
  }
  
  return(results)
}

extract_als_metrics <- function(gedi_buf, als_stack) {
  # Perform zonal statistics
  # 
  zonal_als <- exactextractr::exact_extract(als_stack, gedi_buf, fun = "mean")
  zonal_df <- as.data.frame(zonal_als) %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
  
  # Clean names (remove exact_extract prefixes if any, handle dots)
  names(zonal_df) <- gsub("mean\\.", "mean_", names(zonal_df)) 
  names(zonal_df) <- gsub("\\.", "_", names(zonal_df))
  
  # Clean NAs in PAI profiles
  als_profile_cols <- grep("als_pai_", names(zonal_df), value = TRUE)
  if(length(als_profile_cols) > 0){
    zonal_df[als_profile_cols] <- lapply(zonal_df[als_profile_cols], function(x) {
      x[is.na(x) | is.nan(x)] <- 0
      return(x)
    })
  }
  return(zonal_df)
}

# ==============================================================================
# PART 3: ADVANCED METRICS & SHAPE ANALYSIS
# ==============================================================================

calculate_shape_metrics_vectorized <- function(data) {
  # Define columns
  gedi_cols <- paste0("PAI_", seq(0, 35, 5), "_", seq(5, 40, 5), "m")
  als_cols  <- paste0("mean_als_pai_", seq(0, 35, 5), "_", seq(5, 40, 5), "m")
  
  # Check existence
  if (!all(c(gedi_cols, als_cols) %in% names(data))) {
    warning("Columns for shape metrics missing.")
    return(data)
  }
  
  # Inner function for apply
  calc_row <- function(row_idx) {
    y_gedi <- as.numeric(data[row_idx, gedi_cols])
    y_als  <- as.numeric(data[row_idx, als_cols])
    lai_gedi <- data[row_idx, "lai_gedi"]
    lai_lidar <- data[row_idx, "lai_lidar"]
    
    if (sum(y_gedi, na.rm=T) == 0 || sum(y_als, na.rm=T) == 0) {
      return(c(NA, NA, NA, sum(y_gedi) - sum(y_als)))
    }
    
    # Correlation
    prof_cor <- cor(y_gedi, y_als, method = "pearson")
    
    # EMD (Earth Mover's Distance)
    pdf_gedi <- y_gedi / lai_gedi # sum(y_gedi)
    pdf_als  <- y_als  / lai_lidar # sum(y_als)
    emd_dist <- sum(abs(cumsum(pdf_gedi) - cumsum(pdf_als))) * 5 # 5m bin width
    
    # CoG
    mid_heights <- seq(2.5, 37.5, by = 5)
    cog_gedi <- sum(y_gedi * mid_heights) / lai_gedi
    cog_als  <- sum(y_als * mid_heights) / lai_lidar
    # cog_gedi <- lai_gedi * mid_heights / lai_gedi
    # cog_als  <- lai_lidar * mid_heights / lai_lidar
    
    return(c(prof_cor, emd_dist, cog_gedi - cog_als, 
             # sum(y_gedi) - sum(y_als) * 5)
             lai_gedi - lai_lidar
    )
    )
  }
  
  # Apply
  metrics_mat <- t(sapply(1:nrow(data), calc_row))
  colnames(metrics_mat) <- c("prof_cor", "emd_dist", "delta_cog", "total_pai_bias")
  
  return(cbind(data, as.data.frame(metrics_mat)))
}

# ==============================================================================
# PART 4: MODELING
# ==============================================================================

run_random_forest <- function(data) {
  # Feature Engineering
  df_ml <- data %>%
    mutate(doy = yday(date),
           sin_doy = sin(2 * pi * doy / 365),
           cos_doy = cos(2 * pi * doy / 365)) %>%
    dplyr::select(lai_gedi, mean_lai_s2, mean_lai_als, 
                  mean_h_max_p98, mean_dsm_sd, sin_doy, cos_doy, site) %>%
    na.omit()
  
  # Split
  set.seed(42)
  train_idx <- createDataPartition(df_ml$lai_gedi, p = 0.7, list = FALSE)
  train <- df_ml[train_idx, ]
  test  <- df_ml[-train_idx, ]
  
  # Train
  rf_model <- randomForest(lai_gedi ~ . - site, data = train, ntree = 500, mtry = 2, importance = TRUE)
  
  # Predict
  preds <- predict(rf_model, test)
  
  # Metrics
  metrics <- list(
    rmse = Metrics::rmse(test$lai_gedi, preds),
    r2   = broom::glance(lm(preds ~ test$lai_gedi))$r.squared,
    bias = mean(preds - test$lai_gedi)
  )
  
  return(list(model = rf_model, test_data = cbind(test, Predicted = preds), metrics = metrics))
}

plot_stratified_profiles <- function(data, 
                                     site_name, 
                                     n_per_class = 3, 
                                     stratify_col = "mean_h_max_p95") {
  
  # 1. Filter for specific site
  site_df <- data %>% filter(site == site_name)
  
  if(nrow(site_df) == 0) {
    message("No data found for site: ", site_name)
    return(NULL)
  }
  
  # 2. Stratify the data (Low, Medium, High based on ALS metric)
  # We use ntile to split data into 3 equal-sized groups
  site_df <- site_df %>%
    mutate(strata = ntile(.data[[stratify_col]], 3)) %>%
    mutate(strata_label = case_when(
      strata == 1 ~ "Low Canopy",
      strata == 2 ~ "Medium Canopy",
      strata == 3 ~ "High Canopy"
    ))
  
  # 3. Sample from each stratum
  set.seed(123)
  samples <- site_df %>%
    group_by(strata_label) %>%
    slice_sample(n = n_per_class) %>%
    ungroup() %>%
    # Create a clean ID for the plot title
    mutate(plot_id = paste0(strata_label, "\n(ID: ", row_number(), ")"))
  
  # 4. Reshape to Long Format (The Regex Magic)
  long_df <- samples %>%
    dplyr::select(plot_id, matches("^PAI_\\d+_\\d+m$"), 
                  matches("mean_als_pai_\\d+_\\d+m$")) %>%
    pivot_longer(
      cols = -plot_id,
      names_to = "original_col",
      values_to = "PAI_Value"
    ) %>%
    mutate(
      # Identify Sensor
      Sensor = ifelse(grepl("^mean_als_", original_col), "ALS", "GEDI"),
      
      # Extract numeric bounds
      bin_str = original_col,
      bin_str = str_remove(bin_str, "mean_als_pai_"),
      bin_str = str_remove(bin_str, "^PAI_"),
      bin_str = str_remove(bin_str, "m$")
    ) %>%
    separate(bin_str, into = c("h_min", "h_max"), sep = "_", convert = TRUE) %>%
    mutate(h_mid = (h_min + h_max) / 2)
  
  # 5. Plot
  p <- ggplot(long_df, aes(x = PAI_Value, y = h_mid, color = Sensor)) +
    # Add a fill for the area under curve (optional, looks nice)
    geom_ribbon(aes(xmin=0, xmax=PAI_Value, fill=Sensor), alpha=0.1, color=NA) + 
    geom_path(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_wrap(~ plot_id, scales = "free_x") + 
    scale_y_continuous(breaks = seq(0, 50, 5)) +
    scale_color_manual(values = c("ALS" = "black", "GEDI" = "#E69F00")) +
    scale_fill_manual(values = c("ALS" = "black", "GEDI" = "#E69F00")) +
    labs(
      title = paste0("Vertical Profiles: ", site_name),
      subtitle = paste0("Stratified by ", stratify_col, " (Low/Med/High)"),
      x = bquote("PAI"~(m^2/m^3)),
      y = "Height (m)"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}


plot_average_profiles <- function(data, 
                                  site_name, 
                                  stratify_col = "mean_h_max_p95") {
  
  # 1. Filter for Site
  site_df <- data %>% filter(site == site_name)
  if(nrow(site_df) == 0) return(NULL)
  
  # 2. Stratify (Low, Medium, High)
  # We divide into 3 quantiles based on height
  stratified_df <- site_df %>%
    mutate(strata_id = ntile(.data[[stratify_col]], 3)) %>%
    mutate(strata_label = case_when(
      strata_id == 1 ~ "Low Canopy",
      strata_id == 2 ~ "Medium Canopy",
      strata_id == 3 ~ "High Canopy"
    ))
  
  # 3. Create "All Data" category
  # We duplicate the data, label it "All Data", and bind it to the stratified data
  all_df <- stratified_df
  all_df$strata_label <- "All Data"
  
  combined_df <- bind_rows(stratified_df, all_df)
  
  # 4. Reshape to Long Format
  long_df <- combined_df %>%
    dplyr::select(strata_label, matches("^PAI_\\d+_\\d+m$"), 
                  matches("mean_als_pai_\\d+_\\d+m$")) %>%
    pivot_longer(
      cols = -strata_label,
      names_to = "original_col",
      values_to = "PAI_Value"
    ) %>%
    mutate(
      Sensor = ifelse(grepl("^mean_als_", original_col), "ALS", "GEDI"),
      
      # Clean string to get height
      bin_str = original_col,
      bin_str = str_remove(bin_str, "mean_als_pai_"),
      bin_str = str_remove(bin_str, "^PAI_"),
      bin_str = str_remove(bin_str, "m$")
    ) %>%
    separate(bin_str, into = c("h_min", "h_max"), sep = "_", convert = TRUE) %>%
    mutate(h_mid = (h_min + h_max) / 2)
  
  # 5. Aggregate: Calculate Mean and SD per bin/sensor/strata
  agg_df <- long_df %>%
    group_by(strata_label, Sensor, h_mid) %>%
    summarise(
      mean_pai = mean(PAI_Value, na.rm = TRUE),
      sd_pai   = sd(PAI_Value, na.rm = TRUE),
      n        = n(),
      .groups = "drop"
    ) %>%
    # Calculate Standard Error or confidence interval if preferred (here using Mean +/- 0.5 SD for visibility)
    mutate(
      ymin = pmax(0, mean_pai - 0.5 * sd_pai), # prevent negative PAI in plot
      ymax = mean_pai + 0.5 * sd_pai
    )
  
  # Order the factors so "All Data" appears last or first as preferred
  agg_df$strata_label <- factor(agg_df$strata_label, 
                                levels = c("Low Canopy", "Medium Canopy", "High Canopy", "All Data"))
  
  # 6. Plot
  p <- ggplot(agg_df, aes(x = mean_pai, y = h_mid, color = Sensor, fill = Sensor)) +
    # Add ribbon for variability (Standard Deviation)
    geom_ribbon(aes(xmin = ymin, xmax = ymax), alpha = 0.2, color = NA) +
    # Main average line
    geom_path(linewidth = 1.2) +
    geom_point(size = 1.5) +
    facet_wrap(~ strata_label, nrow = 1) +
    scale_y_continuous(breaks = seq(0, 50, 5)) +
    scale_color_manual(values = c("ALS" = "black", "GEDI" = "#D55E00")) +
    scale_fill_manual(values = c("ALS" = "grey30", "GEDI" = "#D55E00")) +
    labs(
      title = paste0("Average Vertical Profiles: ", site_name),
      subtitle = paste0("Mean Profile ± 0.5 SD | Stratified by ", stratify_col),
      x = bquote("Mean PAI"~(m^2/m^2)),
      y = "Height (m)"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

plot_average_profiles2 <- function(data, 
                                   site_name, 
                                   stratify_col = "mean_h_max_p95") {
  
  # 1. Filter for Site
  site_df <- data %>% filter(site == site_name)
  if(nrow(site_df) == 0) return(NULL)
  
  # 2. Stratify
  stratified_df <- site_df %>%
    mutate(strata_id = ntile(.data[[stratify_col]], 3)) %>%
    mutate(strata_label = case_when(
      strata_id == 1 ~ "Low Canopy",
      strata_id == 2 ~ "Medium Canopy",
      strata_id == 3 ~ "High Canopy"
    ))
  
  # 3. Create "All Data" category
  all_df <- stratified_df
  all_df$strata_label <- "All Data"
  combined_df <- bind_rows(stratified_df, all_df)
  
  # 4. Reshape to Long Format
  long_df <- combined_df %>%
    dplyr::select(strata_label, 
                  starts_with("PAI_"),          # GEDI
                  starts_with("ALS_Org_PAI_"),  # ALS Original
                  starts_with("ALS_Cor_PAI_"))  # ALS Corrected
  
  long_df <- long_df %>%
    pivot_longer(
      cols = -strata_label,
      names_to = "original_col",
      values_to = "PAI_Value"
    ) %>%
    mutate(
      # Identify Sensor/Type
      Sensor = case_when(
        grepl("^PAI_", original_col) ~ "GEDI",
        grepl("ALS_Org_", original_col) ~ "ALS (Original)",
        grepl("ALS_Cor_", original_col) ~ "ALS (Corrected)"
      ),
      
      # Clean string to get height
      bin_str = original_col,
      # Remove prefixes to isolate "0_5m", "5_10m" etc
      bin_str = gsub("ALS_Org_PAI_", "", bin_str),
      bin_str = gsub("ALS_Cor_PAI_", "", bin_str),
      bin_str = gsub("^PAI_", "", bin_str),
      bin_str = str_remove(bin_str, "m$")
    ) %>%
    separate(bin_str, into = c("h_min", "h_max"), sep = "_", convert = TRUE) %>%
    mutate(h_mid = (h_min + h_max) / 2)
  
  # 5. Aggregate
  agg_df <- long_df %>%
    group_by(strata_label, Sensor, h_mid) %>%
    summarise(
      mean_pai = mean(PAI_Value, na.rm = TRUE),
      sd_pai   = sd(PAI_Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(0, mean_pai - 0.5 * sd_pai),
      ymax = mean_pai + 0.5 * sd_pai
    )
  
  # Set Factor Levels
  agg_df$strata_label <- factor(agg_df$strata_label, 
                                levels = c("Low Canopy", "Medium Canopy", "High Canopy", "All Data"))
  
  agg_df$Sensor <- factor(agg_df$Sensor, 
                          levels = c("GEDI", "ALS (Original)", "ALS (Corrected)"))
  
  # 6. Plot
  p <- ggplot(agg_df, aes(x = mean_pai, y = h_mid, color = Sensor, fill = Sensor)) +
    geom_ribbon(aes(xmin = ymin, xmax = ymax), alpha = 0.15, color = NA) +
    geom_path(linewidth = 1.1) +
    geom_point(size = 1.2) +
    facet_wrap(~ strata_label, nrow = 1) +
    
    # Custom Colors: 
    # GEDI = Orange
    # ALS Org = Gray/Black (Reference 1)
    # ALS Cor = Blue (Reference 2 - "Improved")
    scale_color_manual(values = c("GEDI" = "#D55E00", 
                                  "ALS (Original)" = "gray40", 
                                  "ALS (Corrected)" = "blue")) +
    scale_fill_manual(values = c("GEDI" = "#D55E00", 
                                 "ALS (Original)" = "gray40", 
                                 "ALS (Corrected)" = "blue")) +
    
    labs(
      title = paste0("Average Vertical Profiles: ", site_name),
      subtitle = "GEDI vs ALS (Original Footprint) vs ALS (Corrected Footprint)",
      x = bquote("Mean PAI"~(m^2/m^2)),
      y = "Height (m)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

calculate_lad_metrics <- function(las, bin_size = 5, k = 0.5) {
  if (is.empty(las)) {
    print("LAS is NULL")
    return(NULL)
  }
  
  # Height Normalization (if needed)
  # las <- lidR::normalize_height(las, algorithm = tin())
  
  # Define Bins & Calculate LAD
  z_vals <- las$Z[las$Z >= 0]
  if(length(z_vals) < 10) return(NULL) 
  
  lad_df <- tryCatch({
    lidR::LAD(z_vals, dz = bin_size, k = k, z0 = 0)
  }, error = function(e) return(NULL))
  if (is.null(lad_df) || nrow(lad_df) == 0) return(NULL)
  
  max_height <- 40
  desired_centers <- seq(bin_size/2, max_height - bin_size/2, by = bin_size)
  full_profile <- data.frame(z = desired_centers)
  merged_profile <- merge(full_profile, lad_df, by = "z", all.x = TRUE)
  merged_profile$lad[is.na(merged_profile$lad)] <- 0
  lad_profile <- list()
  lad_profile$LAI_lidar <- sum(merged_profile$lad * bin_size)
  
  for (i in 1:nrow(merged_profile)) {
    z_center <- merged_profile$z[i]
    val      <- merged_profile$lad[i]
    
    # Construct exact names: LAD_0_5, LAD_5_10, etc.
    z_bottom <- as.integer(z_center - bin_size/2)
    z_top    <- as.integer(z_center + bin_size/2)
    
    col_name <- paste0("LAD_", z_bottom, "_", z_top)
    lad_profile[[col_name]] <- val
  }
  
  return(lad_profile)
}

add_suffix <- function(metrics_list, suffix) {
  if (is.null(metrics_list)) return(NULL)
  names(metrics_list) <- paste0(names(metrics_list), suffix)
  return(metrics_list)
}

get_scalar_stats <- function(observed, predicted) {
  # Remove NAs
  valid <- complete.cases(observed, predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  if(length(obs) < 2) return(c(r = NA, R2=NA, RMSE=NA, Bias=NA))
  
  # Linear model for R2
  mod <- lm(pred ~ obs)
  
  r    <- round(cor(pred, obs), 2)
  R2   <- round(summary(mod)$r.squared, 2)
  RMSE <- round(sqrt(mean((pred - obs)^2)), 2)
  Bias <- round(mean(pred - obs), 2)
  
  return(c(r = r, R2 = R2, RMSE = RMSE, Bias = Bias))
}