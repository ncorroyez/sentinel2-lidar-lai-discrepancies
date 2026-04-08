# ================================================================
# GEDI extractor - Multi-site version (no manual file selection)
# ================================================================

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(rhdf5)
library(data.table)
library(sf)
library(stringr)
library(lubridate)
library(bit64)
library(dplyr)
library(terra)
library(salpa)

# ---------------------
# Parameters
# ---------------------
# data_dir <- "../../01_DATA"
data_dir <- "/media/corroyez/MyPassport/01_DATA"
results_dir <- "../../03_RESULTS"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Aigoual")
dz <- 5  # bin size for pavd_z profiles

# ---------------------
# Helper functions
# ---------------------

# l2b_h5file <- "/media/corroyez/MyPassport/01_DATA/Aigoual/GEDI/GEDI02_B_2021184062047_O14477_03_T10849_02_003_01_V002.h5"
# l2b_h5_contents <- h5ls(l2b_h5file, recursive = TRUE)
# l2b_beam_subset <- l2b_h5_contents[grepl("BEAM0000", l2b_h5_contents$group), ]
# write.csv(l2b_beam_subset, file = "GEDI_L2B_beam_subset.csv", row.names = F)

# l2a_h5file <- "/media/corroyez/MyPassport/01_DATA/Aigoual/GEDI/GEDI02_A_2021184062047_O14477_03_T10849_02_003_02_V002.h5"
# l2a_h5_contents <- h5ls(l2a_h5file, recursive = TRUE)
# l2a_beam_subset <- l2a_h5_contents[grepl("BEAM0000", l2a_h5_contents$group), ]
# write.csv(l2a_beam_subset, file = "GEDI_L2A_beam_subset.csv", row.names = F)

# Calculate Hmax from PAVD profile
Hmax_padz <- function(padz, dz = 5) {
  if (length(padz) == 0 || all(is.na(padz))) return(0)
  padz <- as.numeric(padz)
  nz <- which(padz > 0)
  if (length(nz) == 0) return(0)
  last_idx <- max(nz)
  (length(padz) - last_idx) * dz
}

# Extract RH95 and RH98 from L2A file
get_l2a_metrics <- function(l2a_file, beam_name) {
  fid_a <- tryCatch(H5Fopen(l2a_file), error = function(e) NULL)
  if (is.null(fid_a)) return(NULL)
  on.exit(H5Fclose(fid_a), add = TRUE)
  
  path_sn <- paste0("/", beam_name, "/geolocation/shot_number")
  path_rh <- paste0("/", beam_name, "/rh")
  
  if (!H5Lexists(fid_a, path_sn) || !H5Lexists(fid_a, path_rh)) return(NULL)
  
  sn <- H5Dread(eval(parse(text=paste0("fid_a&'", path_sn, "'"))), bit64conversion = "bit64")
  rh_subset <- h5read(fid_a, path_rh, index = list(c(96, 99), NULL))
  
  dt_a <- data.table(
    shot_number = as.character(sn),
    RH95 = as.numeric(rh_subset[1, ]), 
    RH98 = as.numeric(rh_subset[2, ])
  )
  dt_a <- unique(dt_a, by = "shot_number")
  return(dt_a)
}

# Main function to extract Beam Data
infos_beam <- function(h5file_path, beam_name, l2a_path = NULL, dz = 5) {
  fid <- H5Fopen(h5file_path)
  base <- paste0("/", beam_name)
  on.exit(H5Fclose(fid), add = TRUE)
  
  safe_read <- function(dataset, bit64conv = FALSE) {
    path <- paste0(base, dataset)
    if (H5Lexists(fid, path)) {
      if (bit64conv) {
        return(H5Dread(eval(parse(text=paste0("fid&'", path, "'"))), bit64conversion = "bit64"))
      } else {
        return(h5read(fid, path))
      }
    } else {
      return(NULL)
    }
  }
  
  # Read flags first for quick filtering
  algorithmrun_flag <- safe_read("/algorithmrun_flag")
  shot_number       <- safe_read("/geolocation/shot_number", bit64conv = TRUE)
  degrade_flag      <- safe_read("/geolocation/degrade_flag")
  elev_bin0         <- safe_read("/geolocation/elevation_bin0")
  elev_dem          <- safe_read("/geolocation/digital_elevation_model")
  if (is.null(algorithmrun_flag) || is.null(elev_bin0) ||
      is.null(elev_dem) || is.null(shot_number)) return(NULL)
  
  diff_alt0 <- elev_bin0 - elev_dem
  subset_idx <- which(
    algorithmrun_flag == 1 &
      diff_alt0 < 100 & diff_alt0 > 0 &
      (if (!is.null(degrade_flag)) degrade_flag == 0 else TRUE)
  )
  # sensitivity >= 0.95 (high waveform sensitivity) & rh100 > 0 (positive only)
  # SALPA paper, + check filters Anouk
  if (length(subset_idx) == 0) return(NULL)
  
  # Free memory
  rm(algorithmrun_flag, elev_bin0, elev_dem, degrade_flag)
  
  # Read data variables
  delta_time        <- safe_read("/geolocation/delta_time")
  lat_lowestmode    <- safe_read("/geolocation/lat_lowestmode")
  lon_lowestmode    <- safe_read("/geolocation/lon_lowestmode")
  elev_lowestmode   <- safe_read("/geolocation/elev_lowestmode")
  rh100             <- safe_read("/rh100")
  pai               <- safe_read("/pai") # Total PAI (Target for normalization)
  cover             <- safe_read("/cover")
  fhd_normal        <- safe_read("/fhd_normal")
  omega             <- safe_read("/omega")
  # pai_z             <- safe_read("/pai_z")
  pavd_z            <- safe_read("/pavd_z") # Density Profile (Shape source)
  
  pai_subset <- as.numeric(pai[subset_idx])
  dt <- data.table(
    shot_number = as.character(shot_number[subset_idx]),
    delta_time  = as.numeric(delta_time[subset_idx]),
    lat         = as.numeric(lat_lowestmode[subset_idx]),
    lon         = as.numeric(lon_lowestmode[subset_idx]),
    elev        = as.numeric(elev_lowestmode[subset_idx]),
    RH100       = as.numeric(rh100[subset_idx]) / 100,
    pai         = pai_subset,
    cover       = as.numeric(cover[subset_idx]),
    fhd_normal  = as.numeric(fhd_normal[subset_idx]),
    omega       = as.numeric(omega[subset_idx]),
    beam        = beam_name
  )
  
  # if (!is.null(pavd_z) && is.matrix(pavd_z)) {
  #   if (ncol(pavd_z) >= max(subset_idx)) {
  #     Hmax_vals <- apply(pavd_z[, subset_idx, drop = FALSE], 2, Hmax_padz, dz = dz)
  #     dt[, Hmax := as.numeric(Hmax_vals)]
  #   } else {
  #     dt[, Hmax := NA_real_]
  #   }
  # } else {
  #   dt[, Hmax := NA_real_]
  # }
  
  # --- PAI Profiles Processing ---.
  nb_bins <- 30 
  bin_names <- paste0("PAI_", seq(0, by=dz, length.out=nb_bins), "_", seq(dz, by=dz, length.out=nb_bins), "m")
  
  if (!is.null(pavd_z) && is.matrix(pavd_z) && ncol(pavd_z) >= max(subset_idx)) {
    
    # Subset matrix immediately
    pavd_sub <- pavd_z[, subset_idx, drop = FALSE]
    rm(pavd_z) # Free the huge original matrix
    
    # Pad matrix if rows missing
    if(nrow(pavd_sub) < nb_bins) {
      missing_rows <- nb_bins - nrow(pavd_sub)
      pavd_sub <- rbind(pavd_sub, matrix(0, nrow=missing_rows, ncol=ncol(pavd_sub)))
    } else {
      pavd_sub <- pavd_sub[1:nb_bins, ] 
    }
    
    # Calculate PAI Area per bin (Density * dz)
    raw_pai_matrix <- pavd_sub * dz
    sums_raw <- colSums(raw_pai_matrix, na.rm = TRUE)
    
    # Normalize so sum(bins) == Total PAI
    factors <- pai_subset / sums_raw
    factors[is.na(factors) | is.infinite(factors) | sums_raw < 0.0001] <- 1
    
    # Apply factor
    final_pai_matrix <- t(t(raw_pai_matrix) * factors)
    rm(raw_pai_matrix)
    
    # Add to DT as columns
    dt_bins <- as.data.table(t(final_pai_matrix))
    colnames(dt_bins) <- bin_names
    dt <- cbind(dt, dt_bins)
    
    dt[, Hmax := apply(pavd_sub, 2, Hmax_padz, dz = dz)]
    rm(pavd_sub, final_pai_matrix)
    
  } else {
    # Fallback if no profile data
    dt[, Hmax := NA_real_]
    dt_bins_na <- as.data.table(matrix(NA_real_, nrow=nrow(dt), ncol=length(bin_names)))
    colnames(dt_bins_na) <- bin_names
    dt <- cbind(dt, dt_bins_na)
  }
  
  # Aggregate values > 40m to PAI_35_40m
  col_target <- bin_names[8]   # "PAI_35_40m"
  cols_high  <- bin_names[9:30] # "PAI_40_45m" ... "PAI_145_150m"
  if (all(cols_high %in% colnames(dt)) && col_target %in% colnames(dt)) {
    sum_high <- rowSums(dt[, ..cols_high], na.rm = TRUE)
    dt[is.na(get(col_target)), (col_target) := 0]
    dt[, (col_target) := get(col_target) + sum_high]
    dt[, (cols_high) := NULL]
  }
  
  # Merge L2A metrics
  if (!is.null(l2a_path) && file.exists(l2a_path)) {
    dt_l2a <- get_l2a_metrics(l2a_file = l2a_path, beam_name = beam_name)
    if (!is.null(dt_l2a)) {
      dt <- merge(dt, dt_l2a, by = "shot_number", all.x = TRUE)
    } else {
      dt[, `:=`(RH95 = NA_real_, RH98 = NA_real_)]
    }
  } else {
    dt[, `:=`(RH95 = NA_real_, RH98 = NA_real_)]
  }
  
  return(dt)
}

# Wrapper to loop through all beams in a file
infos_allbeams <- function(h5file_path, l2a_path = NULL, dz = 5) {
  groups <- h5ls(h5file_path, recursive = FALSE)$name
  groups <- groups[groups != "metadata"]
  result_list <- list()
  
  for (beam in groups) {
    dt_beam <- tryCatch(
      infos_beam(h5file_path, beam, l2a_path = l2a_path, dz = dz),
      error = function(e) {
        message("Error reading beam ", beam, " in ", basename(h5file_path), ": ", e$message)
        NULL
      }
    )
    if (!is.null(dt_beam) && nrow(dt_beam) > 0) result_list[[beam]] <- dt_beam
  }
  
  if (length(result_list) == 0) return(NULL)
  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# ================================================================
# MAIN LOOP OVER SITES
# ================================================================

for (site in sites) {
  message("\n==============================")
  message("Processing site: ", site)
  message("==============================")
  
  gedi_dir <- file.path(data_dir, site, "GEDI")
  shp_path <- file.path(data_dir, site, "Geo_Files", "utm_init.shp")
  dtm_path <- file.path(results_dir, site, "LiDAR/dtm/res_1_m", 
                        "rasterize_terrain.vrt")
  
  if (!dir.exists(gedi_dir)) next
  if (!file.exists(shp_path)) next
  
  # Load Shapefile (Memory Opt: Load once)
  shp_utm <- st_read(shp_path, quiet = TRUE)
  if (st_crs(shp_utm)$epsg != 32631) {
    shp_utm <- st_transform(shp_utm, crs = 32631)
  }
  # Create WGS84 bbox for fast filtering
  shp_wgs84 <- st_transform(shp_utm, crs = 4326)
  bbox_wgs84 <- st_bbox(shp_wgs84)
  
  h5_files <- list.files(gedi_dir, pattern = "GEDI02_B_.*\\.h5$", full.names = TRUE)
  if (length(h5_files) == 0) next
  
  all_filtered_shots <- list()
  
  for (i in seq_along(h5_files)) {
    f_l2b <- h5_files[i]
    fname <- basename(f_l2b)
    
    if (i %% 5 == 0) gc()
    
    # Match with L2A file
    match_id <- str_extract(fname, "(?<=GEDI02_B_)[0-9]+_O[^_]+_[0-9]+_T[^_]+")
    f_l2a <- NULL
    if (!is.na(match_id)) {
      l2a_pattern <- paste0("GEDI02_A_", match_id, ".*\\.h5$")
      found_l2a <- list.files(gedi_dir, pattern = l2a_pattern, full.names = TRUE)
      if (length(found_l2a) > 0) f_l2a <- found_l2a[1]
    }
    
    message(sprintf("Processing file %d / %d : %s %s", 
                    i, length(h5_files), fname, 
                    ifelse(!is.null(f_l2a), paste0("[+ L2A]"), "[No L2A]")))
    
    # Extract Raw Data
    dt_file <- infos_allbeams(f_l2b, l2a_path = f_l2a, dz = dz)
    if (is.null(dt_file) || nrow(dt_file) == 0) next
    
    # Memory Optimization: Spatial Filter
    
    # Fast Bounding Box Check
    file_lat_rng <- range(dt_file$lat, na.rm = TRUE)
    file_lon_rng <- range(dt_file$lon, na.rm = TRUE)
    
    overlap_lat <- (file_lat_rng[1] <= bbox_wgs84$ymax && file_lat_rng[2] >= bbox_wgs84$ymin)
    overlap_lon <- (file_lon_rng[1] <= bbox_wgs84$xmax && file_lon_rng[2] >= bbox_wgs84$xmin)
    
    if (!overlap_lat || !overlap_lon) {
      rm(dt_file)
      next 
    }
    
    # Precise Clipping
    dt_sf <- st_as_sf(dt_file, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    dt_sf_utm <- st_transform(dt_sf, crs = 32631)
    
    dt_clipped <- st_intersection(dt_sf_utm, shp_utm)
    
    if (nrow(dt_clipped) == 0) {
      rm(dt_file, dt_sf, dt_sf_utm, dt_clipped)
      next
    }
    
    # Store only valid shots
    dt_keep <- as.data.table(dt_clipped)
    
    parts  <- str_split(fname, "_")[[1]]
    orbit  <- if (length(parts) >= 4) substr(parts[4], 2, 6) else NA_character_
    dt_keep[, orbit := orbit]
    
    coords <- st_coordinates(dt_clipped)
    dt_keep[, `:=`(X_utm = round(coords[,1], 5), Y_utm = round(coords[,2], 5))]
    dt_keep[, geometry := NULL]
    
    all_filtered_shots[[fname]] <- dt_keep
    
    rm(dt_file, dt_sf, dt_sf_utm, dt_clipped, dt_keep)
  }
  
  if (length(all_filtered_shots) == 0) {
    message("  No shots found inside shapefile for ", site)
    next
  }
  
  # Bind and Finalize
  df_final <- rbindlist(all_filtered_shots, use.names = TRUE, fill = TRUE)
  
  # Time conversion
  gedi_epoch <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
  df_final[, time_utc := gedi_epoch + as.numeric(delta_time)]
  df_final[, `:=`(year = format(time_utc, "%Y"), 
                  month = format(time_utc, "%m"),
                  day = format(time_utc, "%d"),
                  heure = format(time_utc, "%H:%M:%S"))]
  
  # SALPA coordinates correction
  # if (file.exists(dtm_path)) {
  #   message("  Applying SALPA coordinate correction using: ", basename(dtm_path))
  #   
  #   # Convert to SF for SALPA
  #   sf_for_correction <- st_as_sf(df_final, coords = c("X_utm", "Y_utm"), crs = 32631, remove = FALSE)
  #   
  #   # Load Reference Raster (ALS DTM)
  #   ref_dtm <- terra::rast(dtm_path)
  #   
  #   # Run Correction
  #   sf_corrected <- tryCatch({
  #     # Optional: Linear Alignment (if you suspect track shift)
  #     # aligned_pts <- linear_alignment(sf_for_correction, crs_code = 32631) 
  #     
  #     # Positional Correction
  #     positional_correction(
  #       lidar_footprints = sf_for_correction,
  #       input_rast = ref_dtm,
  #       minimizing_method = "euclidean", 
  #       target_variable = "mean", 
  #       lidar_value = "elev", 
  #       buf = 12.5,
  #       lower_bounds = c(-25, -25),
  #       upper_bounds = c(25, 25),
  #       crs_code = 32631,
  #       max_iter = 10,
  #       parallel = F
  #     )
  #   }, error = function(e) {
  #     message("  SALPA Error: ", e$message)
  #     return(NULL)
  #   })
  #   
  #   if (!is.null(sf_corrected)) {
  #     message("  SALPA Correction successful.")
  #     # Update the X_utm/Y_utm in the data.table with the new coordinates
  #     new_coords <- st_coordinates(sf_corrected$position_adjustment)
  #     df_final[, X_utm_corr := round(new_coords[,1], 2)]
  #     df_final[, Y_utm_corr := round(new_coords[,2], 2)]
  #     
  #     # Update Lat/Lon to match new UTM
  #     sf_wgs <- st_transform(sf_corrected$position_adjustment, 4326)
  #     new_wgs <- st_coordinates(sf_wgs)
  #     df_final[, lon_corr := round(new_wgs[,1], 6)]
  #     df_final[, lat_corr := round(new_wgs[,2], 6)]
  #     
  #     df_final[, shift_m := round(sqrt((X_utm_corr - X_utm)^2 + (Y_utm_corr - Y_utm)^2), 2)]
  #     df_final[, salpa_corrected := TRUE]
  #   } else {
  #     # If SALPA returned NULL (but didn't crash), fill new cols with originals or NA
  #     df_final[, `:=`(X_utm_corr = X_utm, Y_utm_corr = Y_utm, 
  #                     lon_corr = lon, lat_corr = lat, 
  #                     shift_m = 0, salpa_corrected = FALSE)]
  #   }
  # } else {
  #   message("  WARNING: DTM not found for SALPA at ", dtm_path, ". Skipping correction.")
  #   df_final[, `:=`(X_utm_corr = X_utm, Y_utm_corr = Y_utm, 
  #                   lon_corr = lon, lat_corr = lat, 
  #                   shift_m = 0, salpa_corrected = FALSE)]
  #   }
  
  # Rounding
  # cols_round <- c("lat", "lon", "elev", "RH100", "RH95", "RH98", "pai", "cover", "Hmax", "fhd_normal", "omega")
  cols_round <- c("lat", "lon", "lat_corr", "lon_corr", "shift_m",
                  "elev", "RH100", "RH95", "RH98",
                  "pai", "cover", "Hmax", "fhd_normal", "omega")
  pai_cols <- grep("^PAI_", colnames(df_final), value = TRUE)
  cols_round <- c(cols_round, pai_cols)
  cols_round <- intersect(cols_round, colnames(df_final))
  
  df_final[, (cols_round) := lapply(.SD, round, 5), .SDcols = cols_round]
  
  # Export
  out_base <- file.path(data_dir, site, "Geo_Files")
  if (!dir.exists(out_base)) dir.create(out_base)
  
  csv_path <- file.path(out_base, paste0("gedi_shots_", site, ".csv"))
  fwrite(df_final, csv_path)
  message("Saved CSV: ", basename(csv_path))
  
  gpkg_path <- file.path(out_base, paste0("gedi_shots_", site, ".gpkg"))
  
  df_sf_final <- st_as_sf(df_final, coords = c("X_utm", "Y_utm"),
                          crs = 32631, remove = FALSE)
  
  # df_sf_original <- st_as_sf(df_final, coords = c("X_utm", "Y_utm"), crs = 32631, remove = FALSE)
  # st_write(df_sf_original, gpkg_path, layer = "original_shots", delete_dsn = TRUE, quiet = TRUE)

  # We check if correction exists; if so, save as "corrected_shots"
  # if ("X_utm_corr" %in% names(df_final)) {
  #   # Note: delete_dsn=FALSE here so we append to the file, not overwrite the first layer
  #   df_sf_corrected <- st_as_sf(df_final, coords = c("X_utm_corr", "Y_utm_corr"), crs = 32631, remove = FALSE)
  #   st_write(df_sf_corrected, gpkg_path, layer = "corrected_shots", delete_dsn = FALSE, quiet = TRUE)
  # }
  
  st_write(df_sf_final, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
  message("Saved GPKG: ", basename(gpkg_path))
  
  rm(all_filtered_shots, df_final, df_sf_final)
  # if (exists("df_sf_corrected")) rm(df_sf_corrected)
  gc()
}

message("\nAll sites processed.")