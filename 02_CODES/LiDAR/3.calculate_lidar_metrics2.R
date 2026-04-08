# ---
# title: "3.calculate_lidar_metrics2.R"
# ---

# ----------------------------- Clear environment -------------------------------------
rm(list=ls(all=TRUE))
gc()
# -------------------------------------------------------------------------------------

library("lidR")
library("terra")

# Set working directory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

# Source external functions
source("../libraries/functions_lidar2.R")
source("../libraries/functions_create_masks.R")

# Base configuration - Blois ONLY
sites <- c("Blois")
data_dir <- "/media/corroyez/MyPassport/01_DATA"
leaf_state <- "leaf_on"
target_resolution <- 20

for (site in sites) {
  print(paste("Processing:", site))
  
  # Load Sentinel-2 raster for strict grid alignment
  s2_rast <- terra::rast('/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31TCN_A031222_20210614T105443_Refl')[[1]]
  s2_res <- terra::res(s2_rast)[1]
  
  if (target_resolution == s2_res) {
    template_rast <- s2_rast
  } else {
    # Aggregate the S2 grid to the target resolution (e.g., from 10m to 20m -> fact = 2)
    fact_val <- target_resolution / s2_res
    template_rast <- terra::aggregate(s2_rast, fact = fact_val)
  }
  
  results_path <- file.path("../../03_RESULTS", site)
  masks_dir <- file.path(results_path, "LiDAR/Heterogeneity_Masks")
  
  metrics_dir <- file.path(results_path, "Metrics")
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # LiDAR point clouds
  las_norm_utm_dir <- file.path(data_dir, site, "LiDAR", leaf_state, "3-las_normalized_utm")
  normalized_ctg <- readLAScatalog(las_norm_utm_dir)
  
  # LiDR processing options
  opt_chunk_size(normalized_ctg) <- 0 
  opt_chunk_buffer(normalized_ctg) <- 10
  opt_select(normalized_ctg) <- "xyz"
  set_lidr_threads(1)
  
  # CHM input (assuming 1m resolution)
  chm_path <- file.path(results_path, "LiDAR", 
                        # leaf_state, 
                        "chm/res_1_m/chm.tif")
  chm <- terra::rast(chm_path)
  
  # Aggregation factor from CHM to target resolution (approx 10m)
  # Actual alignment will be forced by terra::resample later
  # target_res <- terra::res(s2_rast)[1]
  # fact_val <- round(target_res / terra::res(chm)[1])
  
  # ----------------------------------------------------------------------------------
  # 1. Hmax (95th and 98th percentiles)
  # Calculated directly from the point cloud for PERFECT alignment with template_rast
  # ----------------------------------------------------------------------------------
  print("Calculating Hmax...")
  hmax_metrics <- lidR::pixel_metrics(normalized_ctg, 
                                      ~as.list(quantile(Z, probs = c(0.95, 0.98), na.rm = TRUE)), 
                                      res = template_rast)
  
  # Rename layers and save
  names(hmax_metrics) <- c("hmax_p95", "hmax_p98")
  
  for (i in 1:terra::nlyr(hmax_metrics)) {
    p_name <- names(hmax_metrics)[i]
    out_filename <- paste0(p_name, "_res_", target_resolution, "_m.tif")
    
    terra::writeRaster(hmax_metrics[[i]],
                       filename = file.path(raw_dir, out_filename),
                       filetype = "GTiff",
                       overwrite = TRUE)
    apply_and_save_masks(hmax_metrics[[i]], out_filename, masks_dir, metrics_dir)
  }
  
  # ----------------------------------------------------------------------------------
  # 2. Dynamic fCover 
  # Calculated from point cloud (proportion of points >= 30% of local Hmax95)
  # Note: This is a point-based gap fraction, replacing the CHM-based gap fraction.
  # ----------------------------------------------------------------------------------
  print("Calculating fCover...")
  fcover_dynamic_pc <- function(z) {
    if (length(z) == 0) return(list(fCover_dynamic = NA_real_))
    local_hmax95 <- quantile(z, probs = 0.95, na.rm = TRUE)
    threshold <- 0.30 * local_hmax95
    fcover_val <- sum(z >= threshold) / length(z)
    return(list(fCover_dynamic = fcover_val))
  }
  
  fcover_raster <- lidR::pixel_metrics(normalized_ctg, ~fcover_dynamic_pc(Z), res = template_rast)
  
  fcover_file <- paste0("fCover_dynamic_res_", target_resolution, "_m.tif")
  terra::writeRaster(fcover_raster,
                     filename = file.path(raw_dir, fcover_file),
                     filetype = "GTiff",
                     overwrite = TRUE)
  apply_and_save_masks(fcover_raster, fcover_file, masks_dir, metrics_dir)
  
  # ----------------------------------------------------------------------------------
  # 3. LAI (z0 = 1m)
  # ----------------------------------------------------------------------------------
  print("Calculating LAI...")
  lai_raster <- lidR::pixel_metrics(normalized_ctg, ~myPAI(Z, zmin = 1), res = template_rast)
  lai_raster[lai_raster == Inf] <- NA
  lai_file <- paste0("lai_z1_res_", target_resolution, "_m.tif")
  
  terra::writeRaster(lai_raster,
                     filename = file.path(raw_dir, lai_file),
                     filetype = "GTiff",
                     overwrite = TRUE)
  apply_and_save_masks(lai_raster, lai_file, masks_dir, metrics_dir)
  
  # ----------------------------------------------------------------------------------
  # 4. LAD Vertical Profiles (z0 = 1m)
  # ----------------------------------------------------------------------------------
  print("Calculating LAD...")
  pad_rasters <- lidR::pixel_metrics(normalized_ctg, ~myPAD_z1(Z), res = template_rast)
  lad_filename <- paste0("lad_profiles_z1_res_", target_resolution, "_m.tif")
  
  terra::writeRaster(pad_rasters,
                     filename = file.path(raw_dir, lad_filename),
                     filetype = "GTiff",
                     overwrite = TRUE)
  apply_and_save_masks(pad_rasters, lad_filename, masks_dir, metrics_dir)
  
  for (i in 1:terra::nlyr(pad_rasters)) {
    layer <- subset(pad_rasters, i)
    layer_filename <- paste0(names(pad_rasters)[i], "_res_", target_resolution, "_m.tif")
    apply_and_save_masks(layer, layer_filename, masks_dir, metrics_dir)
  }
}