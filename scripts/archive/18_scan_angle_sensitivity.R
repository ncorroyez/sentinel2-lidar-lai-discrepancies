# ---
# title:  18_scan_angle_sensitivity.R
# desc:   Option-b scan angle sensitivity — quantifies the multiplicative
#         correction cos(θ̄) applied to existing LAI_ALS (ladstack.tif) on
#         Aigoual.  Does NOT reprocess LAS point clouds; instead computes
#         cos_theta_mean.tif = pixel_metrics(ctg, mean(cos(|ScanAngle|)),
#         res = 10) and multiplies each ladstack layer.
#
#         Beer-Lambert without angle correction:
#           LAD_obs = LAD_true / cos(θ)  → overestimation
#         Correction:
#           LAD_corrected = LAD_obs × cos(θ̄)
#
#         ScanAngle in LAS 1.4 point format 6 is stored in degrees (lidR
#         reads it directly as degrees — do NOT multiply by 0.006).
#
#         Dry-run mode (default): single tile, fast feedback on timing and
#         angle distribution.  Set DRY_RUN <- FALSE for the full Aigoual
#         catalog.
#
# Output:
#   output/intermediate/scan_angle/Aigoual/
#     cos_theta_mean_dryrun.tif   (dry run, one tile)
#     cos_theta_mean.tif          (full run, all tiles)
#     scan_angle_summary.csv      (distribution stats)
#
# Run from NC_Full/:  source("scripts/18_scan_angle_sensitivity.R")
# ---

library("here")
library("terra")
library("lidR")

source(here::here("R", "paths.R"))

DRY_RUN <- TRUE   # set FALSE after validating dry-run timing

# ── Paths ─────────────────────────────────────────────────────────────────────

smb_base <- "/run/user/$(id -u)/gvfs"   # placeholder; resolved below
uid <- system("id -u", intern = TRUE)
smb_base <- file.path(
  paste0("/run/user/", uid, "/gvfs"),
  "smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse"
)
las_dir <- file.path(
  smb_base,
  "root/_PROJETS/2023_2026_These_Nathan_Corroyez",
  "LiDAR_Points_Clouds/Aigoual/2-las_utm"
)
if (!dir.exists(las_dir))
  stop("SMB not mounted or path wrong:\n  ", las_dir)

ladstack_path <- file.path(
  paths$prosail_lidar, "Aigoual", "LiDAR",
  "PAD_Profiles_Classic", "ladstack.tif"
)
if (!file.exists(ladstack_path))
  stop("ladstack.tif not found:\n  ", ladstack_path)

out_dir <- file.path(paths$output, "intermediate", "scan_angle", "Aigoual")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Compute cos_theta_mean raster ─────────────────────────────────────────────

cos_theta_fun <- function(ScanAngle) {
  list(cos_theta_mean = mean(cos(abs(ScanAngle) * pi / 180), na.rm = TRUE))
}

if (DRY_RUN) {
  # Single tile for timing and angle inspection
  tile_path <- file.path(las_dir, "LAS_737500_6338000.las")
  if (!file.exists(tile_path))
    stop("Dry-run tile not found:\n  ", tile_path)

  cat("DRY RUN — single tile:", basename(tile_path), "\n")
  t0 <- proc.time()

  las_tile <- readLAS(tile_path, select = "a")   # "a" = scan angle only
  cat("Points read:", nrow(las_tile@data), "\n")
  cat("ScanAngle range: [",
      round(min(las_tile@data$ScanAngle), 1), ",",
      round(max(las_tile@data$ScanAngle), 1), "] degrees\n")
  cat("cos(|ScanAngle|) range: [",
      round(min(cos(abs(las_tile@data$ScanAngle) * pi / 180)), 4), ",",
      round(max(cos(abs(las_tile@data$ScanAngle) * pi / 180)), 4), "]\n")

  cos_theta_tile <- pixel_metrics(las_tile, cos_theta_fun, res = 10)
  cos_theta_tile <- rast(cos_theta_tile)
  crs(cos_theta_tile) <- "EPSG:32631"

  out_path <- file.path(out_dir, "cos_theta_mean_dryrun.tif")
  writeRaster(cos_theta_tile, out_path, overwrite = TRUE, gdal = "COMPRESS=LZW")

  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  cat("cos_theta_mean tile written:", out_path, "\n")
  cat("Tile elapsed:", elapsed, "s\n")
  cat("Tile extent:", as.character(ext(cos_theta_tile)), "\n")

  # Summarise angle correction factor for this tile
  vals <- terra::values(cos_theta_tile, na.rm = TRUE)
  cat("\n--- cos(θ̄) distribution (tile) ---\n")
  cat("Min :", round(min(vals), 4), "→ max LAI inflation factor:",
      round(1 / min(vals), 3), "×\n")
  cat("Mean:", round(mean(vals), 4), "→ mean LAI reduction:",
      round((1 - mean(vals)) * 100, 1), "%\n")
  cat("Max :", round(max(vals), 4), "\n")

  # Quick impact on ladstack in the tile footprint
  ladstack <- rast(ladstack_path)
  tile_crop <- crop(ladstack, ext(cos_theta_tile))
  cos_tile_crop <- resample(cos_theta_tile, tile_crop, method = "bilinear")

  lai_orig <- sum(tile_crop, na.rm = TRUE)
  lai_corr <- sum(tile_crop * cos_tile_crop, na.rm = TRUE)

  orig_vals <- terra::values(lai_orig, na.rm = TRUE)
  corr_vals <- terra::values(lai_corr, na.rm = TRUE)

  cat("\n--- LAI_ALS impact (tile footprint) ---\n")
  cat("Original  mean LAI_ALS:", round(mean(orig_vals), 3), "\n")
  cat("Corrected mean LAI_ALS:", round(mean(corr_vals), 3), "\n")
  cat("Mean reduction:", round((1 - mean(corr_vals) / mean(orig_vals)) * 100, 1), "%\n")

  cat("\nDRY RUN complete. Review results above, then set DRY_RUN <- FALSE for full run.\n")
  cat("Estimated full-run time (~115 tiles × ", elapsed, "s): ~",
      round(115 * elapsed / 60, 0), "min\n", sep = "")

} else {
  # Full catalog
  cat("FULL RUN — building LAScatalog for Aigoual (", las_dir, ")\n")
  t0 <- proc.time()

  las_files <- list.files(las_dir, pattern = "^LAS_.*\\.las$",
                           full.names = TRUE, ignore.case = FALSE)
  # Fallback: some filesystems don't pass glob through ls
  if (length(las_files) == 0)
    las_files <- list.files(las_dir, full.names = TRUE)[
      grepl("LAS_[0-9]", list.files(las_dir))]
  cat("Tiles found:", length(las_files), "\n")

  ctg <- readLAScatalog(las_files)
  opt_select(ctg) <- "a"          # scan angle only → fast I/O
  opt_chunk_buffer(ctg) <- 0
  opt_progress(ctg) <- TRUE

  cos_theta_full <- pixel_metrics(ctg, cos_theta_fun, res = 10)
  cos_theta_full <- rast(cos_theta_full)
  crs(cos_theta_full) <- "EPSG:32631"

  out_path <- file.path(out_dir, "cos_theta_mean.tif")
  writeRaster(cos_theta_full, out_path, overwrite = TRUE, gdal = "COMPRESS=LZW")

  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  cat("cos_theta_mean.tif written:", out_path, "\n")
  cat("Full run elapsed:", elapsed, "s\n")

  # ── LAI_ALS corrected ────────────────────────────────────────────────────────

  ladstack <- rast(ladstack_path)
  cos_rast  <- resample(cos_theta_full, ladstack, method = "bilinear")

  lai_orig <- sum(ladstack, na.rm = TRUE)
  lai_corr <- sum(ladstack * cos_rast, na.rm = TRUE)

  out_lai_orig <- file.path(out_dir, "lai_als_original.tif")
  out_lai_corr <- file.path(out_dir, "lai_als_corrected.tif")
  writeRaster(lai_orig, out_lai_orig, overwrite = TRUE, gdal = "COMPRESS=LZW")
  writeRaster(lai_corr, out_lai_corr, overwrite = TRUE, gdal = "COMPRESS=LZW")

  orig_vals <- terra::values(lai_orig, na.rm = TRUE)
  corr_vals <- terra::values(lai_corr, na.rm = TRUE)
  cos_vals  <- terra::values(cos_rast, na.rm = TRUE)

  summary_dt <- data.table::data.table(
    metric = c("cos_theta_mean_mean", "cos_theta_mean_sd",
               "cos_theta_mean_min", "cos_theta_mean_max",
               "lai_als_orig_mean", "lai_als_corr_mean",
               "lai_als_reduction_pct"),
    value  = c(mean(cos_vals, na.rm = TRUE), sd(cos_vals, na.rm = TRUE),
               min(cos_vals, na.rm = TRUE),  max(cos_vals, na.rm = TRUE),
               mean(orig_vals, na.rm = TRUE), mean(corr_vals, na.rm = TRUE),
               (1 - mean(corr_vals, na.rm = TRUE) / mean(orig_vals, na.rm = TRUE)) * 100)
  )
  data.table::fwrite(summary_dt, file.path(out_dir, "scan_angle_summary.csv"))

  cat("\n--- Full Aigoual summary ---\n")
  print(summary_dt)
}
