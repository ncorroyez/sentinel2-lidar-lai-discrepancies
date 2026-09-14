# ---
# invert_shoulder_dates.R — DRY-RUN. Invert Sentinel-2 reflectance to LAI for the
# leaf-out / leaf-off shoulder dates using the prosail-3.0.0-native helpers
# (R/prosail_lut.R, R/prosail_inversion.R, R/get_s2_angles.R) — the same engine as
# scripts/ch2/steps/25_train_apply_prosail_atbd_rasters.R (ATBD, codistribution_lai=TRUE).
#
# Inputs (read-only, external drive):
#   01_DATA/Blois/Sentinel-2/<date>/raster_samples/Blois_001_<date>.tiff   (10-band BOA)
#   01_DATA/Blois/Sentinel-2/<date>/raster_samples/Blois_001_<date>_BIN_v2.tiff (cloud mask)
#   01_DATA/Blois/Sentinel-2/<date>/geom_acq_S2/{saa,sza,vaa,vza}_<date>.tiff
#
# Outputs (scratch only — NOT 01_DATA, NOT Not_Masked yet):
#   output/intermediate/shoulder_inversion/s2lai_<date>_atbd_res_10_m.tif
#
# This run ALSO inverts 2021-06-14 and compares to the existing
# 03_RESULTS/Blois/Metrics/Not_Masked/s2lai_2021-06-14_atbd_res_10_m.tif
# to check that the modern engine reproduces the old summer product (consistency gate)
# before we trust the shoulder LAI.
# ---

suppressPackageStartupMessages({
  library(here); library(terra); library(prosail); library(liquidSVM)
})
source(here::here("R", "prosail_lut.R"))
source(here::here("R", "prosail_inversion.R"))
source(here::here("R", "get_s2_angles.R"))

DRIVE   <- "/media/corroyez/MyPassport/01_DATA/Blois/Sentinel-2"
EXIST   <- here::here("03_RESULTS/Blois/Metrics/Not_Masked")
OUT     <- here::here("output/intermediate/shoulder_inversion")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

NB_SAMPLES    <- 5000L
BAND_SELECT   <- c("B3", "B4", "B8")
BAND_PREFIX   <- c("B03", "B04", "B08")
BAND_NAMES    <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
set.seed(42)

srf  <- prosail::get_srf_sensor("Sentinel_2")
bsel <- match(BAND_SELECT, srf$Spectral_Bands)

invert_one <- function(date, apply_cloud_mask = TRUE) {
  cat(sprintf("\n=== invert %s ===\n", date))
  t0 <- proc.time()
  refl_path <- file.path(DRIVE, date, "raster_samples", sprintf("Blois_001_%s.tiff", date))
  bin_path  <- file.path(DRIVE, date, "raster_samples", sprintf("Blois_001_%s_BIN_v2.tiff", date))
  geom_dir  <- file.path(DRIVE, date, "geom_acq_S2")
  stopifnot(file.exists(refl_path), dir.exists(geom_dir))

  s2 <- terra::rast(refl_path)
  stopifnot(terra::nlyr(s2) == length(BAND_NAMES))
  names(s2) <- BAND_NAMES

  # bbox for angle cropping: derive from the reflectance footprint
  bbox_path <- file.path(OUT, "bbox_tmp.gpkg")
  terra::writeVector(terra::as.polygons(terra::ext(s2), crs = terra::crs(s2)),
                     bbox_path, overwrite = TRUE)
  geom_s2 <- get_s2_angles(path_angles = geom_dir, path_bbox = bbox_path, dateAcq = date)
  cat(sprintf("  geometry: sza %.1f-%.1f  vza %.1f-%.1f  psi %.1f-%.1f\n",
              geom_s2$MinAngle["sza"], geom_s2$MaxAngle["sza"],
              geom_s2$MinAngle["vza"], geom_s2$MaxAngle["vza"],
              geom_s2$MinAngle["psi"], geom_s2$MaxAngle["psi"]))

  geom_acq <- list(
    min = data.frame(tto = geom_s2$MinAngle["vza"], tts = geom_s2$MinAngle["sza"], psi = geom_s2$MinAngle["psi"]),
    max = data.frame(tto = geom_s2$MaxAngle["vza"], tts = geom_s2$MaxAngle["sza"], psi = geom_s2$MaxAngle["psi"]))

  ip  <- prosail::get_atbd_lut_input(nb_samples = NB_SAMPLES, geom_acq = geom_acq,
                                     codistribution_lai = TRUE)
  brf <- build_prosail_lut(ip, srf, bsel)
  svr <- train_svr_ensemble(brf, ip$lai, NB_SAMPLES)
  svr <- lapply(svr, liquidSVM::unserialize.liquidSVM)

  mask_rast <- NULL
  if (apply_cloud_mask && file.exists(bin_path)) {
    m <- terra::rast(bin_path); m[m != 1] <- NA          # 1 = valid
    mask_rast <- m
  }
  lai <- predict_lai_full_raster(svr, s2, BAND_PREFIX, mask_rast = mask_rast)
  out <- file.path(OUT, sprintf("s2lai_%s_atbd_res_10_m.tif", date))
  terra::writeRaster(lai, out, overwrite = TRUE, gdal = "COMPRESS=LZW")

  v <- terra::values(lai); v <- v[is.finite(v)]
  dt <- round((proc.time() - t0)[["elapsed"]], 1)
  cat(sprintf("  LAI: n=%d valid  mean=%.2f  median=%.2f  range=%.2f-%.2f  (%.0fs)\n",
              length(v), mean(v), median(v), min(v), max(v), dt))
  list(rast = lai, out = out, vals = v)
}

# ---- 1) consistency gate: 2021-06-14 modern vs existing old product ----
r0614 <- invert_one("2021-06-14", apply_cloud_mask = FALSE)
old_path <- file.path(EXIST, "s2lai_2021-06-14_atbd_res_10_m.tif")
if (file.exists(old_path)) {
  old <- terra::rast(old_path)
  newr <- terra::resample(r0614$rast, old, method = "bilinear")
  d <- data.frame(old = terra::values(old)[,1], new = terra::values(newr)[,1])
  d <- d[is.finite(d$old) & is.finite(d$new), ]
  cat(sprintf("\n=== CONSISTENCY GATE (06-14 modern vs existing old) ===\n"))
  cat(sprintf("  n=%d  r=%.3f  mean old=%.2f new=%.2f  ratio(new/old)=%.3f  RMSE=%.3f\n",
              nrow(d), cor(d$old, d$new), mean(d$old), mean(d$new),
              mean(d$new)/mean(d$old), sqrt(mean((d$new-d$old)^2))))
} else {
  cat("\n[warn] existing s2lai_2021-06-14 not found at", old_path, "\n")
}

# ---- 2) the real target: leaf-out 2021-04-23 ----
r0423 <- invert_one("2021-04-23", apply_cloud_mask = TRUE)

cat("\nDONE dry-run. Outputs in", OUT, "\n")
