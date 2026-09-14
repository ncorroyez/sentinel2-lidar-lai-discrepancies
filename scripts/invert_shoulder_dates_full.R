# ---
# invert_shoulder_dates_full.R — FULL RUN. Invert all five missing Blois S2 dates
# (leaf-out, senescence, leaf-off + winter anchors) to LAI with the prosail-3.0.0
# ATBD engine (consistency-validated against the old summer product: r=0.997).
# Each LAI raster is ALIGNED to the existing summer Not_Masked grid and written into
#   03_RESULTS/Blois/Metrics/Not_Masked/s2lai_<date>_atbd_res_10_m.tif
# (additive — new dates only, existing summer files untouched). Authorised by user.
# ---

suppressPackageStartupMessages({
  library(here); library(terra); library(prosail); library(liquidSVM)
})
source(here::here("R", "prosail_lut.R"))
source(here::here("R", "prosail_inversion.R"))
source(here::here("R", "get_s2_angles.R"))

DRIVE  <- "/media/corroyez/MyPassport/01_DATA/Blois/Sentinel-2"
NMDIR  <- here::here("03_RESULTS/Blois/Metrics/Not_Masked")
SCRATCH<- here::here("output/intermediate/shoulder_inversion")
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)

# Oct 15 has no extracted reflectance on the drive (no raster_samples, no SAFE) — skipped.
# Feb 24 + Apr 23 already written in the first pass; this pass adds the leaf-off dates.
DATES <- c("2021-11-09", "2021-12-21")

NB_SAMPLES  <- 5000L
BAND_SELECT <- c("B3", "B4", "B8")
BAND_PREFIX <- c("B03", "B04", "B08")
BAND_NAMES  <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
set.seed(42)

srf      <- prosail::get_srf_sensor("Sentinel_2")
bsel     <- match(BAND_SELECT, srf$Spectral_Bands)
template <- terra::rast(file.path(NMDIR, "s2lai_2021-06-14_atbd_res_10_m.tif"))

invert_one <- function(date) {
  cat(sprintf("\n=== %s ===\n", date))
  t0 <- proc.time()
  refl_path <- file.path(DRIVE, date, "raster_samples", sprintf("Blois_001_%s.tiff", date))
  bin_path  <- file.path(DRIVE, date, "raster_samples", sprintf("Blois_001_%s_BIN_v2.tiff", date))
  geom_dir  <- file.path(DRIVE, date, "geom_acq_S2")
  stopifnot(file.exists(refl_path), dir.exists(geom_dir))

  s2 <- terra::rast(refl_path); names(s2) <- BAND_NAMES
  bbox_path <- file.path(SCRATCH, "bbox_tmp.gpkg")
  terra::writeVector(terra::as.polygons(terra::ext(s2), crs = terra::crs(s2)),
                     bbox_path, overwrite = TRUE)
  geom_s2 <- get_s2_angles(path_angles = geom_dir, path_bbox = bbox_path, dateAcq = date)
  cat(sprintf("  geom: sza %.1f vza %.1f psi %.1f\n",
              mean(c(geom_s2$MinAngle["sza"],geom_s2$MaxAngle["sza"])),
              mean(c(geom_s2$MinAngle["vza"],geom_s2$MaxAngle["vza"])),
              mean(c(geom_s2$MinAngle["psi"],geom_s2$MaxAngle["psi"]))))
  geom_acq <- list(
    min = data.frame(tto = geom_s2$MinAngle["vza"], tts = geom_s2$MinAngle["sza"], psi = geom_s2$MinAngle["psi"]),
    max = data.frame(tto = geom_s2$MaxAngle["vza"], tts = geom_s2$MaxAngle["sza"], psi = geom_s2$MaxAngle["psi"]))

  ip  <- prosail::get_atbd_lut_input(nb_samples = NB_SAMPLES, geom_acq = geom_acq,
                                     codistribution_lai = TRUE)
  brf <- build_prosail_lut(ip, srf, bsel)
  svr <- lapply(train_svr_ensemble(brf, ip$lai, NB_SAMPLES), liquidSVM::unserialize.liquidSVM)

  m <- terra::rast(bin_path); m[m != 1] <- NA
  lai <- predict_lai_full_raster(svr, s2, BAND_PREFIX, mask_rast = m)

  # align to existing summer Not_Masked grid (10 m offset) so dates stack cleanly
  lai_al <- terra::resample(lai, template, method = "bilinear")
  out <- file.path(NMDIR, sprintf("s2lai_%s_atbd_res_10_m.tif", date))
  terra::writeRaster(lai_al, out, overwrite = TRUE, gdal = "COMPRESS=LZW")
  terra::writeRaster(lai_al, file.path(SCRATCH, basename(out)), overwrite = TRUE, gdal = "COMPRESS=LZW")

  v <- terra::values(lai_al); v <- v[is.finite(v)]
  cat(sprintf("  LAI aligned: n=%d mean=%.2f median=%.2f range=%.2f-%.2f  (%.0fs) -> %s\n",
              length(v), mean(v), median(v), min(v), max(v),
              (proc.time()-t0)[["elapsed"]], basename(out)))
}

for (d in DATES) invert_one(d)
cat("\nDONE full shoulder inversion. New dates written into Not_Masked:\n  ",
    paste(DATES, collapse = "\n   "), "\n")
