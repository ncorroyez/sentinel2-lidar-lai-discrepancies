# ---
# title:  25_train_apply_prosail_atbd_rasters.R
# desc:   Produces full-raster LAI maps for the two ATBD validation variants:
#           ATBD_T : get_atbd_v3_lut_input(codistribution_lai = TRUE)
#           ATBD_F : get_atbd_v3_lut_input(codistribution_lai = FALSE)
#                    (reuses SVR from 02a if present, else retrains)
#         For each site, applies the per-variant SVR ensemble to the 10-m
#         Sentinel-2 reflectance raster and writes one GeoTIFF per variant.
#         Outputs are consumed by 24_plot_s2_atbd_validation.R.
#
#         Reads (read-only):
#           {ext_results}/{site}/PROSAIL_Optimization/geomAcq_S2/
#           {raw_data}/{site}/Geo_Files/aoi_bbox.GPKG
#           03_RESULTS/{site}/PROSAIL_Optimization/{l2a_id}/Reflectance/{l2a_id}_Refl
#           03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/
#             artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi
#
#         Outputs:
#           revision/output/intermediate/sm6/{site}/
#             s2lai_summer_atbd_T_res_10_m.tif
#             s2lai_summer_atbd_F_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/25_train_apply_prosail_atbd_rasters.R")
# ---

library(here)
library(terra)
library(prosail)
library(liquidSVM)
library(cli)

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "prosail_lut.R"))
source(here::here("revision", "R", "prosail_inversion.R"))
source(here::here("revision", "R", "get_s2_angles.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal")
nb_samples    <- 5000L
band_select   <- c("B3", "B4", "B8")    # training bands (matches 02a)
band_prefixes <- c("B03", "B04", "B08") # raster layer name prefixes

dates <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

l2a_ids <- c(
  Aigoual = "L2A_T31TEJ_A031608_20210711T104217",
  Blois   = "L2A_T31TCN_A031222_20210614T105443",
  Mormal  = "L2A_T31UER_A031222_20210614T105443"
)

out_root <- file.path(paths$output, "intermediate", "sm6")
mdl_root <- file.path(paths$output, "intermediate", "PROSAIL_Models")

set.seed(42)

# ── SRF (computed once) ────────────────────────────────────────────────────────

cli::cli_progress_step("Loading Sentinel-2 SRF")
srf             <- prosail::get_srf_sensor("Sentinel_2")
bands_to_select <- match(band_select, srf$Spectral_Bands)

# ── Helper: train or load an ATBD SVR ─────────────────────────────────────────

train_or_load_atbd_svr <- function(rds_path, geom_s2, codistribution_lai,
                                   nb_samples, srf, bands_to_select,
                                   overwrite = FALSE) {

  variant_label <- if (codistribution_lai) "ATBD_T" else "ATBD_F"
  nb_models     <- round(nb_samples / 100L)

  if (file.exists(rds_path) && !overwrite) {
    cli::cli_progress_step(
      "{variant_label}: loading cached SVR ({nb_models} models) from disk"
    )
    svr <- load_svr_ensemble(rds_path)
    return(svr)
  }

  geom_acq <- list(
    min = data.frame(
      tto = geom_s2$MinAngle["vza"],
      tts = geom_s2$MinAngle["sza"],
      psi = geom_s2$MinAngle["psi"]
    ),
    max = data.frame(
      tto = geom_s2$MaxAngle["vza"],
      tts = geom_s2$MaxAngle["sza"],
      psi = geom_s2$MaxAngle["psi"]
    )
  )

  # Step 1: LUT sampling
  cli::cli_progress_step(
    "{variant_label}: sampling ATBD LUT ({nb_samples} rows, codist={codistribution_lai})"
  )
  ip <- prosail::get_atbd_v3_lut_input(
    nb_samples         = nb_samples,
    geom_acq           = geom_acq,
    codistribution_lai = codistribution_lai
  )

  # Step 2: PROSAIL simulation + SRF convolution + noise
  cli::cli_progress_step(
    "{variant_label}: running PROSAIL simulation + SRF + noise"
  )
  brf <- build_prosail_lut(ip, srf, bands_to_select)

  # Step 3: SVR training (nb_models ensemble members)
  cli::cli_progress_step(
    "{variant_label}: training SVR ensemble ({nb_models} models on {nb_samples} samples)"
  )
  svr_ser <- train_svr_ensemble(brf, ip$lai, nb_samples)

  # Step 4: save
  cli::cli_progress_step("{variant_label}: saving SVR to disk")
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(svr_ser, rds_path)

  lapply(svr_ser, liquidSVM::unserialize.liquidSVM)
}

# ── Main loop: per site ────────────────────────────────────────────────────────

t0 <- proc.time()

for (i_site in seq_along(sites)) {

  site <- sites[[i_site]]
  cli::cli_h1("[{i_site}/{length(sites)}] {site}")

  site_out_dir <- file.path(out_root, site)
  dir.create(site_out_dir, recursive = TRUE, showWarnings = FALSE)
  site_mdl_dir <- file.path(mdl_root, site, "LIDFa_lai_LMA_BROWN", "atbd")

  # Step 1: geometry
  cli::cli_progress_step("{site} [1/6]: loading S2 acquisition geometry")
  path_angles <- file.path(paths$ext_results, site, "PROSAIL_Optimization", "geomAcq_S2")
  path_bbox   <- file.path(paths$raw_data, site, "Geo_Files", "aoi_bbox.GPKG")
  geom_s2     <- get_s2_angles(path_angles = path_angles,
                                dateAcq    = dates[[site]],
                                path_bbox  = path_bbox)

  # Step 2: S2 reflectance raster
  l2a     <- l2a_ids[[site]]
  s2_path <- file.path(paths$ext_results, site, "PROSAIL_Optimization", l2a,
                        "Reflectance", paste0(l2a, "_Refl"))
  if (!file.exists(s2_path)) stop("S2 raster not found:\n  ", s2_path)
  cli::cli_progress_step("{site} [2/6]: loading S2 reflectance raster")
  s2_rast  <- terra::rast(s2_path)
  n_pixels <- terra::ncell(s2_rast)

  # Step 3: mask
  mask_path <- file.path(
    paths$ext_results, site, "LiDAR", "Heterogeneity_Masks",
    "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"
  )
  cli::cli_progress_step("{site} [3/6]: loading deciduous mask")
  mask_rast <- if (file.exists(mask_path)) {
    terra::rast(mask_path)
  } else {
    cli::cli_alert_warning("{site}: mask not found — output will not be masked")
    NULL
  }

  # Step 4: ATBD_T (train or load + predict — masked + unmasked)
  rds_T <- file.path(site_mdl_dir, "LIDFa=1_lai=1_LMA=1_BROWN=1_codistT.rds")
  svr_T <- train_or_load_atbd_svr(rds_T, geom_s2,
                                   codistribution_lai = TRUE,
                                   nb_samples, srf, bands_to_select,
                                   overwrite = TRUE)

  cli::cli_progress_step(
    "{site} [4/7]: ATBD_T predicting LAI on {format(n_pixels, big.mark=',')} pixels"
  )
  lai_T_raw <- predict_lai_full_raster(svr_T, s2_rast, band_prefixes, NULL)
  out_T_raw <- file.path(site_out_dir, "s2lai_summer_atbd_T_raw_res_10_m.tif")
  terra::writeRaster(lai_T_raw, out_T_raw, overwrite = TRUE, gdal = "COMPRESS=LZW")

  lai_T <- if (!is.null(mask_rast)) terra::mask(lai_T_raw, mask_rast) else lai_T_raw
  out_T <- file.path(site_out_dir, "s2lai_summer_atbd_T_res_10_m.tif")
  terra::writeRaster(lai_T, out_T, overwrite = TRUE, gdal = "COMPRESS=LZW")
  rm(svr_T, lai_T, lai_T_raw); gc()

  # Step 5: ATBD_F (train or load + predict — masked + unmasked)
  rds_F <- file.path(site_mdl_dir, "LIDFa=1_lai=1_LMA=1_BROWN=1.rds")
  svr_F <- train_or_load_atbd_svr(rds_F, geom_s2,
                                   codistribution_lai = FALSE,
                                   nb_samples, srf, bands_to_select,
                                   overwrite = TRUE)

  cli::cli_progress_step(
    "{site} [5/7]: ATBD_F predicting LAI on {format(n_pixels, big.mark=',')} pixels"
  )
  lai_F_raw <- predict_lai_full_raster(svr_F, s2_rast, band_prefixes, NULL)
  out_F_raw <- file.path(site_out_dir, "s2lai_summer_atbd_F_raw_res_10_m.tif")
  terra::writeRaster(lai_F_raw, out_F_raw, overwrite = TRUE, gdal = "COMPRESS=LZW")

  lai_F <- if (!is.null(mask_rast)) terra::mask(lai_F_raw, mask_rast) else lai_F_raw
  out_F <- file.path(site_out_dir, "s2lai_summer_atbd_F_res_10_m.tif")
  terra::writeRaster(lai_F, out_F, overwrite = TRUE, gdal = "COMPRESS=LZW")
  rm(svr_F, lai_F, lai_F_raw); gc()

  # Step 6: finalise (forces previous step to print its ✔ with timing)
  cli::cli_progress_step("{site} [6/7]: writing outputs")
  cli::cli_alert_success("  ATBD_T masked : {out_T}")
  cli::cli_alert_success("  ATBD_T raw    : {out_T_raw}")
  cli::cli_alert_success("  ATBD_F masked : {out_F}")
  cli::cli_alert_success("  ATBD_F raw    : {out_F_raw}")
  cli::cli_progress_step("{site} [7/7]: done")
}

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
cli::cli_h1("Done — {elapsed} s total")
cli::cli_bullets(c(
  "v" = "Outputs: revision/output/intermediate/sm6/{{site}}/s2lai_summer_atbd_{{T,F}}_res_10_m.tif",
  ">" = "Next: source('revision/scripts/24_plot_s2_atbd_validation.R')"
))
