# ---
# title:  02_train_prosail_lut.R
# desc:   Orchestration — defines the PROSAIL simulation grid (270 combinations
#         of LIDFa/lai/LMA/BROWN distributions) and trains one SVR ensemble per
#         combination for each of the three study sites (Aigoual, Blois, Mormal).
#
#         Reproduces faithfully the training section of
#         Main_s2_03_PROSAIL_Inversion_bis.R (lines 15-150), refactored with
#         here::here() paths and revision/R/prosail_lut.R. The inversion section
#         of that script (lines 151 onwards) is the scope of SM4.
#
# Prerequisite: Sentinel-2 geometry-of-acquisition rasters must be present at
#   03_RESULTS/{site}/PROSAIL_Optimization/geomAcq_S2/
#   (produced by 01b_download_s2_from_cdse.R via preprocS2::get_s2_raster()).
#
# Prerequisite: LiDAR rasters must be present at
#   01_DATA/{site}/LiDAR/PAD_Profiles_Classic/ladstack.tif   (Phase A + Phase B)
#   01_DATA/{site}/LiDAR/PAD_4m.tif                         (Phase A)
#   01_DATA/{site}/LiDAR/max_res_10_m.tif                   (Phase A mask)
#   01_DATA/{site}/LiDAR/lidar_lai_best_common_depth.tif     (Phase B — used in production)
#   01_DATA/{site}/LiDAR/lidar_lai_3m.tif                   (Phase B — explicit here, was
#                                                             injected manually in legacy workflow)
#
# NOT read (dead code in Main_s2_03_PROSAIL_Inversion_bis.R):
#   01_DATA/{site}/LiDAR/lidar_lai_best_site_depth.tif  (read at line 134 of _03_bis_
#                                                         but overridden by lai4m_values_vec)
#
# Output: one RDS per combination per site, at:
#   revision/output/intermediate/PROSAIL_Models/{site}/{name_strategy}/
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/02_train_prosail_lut.R")
# ---

source(here::here("revision", "R", "prosail_lut.R"))

# TODO: move get_s2_angles into revision/R/ once the S2-angles module is refactored.
source(here::here("PROSAIL-Optimization", "02_CODES", "libraries", "get_s2_angles.R"))

# ── Global parameters ─────────────────────────────────────────────────────────

sites    <- c("Aigoual", "Blois", "Mormal")
dates    <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

parms2test              <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy           <- "LIDFa_lai_LMA_BROWN"
nbSamples_train         <- 1000
S2BandSelect            <- c("B3", "B4", "B8")
nbCPU                   <- 1
nbValuesPerSiteForPool  <- 10000   # values drawn per site in the cross-site pool (Phase A)

# set.seed position mirrors Main_s2_03_PROSAIL_Inversion_bis.R line 18:
# immediately after parameter definition, before any stochastic operation.
set.seed(42)

# ── Simulation grid (same for all sites) ──────────────────────────────────────

strategy_dir <- here::here("revision", "output", "intermediate", "PROSAIL_Models")

simulations <- get_parameter_combinations(
  parms2test    = parms2test,
  output_dir    = strategy_dir,
  name_strategy = name_strategy,
  overwrite     = TRUE
)

cat("Simulation grid:", nrow(simulations$grid_simu), "combinations\n")

combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

# ── Geometry of acquisition (all sites, before training loop) ─────────────────

geom_s2 <- list()
for (site in sites) {
  path_angles <- here::here("03_RESULTS", site, "PROSAIL_Optimization", "geomAcq_S2")
  path_bbox   <- here::here("01_DATA", site, "Geo_Files", "aoi_bbox.GPKG")
  geom_s2[[site]] <- get_s2_angles(
    path_angles = path_angles,
    dateAcq     = dates[[site]],
    path_bbox   = path_bbox
  )
}

# ── Phase A — Cross-site LAI pool ─────────────────────────────────────────────
# Reproduces Main_s2_03_PROSAIL_Inversion_bis.R lines 90-122.
# For each site: read ladstack + PAD_4m, apply coherence mask
# (max_height in [10,40] AND lidar_lai >= 2), draw nbValuesPerSiteForPool values.
# Then concatenate into two cross-site vectors of length 3 * nbValuesPerSiteForPool.
#
# CRITICAL: the order of the two sample() calls inside the loop is preserved
# exactly (lidar_lai_rast first, lidar_lai_4m_rast second) so that the RNG
# state matches _03_bis_ and the trained SVR models are byte-identical to those
# of the production pipeline.

lai_values_all   <- list()
lai4m_values_all <- list()

for (site in sites) {

  lidar_lai_path <- here::here("01_DATA", site, "LiDAR", "PAD_Profiles_Classic",
                                "ladstack.tif")

  # Sum layers of ladstack to get total LAI per pixel
  lidar_lai_rast    <- sum(terra::rast(lidar_lai_path), na.rm = TRUE)

  # 4 m depth raster
  lidar_lai_4m_rast <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "PAD_4m.tif")
  )

  max_rast  <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "max_res_10_m.tif")
  )

  # Coherence mask: canopy height in (10, 40] m AND LAI >= 2
  mask_rast <- (max_rast >= 10 & max_rast <= 40 & lidar_lai_rast >= 2)

  # Apply mask (FALSE cells → NA)
  lidar_lai_rast    <- terra::mask(lidar_lai_rast,    mask_rast, maskvalues = 0)
  lidar_lai_4m_rast <- terra::mask(lidar_lai_4m_rast, mask_rast, maskvalues = 0)
  max_rast          <- terra::mask(max_rast,           mask_rast, maskvalues = 0)

  # Draw random sample of pixel values (order of calls preserved for RNG parity)
  lai_values_all[[site]]   <- sample(terra::values(lidar_lai_rast,    na.rm = TRUE),
                                      nbValuesPerSiteForPool)
  lai4m_values_all[[site]] <- sample(terra::values(lidar_lai_4m_rast, na.rm = TRUE),
                                      nbValuesPerSiteForPool)
}

# Concatenate into flat cross-site vectors
lai_values_vec   <- unlist(lai_values_all,   use.names = FALSE)
lai4m_values_vec <- unlist(lai4m_values_all, use.names = FALSE)

# Sanity check: each vector must contain exactly 3 * nbValuesPerSiteForPool values
stopifnot(length(lai_values_vec)   == length(sites) * nbValuesPerSiteForPool)
stopifnot(length(lai4m_values_vec) == length(sites) * nbValuesPerSiteForPool)

cat("Cross-site pool built:",
    length(lai_values_vec),   "values in lai_values_vec,",
    length(lai4m_values_vec), "in lai4m_values_vec\n")

# ── Phase B — Per-site SVR training ───────────────────────────────────────────
# Reproduces Main_s2_03_PROSAIL_Inversion_bis.R lines 126-149.
# lidar_lai_rast and lidar_lai_site_rast are the cross-site pools built in Phase A
# (same vectors for all three sites — matching production behaviour).
# lidar_lai_common_rast and lidar_lai_3m_rast are site-local SpatRasters.
#
# Dead code from _03_bis_ deliberately omitted:
#   - lidar_lai_rast local re-read  (line 130-132, overridden by lai_values_vec)
#   - lidar_lai_site_rast           (line 134,     overridden by lai4m_values_vec)
#   - lidar_lai_4m_rast re-read     (line 136,     never used after assignment)

for (site in sites) {

  cat("\n── Training PROSAIL for site:", site, "──\n")
  t_start <- Sys.time()

  # Site-local rasters passed as lidar_lai_common_rast and lidar_lai_3m_rast
  lidar_lai_common_rast <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "lidar_lai_best_common_depth.tif")
  )

  # NOTE: this line was missing from the original
  # Main_s2_03_PROSAIL_Inversion_bis.R. In the legacy workflow,
  # lidar_lai_3m_rast was defined interactively in the R session before
  # sourcing the script. This refactor makes the definition explicit for
  # reproducibility.
  lidar_lai_3m_rast <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "lidar_lai_3m.tif")
  )

  # Output directory for SVR model RDS files
  models_dir <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy
  )
  if (!dir.exists(models_dir))
    dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

  filename_svr <- file.path(models_dir, paste0(combination_labels, ".rds"))

  # Train all 270 SVR models for this site
  train_all_simulations(
    simulations           = simulations,
    geom_s2               = geom_s2[[site]],
    lidar_lai_rast        = lai_values_vec,      # cross-site pool (Phase A)
    lidar_lai_site_rast   = lai4m_values_vec,    # cross-site pool PAD_4m (Phase A)
    lidar_lai_common_rast = lidar_lai_common_rast,
    lidar_lai_3m_rast     = lidar_lai_3m_rast,
    nbSamples             = nbSamples_train,
    filename              = filename_svr,
    nbCPU                 = nbCPU,
    S2BandSelect          = S2BandSelect
  )

  cat("Models saved to:", models_dir, "\n")
  cat("Elapsed:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")
}

cat("\n── All sites done ──\n")
