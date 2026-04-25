# ---
# title:  02a_train_prosail_atbd.R
# desc:   Orchestration — trains ONE SVR model per site (the ATBD baseline
#         configuration: LIDFa=1, lai=1, LMA=1, BROWN=1) using existing
#         LiDAR LAI rasters. This first-pass model is used downstream only
#         to determine d_opt (script 06) — it is NOT the final inversion model.
#
#         The final 270-configuration models (per LAI scenario) are produced
#         by script 08, which requires d_opt from script 06 and the LAI_ALS_dopt
#         rasters from script 07.
#
# Prerequisite: LiDAR rasters must be present at
#   PROSAIL-Optimization/01_DATA/{site}/LiDAR/
#   (ladstack.tif, PAD_4m.tif, max_res_10_m.tif,
#    lidar_lai_best_common_depth.tif, lidar_lai_3m.tif)
#
# Output: one RDS per site at:
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/atbd/LIDFa=1_lai=1_LMA=1_BROWN=1.rds
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/02a_train_prosail_atbd.R")
# ---

library("here")
library("terra")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "prosail_lut.R"))
source(file.path(paths$prosail_codes, "libraries", "get_s2_angles.R"))

# ── Global parameters ─────────────────────────────────────────────────────────

sites    <- c("Aigoual", "Blois", "Mormal")
dates    <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

parms2test           <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy        <- "LIDFa_lai_LMA_BROWN"
nbSamples_train      <- 5000
S2BandSelect         <- c("B3", "B4", "B8")
nbValuesPerSiteForPool <- 5000

# k / scan-angle rescaling — applied to all LiDAR LAI pools (computed at
# k_ref = 0.5, theta = 0° in the ladstack / PAD rasters).
# Must match 07, 09, 11.
k_ref        <- 0.5
k_select     <- 0.6
theta_select <- 30   # degrees from nadir
scale_factor <- (k_ref / k_select) * cos(theta_select * pi / 180)

set.seed(42)

# ── Simulation grid — build once, keep only ATBD row ─────────────────────────

strategy_dir <- here::here("revision", "output", "intermediate", "PROSAIL_Models")

simulations <- get_parameter_combinations(
  parms2test    = parms2test,
  output_dir    = strategy_dir,
  name_strategy = name_strategy,
  overwrite     = TRUE
)

# ATBD row: all parameters at index 1
atbd_row    <- simulations$grid_simu[1L, , drop = FALSE]
atbd_label  <- paste(names(atbd_row), atbd_row, sep = "=", collapse = "_")

cat("ATBD label:", atbd_label, "\n")

# ── Geometry of acquisition ───────────────────────────────────────────────────

geom_s2 <- list()
for (site in sites) {
  path_angles <- file.path(paths$ext_results, site, "PROSAIL_Optimization", "geomAcq_S2")
  path_bbox   <- file.path(paths$raw_data, site, "Geo_Files", "aoi_bbox.GPKG")
  geom_s2[[site]] <- get_s2_angles(
    path_angles = path_angles,
    dateAcq     = dates[[site]],
    path_bbox   = path_bbox
  )
}

# ── Phase A — Cross-site LAI pool (full canopy + PAD_4m, random draws) ───────

lai_values_all   <- list()
lai4m_values_all <- list()

for (site in sites) {
  lidar_dir <- file.path(paths$prosail_lidar, site, "LiDAR")

  lidar_lai_rast    <- sum(terra::rast(
    file.path(lidar_dir, "PAD_Profiles_Classic", "ladstack.tif")
  ), na.rm = TRUE)
  lidar_lai_4m_rast <- terra::rast(file.path(lidar_dir, "PAD_4m.tif"))
  max_rast          <- terra::rast(file.path(lidar_dir, "max_res_10_m.tif"))

  mask_rast <- (max_rast >= 10 & max_rast <= 40 & lidar_lai_rast >= 2)
  lidar_lai_rast    <- terra::mask(lidar_lai_rast,    mask_rast, maskvalues = 0)
  lidar_lai_4m_rast <- terra::mask(lidar_lai_4m_rast, mask_rast, maskvalues = 0)

  lai_pool   <- terra::values(lidar_lai_rast,    na.rm = TRUE)
  lai4m_pool <- terra::values(lidar_lai_4m_rast, na.rm = TRUE)

  lai_values_all[[site]]   <- sample(lai_pool,
                                     size    = nbValuesPerSiteForPool,
                                     replace = length(lai_pool)   < nbValuesPerSiteForPool)
  lai4m_values_all[[site]] <- sample(lai4m_pool,
                                     size    = nbValuesPerSiteForPool,
                                     replace = length(lai4m_pool) < nbValuesPerSiteForPool)
}

lai_values_vec   <- unlist(lai_values_all,   use.names = FALSE) * scale_factor
lai4m_values_vec <- unlist(lai4m_values_all, use.names = FALSE) * scale_factor

# ── Phase B — Per-site ATBD SVR training ─────────────────────────────────────

for (site in sites) {

  cat("\n── Training ATBD SVR for site:", site, "──\n")
  t_start <- Sys.time()

  lidar_dir_b <- file.path(paths$prosail_lidar, site, "LiDAR")

  lidar_lai_common_rast <- terra::rast(
    file.path(lidar_dir_b, "lidar_lai_best_common_depth.tif")
  ) * scale_factor
  lidar_lai_3m_rast <- terra::rast(
    file.path(lidar_dir_b, "lidar_lai_3m.tif")
  ) * scale_factor

  models_dir <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy, "atbd"
  )
  if (!dir.exists(models_dir))
    dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

  filename_svr <- file.path(models_dir, paste0(atbd_label, ".rds"))

  atbd_simulations <- list(
    grid_simu   = atbd_row,
    combination = simulations$combination
  )

  train_all_simulations(
    simulations           = atbd_simulations,
    geom_s2               = geom_s2[[site]],
    lidar_lai_rast        = lai_values_vec,
    lidar_lai_site_rast   = lai4m_values_vec,
    lidar_lai_common_rast = lidar_lai_common_rast,
    lidar_lai_3m_rast     = lidar_lai_3m_rast,
    nbSamples             = nbSamples_train,
    filename              = filename_svr,
    nbCPU                 = 1L,
    S2BandSelect          = S2BandSelect
  )

  cat("Model saved to:", models_dir, "\n")
  cat("Elapsed:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "min\n")
}

cat("\n── 02a done — ATBD SVR models ready for 04a ──\n")
