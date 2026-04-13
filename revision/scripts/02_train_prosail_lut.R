# ---
# title:  02_train_prosail_lut.R
# desc:   Orchestration — defines the PROSAIL simulation grid (270 combinations
#         of LIDFa/lai/LMA/BROWN distributions) and trains one SVR ensemble per
#         combination for each of the three study sites (Aigoual, Blois, Mormal).
#
#         Mirrors PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_03_PROSAIL_Inversion.R
#         (training section only, lines 36-113) refactored with here::here()
#         paths and revision/R/prosail_lut.R.
#
# Prerequisite: Sentinel-2 geometry-of-acquisition rasters must be present at
#   03_RESULTS/{site}/PROSAIL_Optimization/geomAcq_S2/
#   (produced by 01b_download_s2_from_cdse.R via preprocS2::get_s2_raster()).
#
# Prerequisite: LiDAR LAI rasters must be present at
#   01_DATA/{site}/LiDAR/PAD_Profiles_Classic/ladstack.tif
#   01_DATA/{site}/LiDAR/lidar_lai_best_site_depth.tif
#   01_DATA/{site}/LiDAR/lidar_lai_best_common_depth.tif   (dead code for 270-combo grid)
#   01_DATA/{site}/LiDAR/lidar_lai_3m.tif                  (dead code for 270-combo grid)
#
# Output: one RDS per combination per site, at:
#   revision/output/intermediate/PROSAIL_Models/{site}/{name_strategy}/
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/02_train_prosail_lut.R")
# ---

library("here")
library("terra")
library("prosail")
library("liquidSVM")
library("truncnorm")
library("progressr")
library("future")
library("future.apply")

source(here::here("revision", "R", "prosail_lut.R"))

# TODO: move get_s2_angles into revision/R/ once the S2-angles module is refactored.
source(here::here("PROSAIL-Optimization", "02_CODES", "libraries", "get_s2_angles.R"))

# ── Global parameters ─────────────────────────────────────────────────────────

set.seed(42)

sites    <- c("Aigoual", "Blois", "Mormal")
dates    <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

parms2test    <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy <- "LIDFa_lai_LMA_BROWN"
nbSamples_train <- 1000
S2BandSelect    <- c("B3", "B4", "B8")
nbCPU           <- 1   # set > 1 for parallel; dry-run with 1 first

# ── Step 1: define simulation grid (same for all sites) ───────────────────────

# Output root for the simulation strategy RDS (site-agnostic here, stored under
# revision/output/intermediate so git ignores large files).
strategy_dir <- here::here("revision", "output", "intermediate", "PROSAIL_Models")

simulations <- get_parameter_combinations(
  parms2test    = parms2test,
  output_dir    = strategy_dir,
  name_strategy = name_strategy,
  overwrite     = TRUE
)

cat("Simulation grid: ", nrow(simulations$grid_simu), "combinations\n")

# Combination labels: used to name individual model RDS files
combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

# ── Step 2: per-site training loop ────────────────────────────────────────────

for (site in sites) {

  cat("\n── Training PROSAIL for site:", site, "──\n")

  # 2a. Output directory for SVR models
  models_dir <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy
  )
  dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

  filename_svr <- file.path(models_dir, paste0(combination_labels, ".rds"))

  # 2b. Load LiDAR LAI rasters
  #
  # lidar_lai_rast / lidar_lai_site_rast: pre-extract values here so the worker
  # function can call sample() on a plain numeric vector (avoids passing a
  # SpatRaster across parallel workers).
  #
  lai_stack_path <- here::here(
    "01_DATA", site, "LiDAR", "PAD_Profiles_Classic", "ladstack.tif"
  )
  lidar_lai_rast <- terra::values(
    sum(terra::rast(lai_stack_path), na.rm = TRUE),
    na.rm = TRUE
  )

  lidar_lai_site_rast <- terra::values(
    terra::rast(
      here::here("01_DATA", site, "LiDAR", "lidar_lai_best_site_depth.tif")
    ),
    na.rm = TRUE
  )

  # These two are dead code for the 270-combination grid (LiDAR_LAI_Best_Common_Depth
  # and LiDAR_LAI_3m? variants are not present in the active get_combination()
  # definitions). Kept as SpatRasters; terra::values() is called inside the worker.
  lidar_lai_common_rast <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "lidar_lai_best_common_depth.tif")
  )
  lidar_lai_3m_rast <- terra::rast(
    here::here("01_DATA", site, "LiDAR", "lidar_lai_3m.tif")
  )

  # 2c. Load Sentinel-2 geometry of acquisition for this site
  path_angles <- here::here(
    "03_RESULTS", site, "PROSAIL_Optimization", "geomAcq_S2"
  )
  path_bbox <- here::here("01_DATA", site, "Geo_Files", "aoi_bbox.GPKG")
  geom_s2   <- get_s2_angles(
    path_angles = path_angles,
    dateAcq     = dates[[site]],
    path_bbox   = path_bbox
  )

  # 2d. Train all 270 SVR models for this site
  train_all_simulations(
    simulations           = simulations,
    geom_s2               = geom_s2,
    lidar_lai_rast        = lidar_lai_rast,
    lidar_lai_site_rast   = lidar_lai_site_rast,
    lidar_lai_common_rast = lidar_lai_common_rast,
    lidar_lai_3m_rast     = lidar_lai_3m_rast,
    nbSamples             = nbSamples_train,
    filename              = filename_svr,
    nbCPU                 = nbCPU,
    S2BandSelect          = S2BandSelect
  )

  cat("Models saved to:", models_dir, "\n")
}

cat("\n── All sites done ──\n")
