# ---
# title:  08_train_prosail_full.R
# desc:   Orchestration — trains all 270 SVR configurations (LIDFa/lai/LMA/BROWN
#         grid) under three LAI scenarios, for each of the three study sites.
#
#         The two scenarios differ only in the LAI source used for the
#         LiDAR_LAI_Best_Site_Depth variant of the simulation grid (lai index 6):
#
#           per_site  : LAI_ALS_dopt per site's own Pareto d_opt
#           common    : LAI_ALS_dopt at the combined-sites Pareto d_opt
#
#         LAI_ALS_dopt rasters are read from script 07 output.
#         The LiDAR_LAI variant (lai index 5) always uses the full-canopy
#         cross-site pool, unchanged across scenarios.
#
# Prerequisites:
#   02a — simulation strategy RDS (written once, shared with 02a)
#   07  — LAI_ALS_dopt rasters at
#         revision/output/intermediate/lai_als_dopt/{site}/LAI_ALS_dopt_{scenario}.tif
#   S2 geometry rasters at 03_RESULTS/{site}/PROSAIL_Optimization/geomAcq_S2/
#
# Output: 270 RDS per site × scenario at:
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/{scenario}/
#   (scenario ∈ {per_site, common})
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/08_train_prosail_full.R")
# ---

library("here")
library("terra")
library("progressr")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "prosail_lut.R"))
source(here::here("revision", "R", "get_s2_angles.R"))

# global = TRUE propagates the handler into nested with_progress() calls
# (e.g. the one inside train_all_simulations in prosail_lut.R)
progressr::handlers(progressr::handler_txtprogressbar(width = 60), global = TRUE)

# ── Global parameters ─────────────────────────────────────────────────────────

sites    <- c("Aigoual", "Blois", "Mormal")
dates    <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

parms2test             <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy          <- "LIDFa_lai_LMA_BROWN"
nbSamples_train        <- 2000
S2BandSelect           <- c("B3", "B4", "B8")
nbCPU                  <- 1L
nbValuesPerSiteForPool <- 5000

lai_scenarios <- c("common")

# k / scan-angle rescaling — must match 07 and 11.
# scale_factor = (k_ref / k_select) × cos(theta_select)
# applied to the full-canopy ladstack sum pool (computed at k_ref=0.5, theta=0°).
# lai_dopt_values_vec from 07 rasters is already scaled — no second scaling needed.
k_ref        <- 0.5
k_select     <- 0.6
theta_select <- 0    # degrees from nadir (scan angle correction handled separately)
scale_factor <- (k_ref / k_select) * cos(theta_select * pi / 180)

set.seed(42)

# ── Simulation grid ───────────────────────────────────────────────────────────
# Strategy RDS was already written by 02a; overwrite = FALSE would stop here.
# Pass overwrite = TRUE so 08 can be re-run standalone after grid changes.

strategy_dir <- file.path(paths$output, "intermediate", "PROSAIL_Models")

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

# ── Phase A — Cross-site full-canopy LAI pool (common to all scenarios) ───────
# Random draws from the full ladstack sum, masked to coherent forest pixels.
# Provides values for the LiDAR_LAI variant (lai index 5).

lai_values_all <- list()

for (site in sites) {
  lidar_dir <- file.path(paths$prosail_lidar, site, "LiDAR")

  lidar_lai_rast <- sum(terra::rast(
    file.path(lidar_dir, "PAD_Profiles_Classic", "ladstack.tif")
  ), na.rm = TRUE)
  max_rast <- terra::rast(file.path(lidar_dir, "max_res_10_m.tif"))

  mask_rast      <- (max_rast >= 10 & max_rast <= 40 & lidar_lai_rast >= 2)
  lidar_lai_rast <- terra::mask(lidar_lai_rast, mask_rast, maskvalues = 0)

  lai_pool <- terra::values(lidar_lai_rast, na.rm = TRUE)
  lai_values_all[[site]] <- sample(
    lai_pool,
    size    = nbValuesPerSiteForPool,
    replace = length(lai_pool) < nbValuesPerSiteForPool
  )
}

lai_values_vec <- unlist(lai_values_all, use.names = FALSE) * scale_factor
stopifnot(length(lai_values_vec) == length(sites) * nbValuesPerSiteForPool)
cat("Full-canopy cross-site pool:", length(lai_values_vec), "values (scale_factor=",
    round(scale_factor, 4), ")\n")

# ── LAI_ALS_dopt raster root ──────────────────────────────────────────────────

dopt_rast_root <- file.path(paths$output, "intermediate", "lai_als_dopt")

# ── Outer loop: LAI scenario → Phase A dopt pool → Phase B training ───────────

for (lai_scenario in lai_scenarios) {

  cat("\n══ LAI scenario:", lai_scenario, "══\n")
  t_scenario <- Sys.time()

  # Phase A (scenario-specific) — cross-site LAI_ALS_dopt pool
  # Provides values for the LiDAR_LAI_Best_Site_Depth variant (lai index 6).
  lai_dopt_values_all <- list()

  for (site in sites) {
    dopt_tif <- file.path(dopt_rast_root, site,
                           paste0("LAI_ALS_dopt_", lai_scenario, ".tif"))
    if (!file.exists(dopt_tif))
      stop("LAI_ALS_dopt raster not found — run script 07 first:\n  ", dopt_tif)

    lidar_dir <- file.path(paths$prosail_lidar, site, "LiDAR")
    max_rast  <- terra::rast(file.path(lidar_dir, "max_res_10_m.tif"))

    lai_dopt_rast <- terra::rast(dopt_tif)

    # Apply same coherence mask as Phase A full-canopy pool
    full_lai_rast <- sum(terra::rast(
      file.path(lidar_dir, "PAD_Profiles_Classic", "ladstack.tif")
    ), na.rm = TRUE)
    mask_rast     <- (max_rast >= 10 & max_rast <= 40 & full_lai_rast >= 2)
    lai_dopt_rast <- terra::mask(lai_dopt_rast, mask_rast, maskvalues = 0)

    dopt_pool <- terra::values(lai_dopt_rast, na.rm = TRUE)
    lai_dopt_values_all[[site]] <- sample(
      dopt_pool,
      size    = nbValuesPerSiteForPool,
      replace = length(dopt_pool) < nbValuesPerSiteForPool
    )
  }

  lai_dopt_values_vec <- unlist(lai_dopt_values_all, use.names = FALSE)
  stopifnot(length(lai_dopt_values_vec) == length(sites) * nbValuesPerSiteForPool)
  cat("LAI_dopt cross-site pool:", length(lai_dopt_values_vec), "values\n")

  # Phase B — per-site training with scenario-specific dopt pool
  n_sites   <- length(sites)
  t_phase_b <- Sys.time()

  for (i_site in seq_along(sites)) {

    site    <- sites[[i_site]]
    t_start <- Sys.time()
    cat(sprintf("\n[%d/%d] Training site: %s | scenario: %s (%d models)\n",
                i_site, n_sites, site, lai_scenario, nrow(simulations$grid_simu)))

    models_dir   <- here::here(
      "revision", "output", "intermediate", "PROSAIL_Models",
      site, name_strategy, lai_scenario
    )
    filename_svr <- file.path(models_dir, paste0(combination_labels, ".rds"))

    if (!dir.exists(models_dir))
      dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

    train_all_simulations(
      simulations         = simulations,
      geom_s2             = geom_s2[[site]],
      lidar_lai_rast      = lai_values_vec,      # full-canopy pool (variant 5)
      lidar_lai_site_rast = lai_dopt_values_vec, # d_opt pool (variant 6)
      nbSamples           = nbSamples_train,
      filename            = filename_svr,
      nbCPU               = nbCPU,
      S2BandSelect        = S2BandSelect
    )

    elapsed_site <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
    elapsed_tot  <- round(difftime(Sys.time(), t_phase_b, units = "mins"), 1)
    cat(sprintf("  Done in %.1f min  (cumul: %.1f min)\n", elapsed_site, elapsed_tot))
  }

  cat("Scenario", lai_scenario, "done —",
      round(difftime(Sys.time(), t_scenario, units = "mins"), 1), "min\n")
}

cat("\n── 08 done — all scenarios trained ──\n")
