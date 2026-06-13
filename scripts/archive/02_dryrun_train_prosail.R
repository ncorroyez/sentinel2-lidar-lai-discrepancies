# ---
# title:  02_dryrun_train_prosail.R
# desc:   Dry run of 02_train_prosail_lut.R — 1 site (Blois), 5 first
#         combinations only. Validates prosail 3.0.0 API and estimates
#         full training time before committing to the 810-model run.
#
# Run from the project root (NC_Full/):
#   source("scripts/02_dryrun_train_prosail.R")
# ---

source(here::here("R", "paths.R"))
source(here::here("R", "prosail_lut.R"))
source(file.path(paths$prosail_codes, "libraries", "get_s2_angles.R"))

# ── Dry run parameters ────────────────────────────────────────────────────────

site_dryrun    <- "Blois"
date_dryrun    <- "2021-06-14"
n_combos       <- 5          # only first 5 of 270
nbSamples_train          <- 5000
S2BandSelect             <- c("B3", "B4", "B8")
nbCPU                    <- 1
nbValuesPerSiteForPool   <- 5000
parms2test               <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy            <- "LIDFa_lai_LMA_BROWN"

set.seed(42)

# ── Simulation grid ───────────────────────────────────────────────────────────

strategy_dir <- file.path(paths$output, "intermediate", "PROSAIL_Models")

simulations <- get_parameter_combinations(
  parms2test    = parms2test,
  output_dir    = strategy_dir,
  name_strategy = name_strategy,
  overwrite     = TRUE
)

cat("Full grid:", nrow(simulations$grid_simu), "combinations\n")
cat("Dry run: first", n_combos, "combinations on site", site_dryrun, "\n\n")

combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

# Subset to first n_combos
simulations_dry <- simulations
simulations_dry$grid_simu <- simulations$grid_simu[seq_len(n_combos), , drop = FALSE]
combination_labels_dry <- combination_labels[seq_len(n_combos)]

cat("Combinations to train:\n")
cat(paste(" -", combination_labels_dry), sep = "\n")
cat("\n")

# ── Geometry of acquisition ───────────────────────────────────────────────────

path_angles <- file.path(paths$ext_results, site_dryrun, "PROSAIL_Optimization", "geomAcq_S2")
path_bbox   <- file.path(paths$raw_data, site_dryrun, "Geo_Files", "aoi_bbox.GPKG")
geom_s2     <- get_s2_angles(
  path_angles = path_angles,
  dateAcq     = date_dryrun,
  path_bbox   = path_bbox
)
cat("Geometry loaded for", site_dryrun, "\n")

# ── Phase A — cross-site LAI pool (Blois only for dry run) ───────────────────

lidar_dir <- file.path(paths$prosail_lidar, site_dryrun, "LiDAR")

lidar_lai_rast    <- sum(terra::rast(
  file.path(lidar_dir, "PAD_Profiles_Classic", "ladstack.tif")
), na.rm = TRUE)
lidar_lai_4m_rast <- terra::rast(file.path(lidar_dir, "PAD_4m.tif"))
max_rast          <- terra::rast(file.path(lidar_dir, "max_res_10_m.tif"))

mask_rast         <- (max_rast >= 10 & max_rast <= 40 & lidar_lai_rast >= 2)
lidar_lai_rast    <- terra::mask(lidar_lai_rast,    mask_rast, maskvalues = 0)
lidar_lai_4m_rast <- terra::mask(lidar_lai_4m_rast, mask_rast, maskvalues = 0)

lai_values_vec   <- sample_lai_stratified_uniform(
  terra::values(lidar_lai_rast,    na.rm = TRUE), n_samples = nbValuesPerSiteForPool
)
lai4m_values_vec <- sample_lai_stratified_uniform(
  terra::values(lidar_lai_4m_rast, na.rm = TRUE), n_samples = nbValuesPerSiteForPool
)

cat("LAI pool (single site, dry run):",
    length(lai_values_vec), "values\n\n")

# ── Phase B — train 5 SVR models ─────────────────────────────────────────────

lidar_lai_common_rast <- terra::rast(
  file.path(lidar_dir, "lidar_lai_best_common_depth.tif")
)
lidar_lai_3m_rast <- terra::rast(
  file.path(lidar_dir, "lidar_lai_3m.tif")
)

models_dir <- here::here(
  "output", "intermediate", "PROSAIL_Models",
  site_dryrun, paste0(name_strategy, "_dryrun")
)
if (!dir.exists(models_dir))
  dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

filename_svr <- file.path(models_dir, paste0(combination_labels_dry, ".rds"))

t_start <- Sys.time()

train_all_simulations(
  simulations           = simulations_dry,
  geom_s2               = geom_s2,
  lidar_lai_rast        = lai_values_vec,
  lidar_lai_site_rast   = lai4m_values_vec,
  lidar_lai_common_rast = lidar_lai_common_rast,
  lidar_lai_3m_rast     = lidar_lai_3m_rast,
  nbSamples             = nbSamples_train,
  filename              = filename_svr,
  nbCPU                 = nbCPU,
  S2BandSelect          = S2BandSelect
)

elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

# ── Time estimate ─────────────────────────────────────────────────────────────

n_total   <- nrow(simulations$grid_simu) * length(c("Aigoual", "Blois", "Mormal"))
time_full <- elapsed_min / n_combos * n_total

cat("\n── Dry run complete ──\n")
cat(sprintf("  %d models in %.1f min (%.1f s/model)\n",
            n_combos, elapsed_min, elapsed_min * 60 / n_combos))
cat(sprintf("  Estimated full run (%d models × 3 sites = %d): %.0f min (~%.1f h)\n",
            nrow(simulations$grid_simu), n_total, time_full, time_full / 60))
cat("\nDry run models saved to:", models_dir, "\n")
cat("Check the .rds files exist and are valid before running 02_train_prosail_lut.R\n")
