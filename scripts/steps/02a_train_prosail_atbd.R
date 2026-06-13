# ---
# title:  02a_train_prosail_atbd.R
# desc:   Orchestration — trains ONE SVR model per site (the ATBD baseline
#         configuration: LIDFa=1, lai=1, LMA=1, BROWN=1). The ATBD row uses
#         the "ATBD" distribution for every parameter, so no LiDAR rasters
#         are sampled — only geometry and the LUT baseline are needed.
#         This first-pass model is used downstream only to determine d_opt
#         (script 06) — it is NOT the final inversion model.
#
# Output: one RDS per site at:
#   output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/atbd/LIDFa=1_lai=1_LMA=1_BROWN=1.rds
#
# Run from the project root (NC_Full/):
#   source("scripts/steps/02a_train_prosail_atbd.R")
# ---

library("here")

source(here::here("R", "paths.R"))
source(here::here("R", "prosail_lut.R"))
source(here::here("R", "get_s2_angles.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites    <- c("Aigoual", "Blois", "Mormal")
dates    <- c(Aigoual = "2021-07-11", Blois = "2021-06-14", Mormal = "2021-06-14")

parms2test      <- c("LIDFa", "lai", "LMA", "BROWN")
name_strategy   <- "LIDFa_lai_LMA_BROWN"
nbSamples_train <- 2000
S2BandSelect    <- c("B3", "B4", "B8")

# ── Simulation grid — keep only the ATBD row (all params at index 1) ──────────

strategy_dir <- file.path(paths$output, "intermediate", "PROSAIL_Models")

simulations <- get_parameter_combinations(
  parms2test    = parms2test,
  output_dir    = strategy_dir,
  name_strategy = name_strategy,
  overwrite     = TRUE
)

atbd_row   <- simulations$grid_simu[1L, , drop = FALSE]
atbd_label <- paste(names(atbd_row), atbd_row, sep = "=", collapse = "_")

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

# ── Per-site ATBD SVR training ────────────────────────────────────────────────
# All parameters are at ATBD index 1 → no LiDAR raster sampling.

for (site in sites) {

  cat("\n── Training ATBD SVR for site:", site, "──\n")
  t0_site <- Sys.time()

  models_dir <- file.path(paths$output, "intermediate", "PROSAIL_Models",
    site, name_strategy, "atbd"
  )
  if (!dir.exists(models_dir))
    dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

  filename_svr <- file.path(models_dir, paste0(atbd_label, ".rds"))

  train_all_simulations(
    simulations  = list(grid_simu = atbd_row, combination = simulations$combination),
    geom_s2      = geom_s2[[site]],
    nbSamples    = nbSamples_train,
    filename     = filename_svr,
    nbCPU        = 1L,
    S2BandSelect = S2BandSelect
  )

  cat("Model saved to:", models_dir, "\n")
  cat("Elapsed:", round(difftime(Sys.time(), t0_site, units = "mins"), 1), "min\n")
}

cat("\n── 02a done — ATBD SVR models ready for 04a ──\n")
