# ---
# title:  01_prepare.R
# desc:   Phase 1 — Preparation: trains the ATBD SVR baseline, samples S2
#         pixels, extracts LiDAR PAD values at sample positions, and applies
#         the ATBD SVR to produce first-pass LAI_S2_atbd estimates.
#
#         Steps (in order):
#           02a  train_prosail_atbd  — one SVR per site (ATBD config)
#           03a  sample_s2_pixels   — stratified-uniform S2 pixel sampling
#           03b  extract_lidar      — extract PAD at sample positions
#           04a  apply_prosail_atbd — LAI_S2_atbd from ATBD SVR
#
#         Outputs:
#           output/intermediate/PROSAIL_Models/{site}/
#             LIDFa_lai_LMA_BROWN/atbd/LIDFa=1_lai=1_LMA=1_BROWN=1.rds
#           03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#             Sampling_*.GPKG, S2_reflectance_*.csv, PAD_*.csv,
#             LAI_estimated_atbd_*.csv
#
# Run from project root (NC_Full/):
#   source("scripts/01_prepare.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[01] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("scripts/steps/02a_train_prosail_atbd.R")
step("scripts/steps/03a_sample_s2_pixels.R")
step("scripts/steps/03b_extract_lidar_at_samples.R")
step("scripts/steps/04a_apply_prosail_atbd.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("01_prepare done — {elapsed_total} min")
