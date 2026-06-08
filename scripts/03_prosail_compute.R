# ---
# title:  03_prosail_compute.R
# desc:   Phase 3 — PROSAIL full inversion: trains 270 SVR configurations,
#         applies them to S2 reflectances, computes metrics, selects the
#         Pareto-optimal configuration, predicts LAI_S2_opt rasters, and
#         produces ATBD_T / ATBD_F validation rasters (prosail 3.0.0).
#
#         Steps (in order):
#           09  train_prosail_full  — 270 SVR per site × scenario
#           10  apply_prosail_full  — LAI_S2 predictions for all configs
#           11  compute_metrics_full — R², RMSE, Bias, Slope for all configs
#           12  select_prosail_opt  — Pareto selection → prosail_opt.csv
#           13  predict_lai_raster  — LAI_S2_opt GeoTIFFs per site × scenario
#           25  train_apply_atbd_rasters — ATBD_T / ATBD_F full rasters
#
#         Outputs:
#           revision/output/intermediate/PROSAIL_Models/{site}/
#             LIDFa_lai_LMA_BROWN/{scenario}/*.rds
#           revision/output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv, prosail_opt.csv
#           revision/output/intermediate/sm6/{site}/
#             s2lai_*_{scenario}_res_10_m.tif, s2lai_*_atbd_*_res_10_m.tif
#
# Run from project root (NC_Full/):
#   source("revision/scripts/03_prosail_compute.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[03] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("revision/scripts/steps/09_train_prosail_full.R")
step("revision/scripts/steps/10_apply_prosail_full.R")
step("revision/scripts/steps/11_sm5_compute_metrics_full.R")
step("revision/scripts/steps/12_sm5_select_prosail_opt.R")
step("revision/scripts/steps/13_sm5_predict_lai_raster.R")
step("revision/scripts/steps/25_train_apply_prosail_atbd_rasters.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("03_prosail_compute done — {elapsed_total} min")
