# ---
# title:  02_dopt_compute.R
# desc:   Phase 2 — d_opt selection: computes ATBD metrics at all depths,
#         selects the optimal canopy integration depth (d_opt) via Pareto
#         multi-criteria method, and produces LAI_ALS_dopt rasters.
#
#         Steps (in order):
#           05a  compute_metrics_atbd — R², RMSE, Bias, Slope per depth
#           05c  k_zmin_sensitivity   — sweep (k, theta) → kt_sensitivity CSV
#                                       (needed by 06 when k_select ≠ 0.5)
#           06   select_dopt          — Pareto d_opt → dopt_reference.csv
#           07   compute_lai_als_dopt — LAI_ALS_dopt GeoTIFFs per site × scenario
#
#         Outputs:
#           output/intermediate/sm5/
#             all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#             dopt_reference.csv, dopt_all.csv
#           output/intermediate/lai_als_dopt/{site}/
#             LAI_ALS_dopt_{per_site,common}.tif
#
# Run from project root (NC_Full/):
#   source("scripts/02_dopt_compute.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[02] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

# step("scripts/steps/05a_sm5_compute_metrics_atbd.R")
step("scripts/steps/05c_sm5_k_zmin_sensitivity_dopt.R")
step("scripts/steps/06_sm5_select_dopt.R")
step("scripts/steps/07_compute_lai_als_dopt.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("02_dopt_compute done — {elapsed_total} min")
