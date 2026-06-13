# ---
# title:  04_het_compute.R
# desc:   Phase 4 — Heterogeneity analysis: computes DSM_sd / CHM_sd metrics,
#         runs the main SM6b regression (d_opt_source = all_sites), and sweeps
#         all configurations to identify the best heterogeneity predictor.
#
#         Steps (in order):
#           14  compute_heterogeneity — DSM_sd, CHM_sd rasters (3×3, 10 m)
#           15  sm6_analysis          — main SM6b regression (all_sites d_opt)
#           16  sm6_analysis_sweep    — sweep all configs (DSM_sd only)
#
#         Outputs:
#           output/intermediate/sm6/{site}/
#             dsm_sd_res_10_m.tif, chm_sd_res_10_m.tif
#           output/intermediate/sm6/
#             sm6_results_all_sites.rds, sm6_sweep_results.rds
#
# Run from project root (NC_Full/):
#   source("scripts/04_het_compute.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[04] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("scripts/steps/14_sm6_compute_heterogeneity.R")
step("scripts/steps/15_sm6_analysis.R")
step("scripts/steps/16_sm6_analysis_sweep.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("04_het_compute done — {elapsed_total} min")
