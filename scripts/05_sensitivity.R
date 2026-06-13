# ---
# title:  05_sensitivity.R
# desc:   Phase 5 — Sensitivity analyses (compute): tests sensitivity of
#         d_opt and LAI metrics to k (extinction coefficient), z_min (minimum
#         height), fCover threshold, h_min (height cutoff), and their joint
#         variation (fCover × h_min).
#
#         Steps (in order):
#           05c  k_zmin_sensitivity_dopt    — d_opt vs k ∈ {0.4,0.5,0.6}, z_min
#           17   k_sensitivity             — LAI_ALS vs LAI_S2 vs k
#           18   fcover_sensitivity        — LAI comparison vs fCover ∈ {80,90,95%}
#           19   h_min_sensitivity         — LAI comparison vs h_min ∈ {2,3,4,5} m
#           20j  fcover_hmin_joint         — joint fCover × h_min sensitivity
#
#         Outputs (output/intermediate/reviewers/):
#           kt_sensitivity_atbd_*.csv
#           k_sensitivity_results.csv
#           fcover_sensitivity_atbd.csv
#           h_min_sensitivity_atbd.csv
#           joint_fcover_hmin_sensitivity_atbd.csv
#
# Run from project root (NC_Full/):
#   source("scripts/05_sensitivity.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[05] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("scripts/steps/05c_sm5_k_zmin_sensitivity_dopt.R")
step("scripts/steps/17_k_sensitivity.R")
step("scripts/steps/18_fcover_sensitivity.R")
step("scripts/steps/19_h_min_sensitivity.R")
step("scripts/steps/20_fcover_hmin_joint_sensitivity.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("05_sensitivity done — {elapsed_total} min")
