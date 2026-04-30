# ---
# title:  05_sensitivity_figures.R
# desc:   Phase 5 — Sensitivity figures: visualises the k/z_min, fCover,
#         h_min, and joint fCover×h_min sensitivity analyses from 05_sensitivity.R.
#
#         Steps (in order):
#           05d  plot_k_zmin_sensitivity — d_opt vs k/z_min facet grid
#           17b  k_sensitivity_plot      — LAI vs k comparison figure
#           18c  fcover_sensitivity_plot — LAI metrics vs fCover threshold
#           20b  fcover_hmin_combined    — 1D fCover + h_min combined panel
#           20c  fcover_hmin_joint_plot  — 2D joint sensitivity heatmap
#
#         Outputs (revision/output/figures/reviewers/):
#           fig_dopt_k_sensitivity_*.pdf, fig_lai_k_sensitivity_*.pdf,
#           fig_fcover_sensitivity_*.pdf, fig_fcover_hmin_combined_*.pdf,
#           fig_fcover_hmin_joint_*.pdf
#
# Run from project root (NC_Full/):
#   source("revision/scripts/05_sensitivity_figures.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[05f] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("revision/scripts/steps/05d_sm5_plot_k_zmin_sensitivity_dopt.R")
step("revision/scripts/steps/17b_k_sensitivity_plot.R")
step("revision/scripts/steps/18c_fcover_sensitivity_plot.R")
step("revision/scripts/steps/20b_fcover_hmin_combined_plot.R")
step("revision/scripts/steps/20c_fcover_hmin_joint_plot.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("05_sensitivity_figures done — {elapsed_total} min")
