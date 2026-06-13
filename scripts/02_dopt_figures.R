# ---
# title:  02_dopt_figures.R
# desc:   Phase 2 — d_opt figures: Pearson r vs depth curves (SM5 Fig.),
#         LAI_ALS vs LAI_S2_ATBD scatter and histogram comparison.
#
#         Steps (in order):
#           05b  plot_dopt_metrics — r vs depth curves, d_opt arrows, SM5 figs
#         NOTE: 08 (scatter_lai_atbd) moved to 03_prosail_figures.R because
#         it depends on s2lai_summer_atbd_T_*.tif produced by 25 (phase 3).
#         Keeping 08 in phase 2 broke full-rerun pipelines.
#
#         Outputs (output/figures/):
#           fig_sm5_r_vs_depth_*.pdf, fig_sm5_scatter_*.pdf, etc.
#
# Run from project root (NC_Full/):
#   source("scripts/02_dopt_figures.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[02f] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("scripts/steps/05b_sm5_plot_dopt_metrics.R")
# 08_sm5_scatter_lai_atbd.R is now in 03_prosail_figures.R (it consumes the
# ATBD_T raster produced by 25 in phase 3).

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("02_dopt_figures done — {elapsed_total} min")
