# ---
# title:  03_prosail_figures.R
# desc:   Phase 3 — PROSAIL figures: Pareto selection plot, parameter
#         distribution histograms, and ATBD_T/ATBD_F validation scatter.
#
#         Steps (in order):
#           12b  plot_prosail_pareto         — Pareto front visualisation
#           22   plot_prosail_param_distrib  — LUT parameter distributions
#           24   plot_s2_atbd_validation     — ATBD_T vs ATBD_F vs LAI_ALS
#
#         Outputs (revision/output/figures/sm5/):
#           fig_sm5_pareto_*.pdf, fig_sm5_param_distrib_*.pdf,
#           fig_sm5_atbd_validation_*.pdf
#
# Run from project root (NC_Full/):
#   source("revision/scripts/03_prosail_figures.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[03f] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("revision/scripts/steps/12b_sm5_plot_prosail_pareto.R")
step("revision/scripts/steps/22_sm5_plot_prosail_param_distributions.R")
step("revision/scripts/steps/24_plot_s2_atbd_validation.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("03_prosail_figures done — {elapsed_total} min")
