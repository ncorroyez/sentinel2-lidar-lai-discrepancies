# ---
# title:  02_dopt_figures.R
# desc:   Phase 2 — d_opt figures: Pearson r vs depth curves (SM5 Fig.),
#         LAI_ALS vs LAI_S2_ATBD scatter and histogram comparison.
#
#         Steps (in order):
#           05b  plot_dopt_metrics — r vs depth curves, d_opt arrows, SM5 figs
#           08   scatter_lai_atbd  — 2×3 scatter + histogram comparison figure
#
#         Outputs (revision/output/figures/sm5/):
#           fig_sm5_r_vs_depth_*.pdf, fig_sm5_scatter_*.pdf, etc.
#
# Run from project root (NC_Full/):
#   source("revision/scripts/02_dopt_figures.R")
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

step("revision/scripts/steps/05b_sm5_plot_dopt_metrics.R")
step("revision/scripts/steps/08_sm5_scatter_lai_atbd.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("02_dopt_figures done — {elapsed_total} min")
