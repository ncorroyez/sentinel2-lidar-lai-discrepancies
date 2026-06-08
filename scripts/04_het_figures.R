# ---
# title:  04_het_figures.R
# desc:   Phase 4 — Heterogeneity figures: SM6 heterogeneity maps and
#         regression panels, density scatter comparison, S2 reflectance
#         histograms, ALS metrics scatter, max height density, and PAD profiles.
#
#         Steps (in order):
#           20_sm6  plot_heterogeneity     — SM6b main figure (het grid)
#           21      sm6_density_scatter    — density scatter LAI_ALS vs LAI_S2
#           23      plot_s2_reflectances   — S2 reflectance histogram grid
#           26      plot_als_metrics       — ALS metrics scatter (LAI, fCover…)
#           27      plot_max_height_density — tree height density by site
#           28      plot_pad_profiles      — PAD profiles per site
#
#         Outputs (revision/output/figures/):
#           sm6/fig_sm6_het_*.pdf, fig_sm6_density_scatter_*.pdf,
#           fig_s2_reflectances_*.pdf, fig_als_metrics_*.pdf,
#           fig_max_height_*.pdf, fig_pad_profiles_*.pdf
#
# Run from project root (NC_Full/):
#   source("revision/scripts/04_het_figures.R")
# ---

library(here)
library(cli)

step <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- basename(path)
  cli::cli_h2("[04f] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[done] {label}  ({elapsed} min)")
  invisible(elapsed)
}

t_start <- proc.time()

step("revision/scripts/steps/20_sm6_plot_heterogeneity.R")
step("revision/scripts/steps/21_sm6_density_scatter.R")
step("revision/scripts/steps/23_plot_s2_reflectance_histograms.R")
step("revision/scripts/steps/26_plot_als_metrics_scatter.R")
step("revision/scripts/steps/27_plot_max_height_density.R")
step("revision/scripts/steps/28_plot_pad_profiles.R")

elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_h1("04_het_figures done — {elapsed_total} min")
