# ---
# title:  00_run_all.R
# desc:   Master pipeline — runs all 5 phases in order.
#         Generates all figures, CSVs, and intermediate rasters for the
#         RSE revision (Corroyez et al.).
#
#         Prerequisites:
#           - config.yml profile pointing to valid data_root and output_root
#             (set profile: local or profile: smb)
#           - R packages: here, terra, prosail, liquidSVM, data.table,
#             cli, ggplot2, patchwork, readr, scales, yaml, sgsR, sf, dplyr
#
#         Exclusions (run manually if needed):
#           - 01a / 01b  : S2 downloads (scripts/archive/)
#           - 18_scan_angle_* : requires SMB-mounted LAS files (archive/)
#
#         Estimated runtime: ~3–6 h (SVR training dominates)
#
# Run from project root:
#   source("scripts/00_run_all.R")
#
# ── Key parameters — where to tune them ──────────────────────────────────────
#
# The pipeline runs with the defaults used in the manuscript. To change a key
# parameter, edit the header of the corresponding step file:
#
#   nbSamples = 5000              scripts/steps/03a_sample_s2_pixels.R
#   h_min sweep = c(10, 15, 20)   scripts/steps/03a_sample_s2_pixels.R
#                                 (dthr in the manuscript)
#   sampling method = stratified_uniform
#                                 scripts/steps/03a_sample_s2_pixels.R
#   fCover threshold = 0.90       R/sentinel2_sampling.R (stratified_sampling_uniform)
#
#   k_ref = 0.5 (PADs precomputed at this k)
#   k_select = 0.65 (k used in this study) — set consistently in:
#                                 scripts/steps/06_sm5_select_dopt.R
#                                 scripts/steps/07_compute_lai_als_dopt.R
#                                 scripts/steps/18_fcover_sensitivity.R
#                                 scripts/steps/19_h_min_sensitivity.R
#                                 scripts/steps/20_fcover_hmin_joint_sensitivity.R
#                                 (rescaled via rescale_lai_for_k() — R/k_sensitivity.R)
#   theta_select = 0              scripts/steps/06_sm5_select_dopt.R
#   d_opt selection method = pareto
#                                 R/sm5_dopt.R (select_dopt) and step 06
#   combined_mode = average       scripts/steps/06_sm5_select_dopt.R
#
#   PROSAIL LAI scenarios         scripts/steps/09_train_prosail_full.R
#   LUT base config = ATBD v2     R/prosail_lut.R (get_atbd_lut_input)
#   soil spectrum = atbd_v2       R/prosail_lut.R (spec_soil_atbd_v2)
#
#   k sweep = {0.4, 0.5, 0.6}     scripts/steps/17_k_sensitivity.R
#   fCover sweep = {0.80…0.95}    scripts/steps/18_fcover_sensitivity.R
#   h_min sweep                   scripts/steps/19_h_min_sensitivity.R
# ---

library(here)
library(cli)

phase <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- sub("\\.R$", "", basename(path))
  cli::cli_rule(left = "[PHASE] {label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("[PHASE done] {label}  ({elapsed} min)")
  cat("\n")
  invisible(elapsed)
}

t_start <- proc.time()

# ── Phase 1 : Preparation (SVR ATBD + S2 sampling + PAD extraction) ──────────
phase("scripts/01_prepare.R")

# ── Phase 2 : d_opt selection + LAI_ALS_dopt rasters ─────────────────────────
phase("scripts/02_dopt_compute.R")
phase("scripts/02_dopt_figures.R")

# ── Phase 3 : PROSAIL full inversion + LAI_S2_opt rasters ────────────────────
phase("scripts/03_prosail_compute.R")
phase("scripts/03_prosail_figures.R")

# ── Phase 4 : Heterogeneity analysis (SM6) ───────────────────────────────────
phase("scripts/04_het_compute.R")
phase("scripts/04_het_figures.R")

# ── Phase 5 : Sensitivity analyses ───────────────────────────────────────────
phase("scripts/05_sensitivity.R")
phase("scripts/05_sensitivity_figures.R")

# ── Summary ───────────────────────────────────────────────────────────────────
elapsed_total <- round((proc.time() - t_start)[["elapsed"]] / 60, 1)
cli::cli_rule()
cli::cli_h1("Pipeline complet — {elapsed_total} min")
cli::cli_bullets(c(
  "v" = "Phase 1 : preparation (ATBD SVR, S2 sampling, PAD extraction)",
  "v" = "Phase 2 : d_opt selection + LAI_ALS_dopt rasters + figures",
  "v" = "Phase 3 : PROSAIL inversion + LAI_S2_opt rasters + figures",
  "v" = "Phase 4 : heterogeneity computation + figures",
  "v" = "Phase 5 : sensitivity analyses + figures",
  ">" = "Exclus  : S2 downloads (01a/01b), scan angle correction (archive/)"
))
