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
#           - 01a / 01b  : S2 downloads (revision/scripts/archive/)
#           - 18_scan_angle_* : requires SMB-mounted LAS files (archive/)
#
#         Estimated runtime: ~3–6 h (SVR training dominates)
#
# Run from project root (NC_Full/):
#   source("revision/scripts/00_run_all.R")
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
phase("revision/scripts/01_prepare.R")

# ── Phase 2 : d_opt selection + LAI_ALS_dopt rasters ─────────────────────────
phase("revision/scripts/02_dopt_compute.R")
phase("revision/scripts/02_dopt_figures.R")

# ── Phase 3 : PROSAIL full inversion + LAI_S2_opt rasters ────────────────────
phase("revision/scripts/03_prosail_compute.R")
phase("revision/scripts/03_prosail_figures.R")

# ── Phase 4 : Heterogeneity analysis (SM6) ───────────────────────────────────
phase("revision/scripts/04_het_compute.R")
phase("revision/scripts/04_het_figures.R")

# ── Phase 5 : Sensitivity analyses ───────────────────────────────────────────
phase("revision/scripts/05_sensitivity.R")
phase("revision/scripts/05_sensitivity_figures.R")

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
