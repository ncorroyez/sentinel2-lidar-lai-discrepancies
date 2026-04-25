#!/usr/bin/env bash
# Pipeline k=0.6, theta=30° — du training PROSAIL ATBD à l'analyse d'hétérogénéité.
# Run from the project root (NC_Full/):  bash run_pipeline.sh
set -euo pipefail
cd "$(dirname "$0")"

scripts=(
  # ── Step 1: ATBD baseline SVR (k=0.6 pool) ───────────────────────────────
  revision/scripts/02a_train_prosail_atbd.R
  # ── Step 2: Apply ATBD SVR → LAI_S2_atbd estimates ────────────────────────
  revision/scripts/04a_apply_prosail_atbd.R
  # ── Step 3: k×theta factorial (LiDAR) vs S2 ATBD → kt_sensitivity CSV ────
  revision/scripts/05c_sm5_k_zmin_sensitivity_dopt.R
  # ── Step 4: Select d_opt at k=0.6, theta=30° (reads 05c CSV) ─────────────
  revision/scripts/06_sm5_select_dopt.R
  # ── Step 5: Compute LAI_ALS_dopt rasters (rescaled to k=0.6, theta=30°) ──
  revision/scripts/07_compute_lai_als_dopt.R
  # ── Step 6: Train 270 SVR configs × 3 sites × 2 scenarios ─────────────────
  revision/scripts/09_train_prosail_full.R
  # ── Step 7: Apply SVRs → LAI_estimated CSVs ───────────────────────────────
  revision/scripts/10_apply_prosail_full.R
  # ── Step 8: Full metrics table (lidar_scale applied) ──────────────────────
  revision/scripts/11_sm5_compute_metrics_full.R
  # ── Step 9: Select optimal PROSAIL config ─────────────────────────────────
  revision/scripts/12_sm5_select_prosail_opt.R
  # ── Step 9b: Pareto figure (R² vs RMSE, all configs) ──────────────────────
  revision/scripts/12b_sm5_plot_prosail_pareto.R
  # ── Step 10: Predict LAI raster ───────────────────────────────────────────
  revision/scripts/13_sm5_predict_lai_raster.R
  # ── Step 11–13: Heterogeneity analysis ────────────────────────────────────
  revision/scripts/14_sm6_compute_heterogeneity.R
  revision/scripts/15_sm6_analysis.R
  revision/scripts/16_sm6_analysis_sweep.R
  # ── Step 14: Heterogeneity figure (SM6) ───────────────────────────────────
  revision/scripts/20_sm6_plot_heterogeneity.R
)

t_total=$SECONDS

for script in "${scripts[@]}"; do
  # Skip comment lines
  [[ "$script" == \#* ]] && continue

  echo "▶ $(date '+%H:%M:%S') — $script"
  t0=$SECONDS
  Rscript "$script"
  elapsed=$(( SECONDS - t0 ))
  echo "  ✓ done in ${elapsed}s"
  echo ""
done

echo "Pipeline complet — $(date '+%H:%M:%S') — total $(( SECONDS - t_total ))s"
