# Pipeline steps — inputs, outputs, and produced figures

This document maps each step script under `scripts/ch2/steps/` to its inputs,
outputs, and (if applicable) the manuscript figure(s) it produces. Steps
are listed in the order they run inside the phase orchestrators
(`scripts/0X_*.R`). All paths are relative to the project root.

A step is normally invoked by its phase orchestrator. Each step also runs
standalone via `source("scripts/ch2/steps/<step>.R")` once its upstream
outputs exist.

Legend:
- **Inputs / Outputs**: paths under `data/` (read) or `output/` (write).
- **Figures**: written to `output/figures/<name>.{pdf,png}`.
- *(none)*: a script that only writes intermediate CSV/raster files.

---

## Phase 1 — Preparation (`scripts/ch2/01_prepare.R`)

### `02a_train_prosail_atbd.R`
PROSAIL ATBD-v3 SVR training (per site, codist = TRUE / FALSE).
- **Inputs**: S2 angle TIFFs in `data/results/{site}/PROSAIL_Optimization/`.
- **Outputs**: `output/intermediate/PROSAIL_Models/{site}/atbd_{T,F}/`.
- **Figures**: (none.)

### `03a_sample_s2_pixels.R`
Stratified-uniform pixel sampling (5000 pixels per site, three h_min).
- **Inputs**: S2 reflectance, ATBD LAI raster, fCover mask, AOI shapefiles.
- **Outputs**: `output/intermediate/sampling/{site}/Sampling_*.GPKG`,
  `S2_reflectance_*.csv`.
- **Figures**: (none.)

### `03b_extract_lidar_at_samples.R`
PAD profile extraction at sample positions for 4 normalisation variants.
- **Inputs**: per-depth PAD rasters in `data/prosail_lidar/{site}/LiDAR/testPADs/`.
- **Outputs**: `output/intermediate/sampling/{site}/PAD_*.csv`.
- **Figures**: (none.)

### `04a_apply_prosail_atbd.R`
Apply trained ATBD SVR to the sampled pixels.
- **Inputs**: ATBD SVR models (02a), sampled S2 reflectances (03a).
- **Outputs**: `output/intermediate/sm5/{site}/LAI_estimated_atbd_*.csv`.
- **Figures**: (none.)

---

## Phase 2 — d<sub>opt</sub> selection (`scripts/ch2/02_dopt_compute.R` / `02_dopt_figures.R`)

### `05a_sm5_compute_metrics_atbd.R`
R, R², RMSE, Bias, Slope between LAI<sub>S2-ATBD</sub> and LAI<sub>ALS</sub>
at depths 1–38 m, for all 4 normalisations × 3 sites × 3 h_min.
- **Inputs**: ATBD SVR predictions (04a), PAD CSVs (03b).
- **Outputs**: `output/intermediate/sm5/all_results_atbd_LIDFa_lai_LMA_BROWN.csv`.
- **Figures**: (none.)

### `05c_sm5_k_zmin_sensitivity_dopt.R`
Sensitivity sweep over (k, theta) for d<sub>opt</sub> robustness.
- **Inputs**: same as 05a; uses analytical rescaling `lai × k_ref / k_new`.
- **Outputs**: `output/intermediate/sm5/kt_sensitivity_atbd_*.csv`,
  `z_min_sensitivity_atbd_*.csv`.
- **Figures**: (none — plotted by 05d.)

### `06_sm5_select_dopt.R`
Pareto multi-criteria d<sub>opt</sub> selection at k=0.6, θ=0
(reads 05c if k≠0.5).
- **Inputs**: 05a or 05c CSV depending on (k_select, theta_select).
- **Outputs**: `output/intermediate/sm5/dopt_reference.csv`,
  `dopt_all.csv`.
- **Figures**: (none — plotted by 05b.)

### `07_compute_lai_als_dopt.R`
LAI<sub>ALS,dopt</sub> rasters per site × scenario (rescaled to k_select).
- **Inputs**: 06 dopt_reference, per-depth PADs (`testPADs/`).
- **Outputs**: `output/intermediate/lai_als_dopt/{site}/LAI_ALS_dopt_*.tif`.
- **Figures**: (none.)

### `05b_sm5_plot_dopt_metrics.R`
Per-metric (R, RMSE, |Bias|, |Slope-1|) and joint-Pareto figures for d<sub>opt</sub>.
- **Inputs**: 05a CSV (or 05c if k≠0.5).
- **Outputs**: (none beyond figures.)
- **Figures**: `dopt_{R,RMSE,absBias,absSlope1,Pareto}_ATBD_hmin10.{pdf,png}`,
  `sm_dopt_4metrics_hmin10{,_combined}.{pdf,png}`,
  `sm_dopt_pareto_all_hmin.{pdf,png}`.

### `05d_sm5_plot_k_zmin_sensitivity_dopt.R`
d<sub>opt</sub> sensitivity to k and z_min, derived from 05c.
- **Figures**: `dopt_k_sensitivity.{pdf,png}`,
  `dopt_zmin_sensitivity.{pdf,png}`.

### `08_sm5_scatter_lai_atbd.R`
Scatter plots LAI<sub>ALS,dopt</sub> vs LAI<sub>S2-ATBD</sub>.
- **Inputs**: 04a SVR predictions, 07 LAI rasters, 06 d<sub>opt</sub>.
- **Figures**: `scatter_lai_atbd_*.{pdf,png}`.

---

## Phase 3 — PROSAIL full inversion (`scripts/ch2/03_prosail_*.R`)

### `09_train_prosail_full.R`
270 PROSAIL configurations × LAI scenarios, SVR ensemble per (site, h_min).
- **Inputs**: 03a samples, 03b PAD CSVs, 07 LAI<sub>ALS,dopt</sub> rasters.
- **Outputs**: `output/intermediate/PROSAIL_Models/{site}/{scenario}/...`.
- **Figures**: (none.)

### `10_apply_prosail_full.R`
Apply the 270 trained SVRs to the sampled pixels.
- **Outputs**: `output/intermediate/sm5/{site}/LAI_estimated_*.csv`.

### `11_sm5_compute_metrics_full.R`
Metrics (R², RMSE, Bias, Slope) per PROSAIL config × site × h_min × scenario.
- **Outputs**: `output/intermediate/sm5/all_results_combined_*.csv`
  (~500 k rows).

### `12_sm5_select_prosail_opt.R`
Pareto selection of the best PROSAIL config per site (per_site, all_sites, fixed).
- **Inputs**: 11 combined CSV, 06 dopt_reference.
- **Outputs**: `output/intermediate/sm5/prosail_opt.csv`.

### `12b_sm5_plot_prosail_pareto.R`
Pareto-front PROSAIL configuration plot + parameter frequency.
- **Figures**: `prosail_pareto_DSM.{pdf,png}`,
  `prosail_front1_metrics_DSM.{pdf,png}`,
  `prosail_front1_params_DSM.{pdf,png}`.

### `13_sm5_predict_lai_raster.R`
Full-tile LAI<sub>S2-opt</sub> raster prediction with the selected config.
- **Outputs**: `output/intermediate/sm6/{site}/s2lai_summer_opt_*.tif`.

### `25_train_apply_prosail_atbd_rasters.R`
Full-tile ATBD<sub>T</sub> / ATBD<sub>F</sub> rasters (for SM6 + validation).
- **Outputs**: `output/intermediate/sm6/{site}/s2lai_summer_atbd_{T,F}_res_10_m.tif`.

### `22_sm5_plot_prosail_param_distributions.R`
PROSAIL input parameter distributions (ATBD vs OPT configs).
- **Figures**: `prosail_param_distributions.{pdf,png}`.

### `24_plot_s2_atbd_validation.R`
ATBD<sub>T</sub> / ATBD<sub>F</sub> validation against SNAP biophysical.
- **Inputs**: 25 ATBD rasters, SNAP `lai.img`, masks.
- **Figures**: `atbd_validation_*.{pdf,png}`.

---

## Phase 4 — Heterogeneity (`scripts/ch2/04_het_*.R`)

### `14_sm6_compute_heterogeneity.R`
DSM-std and CHM-std heterogeneity metrics at 10 m.
- **Inputs**: 1 m DSM, 1 m CHM, fCover mask, ATBD rasters.
- **Outputs**: `output/intermediate/sm6/{site}/heterogeneity_*.tif`.

### `15_sm6_analysis.R`
Per-(site × heterogeneity class) goodness-of-fit metrics for three LAI
combinations: ATBD vs ALS, ATBD vs ALS<sub>dopt</sub>, opt vs ALS<sub>dopt</sub>.
Workflow:
1. Site-level uniform-LAI sampling on the whole forest: 1-unit bins from
   2 to floor(p98<sub>site</sub>), n_target = 50,000 per site
   (`sample_uniform_lai_dt`).
2. DSM_sd heterogeneity thresholds calibrated on the LAI-uniform sample
   (balance score; ≥ 5000 per site × class).
3. Pixels classified Low / Medium / High; metrics computed per class for
   the three combinations.
- **Outputs**: `output/intermediate/sm6/heterogeneity_analysis.csv`.

### `16_sm6_analysis_sweep.R`
Sensitivity sweep across heterogeneity thresholds and configs.
- **Outputs**: `output/intermediate/sm6/DSM_sd_thresholds.csv`.

### `20_sm6_plot_heterogeneity.R`
Final SM6 heterogeneity figure (3 sites × metrics).
- **Figures**: `sm6_heterogeneity.{pdf,png}`.

### `20d_sm6_heter_continuous.R`
Continuous-view alternative to `20_sm6_plot_heterogeneity.R`: hex-bin
LAI residuals vs DSM_sd with per-panel statistics annotations, +
DSM_sd histogram per site, + window-size sweep SI figure.
- **Inputs**: 14 heterogeneity rasters, 13 LAI<sub>S2-opt</sub>,
  07 LAI<sub>ALS,dopt</sub>, 25 ATBD rasters.
- **Figures**: `sm6_heter_continuous_hex.{pdf,png}`,
  `sm6_heter_dsm_sd_hist.{pdf,png}`,
  `sm6_heter_window_sweep.{pdf,png}`.

### `21_sm6_density_scatter.R`
Density scatter LAI<sub>ALS,dopt</sub> vs LAI<sub>S2,opt</sub>.
- **Figures**: `sm6_density_scatter.{pdf,png}`.

### `23_plot_s2_reflectance_histograms.R`
S2 reflectance histograms (sampled vs full tile, per band).
- **Figures**: `s2_reflectance_histograms.{pdf,png}`.

### `26_plot_als_metrics_scatter.R`
ALS-vs-ALS metric cross-comparison (LAI<sub>full</sub> vs LAI<sub>dopt</sub>).
- **Figures**: `als_metrics_scatter.{pdf,png}`.

### `27_plot_max_height_density.R`
Density of max canopy height (`max_res_10_m.tif`) per site.
- **Figures**: `max_height_density.{pdf,png}`.

### `28_plot_pad_profiles.R`
Vertical PAD profile figures (per site, per normalisation).
- **Figures**: `pad_profiles_*.{pdf,png}`.

---

## Phase 5 — Sensitivity analyses (`scripts/ch2/05_sensitivity*.R`)

Reviewer-driven additional analyses. All operate post-hoc on precomputed
PAD stacks and Raw ladstacks; no need to rerun the LiDAR pipeline.

### `17_k_sensitivity.R`
Sweep k ∈ {0.3, 0.4, 0.5, 0.6, 0.7, 0.8} (analytical rescaling).
- **Inputs**: `data/results/{site}/Metrics/Deciduous_Only/ladstack_classic.tif`.
- **Outputs**: `output/intermediate/reviewers/k_sensitivity_metrics.csv`.

### `17b_k_sensitivity_plot.R`
- **Figures**: `k_sensitivity_metrics.{pdf,png}`.

### `18_fcover_sensitivity.R`
Sweep fCover threshold ∈ {0.80, 0.85, 0.90, 0.95} at k = 0.6.
- **Inputs**: Raw ladstack, Raw ATBD raster, mask components.
- **Outputs**: `output/intermediate/reviewers/fcover_sensitivity_atbd.csv`.

### `18c_fcover_sensitivity_plot.R`
- **Figures**: `fcover_sensitivity_atbd.{pdf,png}`.

### `19_h_min_sensitivity.R`
Sweep h<sub>min</sub> ∈ {2, 3, 4, 5} m at k = 0.6.
- **Outputs**: `output/intermediate/reviewers/h_min_sensitivity.csv`.

### `20_fcover_hmin_joint_sensitivity.R`
Joint sweep h<sub>min</sub> × fCover at k = 0.6.
- **Outputs**: `output/intermediate/reviewers/joint_sensitivity_atbd.csv`.

### `20b_fcover_hmin_combined_plot.R`
Combined view: h<sub>min</sub> rows × fCover columns.
- **Figures**: `fcover_hmin_combined.{pdf,png}`.

### `20c_fcover_hmin_joint_plot.R`
Heatmap-style joint figure.
- **Figures**: `fcover_hmin_joint.{pdf,png}`.

---

## Optional / out-of-pipeline steps

Not called by `00_run_all.R`. Source them manually after their upstream
outputs are available.

### `08b_sm5_scatter_lai_atbd_lowcanopy.R`
Companion to `08_sm5_scatter_lai_atbd.R` restricted to low-canopy
pixels (max canopy height ≤ 10 m).
- **Inputs**: `data/results/{site}/Metrics/Deciduous_Only/lidarlai_res_10_m.tif`,
  `max_res_10_m.tif`, 25 ATBD raster.
- **Figures**: `scatter_LAI_ALS_vs_S2_ATBD_h_le_10.{pdf,png}`.

### `20e_sm5_heter_classes.R`
Apples-to-apples discrete heterogeneity companion to `20d` (re-classifies
SM5 samples into Low / Medium / High DSM_sd bins).
- **Outputs**: `output/intermediate/sm5/sm5_heter_classes.csv`.

### `20e_plot_sm5_heter_classes.R`
Plots the discrete heterogeneity analysis produced by `20e_sm5_heter_classes`.
- **Figures**: `sm5_heter_classes.{pdf,png}`.

### `20f_sm6a_heter_continuous.R`
Continuous heterogeneity companion to `20d`, restricted to the SM6a
sample pool (versus the full deciduous mask in `20d`).
- **Figures**: `sm6a_heter_continuous_hex.{pdf,png}`.

### `29_prosail_loso_transfer.R`
Leave-one-site-out transfer of the PROSAIL LUT configuration, and the
partial attribution of the improvement to the LiDAR-derived LAI prior.
Re-scores the existing 225-configuration ensemble; no PROSAIL training
and no inversion are re-run. For each site in turn the configuration is
selected on the two other sites among the 135 configurations whose LAI
prior is not derived from LiDAR, then applied unchanged at the held-out
site. Two depth-transfer rules are reported: the mean of the two training
sites' *d*_opt, and a joint re-run of the depth selection on the two
training sites.
- **Inputs**: `output/intermediate/sm5/all_results_combined_LIDFa_lai_LMA_BROWN.csv`,
  `all_results_atbd_LIDFa_lai_LMA_BROWN.csv`, `dopt_reference.csv`,
  `prosail_opt.csv`.
- **Outputs**: `output/tables/Table_A16_partial_attribution_check.csv`,
  `Table_A17_loso_transfer.csv`, `Table_A17b_loso_summary.csv`,
  `Table_A18_rmse_dispersion.csv`.

### `29b_prosail_loso_figure.R`
Figure A.19: dispersion of the inter-sensor RMSE over the
225-configuration ensemble at each site's Pareto *d*_opt, with the ATBD
baseline, the site-specific Pareto optimum, the best LiDAR-free
configuration and the transferred configuration marked on it. Requires
`29_prosail_loso_transfer.R`.
- **Figures**: `Fig_A19_loso_rmse_spread.{pdf,png}`.
