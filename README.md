# LAI consistency between Sentinel-2 and airborne LiDAR in temperate deciduous forests

Revision pipeline for the paper **"Assessing and Reducing Discrepancies Between Sentinel-2 and LiDAR-Derived LAI in Temperate Deciduous Forests"** (*Remote Sensing of Environment*, in revision).

This repository contains the refactored analysis code supporting the major revision of manuscript RSE-D-25-04417. It re-implements the submitted analysis with a modular structure, adds the sensitivity analyses requested by reviewers, and produces machine-readable CSV outputs.

## Study sites

Three French temperate deciduous forests:

- **Aigoual** (Gard, Cévennes uplands — beech-dominated)
- **Blois** (Loir-et-Cher — oak-dominated)
- **Mormal** (Nord — mixed broadleaved)

## Repository structure

- `R/` — analysis modules (functions, with roxygen2 documentation)
- `scripts/` — orchestration scripts (one per analysis step)
- `.gitignore` — excludes data, intermediate outputs, and binary files

Source data (LiDAR point clouds, Sentinel-2 tiles, pre-computed rasters) are **not** distributed in this repository due to their size and licensing. Paths and expected inputs are documented in each script header.

## Pipeline

Scripts are numbered by execution order. The main pipeline runs sequentially; sensitivity analyses (10–12) are independent and can run in parallel.

| Step | Script | Outputs |
|------|--------|---------|
| 01   | `01a_download_s2_from_safe.R` or `01b_download_s2_from_cdse.R` | Sentinel-2 reflectance |
| 02   | `02_train_prosail_lut.R` | PROSAIL LUT + SVR ensemble per configuration |
| 03a  | `03a_sample_s2_pixels.R` | Stratified S2 pixel samples (AOI) |
| 03b  | `03b_extract_lidar_at_samples.R` | LiDAR LAI/PAD at sampled locations |
| 04   | `04_apply_prosail_inversion.R` | S2-derived LAI estimates |
| 05   | `05_sm5_compute_metrics.R` | Per-configuration comparison metrics (SM5) |
| 06   | `06_sm5_select_dopt.R` | Optimal integration depth d_opt (Pearson, RMSE, Bias, Pareto) |
| 07   | `07_sm6_compute_heterogeneity.R` | Canopy heterogeneity rasters (DSM, CHM; block + focal; multiple scales) |
| 08   | `08_sm6_analysis.R` | SM6b heterogeneity-stratified metrics (10 m block) |
| 09   | `09_sm6_analysis_sweep.R` | SM6b sweep across spatial windows |
| 10   | `10_k_sensitivity.R` | LAI_ALS sensitivity to extinction coefficient k ∈ {0.4, 0.5, 0.6} |
| 11   | `11_fcover_sensitivity.R` | LAI comparison sensitivity to fCover threshold ∈ {0.80, 0.90, 0.95} |
| 12   | `12_h_min_sensitivity.R` | LAI_ALS sensitivity to vegetation height threshold h_min ∈ {2, 3, 4, 5} m |

Each script starts with a pre-flight check that verifies all required input files exist and lists missing files in a single error message.

## Reviewer comments addressed by analysis

| Comment ID | Analysis / Script |
|------------|--------------------|
| R2.8 (CHM_std vs DSM_std) | `07` + `09` |
| R3.MAJOR (k sensitivity) | `10` |
| R3.minor.4 (fCover + h_min sensitivity) | `11` + `12` |
| R3.minor.7 (R² alongside r) | All metric scripts |
| R3.minor.12 (effective spatial window) | `07` + `09` |
| R3.minor.13, R4.spec.2 (multi-criteria d_opt) | `06` |

Other reviewer comments are addressed in the response letter through methodological discussion or editorial revisions.

## Conventions

- **Regression convention**: `lm(s2 ~ lidar)` throughout (Sentinel-2 as response, LiDAR as predictor). Slope = 1 tests the 1:1 line.
- **Bias**: computed as `mean(LAI_S2 - LAI_ALS)`. Positive Bias indicates S2 overestimates relative to ALS.
- **R²**: computed as Pearson r² (not `summary(lm())$r.squared`), consistent across all modules.
- **Heterogeneity class thresholds**: dynamically calibrated per variant to maintain ≥ 5 000 pixels per (site × class) pool.
- **Extinction coefficient k**: fixed at 0.5 (Bouvier et al. 2015 convention for temperate broadleaved forests).
- **Vegetation height threshold h_min**: 2 m by default (manuscript baseline).

## Dependencies

R ≥ 4.3 with: `terra`, `data.table`, `here`, `cli`, `lidR`, `fda`, `stats`. PROSAIL training additionally requires `prosail`, `caret`, `kernlab`. Sentinel-2 download requires credentials for Copernicus Data Space Ecosystem (CDSE).

## License

This work is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) (GPL-3).

## Contact

Please don't hesitate to initiate contact with me or one of my supervisors for any questions, remarks, or advice about this work.

- **nathan.corroyez@inrae.fr** (final-year PhD student)
- **nathan.corroyez14@gmail.com** (non-academic email address)
- **sylvie.durrieu@inrae.fr** (director of research)
- **jean-baptiste.feret@inrae.fr** (co-director of research)
- **jerome.ogee@inrae.fr** (co-director of research)

## Acknowledgments

We would like to thank the CNES agency for the TOSCA Grant N°00007689 along with INRAE's MathNum Department and Agence Nationale de la Recherche via the MaCCMic ANR Project (ANR-21-CE32-0012) for the cofunding of the PhD thesis of Nathan Corroyez, and the IMPRINT ANR Project (ANR-19-CE32-0005) that funded the LiDAR acquisitions of the Mormal, Blois, and Aigoual sites.

## References

- Bouvier, M., Durrieu, S., Fournier, R. A., & Renaud, J. P. (2015). Generalizing relationships between canopy structure and light transmission: the case of European temperate forests for a radiative transfer model parametrization. *Remote Sensing of Environment*, 171, 158–171.
- Féret, J. B., & de Boissieu, F. (2024). *prosail: PROSAIL leaf and canopy radiative transfer model and inversion routines*. R package version 2.4.1. https://gitlab.com/jbferet/prosail
- Roussel, J. R., Auty, D., Coops, N. C., Tompalski, P., Goodbody, T. R. H., Sánchez Meador, A., Bourdon, J. F., De Boissieu, F., & Achim, A. (2021). lidR: An R package for analysis of Airborne Laser Scanning (ALS) data. *Remote Sensing of Environment*, 251, 112061. doi:10.1016/j.rse.2020.112061

## Citation

Citation and Zenodo DOI will be added upon acceptance of the manuscript.
