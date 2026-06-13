# ---
# title:  20e_sm5_heter_classes.R
# desc:   Apples-to-apples discrete companion to 20d_sm6_heter_continuous.R.
#         Uses the SAME SM5 stratified-uniform sample as 20d (5,000 points
#         per site, h_min=10, LAI-uniform, fCover > 0.9) but classifies
#         DSM_sd into Low / Medium / High (+ Total) instead of windowing.
#
#         Purpose: defend the continuous view in 20d by showing that the
#         qualitative story holds when the same sample is summarised by
#         discrete heterogeneity classes.
#
#         Thresholds: fixed at DSM_sd 2.5 m (low/medium) and 5.0 m (medium/
#         high), consistent with the SM6a discrete analysis in script 15.
#
#         Reads (via R/sm5_sample_loader.R):
#           output/intermediate/sampling/{site}/Sampling_*.GPKG
#           output/intermediate/sm6/{site}/dsm_sd_*.tif
#           output/intermediate/PROSAIL_Models/{site}/.../atbd/LAI_*.csv
#           output/intermediate/PROSAIL_Models/{site}/.../common/LAI_*.csv
#           output/intermediate/sampling/{site}/PAD_DSM_Depth_*_Samples_*.csv
#           output/intermediate/sm5/prosail_opt.csv   (d_opt + Column_opt)
#
#         Writes:
#           output/intermediate/sm6/heterogeneity_analysis_sm5_per_site.csv
#           output/intermediate/sm6/heterogeneity_analysis_sm5_all_sites.csv
#
#         Schema matches heterogeneity_analysis_*.csv produced by script 15
#         so that scripts/steps/20_sm6_plot_heterogeneity.R can plot it
#         (by pointing its loader at the *_sm5_* files).
#
# Run from project root (NC_Full/):
#   source("scripts/steps/20e_sm5_heter_classes.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_sample_loader.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites              <- c("Aigoual", "Blois", "Mormal")
norm_ref           <- "DSM_keepTrees"
metric_source_lbl  <- "DSM"
dopt_modes_to_run  <- c("per_site", "all_sites")
dsm_sd_filename    <- "dsm_sd_res_10_m.tif"

# Fixed thresholds on DSM_sd — same as script 15 fixed_thresholds list
low_thr            <- 2.5
high_thr           <- 5.0

# k rescaling (PADs precomputed at k_ref, study uses k_select)
k_ref              <- 0.5
k_select           <- 0.65

# Sampling file naming (same as 20d)
h_min              <- 10L
nb_samples         <- 5000L

# LAI min filter — matches compute_metrics_by_class in script 15
lai_min            <- 2.0

# Combinations — same names as script 15 / 20 so plot script is reusable
combinations <- list(
  list(name = "ATBD_vs_ALS",      lidar_col = "lai_als",      s2_col = "lai_s2_atbd"),
  list(name = "ATBD_vs_ALS_dopt", lidar_col = "lai_als_dopt", s2_col = "lai_s2_atbd"),
  list(name = "opt_vs_ALS_dopt",  lidar_col = "lai_als_dopt", s2_col = "lai_s2_opt")
)

# ── Helper: metrics for one (lidar, s2) pair (s2 ~ lidar convention) ──────────

compute_pair_metrics <- function(lidar_v, s2_v) {
  ok <- !is.na(lidar_v) & !is.na(s2_v)
  lidar_v <- lidar_v[ok]; s2_v <- s2_v[ok]
  n_obs <- length(lidar_v)
  if (n_obs < 10L) {
    return(list(n = n_obs, R = NA_real_, R2 = NA_real_, RMSE = NA_real_,
                Bias = NA_real_, Bias_pvalue = NA_real_,
                Slope = NA_real_, Slope_pvalue = NA_real_))
  }
  diff <- s2_v - lidar_v
  r_val    <- stats::cor(lidar_v, s2_v)
  rmse_val <- sqrt(mean(diff^2))
  bias_val <- mean(diff)
  bias_p   <- tryCatch(stats::t.test(diff)$p.value, error = function(e) NA_real_)
  fit      <- stats::lm(s2_v ~ lidar_v)
  cs       <- summary(fit)$coefficients
  slope_val <- cs[2L, 1L]
  slope_p   <- if (nrow(cs) >= 2L) cs[2L, 4L] else NA_real_
  list(n = n_obs, R = r_val, R2 = r_val^2, RMSE = rmse_val,
       Bias = bias_val, Bias_pvalue = bias_p,
       Slope = slope_val, Slope_pvalue = slope_p)
}

# ── Load PROSAIL opt table (d_opt + Column_opt per site/mode) ─────────────────

opt_csv <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv")
if (!file.exists(opt_csv))
  stop("Missing: ", opt_csv, "\nRun step 12 first.")
opt_dt <- data.table::fread(opt_csv)

# ── Loop over dopt mode (per_site, all_sites) ─────────────────────────────────

for (current_mode in dopt_modes_to_run) {

  cli::cli_h1("Mode: {current_mode}")
  cli::cli_alert_info(
    "DSM_sd thresholds: low = {low_thr} m, high = {high_thr} m (fixed)"
  )

  # ── d_opt + Column_opt per site ────────────────────────────────────────────
  opt_ref <- opt_dt[Norm == norm_ref & Site %in% sites &
                     d_opt_source == current_mode]
  if (nrow(opt_ref) != length(sites)) {
    cli::cli_warn("prosail_opt rows missing for mode '{current_mode}' — skipping.")
    next
  }

  # ── Load per-site SM5 sample (DSM_sd + 4 LAI columns) ──────────────────────
  per_site <- vector("list", length(sites))
  for (i in seq_along(sites)) {
    s   <- sites[i]
    row <- opt_ref[Site == s]
    cli::cli_alert_info(
      "{s}: d_opt = {row$d_opt}  |  Column_opt = {row$Column_opt}"
    )
    per_site[[i]] <- build_sm5_site_dt(
      site            = s,
      d_opt_per_site  = row$d_opt,
      column_opt      = row$Column_opt,
      k_select        = k_select,
      k_ref           = k_ref,
      h_min           = h_min,
      nb_samples      = nb_samples,
      paths_obj       = paths,
      ext_results     = paths$ext_results,
      dsm_sd_filename = dsm_sd_filename
    )
  }
  dt_full <- data.table::rbindlist(per_site)

  # ── LAI ≥ lai_min filter (same as compute_metrics_by_class) ────────────────
  dt_full <- dt_full[lai_als >= lai_min & !is.na(dsm_sd)]
  cli::cli_alert_info(
    "n per site (after LAI ≥ {lai_min}): ",
    paste(table(dt_full$site), collapse = " / ")
  )

  # ── Classify DSM_sd into Low / Medium / High ───────────────────────────────
  dt_full[, het_class := data.table::fcase(
    dsm_sd <  low_thr,                     "Low",
    dsm_sd >= low_thr & dsm_sd < high_thr, "Medium",
    dsm_sd >= high_thr,                    "High"
  )]
  dt_full[, het_class := factor(het_class, levels = c("Low", "Medium", "High"))]

  cli::cli_text("Class counts (site x class):")
  print(table(dt_full$site, dt_full$het_class))

  # ── Compute metrics per (site x combination x class), + Total per site ────
  rows_out <- list()
  for (s in sites) {
    dt_s <- dt_full[site == s]
    for (combo in combinations) {

      for (cls in c("Low", "Medium", "High")) {
        sub <- dt_s[het_class == cls]
        m   <- compute_pair_metrics(sub[[combo$lidar_col]], sub[[combo$s2_col]])
        rows_out[[length(rows_out) + 1L]] <- data.table::data.table(
          site          = s,
          metric_source = metric_source_lbl,
          combination   = combo$name,
          het_class     = cls,
          n             = m$n,
          R             = m$R,
          R2            = m$R2,
          RMSE          = m$RMSE,
          Bias          = m$Bias,
          Bias_pvalue   = m$Bias_pvalue,
          Slope         = m$Slope,
          Slope_pvalue  = m$Slope_pvalue,
          d_opt_source  = current_mode
        )
      }

      # Total row: pooled across the 3 classes (= entire valid sample for site)
      m_tot <- compute_pair_metrics(dt_s[[combo$lidar_col]], dt_s[[combo$s2_col]])
      rows_out[[length(rows_out) + 1L]] <- data.table::data.table(
        site          = s,
        metric_source = metric_source_lbl,
        combination   = combo$name,
        het_class     = "Total",
        n             = m_tot$n,
        R             = m_tot$R,
        R2            = m_tot$R2,
        RMSE          = m_tot$RMSE,
        Bias          = m_tot$Bias,
        Bias_pvalue   = m_tot$Bias_pvalue,
        Slope         = m_tot$Slope,
        Slope_pvalue  = m_tot$Slope_pvalue,
        d_opt_source  = current_mode
      )
    }
  }

  out_dt   <- data.table::rbindlist(rows_out, use.names = TRUE)
  out_path <- file.path(paths$output, "intermediate", "sm6",
                        paste0("heterogeneity_analysis_sm5_",
                               current_mode, ".csv"))
  data.table::fwrite(out_dt, out_path)
  cli::cli_alert_success("Written: {out_path} ({nrow(out_dt)} rows)")
}

cli::cli_h1("Done — SM5-based discrete heterogeneity tables")
cli::cli_text(
  "To plot via the existing 20_sm6_plot_heterogeneity.R, edit its loader ",
  "to read heterogeneity_analysis_sm5_<mode>.csv (or duplicate the script ",
  "and swap the input filename pattern)."
)
