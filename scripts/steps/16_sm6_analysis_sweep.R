# ---
# title:  09_sm6_analysis_sweep.R
# desc:   Orchestration — SM6b passe 2 multi-scale heterogeneity sweep.
#         Extends SM6b passe 1 (08_sm6_analysis.R) to consume all 6
#         heterogeneity variants per metric produced by SM6a passe 2:
#           block_5_m, block_10_m, block_20_m, block_50_m,
#           focal_30_m, focal_50_m
#         for each metric ∈ c("DSM", "CHM").
#
#         For each (metric_source × window_method × window_size_m) variant,
#         thresholds are recalibrated independently via calibrate_thresholds()
#         (target_n_per_class = 5000L). A class "Low" in block_10_m is not
#         the same physical region as "Low" in focal_50_m; the calibration
#         guarantees uniform statistical representation (≥5000 pixels/class/
#         site) within each variant, not physical comparability between them.
#
#         Outputs:
#           output/intermediate/sm6/thresholds_sweep.csv
#             12 rows (2 metrics × 6 variants), columns:
#             metric_source, window_method, window_size_m,
#             low_threshold, high_threshold, balance_score,
#             calibration_status,   # "ok" or "failed"
#             count_aigoual_{low,medium,high},
#             count_blois_{low,medium,high},
#             count_mormal_{low,medium,high}
#
#           output/intermediate/sm6/heterogeneity_analysis_sweep.csv
#             ≤378 rows (3 sites × 2 metrics × 6 variants × 3 combos × 3 classes);
#             fewer if calibration fails for some variants.
#             Columns: site, metric_source, window_method, window_size_m,
#             combination, het_class, n,
#             R, R2, RMSE, Bias, Bias_pvalue, Slope, Slope_pvalue
#
#         SM6b passe 1 outputs (heterogeneity_analysis.csv,
#         DSM_sd_thresholds.csv, CHM_sd_thresholds.csv) are not modified.
#
#         Prerequisites:
#           SM6a passe 2 outputs (07_sm6_compute_heterogeneity.R):
#             output/intermediate/sm6/{site}/
#               {dsm,chm}_sd_{block,focal}_{5,10,20,50,30,50}_m.tif
#           External rasters (same as SM6b passe 1 — 08_sm6_analysis.R):
#             ladstack_classic.tif, PAD_36.5_40.tif, max_res_10_m.tif,
#             s2lai_summer_atbd_res_10_m.tif,
#             s2lai_summer_best_indiv_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("scripts/09_sm6_analysis_sweep.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_load_rasters_sweep.R"))
source(here::here("R", "sm6_thresholds.R"))
source(here::here("R", "sm6_classify.R"))
source(here::here("R", "sm6_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites          <- c("Aigoual", "Blois", "Mormal")
canopy_max_m   <- 40         # fixed canopy height ceiling used in PAD naming
metric_sources <- c("DSM", "CHM")

# d_opt for the sweep: load the common (all-sites) value from prosail_opt.csv.
# The sweep uses a single d_opt across sites (common mode).
norm_ref            <- "DSM_keepTrees"
prosail_opt_csv_sw  <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv")
if (!file.exists(prosail_opt_csv_sw))
  stop("prosail_opt.csv not found — run step 12_sm5_select_prosail_opt.R first:\n  ",
       prosail_opt_csv_sw, call. = FALSE)
prosail_opt_sw <- data.table::fread(prosail_opt_csv_sw)
opt_common <- prosail_opt_sw[Norm == norm_ref & d_opt_source == "all_sites"]
if (nrow(opt_common) == 0L)
  stop("prosail_opt.csv has no all_sites row for Norm=", norm_ref,
       ". Re-run step 12 with d_opt_source 'all_sites' enabled.",
       call. = FALSE)
dopt_value <- opt_common$d_opt[[1L]]
cli::cli_alert_info("Sweep d_opt (common, norm = {norm_ref}) = {dopt_value} m")

sweep_variants <- list(
  list(method = "block", size_m = 5L),
  list(method = "block", size_m = 10L),
  list(method = "block", size_m = 20L),
  list(method = "block", size_m = 50L),
  list(method = "focal", size_m = 30L),
  list(method = "focal", size_m = 50L)
)

# Site-level uniform-LAI sampling done BEFORE threshold calibration.
n_target_per_site <- 50000L
lai_min           <- 2.0
lai_max_quantile  <- 0.98
bin_width         <- 1.0

# Heterogeneity thresholds. NULL → dynamic calibration per variant.
# Otherwise: c(low, high) applied to every (metric_source × variant).
fixed_thresholds_sweep <- c(low = 2.5, high = 5.0)
# fixed_thresholds_sweep <- NULL    # uncomment to re-enable per-variant calib.

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")
ext_dir  <- paths$ext_results
out_dir  <- file.path(paths$output, "intermediate", "sm6")

# ── Pre-flight: verify all required input files ─────────────────────────────────

pad_fn   <- pad_filename_sweep(dopt_value, canopy_max_m)
metrics  <- c("DSM", "CHM")

als_dopt_07_dir <- file.path(paths$output, "intermediate", "lai_als_dopt")

required_files <- list()

for (site in sites) {
  base_ext  <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")
  base_sm6a <- file.path(sm6a_dir, site)

  # Fixed LAI rasters (same as SM6b passe 1)
  required_files[[paste0(site, "_ladstack")]] <- list(
    path     = file.path(base_ext, "ladstack_classic.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  # LAI_ALS_dopt: prefer script-07 raster, fall back to raw PAD file.
  als_dopt_07  <- file.path(als_dopt_07_dir, site, "LAI_ALS_dopt_per_site.tif")
  pad_path_raw <- file.path(base_ext, "PAD_Profiles_dsm_keepTrees", pad_fn)
  required_files[[paste0(site, "_lai_als_dopt")]] <- list(
    path     = if (file.exists(als_dopt_07)) als_dopt_07 else pad_path_raw,
    producer = if (file.exists(als_dopt_07))
      paste0("produced by step 07 (LAI_ALS_dopt_per_site.tif, d_opt=", dopt_value, ")")
    else
      paste0("PAD_Profiles_dsm_keepTrees/", pad_fn, " (d_opt=", dopt_value, ")")
  )
  required_files[[paste0(site, "_max")]] <- list(
    path     = file.path(base_ext, "max_res_10_m.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = file.path(sm6a_dir, site, "s2lai_summer_atbd_T_res_10_m.tif"),
    producer = "produced by script 25 (prosail 3.0.0, ATBD_T)"
  )
  required_files[[paste0(site, "_s2_opt")]] <- list(
    path     = file.path(sm6a_dir, site, "s2lai_summer_opt_common_res_10_m.tif"),
    producer = "produced by script 13 (s2lai_summer_opt_common)"
  )

  # SM6a passe 2 het rasters
  for (metric in metrics) {
    for (variant in sweep_variants) {
      fn  <- paste0(
        tolower(metric), "_sd_", variant$method, "_", variant$size_m, "_m.tif"
      )
      key <- paste0(site, "_", tolower(metric), "_", variant$method,
                    "_", variant$size_m)
      required_files[[key]] <- list(
        path     = file.path(base_sm6a, fn),
        producer = "produced by SM6a passe 2 (07_sm6_compute_heterogeneity.R)"
      )
    }
  }
}

missing_msgs <- character(0)
for (key in names(required_files)) {
  entry <- required_files[[key]]
  if (!file.exists(entry$path)) {
    missing_msgs <- c(
      missing_msgs,
      paste0("  MISSING [", entry$producer, "]\n    ", entry$path)
    )
  }
}

if (length(missing_msgs) > 0L) {
  stop(
    "SM6b sweep pre-flight failed — ", length(missing_msgs),
    " file(s) missing:\n",
    paste(missing_msgs, collapse = "\n"),
    "\nRun SM6a passe 2 (07_sm6_compute_heterogeneity.R) first for het rasters.",
    call. = FALSE
  )
}

cli::cli_alert_success(
  "Pre-flight passed — all {length(required_files)} required files found."
)

# ── Output directory ───────────────────────────────────────────────────────────

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── Build multi-site wide data.table ──────────────────────────────────────────

cli::cli_h2("Building multi-site sweep data.table")

lai_als_dopt_paths <- setNames(
  file.path(paths$output, "intermediate", "lai_als_dopt", sites,
             "LAI_ALS_dopt_per_site.tif"),
  sites
)

dt_sweep <- build_multisite_dt_sweep(
  sites              = sites,
  dopt_value         = dopt_value,
  sweep_variants     = sweep_variants,
  sm6a_dir           = sm6a_dir,
  ext_dir            = ext_dir,
  lai_als_dopt_paths = lai_als_dopt_paths
)

cli::cli_alert_info(
  "Total pixels loaded: {nrow(dt_sweep)} across {length(sites)} sites, ",
  "{ncol(dt_sweep)} columns"
)

# Keep a reference to the raw dt for threshold calibration (manuscript order).
dt_sweep_raw <- dt_sweep

# Site-level uniform-LAI sampling (matches 15_sm6_analysis.R workflow)
cli::cli_h2("Site-level uniform-LAI sampling")
dt_sweep <- sample_uniform_lai_dt(
  dt                = dt_sweep,
  n_target_per_site = n_target_per_site,
  lai_min           = lai_min,
  lai_max_quantile  = lai_max_quantile,
  bin_width         = bin_width
)
cli::cli_alert_success(
  "Sampled: {paste(table(dt_sweep$site), collapse=' / ')} pixels"
)

# ── Combinations (unchanged from SM6b passe 1) ────────────────────────────────

combinations <- get_lai_combinations()

# ── Main sweep loop ────────────────────────────────────────────────────────────

t0             <- proc.time()
all_metrics    <- list()
all_thresholds <- list()
thr_idx        <- 0L

n_variants <- length(metric_sources) * length(sweep_variants)
cli::cli_h2(
  "SM6b sweep — {n_variants} calibrations ",
  "({length(metric_sources)} metrics × {length(sweep_variants)} variants)"
)

for (ms in metric_sources) {
  for (variant in sweep_variants) {
    thr_idx <- thr_idx + 1L
    method  <- variant$method
    size_m  <- variant$size_m
    het_col <- paste0(tolower(ms), "_sd_", method, "_", size_m)

    cli::cli_h3(
      "[{thr_idx}/{n_variants}]  metric={ms}  window={method}_{size_m}_m"
    )

    # ── Thresholds: fixed or dynamic balance on the RAW distribution ────────
    if (!is.null(fixed_thresholds_sweep)) {
      low_v  <- fixed_thresholds_sweep[["low"]]
      high_v <- fixed_thresholds_sweep[["high"]]
      vals   <- dt_sweep_raw[[het_col]]
      counts <- vapply(sites, function(s) {
        v <- vals[dt_sweep_raw$site == s]
        c(Low    = sum(v <  low_v,                   na.rm = TRUE),
          Medium = sum(v >= low_v  & v < high_v,     na.rm = TRUE),
          High   = sum(v >= high_v,                  na.rm = TRUE))
      }, numeric(3))
      counts <- t(counts); rownames(counts) <- sites
      thr <- list(low = low_v, high = high_v,
                  balance_score = NA_real_, counts = counts)
      cli::cli_alert_info(
        "Fixed thresholds: low = {low_v}, high = {high_v}"
      )
    } else {
      cli::cli_alert_info("Calibrating thresholds on raw distribution ...")
      thr <- tryCatch(
        calibrate_thresholds(
          df                 = dt_sweep_raw,
          metric_col         = het_col,
          sites              = sites,
          target_n_per_class = 5000L
        ),
        error = function(e) {
          cli::cli_alert_warning(
            "Calibration failed for {ms}/{method}_{size_m}_m: ",
            "{conditionMessage(e)}"
          )
          NULL
        }
      )
    }

    # ── Build threshold row ────────────────────────────────────────────────────
    if (is.null(thr)) {
      thr_row <- data.table::data.table(
        metric_source        = ms,
        window_method        = method,
        window_size_m        = size_m,
        low_threshold        = NA_real_,
        high_threshold       = NA_real_,
        balance_score        = NA_real_,
        calibration_status   = "failed",
        count_aigoual_low    = 0L,
        count_aigoual_medium = 0L,
        count_aigoual_high   = 0L,
        count_blois_low      = 0L,
        count_blois_medium   = 0L,
        count_blois_high     = 0L,
        count_mormal_low     = 0L,
        count_mormal_medium  = 0L,
        count_mormal_high    = 0L
      )
    } else {
      cli::cli_alert_success(
        "Thresholds: low = {round(thr$low, 3)},  ",
        "high = {round(thr$high, 3)},  ",
        "balance_score = {round(thr$balance_score, 4)}"
      )
      cli::cli_text("Count matrix (sites \u00d7 classes):")
      print(thr$counts)

      thr_row <- data.table::data.table(
        metric_source        = ms,
        window_method        = method,
        window_size_m        = size_m,
        low_threshold        = thr$low,
        high_threshold       = thr$high,
        balance_score        = thr$balance_score,
        calibration_status   = "ok",
        count_aigoual_low    = thr$counts["Aigoual", "Low"],
        count_aigoual_medium = thr$counts["Aigoual", "Medium"],
        count_aigoual_high   = thr$counts["Aigoual", "High"],
        count_blois_low      = thr$counts["Blois",   "Low"],
        count_blois_medium   = thr$counts["Blois",   "Medium"],
        count_blois_high     = thr$counts["Blois",   "High"],
        count_mormal_low     = thr$counts["Mormal",  "Low"],
        count_mormal_medium  = thr$counts["Mormal",  "Medium"],
        count_mormal_high    = thr$counts["Mormal",  "High"]
      )
    }
    all_thresholds[[thr_idx]] <- thr_row

    # ── Skip metric computation if calibration failed ─────────────────────────
    if (is.null(thr)) next

    # ── Classify (copy to avoid mutating dt_sweep across iterations) ──────────
    dt_var <- data.table::copy(dt_sweep)
    classify_heterogeneity(
      dt             = dt_var,
      metric_col     = het_col,
      low_threshold  = thr$low,
      high_threshold = thr$high
    )

    # ── Compute metrics ────────────────────────────────────────────────────────
    cli::cli_alert_info("Computing metrics by class ...")
    metrics_var <- compute_metrics_by_class(
      dt               = dt_var,
      combinations     = combinations,
      metric_source    = ms,
      uniform_sample_n = NULL
    )

    # Add sweep identifier columns
    data.table::set(metrics_var, j = "window_method", value = method)
    data.table::set(metrics_var, j = "window_size_m",  value = size_m)

    all_metrics[[length(all_metrics) + 1L]] <- metrics_var
    cli::cli_alert_success("  {nrow(metrics_var)} rows computed")
  }
}

# ── Combine and save results ───────────────────────────────────────────────────

thresholds_dt <- data.table::rbindlist(all_thresholds)
final_dt      <- data.table::rbindlist(all_metrics)

# Reorder columns to spec
data.table::setcolorder(final_dt, c(
  "site", "metric_source", "window_method", "window_size_m",
  "combination", "het_class", "n",
  "R", "R2", "RMSE", "Bias", "Bias_pvalue", "Slope", "Slope_pvalue"
))

thr_csv  <- file.path(out_dir, "thresholds_sweep.csv")
anal_csv <- file.path(out_dir, "heterogeneity_analysis_sweep.csv")

data.table::fwrite(thresholds_dt, thr_csv)
data.table::fwrite(final_dt,      anal_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

n_ok     <- sum(thresholds_dt[["calibration_status"]] == "ok")
n_failed <- sum(thresholds_dt[["calibration_status"]] == "failed")

cli::cli_h1("SM6b passe 2 sweep — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "Variantes calibrées avec succès : {n_ok}/{n_variants}",
  "v" = "Variantes en échec              : {n_failed}/{n_variants}",
  "v" = "Analyse CSV : {anal_csv}  ({nrow(final_dt)} lignes)",
  "v" = "Seuils CSV  : {thr_csv}  ({nrow(thresholds_dt)} lignes)"
))
cli::cli_text("Lignes par (metric_source, window_method, window_size_m, het_class) :")
print(final_dt[, .N, by = .(metric_source, window_method, window_size_m, het_class)])
