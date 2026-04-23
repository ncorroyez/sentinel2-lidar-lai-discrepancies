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
#           revision/output/intermediate/sm6/thresholds_sweep.csv
#             12 rows (2 metrics × 6 variants), columns:
#             metric_source, window_method, window_size_m,
#             low_threshold, high_threshold, balance_score,
#             calibration_status,   # "ok" or "failed"
#             count_aigoual_{low,medium,high},
#             count_blois_{low,medium,high},
#             count_mormal_{low,medium,high}
#
#           revision/output/intermediate/sm6/heterogeneity_analysis_sweep.csv
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
#             revision/output/intermediate/sm6/{site}/
#               {dsm,chm}_sd_{block,focal}_{5,10,20,50,30,50}_m.tif
#           External rasters (same as SM6b passe 1 — 08_sm6_analysis.R):
#             ladstack_classic.tif, PAD_36.5_40.tif, max_res_10_m.tif,
#             s2lai_summer_atbd_res_10_m.tif,
#             s2lai_summer_best_indiv_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/09_sm6_analysis_sweep.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "sm6_load_rasters_sweep.R"))
source(here::here("revision", "R", "sm6_thresholds.R"))
source(here::here("revision", "R", "sm6_classify.R"))
source(here::here("revision", "R", "sm6_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites          <- c("Aigoual", "Blois", "Mormal")
dopt_value     <- 4          # d_opt in metres; see pad_filename_sweep()
canopy_max_m   <- 40         # fixed canopy height ceiling used in PAD naming
metric_sources <- c("DSM", "CHM")

sweep_variants <- list(
  list(method = "block", size_m = 5L),
  list(method = "block", size_m = 10L),
  list(method = "block", size_m = 20L),
  list(method = "block", size_m = 50L),
  list(method = "focal", size_m = 30L),
  list(method = "focal", size_m = 50L)
)

# downsample_n: max pixels per (site × combination × class) group.
# NULL  → use all pixels (default, recommended).
# 5000  → replicates legacy slice_sample(n = 5000) for direct comparability.
downsample_n <- NULL

sm6a_dir <- here::here("revision", "output", "intermediate", "sm6")
ext_dir  <- here::here("03_RESULTS")
out_dir  <- here::here("revision", "output", "intermediate", "sm6")

# ── Pre-flight: verify all required input files ─────────────────────────────────

pad_fn   <- pad_filename_sweep(dopt_value, canopy_max_m)
metrics  <- c("DSM", "CHM")

required_files <- list()

for (site in sites) {
  base_ext  <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")
  base_sm6a <- file.path(sm6a_dir, site)

  # Fixed LAI rasters (same as SM6b passe 1)
  required_files[[paste0(site, "_ladstack")]] <- list(
    path     = file.path(base_ext, "ladstack_classic.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  required_files[[paste0(site, "_pad")]] <- list(
    path     = file.path(base_ext, "PAD_Profiles_dsm_keepTrees", pad_fn),
    producer = "produced by 2.calculate_25m_metrics.R (PAD pipeline)"
  )
  required_files[[paste0(site, "_max")]] <- list(
    path     = file.path(base_ext, "max_res_10_m.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = file.path(base_ext, "s2lai_summer_atbd_res_10_m.tif"),
    producer = "produced by 3_train_predict_prosail.R"
  )
  required_files[[paste0(site, "_s2_opt")]] <- list(
    path     = file.path(base_ext, "s2lai_summer_best_indiv_res_10_m.tif"),
    producer = "produced by 3_train_predict_prosail.R"
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

dt_sweep <- build_multisite_dt_sweep(
  sites          = sites,
  dopt_value     = dopt_value,
  sweep_variants = sweep_variants,
  sm6a_dir       = sm6a_dir,
  ext_dir        = ext_dir
)

cli::cli_alert_info(
  "Total pixels loaded: {nrow(dt_sweep)} across {length(sites)} sites, ",
  "{ncol(dt_sweep)} columns"
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

    # ── Calibrate thresholds (tryCatch: only failed calibrations are skipped) ─
    cli::cli_alert_info("Calibrating thresholds ...")
    thr <- tryCatch(
      calibrate_thresholds(
        df                 = dt_sweep,
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
      dt            = dt_var,
      combinations  = combinations,
      metric_source = ms,
      downsample_n  = downsample_n
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
