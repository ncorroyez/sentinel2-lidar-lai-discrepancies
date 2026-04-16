# ---
# title:  08_sm6_analysis.R
# desc:   Orchestration — SM6b heterogeneity class analysis.
#         Loads per-site rasters, builds a multi-site pixel data.table,
#         calibrates dynamic heterogeneity thresholds, classifies pixels,
#         and computes LAI comparison metrics per (site, class, combination,
#         metric_source).
#
#         Outputs:
#           revision/output/intermediate/sm6/DSM_sd_thresholds.csv
#           revision/output/intermediate/sm6/CHM_sd_thresholds.csv
#           revision/output/intermediate/sm6/heterogeneity_analysis.csv
#
#         Prerequisites — must exist before running:
#           SM6a outputs (produced by 07_sm6_compute_heterogeneity.R):
#             revision/output/intermediate/sm6/{site}/dsm_sd_res_10_m.tif
#             revision/output/intermediate/sm6/{site}/chm_sd_res_10_m.tif
#           External rasters (produced by 2.calculate_25m_metrics.R):
#             03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_classic.tif
#             03_RESULTS/{site}/Metrics/Deciduous_Only/max_res_10_m.tif
#             03_RESULTS/{site}/Metrics/Deciduous_Only/dsm_sd_res_10_m.tif
#                                                       (legacy ref, not read)
#           External rasters (produced by 3_train_predict_prosail.R):
#             03_RESULTS/{site}/Metrics/Deciduous_Only/
#               s2lai_summer_atbd_res_10_m.tif
#               s2lai_summer_best_indiv_res_10_m.tif
#           PAD raster (produced by 2.calculate_25m_metrics.R / PAD pipeline):
#             03_RESULTS/{site}/Metrics/Deciduous_Only/
#               PAD_Profiles_dsm_keepTrees/PAD_{X}_{Y}.tif
#             where X = canopy_max_m - dopt_value + 0.5, Y = canopy_max_m.
#             For dopt_value = 4: PAD_36.5_40.tif.
#
#         Note on lai_s2_common (s2lai_summer_depth_study_common_res_10_m.tif):
#           Not required here — it is only used in the 'Sites combined'
#           analysis (legacy lines 316-328), which is out of scope for SM6b
#           passe 1.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/08_sm6_analysis.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "sm6_load_rasters.R"))
source(here::here("revision", "R", "sm6_dataframe.R"))
source(here::here("revision", "R", "sm6_thresholds.R"))
source(here::here("revision", "R", "sm6_classify.R"))
source(here::here("revision", "R", "sm6_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites          <- c("Aigoual", "Blois", "Mormal")
dopt_value     <- 4          # d_opt in metres; see pad_filename() for convention
canopy_max_m   <- 40         # fixed canopy height ceiling used in PAD naming
metric_sources <- c("DSM", "CHM")

# downsample_n: max pixels per (site x combination x class) group.
# NULL  → use all pixels (default, recommended for the refactored pipeline).
# 5000  → replicates the legacy slice_sample(n = 5000) strategy
#          (three_factors_analysis_final.R lines 229-234) for direct
#          comparability with the submitted manuscript results.
downsample_n <- NULL

sm6a_dir <- here::here("revision", "output", "intermediate", "sm6")
ext_dir  <- here::here("03_RESULTS")
out_dir  <- here::here("revision", "output", "intermediate", "sm6")

# ── Pre-flight: verify all required input files ─────────────────────────────────

pad_fn   <- pad_filename(dopt_value, canopy_max_m)
dec_only <- file.path("03_RESULTS", "{site}", "Metrics", "Deciduous_Only")

required_files <- list()

for (site in sites) {
  base_ext  <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")
  base_sm6a <- file.path(sm6a_dir, site)

  # SM6a outputs (produced by 07_sm6_compute_heterogeneity.R)
  required_files[[paste0(site, "_dsm_sd_sm6a")]] <- list(
    path     = file.path(base_sm6a, "dsm_sd_res_10_m.tif"),
    producer = "produced by SM6a (07_sm6_compute_heterogeneity.R)"
  )
  required_files[[paste0(site, "_chm_sd_sm6a")]] <- list(
    path     = file.path(base_sm6a, "chm_sd_res_10_m.tif"),
    producer = "produced by SM6a (07_sm6_compute_heterogeneity.R)"
  )

  # LiDAR rasters (produced by 2.calculate_25m_metrics.R and PAD pipeline)
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

  # S2 rasters (produced by 3_train_predict_prosail.R)
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = file.path(base_ext, "s2lai_summer_atbd_res_10_m.tif"),
    producer = "produced by 3_train_predict_prosail.R"
  )
  required_files[[paste0(site, "_s2_opt")]] <- list(
    path     = file.path(base_ext, "s2lai_summer_best_indiv_res_10_m.tif"),
    producer = "produced by 3_train_predict_prosail.R"
  )
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
    "SM6b pre-flight failed — ", length(missing_msgs),
    " file(s) missing:\n",
    paste(missing_msgs, collapse = "\n"),
    "\nRun SM6a (07_sm6_compute_heterogeneity.R) first for SM6a outputs.",
    call. = FALSE
  )
}

cli::cli_alert_success("Pre-flight passed — all {length(required_files)} required files found.")

# ── Output directory ───────────────────────────────────────────────────────────

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# ── Build multi-site data.table ────────────────────────────────────────────────

cli::cli_h2("Building multi-site pixel data.table")

dt_full <- build_multisite_dt(
  sites      = sites,
  dopt_value = dopt_value,
  sm6a_dir   = sm6a_dir,
  ext_dir    = ext_dir
)

cli::cli_alert_info(
  "Total pixels loaded: {nrow(dt_full)} across {length(sites)} sites"
)

# ── Combinations ───────────────────────────────────────────────────────────────

combinations <- get_lai_combinations()

# ── Loop over metric sources ───────────────────────────────────────────────────

t0              <- proc.time()
all_metrics     <- vector("list", length(metric_sources))
threshold_paths <- character(length(metric_sources))

for (ms_idx in seq_along(metric_sources)) {
  ms       <- metric_sources[ms_idx]
  ms_col   <- tolower(ms)  # "dsm" → "dsm_sd", "chm" → "chm_sd"
  het_col  <- paste0(ms_col, "_sd")

  cli::cli_h2("metric_source = {ms}  (column: {het_col})")

  # ── Calibrate thresholds ───────────────────────────────────────────────────
  cli::cli_alert_info("Calibrating thresholds on pooled distribution ...")
  thr <- calibrate_thresholds(
    df                 = dt_full,
    metric_col         = het_col,
    sites              = sites,
    target_n_per_class = 5000L
  )
  cli::cli_alert_success(
    "Thresholds: low = {round(thr$low, 3)},  high = {round(thr$high, 3)},",
    "  balance_score = {round(thr$balance_score, 4)}"
  )
  cli::cli_text("Count matrix (sites \u00d7 classes):")
  print(thr$counts)

  # ── Save thresholds CSV ────────────────────────────────────────────────────
  thr_dt <- data.table::data.table(
    metric_source            = ms,
    low_threshold            = thr$low,
    high_threshold           = thr$high,
    balance_score            = thr$balance_score,
    count_aigoual_low        = thr$counts["Aigoual", "Low"],
    count_aigoual_medium     = thr$counts["Aigoual", "Medium"],
    count_aigoual_high       = thr$counts["Aigoual", "High"],
    count_blois_low          = thr$counts["Blois",   "Low"],
    count_blois_medium       = thr$counts["Blois",   "Medium"],
    count_blois_high         = thr$counts["Blois",   "High"],
    count_mormal_low         = thr$counts["Mormal",  "Low"],
    count_mormal_medium      = thr$counts["Mormal",  "Medium"],
    count_mormal_high        = thr$counts["Mormal",  "High"]
  )

  thr_path <- file.path(out_dir, paste0(ms, "_sd_thresholds.csv"))
  data.table::fwrite(thr_dt, thr_path)
  threshold_paths[ms_idx] <- thr_path
  cli::cli_alert_success("Thresholds written: {thr_path}")

  # ── Classify (copy to avoid mutating dt_full across iterations) ────────────
  dt_ms <- data.table::copy(dt_full)
  classify_heterogeneity(
    dt             = dt_ms,
    metric_col     = het_col,
    low_threshold  = thr$low,
    high_threshold = thr$high
  )

  # ── Compute metrics ────────────────────────────────────────────────────────
  cli::cli_alert_info("Computing metrics by class ...")
  metrics_ms <- compute_metrics_by_class(dt_ms, combinations,
                                         metric_source = ms,
                                         downsample_n  = downsample_n)
  all_metrics[[ms_idx]] <- metrics_ms
  cli::cli_alert_success("  {nrow(metrics_ms)} rows computed for {ms}")
}

# ── Combine and save final results ─────────────────────────────────────────────

final_dt <- data.table::rbindlist(all_metrics)

# Expected 3 sites × 2 metric_sources × 3 combos × 3 classes = 54 rows
# (rows with n < 10 are present but have NA metrics)
cli::cli_alert_info("Final table: {nrow(final_dt)} rows")

out_csv <- file.path(out_dir, "heterogeneity_analysis.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("SM6b — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "Analyse CSV : {out_csv}  ({nrow(final_dt)} lignes)",
  "v" = "Seuils DSM  : {threshold_paths[1]}",
  "v" = "Seuils CHM  : {threshold_paths[2]}"
))
cli::cli_text("Lignes par (metric_source, het_class) :")
print(final_dt[, .N, by = .(metric_source, het_class)])
