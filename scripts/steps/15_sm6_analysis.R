# ---
# title:  08_sm6_analysis.R
# desc:   Orchestration — SM6b heterogeneity class analysis.
#         Loads per-site rasters, builds a multi-site pixel data.table,
#         calibrates dynamic heterogeneity thresholds, classifies pixels,
#         and computes LAI comparison metrics per (site, class, combination,
#         metric_source).
#
#         Outputs:
#           output/intermediate/sm6/DSM_sd_thresholds.csv
#           output/intermediate/sm6/CHM_sd_thresholds.csv
#           output/intermediate/sm6/heterogeneity_analysis.csv
#
#         Prerequisites — must exist before running:
#           SM6a outputs (produced by 07_sm6_compute_heterogeneity.R):
#             output/intermediate/sm6/{site}/dsm_sd_res_10_m.tif
#             output/intermediate/sm6/{site}/chm_sd_res_10_m.tif
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
#   source("scripts/08_sm6_analysis.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_load_rasters.R"))
source(here::here("R", "sm6_dataframe.R"))
source(here::here("R", "sm6_thresholds.R"))
source(here::here("R", "sm6_classify.R"))
source(here::here("R", "sm6_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites                <- c("Aigoual", "Blois", "Mormal")
norm_ref             <- "DSM_keepTrees"
canopy_max_m         <- 40
metric_sources       <- c("DSM")
dopt_modes           <- c("per_site", "all_sites")  # produces one set of outputs per mode

# DSM_sd / CHM_sd rasters used for heterogeneity classification.
# Available variants in output/intermediate/sm6/{site}/ :
#   dsm_sd_res_10_m.tif       (10 m default → thresholds ~ 1.5 / 3.0)
#   dsm_sd_block_20_m.tif     (20 m         → thresholds ~ 2.5 / 4.5)
#   dsm_sd_block_50_m.tif     (50 m         → thresholds ~ 4.0 / 6.0)
#   dsm_sd_focal_30_m.tif / dsm_sd_focal_50_m.tif (focal variants)
# dsm_sd_filename      <- "dsm_sd_block_20_m.tif"
dsm_sd_filename      <- "dsm_sd_res_10_m.tif"
chm_sd_filename      <- "chm_sd_res_10_m.tif"

# Heterogeneity thresholds.
# fixed_thresholds = NULL → dynamic calibration via balance score (default).
# fixed_thresholds = named list, e.g.
#     list(DSM = c(low = 2.5, high = 5.0),
#          CHM = c(low = 2.5, high = 5.0))
#     → use these values directly, skip calibration.
fixed_thresholds <- list(
  DSM = c(low = 2.5, high = 5.0),
  CHM = c(low = 2.5, high = 5.0)
)
# fixed_thresholds <- NULL   # uncomment to re-enable dynamic calibration

# Per-mode helpers used inside the dopt_modes loop
s2_opt_fn_for_mode <- function(mode) {
  switch(mode,
    per_site  = "s2lai_summer_opt_per_site_res_10_m.tif",
    all_sites = "s2lai_summer_opt_common_res_10_m.tif",
    fixed_4   = "s2lai_summer_opt_fixed_4_res_10_m.tif",
    stop("Unknown dopt mode: ", mode))
}
als_dopt_fn_for_mode <- function(mode) {
  switch(mode,
    per_site  = "LAI_ALS_dopt_per_site.tif",
    all_sites = "LAI_ALS_dopt_common.tif",
    fixed_4   = "LAI_ALS_dopt_fixed_4.tif",
    stop("Unknown dopt mode: ", mode))
}

set.seed(42)

# ── Site-specific d_opt from SM5 passe 2 (prosail_opt.csv) ────────────────────
prosail_opt_csv <- here::here(
  "output", "intermediate", "sm5", "prosail_opt.csv"
)
if (!file.exists(prosail_opt_csv))
  stop("prosail_opt.csv not found — run step 12_sm5_select_prosail_opt.R first:\n  ",
       prosail_opt_csv, call. = FALSE)

prosail_opt <- data.table::fread(prosail_opt_csv)

# Site-level uniform-LAI sampling, drawn BEFORE heterogeneity thresholding.
# - n_target_per_site = 50000 — gives enough headroom for balance scoring to
#   find a (low, high) pair satisfying target_n_per_class = 5000 across all
#   (site × class) combinations.
# - lai_min = 2, upper = floor(p98_lai), 1-unit bins
n_target_per_site <- 100000L
lai_min           <- 2.0
lai_max_quantile  <- 0.98
bin_width         <- 1.0

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")
ext_dir  <- paths$ext_results
out_dir  <- file.path(paths$output, "intermediate", "sm6")

# ════════════════════════════════════════════════════════════════════════════
# Outer loop : one analysis per d_opt mode (per_site, all_sites, ...)
# Each iteration writes mode-tagged CSVs:
#   {DSM,CHM}_sd_thresholds_{mode}.csv
#   heterogeneity_analysis_{mode}.csv
# ════════════════════════════════════════════════════════════════════════════

mode_summary <- list()

for (current_mode in dopt_modes) {

cli::cli_rule(left = paste0("MODE: ", current_mode))

opt_ref <- prosail_opt[Norm == norm_ref & Site %in% sites &
                        d_opt_source == current_mode]
if (nrow(opt_ref) != length(sites))
  stop("prosail_opt.csv missing rows for some sites (norm=", norm_ref,
       ", source=", current_mode, "). Got: ",
       paste(opt_ref$Site, collapse = ", "), call. = FALSE)

dopt_by_site <- setNames(opt_ref$d_opt, opt_ref$Site)[sites]
cli::cli_alert_info(
  "d_opt ({current_mode}): {paste(names(dopt_by_site), dopt_by_site, sep='=', collapse=', ')}"
)

s2_opt_fn <- s2_opt_fn_for_mode(current_mode)
als_dopt_fn <- als_dopt_fn_for_mode(current_mode)

# ── Pre-flight: verify all required input files ─────────────────────────────

required_files <- list()

lai_als_dopt_paths <- setNames(
  file.path(paths$output, "intermediate", "lai_als_dopt", sites, als_dopt_fn),
  sites
)

for (site in sites) {
  d_val     <- dopt_by_site[[site]]
  pad_fn    <- pad_filename(d_val, canopy_max_m)
  base_ext  <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")
  base_sm6a <- file.path(sm6a_dir, site)

  # SM6a outputs (produced by 07_sm6_compute_heterogeneity.R)
  required_files[[paste0(site, "_dsm_sd_sm6a")]] <- list(
    path     = file.path(base_sm6a, dsm_sd_filename),
    producer = "produced by SM6a (07_sm6_compute_heterogeneity.R)"
  )
  required_files[[paste0(site, "_chm_sd_sm6a")]] <- list(
    path     = file.path(base_sm6a, chm_sd_filename),
    producer = "produced by SM6a (07_sm6_compute_heterogeneity.R)"
  )

  # LiDAR rasters (produced by 2.calculate_25m_metrics.R and PAD pipeline)
  required_files[[paste0(site, "_ladstack")]] <- list(
    path     = file.path(base_ext, "ladstack_classic.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  # LAI_ALS_dopt: prefer the script-07 raster (k_select rescaled);
  # fall back to the raw PAD file when present.
  als_dopt_07  <- lai_als_dopt_paths[[site]]
  pad_path_raw <- file.path(base_ext, "PAD_Profiles_dsm_keepTrees", pad_fn)
  required_files[[paste0(site, "_lai_als_dopt")]] <- list(
    path     = if (file.exists(als_dopt_07)) als_dopt_07 else pad_path_raw,
    producer = if (file.exists(als_dopt_07))
      paste0("produced by step 07 (LAI_ALS_dopt_per_site.tif, d_opt=", d_val, ")")
    else
      paste0("PAD_Profiles_dsm_keepTrees/", pad_fn, " (d_opt=", d_val, ")")
  )
  required_files[[paste0(site, "_max")]] <- list(
    path     = file.path(base_ext, "max_res_10_m.tif"),
    producer = "produced by 2.calculate_25m_metrics.R"
  )

  # S2 ATBD raster: prosail 3.0.0 ATBD_T from script 25 — no legacy fallback.
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = file.path(base_sm6a, "s2lai_summer_atbd_T_res_10_m.tif"),
    producer = "produced by script 25 (prosail 3.0.0, ATBD_T)"
  )

  # S2 opt raster: SM5 passe 3 preferred, legacy fallback tolerated
  opt_path    <- file.path(base_sm6a, s2_opt_fn)
  legacy_path <- file.path(base_ext, "s2lai_summer_best_indiv_res_10_m.tif")
  required_files[[paste0(site, "_s2_opt")]] <- list(
    path     = if (file.exists(opt_path)) opt_path else legacy_path,
    producer = if (file.exists(opt_path))
      paste0("produced by script 13 (", s2_opt_fn, ")")
    else
      "produced by 3_train_predict_prosail.R (legacy fallback)"
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

lai_s2_opt_paths <- setNames(file.path(sm6a_dir, sites, s2_opt_fn), sites)

dt_full <- build_multisite_dt(
  sites              = sites,
  dopt_by_site       = dopt_by_site,
  sm6a_dir           = sm6a_dir,
  ext_dir            = ext_dir,
  lai_s2_opt_paths   = lai_s2_opt_paths,
  lai_als_dopt_paths = lai_als_dopt_paths,
  dsm_sd_filename    = dsm_sd_filename,
  chm_sd_filename    = chm_sd_filename
)

cli::cli_alert_info(
  "Total pixels loaded: {nrow(dt_full)} across {length(sites)} sites"
)

# ── Heterogeneity thresholds on the raw distribution ──────────────────────────
# Either dynamic balance calibration (fixed_thresholds = NULL) or fixed
# user-supplied values. In both cases counts are computed on the raw dt_full
# (manuscript order).
if (is.null(fixed_thresholds)) {
  cli::cli_h2("Calibrating heterogeneity thresholds on raw distribution")
} else {
  cli::cli_h2("Using user-fixed heterogeneity thresholds")
}
het_thresholds_raw <- list()
for (ms in metric_sources) {
  het_col <- paste0(tolower(ms), "_sd")
  if (is.null(fixed_thresholds)) {
    thr <- calibrate_thresholds(
      df                 = dt_full,
      metric_col         = het_col,
      sites              = sites,
      target_n_per_class = 5000L,
      round_step         = 0.5
    )
    cli::cli_alert_success(
      "{ms}_sd: low = {thr$low}  high = {thr$high} ",
      "(balance = {round(thr$balance_score, 3)})"
    )
  } else {
    if (is.null(fixed_thresholds[[ms]]))
      stop("fixed_thresholds has no entry for metric_source '", ms, "'.")
    low_v  <- fixed_thresholds[[ms]][["low"]]
    high_v <- fixed_thresholds[[ms]][["high"]]
    vals   <- dt_full[[het_col]]
    counts <- vapply(sites, function(s) {
      v <- vals[dt_full$site == s]
      c(Low    = sum(v <  low_v,                   na.rm = TRUE),
        Medium = sum(v >= low_v  & v < high_v,     na.rm = TRUE),
        High   = sum(v >= high_v,                  na.rm = TRUE))
    }, numeric(3))
    counts <- t(counts); rownames(counts) <- sites
    thr <- list(low = low_v, high = high_v,
                balance_score = NA_real_, counts = counts)
    cli::cli_alert_success(
      "{ms}_sd: low = {thr$low}  high = {thr$high}  (fixed)"
    )
  }
  het_thresholds_raw[[ms]] <- thr
}

# ── LAI sampling strategy (always per-class) ──────────────────────────────────
# uniform_sample_n_per_class controls the per-(site × class) sampling that
# happens INSIDE compute_metrics_by_class:
#   NULL      → use all pixels in each (site × class)
#   integer   → uniform-LAI sample of that size per (site × class) in
#               [lai_min, floor(p{lai_max_quantile*100}_lai)], bin_width unit
# This applies the same way regardless of whether thresholds are dynamic
# or fixed (the threshold mode only affects WHERE the class boundaries lie).
# uniform_sample_n_per_class <- 2000L
uniform_sample_n_per_class <- NULL   # NULL = use all pixels per (site × class)

uniform_sample_n_metrics <- uniform_sample_n_per_class
if (is.null(uniform_sample_n_metrics)) {
  cli::cli_h2("LAI sampling DISABLED — using all pixels per (site × class)")
} else {
  cli::cli_h2("Per-class uniform-LAI sampling")
  cli::cli_alert_info(
    "Per (site × class): up to {uniform_sample_n_per_class} pixels, ",
    "uniform on LAI in [{lai_min}, floor(p{round(lai_max_quantile*100)}_lai)], ",
    "bin_width = {bin_width}."
  )
}

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

  # ── Use pre-calibrated thresholds (computed on raw distribution above) ────
  thr <- het_thresholds_raw[[ms]]
  cli::cli_alert_info(
    "Applying raw-calibrated thresholds: low = {thr$low} m, high = {thr$high} m"
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

  thr_path <- file.path(out_dir,
                         paste0(ms, "_sd_thresholds_", current_mode, ".csv"))
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
  metrics_ms <- compute_metrics_by_class(
    dt_ms, combinations,
    metric_source    = ms,
    uniform_sample_n = uniform_sample_n_metrics,
    lai_min          = lai_min,
    lai_max_quantile = lai_max_quantile,
    bin_width        = bin_width
  )
  all_metrics[[ms_idx]] <- metrics_ms
  cli::cli_alert_success("  {nrow(metrics_ms)} rows computed for {ms}")
}

# ── Combine and save final results ─────────────────────────────────────────────

final_dt <- data.table::rbindlist(all_metrics)
final_dt[, d_opt_source := current_mode]

# Expected 3 sites × |metric_sources| × 3 combos × (3 classes + 1 Total) rows
cli::cli_alert_info("Final table ({current_mode}): {nrow(final_dt)} rows")

out_csv <- file.path(out_dir,
                      paste0("heterogeneity_analysis_", current_mode, ".csv"))
data.table::fwrite(final_dt, out_csv)

mode_summary[[current_mode]] <- list(
  out_csv  = out_csv,
  thr_csvs = threshold_paths,
  n_rows   = nrow(final_dt)
)

}  # end for current_mode

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("SM6b — résumé (tous modes)")
cli::cli_alert_success("Durée totale : {elapsed} s")
for (m in names(mode_summary)) {
  s <- mode_summary[[m]]
  cli::cli_bullets(c(
    "v" = "[{m}] {s$out_csv}  ({s$n_rows} lignes)",
    "v" = "       seuils: {paste(s$thr_csvs, collapse = ' | ')}"
  ))
}
