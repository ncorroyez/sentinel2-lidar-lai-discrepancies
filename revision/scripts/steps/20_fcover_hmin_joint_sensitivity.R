# ---
# title:  20_fcover_hmin_joint_sensitivity.R
# desc:   Orchestration — joint fCover × h_min sensitivity analysis.
#         Crosses fCover threshold ∈ {0.80, 0.90, 0.95} with minimum
#         vegetation height h_min ∈ {2, 3, 4, 5} m to produce a 3×4 grid
#         of LAI_ALS vs LAI_S2_ATBD comparison metrics per site.
#
#         Input rasters (per site):
#           03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/ (4 component ENVIs)
#           03_RESULTS/{site}/Metrics/Raw/ladstack.tif
#           03_RESULTS/{site}/Metrics/Raw/s2lai_summer_atbd_res_10_m.tif
#
#         Method:
#           1. Load Raw/ladstack.tif; project once to the mask grid (built at
#              the reference threshold p = 0.90).
#           2. For each h_min: sum layers with bin centre >= h_min + 0.5 m
#              (partial-band approximation; see h_min_sensitivity.R header).
#           3. For each fCover threshold: build the 4-component mask.
#           4. For each (h_min, fCover) pair: apply mask, compute metrics.
#
#         Approximation note: h_min is applied via partial-band summation on
#         the projected Raw ladstack, not by relaunching the LiDAR pipeline.
#         This is the same upper-bound approximation as in 19_h_min_sensitivity.R.
#
#         Addresses reviewer comment R3.minor.4 (joint fCover/h_min
#         sensitivity as a complement to the two independent 1D analyses).
#
#         Output:
#           revision/output/intermediate/reviewers/
#             joint_fcover_hmin_sensitivity_atbd.csv
#           36 rows (3 sites × 3 fCover × 4 h_min), 17 columns:
#             site, fcover_threshold, h_min_value, n_pixels,
#             lai_als_{mean,median,p95,max},
#             lai_s2_atbd_{mean,median,p95,max},
#             R, R2, RMSE, Bias, Slope
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/20_fcover_hmin_joint_sensitivity.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "fcover_sensitivity.R"))
source(here::here("revision", "R", "h_min_sensitivity.R"))
source(here::here("revision", "R", "joint_sensitivity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites             <- c("Aigoual", "Blois", "Mormal")
fcover_thresholds <- c(0.80, 0.90, 0.95)
h_min_values      <- c(2, 3, 4, 5)

ext_dir <- paths$ext_results
out_dir <- file.path(paths$output, "intermediate", "reviewers")

# ── Input path helpers ─────────────────────────────────────────────────────────

masks_dir_path <- function(site) {
  file.path(ext_dir, site, "LiDAR", "Heterogeneity_Masks")
}

ladstack_raw_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Raw", "ladstack.tif")
}

s2_atbd_raw_path <- function(site) {
  file.path(paths$output, "intermediate", "sm6",
            site, "s2lai_summer_atbd_T_raw_res_10_m.tif")
}

# ── Pre-flight: verify all required input files ─────────────────────────────────

mask_components <- c(
  "cloud_mask_res_10_m.envi",
  "site_edges_mask_res_10_m.envi",
  "deciduous_only_mask_res_10_m.envi",
  "chm_1_m_thresholded_average_to_res_10_m.envi"
)

required_files <- list()

for (site in sites) {
  mdir <- masks_dir_path(site)
  for (comp in mask_components) {
    key <- paste0(site, "_mask_", gsub("\\.envi$", "", comp))
    required_files[[key]] <- list(
      path     = file.path(mdir, comp),
      producer = "produced by legacy mask pipeline"
    )
  }
  required_files[[paste0(site, "_ladstack_raw")]] <- list(
    path     = ladstack_raw_path(site),
    producer = "produced by 2.calculate_25m_metrics.R (Metrics/Raw output)"
  )
  required_files[[paste0(site, "_s2_atbd_raw")]] <- list(
    path     = s2_atbd_raw_path(site),
    producer = "produced by 3_train_predict_prosail.R (Metrics/Raw output)"
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
    "joint_sensitivity pre-flight failed — ", length(missing_msgs),
    " file(s) missing:\n",
    paste(missing_msgs, collapse = "\n"),
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

# ── Main loop: site × h_min × fcover_threshold ────────────────────────────────
#
# Order: site outer, h_min middle, fCover inner.
# The Raw ladstack is projected once per site (at reference p = 0.90);
# then each h_min sum is computed once per (site, h_min) pair;
# then the fCover mask is applied for each (site, h_min, fCover) combination.

t0      <- proc.time()
n_total <- length(sites) * length(h_min_values) * length(fcover_thresholds)
all_rows <- vector("list", n_total)
row_idx  <- 0L

cli::cli_h2(
  "Joint fCover × h_min sensitivity — {n_total} combinations ",
  "({length(sites)} sites × {length(h_min_values)} h_min × ",
  "{length(fcover_thresholds)} fCover)"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")

  mdir <- masks_dir_path(site)

  # Reference mask at p = 0.90 defines the projection target grid.
  # All fCover masks share this grid (same 4 component files).
  template_mask <- build_fcover_mask(0.90, mdir)

  cli::cli_alert_info("  Loading and projecting Raw LAI rasters ...")

  # Project the full ladstack once; h_min selection operates on the
  # projected multi-layer raster to avoid repeated disk I/O.
  ladstack_raw  <- terra::rast(ladstack_raw_path(site))
  ladstack_proj <- terra::project(ladstack_raw, template_mask)

  s2_raw     <- terra::rast(s2_atbd_raw_path(site))
  lai_s2_raw <- terra::project(s2_raw, template_mask)

  cli::cli_alert_success(
    "  Projected: {terra::nlyr(ladstack_proj)} LAD layers, ",
    "{terra::ncell(template_mask)} cells"
  )

  for (h_min in h_min_values) {
    cli::cli_alert_info("  h_min = {h_min} m")

    lai_als_hmin <- sum_lai_als_at_h_min_from_raster(ladstack_proj, h_min)

    for (threshold in fcover_thresholds) {
      row_idx  <- row_idx + 1L
      mask_p   <- build_fcover_mask(threshold, mdir)

      all_rows[[row_idx]] <- compute_joint_sensitivity_metrics(
        lai_als_rast     = lai_als_hmin,
        lai_s2_atbd_rast = lai_s2_raw,
        mask_rast        = mask_p,
        site             = site,
        fcover_threshold = threshold,
        h_min_value      = h_min
      )

      r <- all_rows[[row_idx]]
      cli::cli_alert_success(
        "    fCover={threshold}  n={r$n_pixels}  ",
        "LAI_ALS={round(r$lai_als_mean, 3)}  R2={round(r$R2, 3)}"
      )
    }
  }
}

# ── Combine, reorder columns, write ───────────────────────────────────────────

final_dt <- data.table::rbindlist(all_rows)

data.table::setcolorder(final_dt, c(
  "site", "fcover_threshold", "h_min_value", "n_pixels",
  "lai_als_mean", "lai_als_median", "lai_als_p95", "lai_als_max",
  "lai_s2_atbd_mean", "lai_s2_atbd_median", "lai_s2_atbd_p95", "lai_s2_atbd_max",
  "R", "R2", "RMSE", "Bias", "Slope"
))

out_csv <- file.path(out_dir, "joint_fcover_hmin_sensitivity_atbd.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("Joint fCover × h_min sensitivity — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "CSV   : {out_csv}  ({nrow(final_dt)} lignes)"
))
cli::cli_text("Résultats (site / h_min / fCover) :")
print(
  final_dt[, .(site, fcover_threshold, h_min_value, n_pixels,
               lai_als_mean = round(lai_als_mean, 3),
               R2           = round(R2, 3),
               RMSE         = round(RMSE, 3),
               Bias         = round(Bias, 3),
               Slope        = round(Slope, 3))]
)
