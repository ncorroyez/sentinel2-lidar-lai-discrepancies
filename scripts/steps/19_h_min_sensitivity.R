# ---
# title:  12_h_min_sensitivity.R
# desc:   Orchestration — h_min threshold sensitivity analysis.
#         Computes how LAI_ALS descriptive statistics and LAI_S2_ATBD vs
#         LAI_ALS comparison metrics vary with the minimum vegetation height
#         threshold h_min ∈ {2, 3, 5} m, via analytical partial-band
#         summation of the existing LAD profile stack without relaunching
#         the LiDAR pipeline.
#
#         Method (see h_min_sensitivity.R for full description):
#           LAI_ALS(h_min) ≈ sum of LAD layers with bin centre ≥ h_min + 0.5 m
#           This is an upper-bound approximation; the exact sensitivity would
#           require relaunching gap-fraction computation with a different h_min.
#           The approximation is flagged explicitly in the response letter.
#
#         Special case verified by equivalence check:
#           compute_lai_als_at_h_min(blois, h_min = 2) == sum(ladstack)
#           (all 38 bins included when h_min = 2, identical to reference LAI)
#
#         Addresses reviewer comment R3.minor.4 (h_min threshold sensitivity;
#         fCover sensitivity is handled in 11_fcover_sensitivity.R).
#
#         Output:
#           output/intermediate/reviewers/h_min_sensitivity_atbd.csv
#           9 rows (3 sites × 3 h_min values), 16 columns:
#             site, h_min_value, n_pixels,
#             lai_als_{mean,median,p95,max},
#             lai_s2_atbd_{mean,median,p95,max},
#             R, R2, RMSE, Bias, Slope
#
#         Input rasters (Deciduous_Only, produced by legacy pipeline):
#           03_RESULTS/{site}/Metrics/Deciduous_Only/ladstack_classic.tif
#           03_RESULTS/{site}/Metrics/Deciduous_Only/s2lai_summer_atbd_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("scripts/12_h_min_sensitivity.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "h_min_sensitivity.R"))
source(here::here("R", "k_sensitivity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites        <- c("Aigoual", "Blois", "Mormal")
h_min_values <- c(2, 3, 4, 5)    # vegetation height thresholds to test (metres)

# k rescaling — PADs were precomputed with k_ref; analysis uses k_select.
k_ref        <- 0.5
k_select     <- 0.65

ext_dir <- paths$ext_results
out_dir <- file.path(paths$output, "intermediate", "reviewers")

# ── Input path helpers ─────────────────────────────────────────────────────────

ladstack_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Deciduous_Only", "ladstack_classic.tif")
}

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")

s2_atbd_path <- function(site) {
  file.path(sm6a_dir, site, "s2lai_summer_atbd_T_res_10_m.tif")
}

# ── Pre-flight: verify all required input files ─────────────────────────────────
#
# 2 rasters × 3 sites = 6 files (same inputs as 10_k_sensitivity.R)

required_files <- list()

for (site in sites) {
  required_files[[paste0(site, "_ladstack")]] <- list(
    path     = ladstack_path(site),
    producer = "produced by 2.calculate_25m_metrics.R"
  )
  required_files[[paste0(site, "_s2_atbd")]] <- list(
    path     = s2_atbd_path(site),
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
    "h_min_sensitivity pre-flight failed — ", length(missing_msgs),
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

# ── Equivalence check: compute_lai_als_at_h_min(blois, h_min=2) == sum(ladstack)
#
# At h_min = 2 m all 38 bins (centres 2.5–39.5 m) satisfy the keep condition,
# so the partial-band sum must equal the full sum. Both rasters must have the
# same dimensions and identical non-NA pixel values.

cli::cli_h2(
  "Equivalence check — compute_lai_als_at_h_min(Blois, h_min=2) vs sum(ladstack)"
)

blois_ladstack_path <- ladstack_path("Blois")
blois_ladstack      <- terra::rast(blois_ladstack_path)
lai_full            <- sum(blois_ladstack, na.rm = TRUE)
lai_h2              <- compute_lai_als_at_h_min(blois_ladstack_path, h_min = 2)

full_vals <- as.numeric(terra::values(lai_full))
h2_vals   <- as.numeric(terra::values(lai_h2))

max_abs_diff <- max(abs(full_vals - h2_vals), na.rm = TRUE)
equiv_ok     <- isTRUE(all.equal(full_vals, h2_vals, tolerance = 1e-6))

if (!equiv_ok) {
  stop(
    "Equivalence check FAILED — compute_lai_als_at_h_min(h_min=2) ",
    "does not equal sum(ladstack).\n",
    "  Max absolute difference: ", max_abs_diff, "\n",
    "  Check get_layer_heights() band-name parsing for Blois ladstack.",
    call. = FALSE
  )
}

cli::cli_alert_success(
  "Equivalence check PASSED — h_min=2 matches sum(ladstack) ",
  "(max|diff| = {format(max_abs_diff, scientific = TRUE, digits = 3)})."
)

# ── Main loop: site × h_min ────────────────────────────────────────────────────

t0       <- proc.time()
n_total  <- length(sites) * length(h_min_values)
all_rows <- vector("list", n_total)
row_idx  <- 0L

cli::cli_h2(
  "h_min sensitivity — {n_total} combinations ",
  "({length(sites)} sites × {length(h_min_values)} h_min values)"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")

  # LAI_S2_ATBD is constant across h_min — load once per site
  lai_s2_atbd_rast <- terra::rast(s2_atbd_path(site))

  for (h_min in h_min_values) {
    row_idx <- row_idx + 1L
    cli::cli_alert_info("  h_min = {h_min} m")

    lai_als_rast <- compute_lai_als_at_h_min(
      ladstack_path = ladstack_path(site),
      h_min         = h_min
    )
    lai_als_rast <- rescale_lai_for_k(lai_als_rast, k_ref = k_ref,
                                       k_new = k_select)

    all_rows[[row_idx]] <- compute_h_min_sensitivity_metrics(
      lai_als_rast      = lai_als_rast,
      lai_s2_atbd_rast  = lai_s2_atbd_rast,
      site              = site,
      h_min_value       = h_min
    )

    cli::cli_alert_success(
      "    n={all_rows[[row_idx]]$n_pixels}  ",
      "LAI_ALS_mean={round(all_rows[[row_idx]]$lai_als_mean, 3)}  ",
      "R2={round(all_rows[[row_idx]]$R2, 3)}"
    )
  }
}

# ── Combine, reorder columns, write ───────────────────────────────────────────

final_dt <- data.table::rbindlist(all_rows)

data.table::setcolorder(final_dt, c(
  "site", "h_min_value", "n_pixels",
  "lai_als_mean", "lai_als_median", "lai_als_p95", "lai_als_max",
  "lai_s2_atbd_mean", "lai_s2_atbd_median", "lai_s2_atbd_p95", "lai_s2_atbd_max",
  "R", "R2", "RMSE", "Bias", "Slope"
))

out_csv <- file.path(out_dir, "h_min_sensitivity_atbd.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("h_min sensitivity — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "CSV   : {out_csv}  ({nrow(final_dt)} lignes)"
))
cli::cli_text("Résultats (pivot site × h_min_value) :")
print(
  final_dt[, .(site, h_min_value, n_pixels,
               lai_als_mean = round(lai_als_mean, 3),
               R2           = round(R2, 3),
               RMSE         = round(RMSE, 3),
               Bias         = round(Bias, 3),
               Slope        = round(Slope, 3))]
)
