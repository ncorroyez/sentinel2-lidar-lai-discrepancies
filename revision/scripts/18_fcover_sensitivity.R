# ---
# title:  11_fcover_sensitivity.R
# desc:   Orchestration — fCover threshold sensitivity analysis.
#         Computes how LAI_ALS descriptive statistics and LAI_S2_ATBD vs
#         LAI_ALS comparison metrics vary with fCover threshold p ∈
#         {0.80, 0.90, 0.95}, by reconstructing the deciduous-forest mask at
#         each threshold from its four component rasters without relaunching
#         the LiDAR pipeline.
#
#         Mask formula (from functions_create_masks.R lines 753-777):
#           mask_p = cloud × edges × deciduous × I(fcover_continuous > p)
#           mask_p[mask_p == 0] ← NA
#
#         LAI_ALS source: Metrics/Raw/ladstack.tif — the unmasked LAD profile
#         stack produced by the LiDAR pipeline. Same underlying data as
#         Deciduous_Only/ladstack_classic.tif (identical 38 LAD layers, band
#         order reversed but irrelevant for the sum). Projected to the mask
#         grid via terra::project() to handle minor extent differences
#         (notable for Mormal: Raw has 1270 columns, mask grid has 1260).
#
#         Equivalence check (Blois, p = 0.90): verifies that
#         build_fcover_mask(0.90, blois_masks_dir) reproduces the legacy mask
#         artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi
#         pixel-for-pixel. The script stops with details if the check fails.
#
#         Addresses reviewer comment R3.minor.4 (fCover threshold sensitivity;
#         h_min sensitivity is handled qualitatively in the response letter).
#
#         Output:
#           revision/output/intermediate/reviewers/fcover_sensitivity_atbd.csv
#           9 rows (3 sites × 3 thresholds), 16 columns:
#             site, fcover_threshold, n_pixels,
#             lai_als_{mean,median,p95,max},
#             lai_s2_atbd_{mean,median,p95,max},
#             R, R2, RMSE, Bias, Slope
#
#         Input rasters (per site):
#           03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/ (4 component ENVIs)
#           03_RESULTS/{site}/Metrics/Raw/ladstack.tif
#           03_RESULTS/{site}/Metrics/Raw/s2lai_summer_atbd_res_10_m.tif
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/11_fcover_sensitivity.R")
# ---

library(here)
library(terra)
library(data.table)
library(cli)

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "fcover_sensitivity.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites             <- c("Aigoual", "Blois", "Mormal")
fcover_thresholds <- c(0.80, 0.90, 0.95)

ext_dir <- paths$ext_results
out_dir <- here::here("revision", "output", "intermediate", "reviewers")

# ── Input path helpers ─────────────────────────────────────────────────────────

masks_dir_path <- function(site) {
  file.path(ext_dir, site, "LiDAR", "Heterogeneity_Masks")
}

ladstack_raw_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Raw", "ladstack.tif")
}

s2_atbd_raw_path <- function(site) {
  file.path(ext_dir, site, "Metrics", "Raw", "s2lai_summer_atbd_res_10_m.tif")
}

legacy_mask_path <- function(site) {
  file.path(
    masks_dir_path(site),
    "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"
  )
}

# ── Pre-flight: verify all required input files ─────────────────────────────────
#
# 4 mask components × 3 sites = 12 files
# 2 LAI rasters × 3 sites     =  6 files
# 1 legacy mask (Blois only)   =  1 file
# Total: 19 files

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

required_files[["blois_legacy_mask"]] <- list(
  path     = legacy_mask_path("Blois"),
  producer = "produced by legacy mask pipeline at p = 0.90"
)

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
    "fcover_sensitivity pre-flight failed — ", length(missing_msgs),
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

# ── Equivalence check: build_fcover_mask(0.9) vs legacy mask (Blois) ──────────
#
# Constructs the mask at p = 0.90 from the 4 component rasters and compares
# it pixel-for-pixel with the pre-computed legacy mask. Two sequential checks:
#   (1) grid geometry must be identical
#   (2) NA positions must be identical (same pixels masked)
#
# If either check fails, the script stops with a detailed error message.

cli::cli_h2(
  "Equivalence check — build_fcover_mask(0.90, Blois) vs legacy mask"
)

blois_masks_dir <- masks_dir_path("Blois")
mask_blois_90   <- build_fcover_mask(0.90, blois_masks_dir)
legacy_mask     <- terra::rast(legacy_mask_path("Blois"))

# Check 1: grid geometry (extent, resolution, CRS, dimensions)
geom_ok <- tryCatch(
  terra::compareGeom(mask_blois_90, legacy_mask,
                     lyrs = FALSE, crs = TRUE, ext = TRUE,
                     rowcol = TRUE, res = TRUE),
  error = function(e) FALSE
)

if (!isTRUE(geom_ok)) {
  stop(
    "Equivalence check FAILED — geometry mismatch:\n",
    "  build_fcover_mask(0.90) : ",
    terra::nrow(mask_blois_90), " rows × ", terra::ncol(mask_blois_90), " cols\n",
    "  legacy mask             : ",
    terra::nrow(legacy_mask), " rows × ", terra::ncol(legacy_mask), " cols\n",
    "The reconstructed mask does not share the grid of the legacy mask. ",
    "Investigate mask component grids.",
    call. = FALSE
  )
}

# Check 2: NA positions (which pixels are masked)
rec_vals <- as.numeric(terra::values(mask_blois_90))
leg_vals <- as.numeric(terra::values(legacy_mask))

rec_na <- is.na(rec_vals)
leg_na <- is.na(leg_vals)

n_pos_diff      <- sum(rec_na != leg_na)
n_reconstructed <- sum(!rec_na)
n_legacy        <- sum(!leg_na)

if (n_pos_diff > 0L) {
  stop(
    "Equivalence check FAILED — NA positions differ:\n",
    "  build_fcover_mask(0.90) : ", n_reconstructed, " non-NA pixels\n",
    "  legacy mask             : ", n_legacy, " non-NA pixels\n",
    "  pixels with differing NA status : ", n_pos_diff, "\n",
    "The 4-component formula does not exactly reproduce the legacy mask. ",
    "The fCover continuous raster or another component may differ from ",
    "the one used to build the legacy mask. Investigate before continuing.",
    call. = FALSE
  )
}

cli::cli_alert_success(
  "Equivalence check PASSED — {n_legacy} non-NA pixels, ",
  "NA positions identical between reconstructed and legacy mask."
)

# ── Main loop: site × fcover_threshold ────────────────────────────────────────

t0       <- proc.time()
n_total  <- length(sites) * length(fcover_thresholds)
all_rows <- vector("list", n_total)
row_idx  <- 0L

cli::cli_h2(
  "fCover sensitivity — {n_total} combinations ",
  "({length(sites)} sites × {length(fcover_thresholds)} thresholds)"
)

for (site in sites) {
  cli::cli_h3("Site: {site}")

  mdir <- masks_dir_path(site)

  # Build a reference mask (p = 0.90) to obtain the target grid for projection.
  # All three threshold masks share the same grid (same 4 component files).
  template_mask <- build_fcover_mask(0.90, mdir)

  # Load raw LAI rasters once per site; project to the mask grid.
  # terra::project() handles:
  #   - Extent differences (e.g. Mormal Raw ladstack: 1270 cols → 1260 cols)
  #   - S2 tile extent (e.g. Mormal Raw s2_atbd covers full tile 5935×7866)
  # Both rasters are already in UTM 31N — no CRS transformation involved.
  cli::cli_alert_info("  Loading and projecting Raw LAI rasters ...")

  ladstack_raw  <- terra::rast(ladstack_raw_path(site))
  ladstack_proj <- terra::project(ladstack_raw, template_mask)
  lai_als_raw   <- sum(ladstack_proj, na.rm = TRUE)

  s2_raw      <- terra::rast(s2_atbd_raw_path(site))
  lai_s2_raw  <- terra::project(s2_raw, template_mask)

  cli::cli_alert_success("  LAI rasters projected: {terra::ncell(template_mask)} cells")

  for (threshold in fcover_thresholds) {
    row_idx <- row_idx + 1L
    cli::cli_alert_info("  p = {threshold}")

    mask_p <- build_fcover_mask(threshold, mdir)

    all_rows[[row_idx]] <- compute_fcover_sensitivity_metrics(
      lai_als_rast     = lai_als_raw,
      lai_s2_atbd_rast = lai_s2_raw,
      mask_rast        = mask_p,
      site             = site,
      fcover_threshold = threshold
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
  "site", "fcover_threshold", "n_pixels",
  "lai_als_mean", "lai_als_median", "lai_als_p95", "lai_als_max",
  "lai_s2_atbd_mean", "lai_s2_atbd_median", "lai_s2_atbd_p95", "lai_s2_atbd_max",
  "R", "R2", "RMSE", "Bias", "Slope"
))

out_csv <- file.path(out_dir, "fcover_sensitivity_atbd.csv")
data.table::fwrite(final_dt, out_csv)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("fCover sensitivity — résumé")
cli::cli_bullets(c(
  "v" = "Durée : {elapsed} s",
  "v" = "CSV   : {out_csv}  ({nrow(final_dt)} lignes)"
))
cli::cli_text("Résultats (pivot site × fcover_threshold) :")
print(
  final_dt[, .(site, fcover_threshold, n_pixels,
               lai_als_mean = round(lai_als_mean, 3),
               R2           = round(R2, 3),
               RMSE         = round(RMSE, 3),
               Bias         = round(Bias, 3),
               Slope        = round(Slope, 3))]
)
