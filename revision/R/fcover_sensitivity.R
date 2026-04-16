# ---
# title:  fcover_sensitivity.R
# desc:   Analytical fCover mask reconstruction and comparison metrics for
#         fCover threshold sensitivity analysis.
#         Exposes build_fcover_mask() and compute_fcover_sensitivity_metrics().
#
#         Scientific basis: the deciduous-forest analysis mask applied in the
#         submitted manuscript (p = 0.90) combines four binary / continuous
#         components via:
#
#             mask_p = cloud × edges × deciduous × I(fcover > p)
#             mask_p[mask_p == 0] ← NA
#
#         Changing p while holding the other three components fixed tests the
#         sensitivity of the LAI_ALS vs LAI_S2_ATBD comparison to the fCover
#         threshold. This module re-implements the mask construction directly
#         (without sourcing the legacy pipeline) and computes the same
#         descriptive statistics and comparison metrics as k_sensitivity.R.
#
#         Addresses reviewer comment R3.minor.4 (fCover threshold sensitivity;
#         h_min sensitivity is handled qualitatively in the response letter).
# ---


# ── build_fcover_mask ───────────────────────────────────────────────────────────

#' @title Build fCover sensitivity mask for one threshold
#'
#' @description
#' Constructs the binary deciduous-forest analysis mask at a given fCover
#' threshold \eqn{p} by multiplying four component rasters from
#' \code{Heterogeneity_Masks/}:
#'
#' \deqn{\text{mask}_p = \text{cloud} \times \text{edges} \times \text{decidu}
#'       \times \mathbf{1}[\text{fcover} > p]}
#'
#' Pixels where the product equals 0 are set to \code{NA}. This formula exactly
#' reproduces the legacy mask at \eqn{p = 0.90} (file:
#' \code{artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi})
#' and generalises to arbitrary \eqn{p} without relaunching the LiDAR pipeline.
#'
#' The four components (all at 10 m resolution, in \code{masks_dir}):
#' \itemize{
#'   \item \code{cloud_mask_res_10_m.envi} — cloud/shadow mask derived from S2
#'     reflectance (binary, 1 = valid).
#'   \item \code{site_edges_mask_res_10_m.envi} — LiDAR footprint edge mask
#'     (binary, 1 = within footprint).
#'   \item \code{deciduous_only_mask_res_10_m.envi} — BDForêt deciduous-forest
#'     composition mask, resampled to 10 m (binary, 1 = deciduous).
#'   \item \code{chm_1_m_thresholded_average_to_res_10_m.envi} — continuous
#'     fCover: proportion of 1 m CHM pixels exceeding 2 m within each 10 m
#'     pixel (values in \eqn{[0, 1]}; written as \code{average_chm_thresh} in
#'     the legacy pipeline).
#' }
#'
#' @param threshold Numeric. fCover threshold in \eqn{(0, 1)}. Pixels with
#'   \code{fcover_continuous <= threshold} are excluded from the analysis.
#'   Values tested in the manuscript revision: \code{0.80}, \code{0.90}
#'   (reference), \code{0.95}.
#' @param masks_dir Character. Full path to the \code{Heterogeneity_Masks/}
#'   directory for the site, containing the four component ENVI rasters.
#'
#' @return A single-band \code{SpatRaster} at 10 m resolution. Pixels
#'   satisfying all four conditions have value 1; all others are \code{NA}.
#'
#' @references
#' Legacy mask construction: \code{02_CODES/libraries/functions_create_masks.R},
#' lines 651-657 (\code{mask_majority_project_chm_thresh()}) and lines 753-777
#' (final mask formula and zero-to-NA conversion). Reproduced here without
#' sourcing the legacy module.
#'
#' @examples
#' \dontrun{
#' masks_dir <- file.path("03_RESULTS", "Blois", "LiDAR", "Heterogeneity_Masks")
#' mask_90 <- build_fcover_mask(0.90, masks_dir)
#' mask_80 <- build_fcover_mask(0.80, masks_dir)
#' # mask_80 should include more pixels than mask_90
#' sum(!is.na(terra::values(mask_80))) >= sum(!is.na(terra::values(mask_90)))
#' }
#'
#' @export
build_fcover_mask <- function(threshold, masks_dir) {
  cloud  <- terra::rast(file.path(masks_dir, "cloud_mask_res_10_m.envi"))
  edges  <- terra::rast(file.path(masks_dir, "site_edges_mask_res_10_m.envi"))
  decidu <- terra::rast(file.path(masks_dir, "deciduous_only_mask_res_10_m.envi"))
  fcover <- terra::rast(
    file.path(masks_dir, "chm_1_m_thresholded_average_to_res_10_m.envi")
  )

  mask_p             <- cloud * edges * decidu * (fcover > threshold)
  mask_p[mask_p == 0] <- NA
  mask_p
}


# ── compute_fcover_sensitivity_metrics ─────────────────────────────────────────

#' @title Compute descriptive statistics and comparison metrics for one
#'   (site, fCover threshold)
#'
#' @description
#' Given a LAI_ALS raster, a LAI_S2_ATBD raster, and an fCover binary mask,
#' applies the mask to both rasters, extracts pixel-level values, applies
#' pair-wise NA filtering, and returns a single-row \code{data.table} with:
#'
#' \enumerate{
#'   \item Descriptive statistics for LAI_ALS: mean, median, 95th percentile,
#'     maximum.
#'   \item Same statistics for LAI_S2_ATBD (constant across thresholds for a
#'     given site, repeated for output readability).
#'   \item Comparison metrics: Pearson \eqn{r}, \eqn{r^2}, RMSE, Bias (mean
#'     signed difference s2 − lidar; positive = S2 overestimates LiDAR), and
#'     the slope of \code{lm(s2 ~ lidar)}.
#' }
#'
#' \strong{Regression convention}: \code{lm(s2 ~ lidar)}, harmonized with
#' the SM5 convention as of 2026-04-15.
#'
#' \strong{No p-values}: not requested for this analysis (R3.minor.4 focuses
#' on the magnitude of fCover sensitivity, not significance testing).
#'
#' @param lai_als_rast      \code{SpatRaster} (single-band). Integrated LAI_ALS
#'   from \code{Raw/ladstack.tif} projected to the mask grid and summed over
#'   height bins. Must be on the same grid as \code{mask_rast}.
#' @param lai_s2_atbd_rast  \code{SpatRaster} (single-band). LAI_S2_ATBD from
#'   \code{Raw/s2lai_summer_atbd_res_10_m.tif} projected to the mask grid.
#'   Must be on the same grid as \code{mask_rast}.
#' @param mask_rast         \code{SpatRaster} (single-band). Binary mask
#'   produced by \code{build_fcover_mask()}: value 1 for valid pixels, \code{NA}
#'   elsewhere. Must share the grid of \code{lai_als_rast} and
#'   \code{lai_s2_atbd_rast}.
#' @param site              Character. Site name, used to populate the
#'   \code{site} column in the output.
#' @param fcover_threshold  Numeric. fCover threshold used to build
#'   \code{mask_rast}, used to populate the \code{fcover_threshold} column.
#'
#' @return A \code{data.table} with one row and 16 columns:
#'   \code{site}, \code{fcover_threshold}, \code{n_pixels},
#'   \code{lai_als_mean}, \code{lai_als_median}, \code{lai_als_p95},
#'   \code{lai_als_max},
#'   \code{lai_s2_atbd_mean}, \code{lai_s2_atbd_median},
#'   \code{lai_s2_atbd_p95}, \code{lai_s2_atbd_max},
#'   \code{R}, \code{R2}, \code{RMSE}, \code{Bias}, \code{Slope}.
#'
#' @export
compute_fcover_sensitivity_metrics <- function(lai_als_rast, lai_s2_atbd_rast,
                                               mask_rast, site,
                                               fcover_threshold) {
  # Apply mask: valid pixels stay, masked pixels → NA
  lai_als_masked <- terra::mask(lai_als_rast,     mask_rast)
  lai_s2_masked  <- terra::mask(lai_s2_atbd_rast, mask_rast)

  lidar_raw <- as.numeric(terra::values(lai_als_masked))
  s2_raw    <- as.numeric(terra::values(lai_s2_masked))

  # Pair-wise complete cases
  complete  <- !is.na(lidar_raw) & !is.na(s2_raw)
  lidar_vec <- lidar_raw[complete]
  s2_vec    <- s2_raw[complete]
  n_pixels  <- length(lidar_vec)

  # ── LAI_ALS descriptive statistics ─────────────────────────────────────────
  als_mean   <- mean(lidar_vec)
  als_median <- stats::median(lidar_vec)
  als_p95    <- stats::quantile(lidar_vec, 0.95, names = FALSE)
  als_max    <- max(lidar_vec)

  # ── LAI_S2_ATBD descriptive statistics ──────────────────────────────────────
  s2_mean    <- mean(s2_vec)
  s2_median  <- stats::median(s2_vec)
  s2_p95     <- stats::quantile(s2_vec, 0.95, names = FALSE)
  s2_max     <- max(s2_vec)

  # ── Comparison metrics ───────────────────────────────────────────────────────
  r_val    <- stats::cor(lidar_vec, s2_vec, use = "complete.obs")
  r2_val   <- r_val ^ 2
  rmse_val <- sqrt(mean((lidar_vec - s2_vec) ^ 2))
  # Bias = mean(s2 - lidar). Harmonized with lm(s2 ~ lidar) convention as of 2026-04-16.
  bias_val <- mean(s2_vec - lidar_vec)

  # lm(s2 ~ lidar): LAI_S2_ATBD as response, LAI_ALS as predictor.
  # Harmonized with SM5 convention as of 2026-04-15.
  lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
  slope_val <- stats::coef(lm_fit)[2L]

  data.table::data.table(
    site               = site,
    fcover_threshold   = fcover_threshold,
    n_pixels           = n_pixels,
    lai_als_mean       = als_mean,
    lai_als_median     = als_median,
    lai_als_p95        = als_p95,
    lai_als_max        = als_max,
    lai_s2_atbd_mean   = s2_mean,
    lai_s2_atbd_median = s2_median,
    lai_s2_atbd_p95    = s2_p95,
    lai_s2_atbd_max    = s2_max,
    R                  = r_val,
    R2                 = r2_val,
    RMSE               = rmse_val,
    Bias               = bias_val,
    Slope              = slope_val
  )
}
