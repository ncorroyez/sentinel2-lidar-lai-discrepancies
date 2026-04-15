# ---
# title:  h_min_sensitivity.R
# desc:   Analytical h_min sensitivity for LAI_ALS via partial-band summation.
#         Exposes get_layer_heights(), compute_lai_als_at_h_min(), and
#         compute_h_min_sensitivity_metrics().
#
#         Scientific basis: the LAD profile stack (ladstack_classic.tif)
#         stores one LAD layer per 1 m height bin, with band names of the
#         form "LAD_Layer_X.5" where X.5 is the bin centre in metres. The
#         vegetation height threshold h_min (default 2 m in the submitted
#         manuscript) excludes returns below h_min from the LAD integration.
#         Changing h_min while keeping all other pipeline parameters fixed
#         is approximated by summing only the bins whose centre ≥ h_min + 0.5.
#
#         Approximation caveat: this is an upper bound of true sensitivity.
#         The legacy pipeline applies h_min upstream during return
#         classification, which also affects the gap-fraction profile used
#         to compute LAD. The partial-band approach ignores that upstream
#         effect and therefore overestimates the impact of changing h_min.
#         The approximation is explicitly flagged in the response letter.
#
#         Addresses reviewer comment R3.minor.4 (h_min threshold sensitivity;
#         fCover sensitivity is handled in fcover_sensitivity.R).
# ---


# ── get_layer_heights ──────────────────────────────────────────────────────────

#' @title Parse band names of a LAD profile stack into numeric heights
#'
#' @description
#' Reads the band names of a multi-layer LAD \code{SpatRaster} produced by the
#' LiDAR pipeline and returns a numeric vector of bin-centre heights in metres.
#' Band names are expected to match the pattern \code{"LAD_Layer_X.5"} where
#' \code{X.5} is the bin centre (e.g. \code{"LAD_Layer_2.5"},
#' \code{"LAD_Layer_39.5"}). The function stops with an informative error if
#' any band name fails to match, preventing silent mis-assignment of heights.
#'
#' Band order in the stack (ascending vs. descending) is irrelevant: heights
#' are derived from the names, not from band position.
#'
#' @param ladstack A multi-layer \code{SpatRaster}. Band names must follow the
#'   \code{"LAD_Layer_X.5"} convention produced by
#'   \code{02_CODES/LiDAR/2.calculate_25m_metrics.R}.
#'
#' @return Numeric vector of length \code{terra::nlyr(ladstack)}, giving the
#'   bin-centre height in metres for each layer, in the same order as the
#'   bands of \code{ladstack}.
#'
#' @references
#' Band naming convention: \code{02_CODES/LiDAR/2.calculate_25m_metrics.R}
#' (LAD profile stack output); \code{02_CODES/libraries/functions_lidar.R}
#' (LAD computation using Bouvier et al. 2015).
#'
#' @examples
#' \dontrun{
#' ls <- terra::rast("03_RESULTS/Blois/Metrics/Deciduous_Only/ladstack_classic.tif")
#' heights <- get_layer_heights(ls)
#' range(heights)   # should be c(2.5, 39.5) for a 38-band stack
#' }
#'
#' @export
get_layer_heights <- function(ladstack) {
  nms <- terra::names(ladstack)

  heights <- suppressWarnings(
    as.numeric(sub("^LAD_Layer_", "", nms))
  )

  bad <- is.na(heights)
  if (any(bad)) {
    stop(
      "get_layer_heights: ", sum(bad), " band name(s) do not match ",
      "'LAD_Layer_X.5' pattern:\n  ",
      paste(nms[bad], collapse = ", "),
      call. = FALSE
    )
  }

  heights
}


# ── compute_lai_als_at_h_min ──────────────────────────────────────────────────

#' @title Compute LAI_ALS at a given minimum vegetation height via partial-band
#'   summation
#'
#' @description
#' Loads the multi-layer LAD profile stack (\code{ladstack_classic.tif}),
#' identifies the layers whose bin centre is \eqn{\geq h\_min + 0.5} m
#' (i.e., bins that start at or above \code{h_min}), and sums those layers
#' to obtain an approximation of integrated LAI at the given height threshold.
#'
#' \strong{Special case}: when \code{h_min = 2} all layers of the standard
#' 38-band stack (bin centres 2.5 m – 39.5 m) satisfy the condition, so the
#' result equals \code{sum(ladstack)} — the reference LAI at k = 0.5.
#'
#' \strong{Approximation}: this is an analytical upper bound. The legacy
#' pipeline applies \code{h_min} upstream during return classification, which
#' modifies the gap-fraction profile before LAD integration. Partial-band
#' summation ignores that upstream effect (see module header for details).
#'
#' @param ladstack_path Character. Full path to \code{ladstack_classic.tif}.
#'   Typically
#'   \code{file.path(ext_dir, site, "Metrics", "Deciduous_Only",
#'   "ladstack_classic.tif")}.
#' @param h_min Numeric. Minimum vegetation height threshold in metres.
#'   Values tested in the revision: \code{2} (reference), \code{3}, \code{5}.
#'   Must be a positive number less than the maximum bin height.
#'
#' @return A single-band \code{SpatRaster} at the same grid as the ladstack,
#'   containing integrated LAI from layers with bin centre \eqn{\geq h\_min
#'   + 0.5} m.
#'
#' @references
#' h_min convention: \code{02_CODES/libraries/functions_lidar.R}, hardcoded
#' at 2 m in the submitted manuscript. Bouvier et al. (2015) Eq. 1 for the
#' LAD formulation.
#'
#' @examples
#' \dontrun{
#' lp <- file.path("03_RESULTS", "Blois", "Metrics", "Deciduous_Only",
#'                 "ladstack_classic.tif")
#' lai_2 <- compute_lai_als_at_h_min(lp, h_min = 2)
#' lai_5 <- compute_lai_als_at_h_min(lp, h_min = 5)
#' # lai_5 <= lai_2 everywhere
#' }
#'
#' @export
compute_lai_als_at_h_min <- function(ladstack_path, h_min) {
  ladstack <- terra::rast(ladstack_path)

  if (terra::nlyr(ladstack) == 1L) {
    stop(
      "compute_lai_als_at_h_min: '", basename(ladstack_path),
      "' has only 1 layer. Expected a multi-layer LAD profile stack ",
      "(one layer per 1 m height bin). Check that the correct file is used.",
      call. = FALSE
    )
  }

  heights     <- get_layer_heights(ladstack)
  keep        <- heights >= (h_min + 0.5)
  n_kept      <- sum(keep)

  if (n_kept == 0L) {
    stop(
      "compute_lai_als_at_h_min: h_min = ", h_min,
      " excludes all ", terra::nlyr(ladstack), " layers ",
      "(max bin centre = ", max(heights), " m). Reduce h_min.",
      call. = FALSE
    )
  }

  sum(terra::subset(ladstack, which(keep)), na.rm = TRUE)
}


# ── compute_h_min_sensitivity_metrics ─────────────────────────────────────────

#' @title Compute descriptive statistics and comparison metrics for one
#'   (site, h_min)
#'
#' @description
#' Given a LAI_ALS raster (summed from layers at or above \code{h_min}) and a
#' LAI_S2_ATBD raster, extracts pixel-level values, applies pair-wise NA
#' filtering, and returns a single-row \code{data.table} with:
#'
#' \enumerate{
#'   \item Descriptive statistics for LAI_ALS: mean, median, 95th percentile,
#'     maximum.
#'   \item Same statistics for LAI_S2_ATBD (constant across \code{h_min}
#'     values for a given site, repeated for output readability).
#'   \item Comparison metrics: Pearson \eqn{r}, \eqn{r^2}, RMSE, Bias (mean
#'     signed difference lidar − S2), and the slope of
#'     \code{lm(s2 ~ lidar)}.
#' }
#'
#' \strong{Regression convention}: \code{lm(s2 ~ lidar)}, S2 as response and
#' LiDAR as predictor. Harmonized with SM5 convention as of 2026-04-15.
#'
#' \strong{No p-values}: not requested for this analysis (R3.minor.4 focuses
#' on the magnitude of h_min sensitivity, not significance testing).
#'
#' @param lai_als_rast     \code{SpatRaster} (single-band). Integrated LAI_ALS
#'   from \code{compute_lai_als_at_h_min()}, already on the Deciduous_Only
#'   grid. Must share the grid of \code{lai_s2_atbd_rast}.
#' @param lai_s2_atbd_rast \code{SpatRaster} (single-band). LAI_S2_ATBD from
#'   \code{Deciduous_Only/s2lai_summer_atbd_res_10_m.tif}. Must share the
#'   grid of \code{lai_als_rast}.
#' @param site             Character. Site name, used to populate the
#'   \code{site} column in the output.
#' @param h_min_value      Numeric. h_min threshold used to build
#'   \code{lai_als_rast}, used to populate the \code{h_min_value} column.
#'
#' @return A \code{data.table} with one row and 16 columns:
#'   \code{site}, \code{h_min_value}, \code{n_pixels},
#'   \code{lai_als_mean}, \code{lai_als_median}, \code{lai_als_p95},
#'   \code{lai_als_max},
#'   \code{lai_s2_atbd_mean}, \code{lai_s2_atbd_median},
#'   \code{lai_s2_atbd_p95}, \code{lai_s2_atbd_max},
#'   \code{R}, \code{R2}, \code{RMSE}, \code{Bias}, \code{Slope}.
#'
#' @export
compute_h_min_sensitivity_metrics <- function(lai_als_rast, lai_s2_atbd_rast,
                                              site, h_min_value) {
  lidar_raw <- as.numeric(terra::values(lai_als_rast))
  s2_raw    <- as.numeric(terra::values(lai_s2_atbd_rast))

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
  bias_val <- mean(lidar_vec - s2_vec)

  # lm(s2 ~ lidar): LAI_S2_ATBD as response, LAI_ALS as predictor.
  # Harmonized with SM5 convention as of 2026-04-15.
  lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
  slope_val <- stats::coef(lm_fit)[2L]

  data.table::data.table(
    site               = site,
    h_min_value        = h_min_value,
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
