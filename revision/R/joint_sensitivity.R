# ---
# title:  joint_sensitivity.R
# desc:   Joint fCover × h_min sensitivity analysis for LAI_ALS vs LAI_S2_ATBD.
#         Exposes sum_lai_als_at_h_min_from_raster() and
#         compute_joint_sensitivity_metrics().
#
#         Depends on get_layer_heights() from h_min_sensitivity.R and
#         build_fcover_mask() from fcover_sensitivity.R — source both before
#         loading this module.
#
#         Scientific basis: the joint analysis crosses two threshold parameters
#         that each independently affect the LAI_ALS vs LAI_S2 comparison:
#
#           h_min  — minimum vegetation height for LAD integration (approx.
#                    via partial-band summation; see h_min_sensitivity.R header)
#           fCover — proportion of CHM pixels > 2 m within a 10 m pixel,
#                    used to filter out gap/edge pixels
#
#         Unlike the two 1D analyses (18_ and 19_), this module operates on
#         Raw/ladstack.tif (not Deciduous_Only/ladstack_classic.tif) so that
#         the fCover mask can be reconstructed at any threshold from the four
#         component rasters. The stack is projected once per site to the mask
#         grid; h_min selection and fCover masking are then applied on the
#         already-projected raster.
# ---


# ── sum_lai_als_at_h_min_from_raster ──────────────────────────────────────────

#' @title Sum LAD layers above h_min from an already-loaded SpatRaster
#'
#' @description
#' Selects layers whose bin-centre height (parsed from band names via
#' \code{get_layer_heights()}) satisfies \eqn{\geq h\_min + 0.5} m and returns
#' their pixel-wise sum. Unlike \code{compute_lai_als_at_h_min()}, this
#' function takes an already-loaded and projected \code{SpatRaster} rather than
#' a file path, avoiding repeated disk reads when iterating over multiple
#' \code{h_min} values for the same site.
#'
#' @param ladstack_rast A multi-layer \code{SpatRaster}. Band names must follow
#'   the \code{"LAD_Layer_X.5"} convention. The raster must already be
#'   projected to the target grid.
#' @param h_min Numeric. Minimum vegetation height threshold in metres.
#'   Layers with bin centre \eqn{< h\_min + 0.5} m are excluded.
#'
#' @return A single-band \code{SpatRaster} at the same grid as
#'   \code{ladstack_rast}, containing integrated LAI from the selected layers.
#'
#' @seealso \code{get_layer_heights()}, \code{compute_lai_als_at_h_min()}
#'
#' @export
sum_lai_als_at_h_min_from_raster <- function(ladstack_rast, h_min) {
  heights <- get_layer_heights(ladstack_rast)
  keep    <- heights >= (h_min + 0.5)

  if (!any(keep)) {
    stop(
      "sum_lai_als_at_h_min_from_raster: h_min = ", h_min,
      " excludes all ", terra::nlyr(ladstack_rast), " layers ",
      "(max bin centre = ", max(heights), " m). Reduce h_min.",
      call. = FALSE
    )
  }

  sum(terra::subset(ladstack_rast, which(keep)), na.rm = TRUE)
}


# ── compute_joint_sensitivity_metrics ─────────────────────────────────────────

#' @title Compute comparison metrics for one (site, fCover threshold, h_min)
#'
#' @description
#' Applies the fCover binary mask to both LAI rasters, extracts pair-wise
#' complete pixel values, and returns a single-row \code{data.table} with
#' descriptive statistics and comparison metrics, identically to
#' \code{compute_fcover_sensitivity_metrics()} but adding the \code{h_min_value}
#' column to identify the h_min level.
#'
#' Comparison metrics: Pearson \eqn{r}, \eqn{r^2}, RMSE, Bias
#' (\eqn{\overline{s2 - \text{lidar}}}; positive = S2 overestimates),
#' and the slope of \code{lm(s2 ~ lidar)}.
#'
#' @param lai_als_rast      \code{SpatRaster} (single-band). LAI_ALS summed via
#'   \code{sum_lai_als_at_h_min_from_raster()}, on the mask grid.
#' @param lai_s2_atbd_rast  \code{SpatRaster} (single-band). LAI_S2_ATBD
#'   projected to the mask grid.
#' @param mask_rast         \code{SpatRaster} (single-band). Binary fCover mask
#'   from \code{build_fcover_mask()}: 1 = valid, \code{NA} = excluded.
#' @param site              Character. Site name.
#' @param fcover_threshold  Numeric. fCover threshold used to build
#'   \code{mask_rast}.
#' @param h_min_value       Numeric. h_min threshold used to build
#'   \code{lai_als_rast}.
#'
#' @return A \code{data.table} with one row and 17 columns:
#'   \code{site}, \code{fcover_threshold}, \code{h_min_value}, \code{n_pixels},
#'   \code{lai_als_{mean,median,p95,max}},
#'   \code{lai_s2_atbd_{mean,median,p95,max}},
#'   \code{R}, \code{R2}, \code{RMSE}, \code{Bias}, \code{Slope}.
#'
#' @export
compute_joint_sensitivity_metrics <- function(lai_als_rast, lai_s2_atbd_rast,
                                              mask_rast, site,
                                              fcover_threshold, h_min_value) {
  lai_als_masked <- terra::mask(lai_als_rast,     mask_rast)
  lai_s2_masked  <- terra::mask(lai_s2_atbd_rast, mask_rast)

  lidar_raw <- as.numeric(terra::values(lai_als_masked))
  s2_raw    <- as.numeric(terra::values(lai_s2_masked))

  complete  <- !is.na(lidar_raw) & !is.na(s2_raw)
  lidar_vec <- lidar_raw[complete]
  s2_vec    <- s2_raw[complete]
  n_pixels  <- length(lidar_vec)

  als_mean   <- mean(lidar_vec)
  als_median <- stats::median(lidar_vec)
  als_p95    <- stats::quantile(lidar_vec, 0.95, names = FALSE)
  als_max    <- max(lidar_vec)

  s2_mean    <- mean(s2_vec)
  s2_median  <- stats::median(s2_vec)
  s2_p95     <- stats::quantile(s2_vec, 0.95, names = FALSE)
  s2_max     <- max(s2_vec)

  r_val     <- stats::cor(lidar_vec, s2_vec, use = "complete.obs")
  r2_val    <- r_val ^ 2
  rmse_val  <- sqrt(mean((lidar_vec - s2_vec) ^ 2))
  bias_val  <- mean(s2_vec - lidar_vec)
  lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
  slope_val <- stats::coef(lm_fit)[2L]

  data.table::data.table(
    site               = site,
    fcover_threshold   = fcover_threshold,
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
