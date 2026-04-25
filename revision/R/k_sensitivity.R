# ---
# title:  k_sensitivity.R
# desc:   Analytical rescaling of LAI_ALS for extinction coefficient sensitivity.
#         Exposes rescale_lai_for_k(), compute_lai_als_at_k(), and
#         compute_k_sensitivity_metrics().
#
#         Scientific basis: Bouvier et al. (2015) formulation of LAD from LiDAR
#         gap fraction gives LAD(z) = -1/k × ln(P_gap(z)). Because P_gap does
#         not depend on k, the ratio between LAD computed at two different k
#         values is constant: LAD(z, k_new) = LAD(z, k_ref) × k_ref / k_new.
#         Integrating over depth gives:
#             LAI(k_new) = LAI(k_ref) × k_ref / k_new
#         This is a linear, exact rescaling with no approximation. Example:
#           k_ref = 0.5, k_new = 0.4  →  LAI(0.4) = LAI(0.5) × 1.25
#           k_ref = 0.5, k_new = 0.6  →  LAI(0.6) = LAI(0.5) × 0.833...
#
#         Addresses reviewer comment R3.MAJOR (sensitivity of LAI_ALS to k).
#         Analysis limited to LAI_S2_ATBD vs LAI_ALS (integrated, not d_opt).
#         No pipeline rerun; no PROSAIL LUT retraining required.
# ---

# ── rescale_lai_for_k ──────────────────────────────────────────────────────────

#' @title Rescale LAI from a reference extinction coefficient to a new value
#'
#' @description
#' Applies the analytical rescaling derived from the Bouvier et al. (2015)
#' LAD formulation:
#'
#' \deqn{\text{LAI}(k_\text{new}) = \text{LAI}(k_\text{ref})
#'       \times \frac{k_\text{ref}}{k_\text{new}}}
#'
#' The rescaling is exact (no approximation): because \eqn{P_\text{gap}(z)}
#' does not depend on \eqn{k}, the ratio between two LAD profiles computed with
#' different \eqn{k} values is a constant equal to \eqn{k_\text{ref} / k_\text{new}},
#' and this constant propagates linearly through the depth integral.
#'
#' The function operates on any numeric input (\code{numeric} vector or
#' \code{SpatRaster}): \code{terra} handles raster arithmetic natively, so
#' passing a \code{SpatRaster} produces a rescaled \code{SpatRaster}.
#'
#' @param lai_at_k_ref Numeric vector or \code{SpatRaster}. LAI values computed
#'   at \code{k_ref} (the reference extinction coefficient).
#' @param k_ref        Numeric. Reference extinction coefficient used when
#'   \code{lai_at_k_ref} was computed. Default \code{0.5} (value used in the
#'   submitted manuscript, Bouvier et al. 2015 convention for broadleaved
#'   forests).
#' @param k_new        Numeric. Target extinction coefficient.
#'
#' @return Same type and dimensions as \code{lai_at_k_ref}, with values
#'   multiplied by \code{k_ref / k_new}.
#'
#' @references
#' Bouvier, M., Durrieu, S., Fournier, R. A., & Renaud, J. P. (2015).
#' Generalizing relationships between canopy structure and light transmission:
#' the case of European temperate forests for a radiative transfer model
#' parametrization. Remote Sensing of Environment, 171, 158-171.
#' Eq. 1: LAD(z) = -1/k × ln(P_gap(z)).
#'
#' @examples
#' rescale_lai_for_k(10, k_ref = 0.5, k_new = 0.4)  # 12.5
#' rescale_lai_for_k(10, k_ref = 0.5, k_new = 0.6)  # 8.333...
#'
#' @export
rescale_lai_for_k <- function(lai_at_k_ref, k_ref = 0.5, k_new) {
  lai_at_k_ref * (k_ref / k_new)
}


# ── rescale_lai_for_theta ──────────────────────────────────────────────────────

#' @title Rescale LAI_ALS for a given scan angle via cosine correction
#'
#' @description
#' The Beer-Lambert formulation with scan angle \eqn{\theta} gives:
#'
#' \deqn{\text{LAI}(\theta) = \text{LAI}(\theta=0) \times \cos(\theta)}
#'
#' because the path length through each voxel increases as \eqn{1/\cos(\theta)},
#' so integrating LAD over slant path and projecting back to vertical gives a
#' factor of \eqn{\cos(\theta)}. This is an exact rescaling when \eqn{\theta}
#' is uniform across the scene (i.e. the PAD stack was computed at nadir).
#'
#' @param lai_at_nadir Numeric vector or \code{SpatRaster}. LAI computed at
#'   nadir (\eqn{\theta = 0}).
#' @param theta_deg    Numeric. Scan angle from nadir in degrees.
#'
#' @return Same type as \code{lai_at_nadir}, multiplied by \eqn{\cos(\theta)}.
#'
#' @examples
#' rescale_lai_for_theta(10, theta_deg = 15)  # 10 * cos(15° in rad) ≈ 9.66
#'
#' @export
rescale_lai_for_theta <- function(lai_at_nadir, theta_deg) {
  lai_at_nadir * cos(theta_deg * pi / 180)
}


# ── compute_lai_als_at_k ───────────────────────────────────────────────────────

#' @title Compute LAI_ALS at a given extinction coefficient via rescaling
#'
#' @description
#' Loads the multi-layer LAD profile stack (\code{ladstack_classic.tif}),
#' sums across layers to obtain integrated LAI at the reference extinction
#' coefficient (\code{k_ref = 0.5}), then applies the analytical rescaling
#' from \code{rescale_lai_for_k()} to obtain LAI at \code{k_new}.
#'
#' The ladstack file stores one LAD layer per 1 m height bin; summing across
#' all layers gives integrated LAI at the reference \eqn{k}. This convention
#' is identical to the loading logic in \code{sm6_load_rasters.R} (SM6b passe
#' 1). A \code{stop()} is raised if the ladstack has only one layer, indicating
#' an unexpected file (should be multi-layer).
#'
#' @param ladstack_path Character. Full path to \code{ladstack_classic.tif}.
#'   Typically
#'   \code{file.path(ext_dir, site, "Metrics", "Deciduous_Only",
#'   "ladstack_classic.tif")}.
#' @param k_new        Numeric. Target extinction coefficient.
#' @param k_ref        Numeric. Extinction coefficient embedded in the
#'   ladstack file. Default \code{0.5} (Bouvier et al. 2015).
#'
#' @return A single-band \code{SpatRaster} at the same grid as the ladstack,
#'   containing integrated LAI rescaled to \code{k_new}.
#'
#' @references
#' Convention for ladstack loading: \code{revision/R/sm6_load_rasters.R},
#' \code{load_site_rasters()} — multi-layer sum giving LAI_ALS at k_ref.
#'
#' Bouvier et al. (2015) — see \code{rescale_lai_for_k()}.
#'
#' @export
compute_lai_als_at_k <- function(ladstack_path, k_new, k_ref = 0.5) {
  ladstack <- terra::rast(ladstack_path)

  if (terra::nlyr(ladstack) == 1L) {
    stop(
      "compute_lai_als_at_k: '", basename(ladstack_path), "' has only 1 layer.",
      " Expected a multi-layer LAD profile stack (one layer per 1 m height bin).",
      " Check that the correct ladstack file is referenced."
    )
  }

  lai_ref <- sum(ladstack, na.rm = TRUE)
  rescale_lai_for_k(lai_ref, k_ref = k_ref, k_new = k_new)
}


# ── compute_k_sensitivity_metrics ──────────────────────────────────────────────

#' @title Compute descriptive statistics and comparison metrics for one (site, k)
#'
#' @description
#' Given a LAI_ALS raster (rescaled to a specific \eqn{k}) and a LAI_S2_ATBD
#' raster, extracts pixel-level values, applies pair-wise NA filtering, and
#' returns a single-row \code{data.table} with:
#'
#' \enumerate{
#'   \item Descriptive statistics for LAI_ALS: mean, median, 95th percentile,
#'     maximum.
#'   \item Same statistics for LAI_S2_ATBD (constant across \eqn{k} values,
#'     repeated for output readability).
#'   \item Comparison metrics: Pearson \eqn{r}, \eqn{r^2}, RMSE, Bias (mean
#'     signed difference s2 − lidar; positive = S2 overestimates LiDAR), and
#'     the slope of \code{lm(s2 ~ lidar)} (LAI_S2_ATBD as response, LAI_ALS as predictor).
#' }
#'
#' \strong{Regression convention}: \code{lm(s2 ~ lidar)}, harmonized with
#' the SM5 convention as of 2026-04-15.
#'
#' \strong{No p-values}: not requested for this reviewer analysis (R3.MAJOR
#' focuses on the magnitude of k sensitivity, not significance testing).
#'
#' @param lai_als_rast      \code{SpatRaster} (single-band). LAI_ALS rescaled
#'   to the target \eqn{k} via \code{compute_lai_als_at_k()}.
#' @param lai_s2_atbd_rast  \code{SpatRaster} (single-band). LAI_S2_ATBD.
#' @param site              Character. Site name (used only to populate the
#'   \code{site} column in the output).
#' @param k_value           Numeric. Extinction coefficient used for
#'   \code{lai_als_rast} (used to populate the \code{k_value} column).
#'
#' @return A \code{data.table} with one row and 15 columns:
#'   \code{site}, \code{k_value}, \code{n_pixels},
#'   \code{lai_als_mean}, \code{lai_als_median}, \code{lai_als_p95},
#'   \code{lai_als_max},
#'   \code{lai_s2_atbd_mean}, \code{lai_s2_atbd_median},
#'   \code{lai_s2_atbd_p95}, \code{lai_s2_atbd_max},
#'   \code{R}, \code{R2}, \code{RMSE}, \code{Bias}, \code{Slope}.
#'
#' @export
compute_k_sensitivity_metrics <- function(lai_als_rast, lai_s2_atbd_rast,
                                          site, k_value, theta_deg = 0) {
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
  # Bias = mean(s2 - lidar). Harmonized with lm(s2 ~ lidar) convention as of 2026-04-16.
  bias_val <- mean(s2_vec - lidar_vec)

  # lm(s2 ~ lidar): LAI_S2_ATBD as response, LAI_ALS as predictor.
  # Harmonized with SM5 convention as of 2026-04-15.
  lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
  slope_val <- stats::coef(lm_fit)[2L]

  data.table::data.table(
    site             = site,
    k_value          = k_value,
    theta_deg        = theta_deg,
    n_pixels         = n_pixels,
    lai_als_mean     = als_mean,
    lai_als_median   = als_median,
    lai_als_p95      = als_p95,
    lai_als_max      = als_max,
    lai_s2_atbd_mean   = s2_mean,
    lai_s2_atbd_median = s2_median,
    lai_s2_atbd_p95    = s2_p95,
    lai_s2_atbd_max    = s2_max,
    R                = r_val,
    R2               = r2_val,
    RMSE             = rmse_val,
    Bias             = bias_val,
    Slope            = slope_val
  )
}
