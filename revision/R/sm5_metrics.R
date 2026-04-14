# ---
# title:  sm5_metrics.R
# desc:   Core metric computation for the d_opt analysis (SM5).
#         Exposes compute_metrics_s2_lidar(), which computes five agreement
#         statistics between a Sentinel-2 LAI estimate vector and a LiDAR
#         PAD-based LAI reference vector.
#
#         Mirrors the inner compute_metrics() function defined locally in
#         PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_04C_analysis_depth.R
#         (lines 33–55), with one extension: R2 column added, computed as r^2
#         (Pearson squared).
# ---

# ── compute_metrics_s2_lidar ───────────────────────────────────────────────────

#' @title Compute agreement metrics between S2 LAI and LiDAR LAI
#'
#' @description
#' Computes five agreement statistics between a vector of Sentinel-2 LAI
#' estimates and a reference vector of LiDAR-derived PAD-based LAI values
#' (Bouvier et al. 2015). NA pairs are removed before any computation via
#' \code{complete.cases}, replicating the behaviour of the legacy inner
#' function (04C lines 35-37).
#'
#' The regression convention is \code{lm(s2 ~ lidar)}: S2 in Y, LiDAR in X.
#' This means the slope quantifies how much S2 LAI changes per unit of
#' LiDAR LAI, and reflects the "evaluation of S2 against a structural
#' reference" framing of the paper.
#'
#' R2 is computed as \code{r^2} (Pearson correlation squared).
#'
#' @param s2_vals     Numeric vector. Sentinel-2 LAI estimates (one column of
#'                    the inversion output CSV produced by SM4).
#' @param lidar_vals  Numeric vector. LiDAR PAD-based LAI values (column
#'                    \code{lidar_values} from the SM3 sampling CSV). Must be
#'                    the same length as \code{s2_vals}.
#'
#' @return A named list with five numeric scalars:
#' \describe{
#'   \item{R}{Pearson correlation (\code{cor(..., use = "complete.obs")}).
#'            Mirrors 04C line 43.}
#'   \item{R2}{Pearson correlation squared (\code{r^2}).
#'             Not in legacy 04C — added for reviewer R3.minor.7.}
#'   \item{RMSE}{Root mean squared error in LAI units.
#'               Mirrors 04C line 44: \code{sqrt(mean((s2 - lidar)^2))}.}
#'   \item{Bias}{Mean signed error \code{mean(s2 - lidar)}.
#'               Mirrors 04C line 45.}
#'   \item{Slope}{Regression slope from \code{lm(s2 ~ lidar)}.
#'                Mirrors 04C line 46: \code{coef(lm(s2_vals ~ lidar_vals))[2]}.}
#' }
#' Returns a list of five \code{NA_real_} when fewer than two complete pairs
#' are available.
#'
#' @references
#' Legacy source: Main_s2_04C_analysis_depth.R lines 33–55
#' (PROSAIL-Optimization/02_CODES/Sentinel2/, December 2025 version).
#'
#' @examples
#' set.seed(42)
#' lidar <- runif(100, 1, 8)
#' s2    <- lidar * 0.9 + rnorm(100, 0, 0.5)
#' compute_metrics_s2_lidar(s2, lidar)
#'
#' @export
compute_metrics_s2_lidar <- function(s2_vals, lidar_vals) {
  stopifnot(length(s2_vals) == length(lidar_vals))

  # Remove NA pairs (mirrors 04C lines 35-37)
  valid   <- stats::complete.cases(s2_vals, lidar_vals)
  s2_v    <- s2_vals[valid]
  lidar_v <- lidar_vals[valid]

  if (length(s2_v) < 2L) {
    return(list(
      R = NA_real_, R2 = NA_real_, RMSE = NA_real_,
      Bias = NA_real_, Slope = NA_real_
    ))
  }

  # Pearson r — mirrors 04C line 43
  r <- stats::cor(s2_v, lidar_v, use = "complete.obs")

  # R2 = r^2 (Pearson squared). Slope from lm(s2 ~ lidar) — mirrors 04C line 46.
  r2    <- r^2
  slope <- stats::coef(stats::lm(s2_v ~ lidar_v))[[2L]]

  # RMSE — mirrors 04C line 44: sqrt(mean((s2_vals - lidar_vals)^2, na.rm=TRUE))
  rmse <- sqrt(mean((s2_v - lidar_v)^2, na.rm = TRUE))

  # Bias = mean(S2 - LiDAR) — mirrors 04C line 45
  bias <- mean(s2_v - lidar_v, na.rm = TRUE)

  list(R = r, R2 = r2, RMSE = rmse, Bias = bias, Slope = slope)
}
