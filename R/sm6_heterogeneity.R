# ---
# title:  sm6_heterogeneity.R
# desc:   Computation of canopy horizontal heterogeneity metrics for SM6.
#         Exposes build_output_filename(), compute_heterogeneity(), and
#         save_heterogeneity_raster().
#
#         Passe 1 — faithful refactor of 02_CODES/LiDAR/2.calculate_25m_metrics.R
#         (active block lines 901-916) with CHM extension for reviewer R2.8.
#
#         Passe 2 — multi-scale sweep:
#           method = "block": terra::aggregate(fact, fun="sd") at fact ∈ {5,10,20,50}
#             → output rasters at 5, 10, 20, 50 m window.
#           method = "focal": terra::focal(w=K, fun="sd") on the block_10_m base,
#             K ∈ {3, 5} → effective sliding windows of 30, 50 m.
#             Note: K=1 is excluded because terra::focal rejects w=1 as
#             non-meaningful. The block_10_m raster (aggregate fact=10, fun='sd')
#             is its semantic equivalent — a 1×1 focal SD on the 10 m grid is
#             identically the per-pixel value, which block_10_m already represents.
#             focal_Ks = c(3L, 5L) in the orchestration.
#
#         Output filename convention (passe 2):
#           block: {metric}_sd_block_{size_m}_m.tif    (e.g. dsm_sd_block_10_m.tif)
#           focal: {metric}_sd_focal_{size_m}_m.tif    (e.g. dsm_sd_focal_30_m.tif)
#           legacy: {metric}_sd_res_10_m.tif           (= block_10_m, rétrocompat SM6b)
# ---

# ── build_output_filename ──────────────────────────────────────────────────────

#' @title Build output filename for a heterogeneity raster
#'
#' @description
#' Constructs the standard filename for a heterogeneity raster produced by
#' \code{compute_heterogeneity()}, following the SM6 passe-2 convention:
#' \code{{metric}_sd_{method}_{size_m}_m.tif}.
#'
#' @param metric  Character. One of \code{"DSM"} or \code{"CHM"} (case-
#'   insensitive — lowercased internally).
#' @param method  Character. \code{"block"} or \code{"focal"}.
#' @param size_m  Integer or numeric. Effective window size in metres.
#'   For block: equals the aggregation factor (1 m source → \code{fact} m
#'   output). For focal: \code{K * 10} (window of K pixels on the 10 m grid).
#'
#' @return Character. Filename, e.g. \code{"dsm_sd_block_10_m.tif"} or
#'   \code{"chm_sd_focal_30_m.tif"}.
#'
#' @examples
#' build_output_filename("DSM", "block", 10)   # "dsm_sd_block_10_m.tif"
#' build_output_filename("CHM", "focal", 30)   # "chm_sd_focal_30_m.tif"
#'
build_output_filename <- function(metric, method, size_m) {
  paste0(tolower(metric), "_sd_", method, "_", size_m, "_m.tif")
}


# ── compute_heterogeneity ──────────────────────────────────────────────────────

#' @title Compute canopy height heterogeneity at a given spatial scale
#'
#' @description
#' Derives a canopy height standard-deviation raster at 10 m resolution from
#' a 1 m source raster (DSM or CHM) using one of two aggregation strategies:
#'
#' \strong{Block} (\code{method = "block"}): Non-overlapping aggregation via
#' \code{terra::aggregate(fact, fun = "sd", na.rm = TRUE)} followed by
#' \code{terra::project} onto the 10 m mask grid and \code{terra::mask}.
#' Each output pixel represents exactly one \code{fact × fact} block of the
#' 1 m source. Setting \code{fact = 10L} reproduces the legacy pipeline from
#' \code{02_CODES/LiDAR/2.calculate_25m_metrics.R} lines 901-916.
#'
#' \strong{Focal} (\code{method = "focal"}): Sliding-window SD computed by
#' \code{terra::focal(w = K, fun = "sd", na.rm = TRUE)} on the block-10 m
#' base raster (either supplied via \code{base_raster} or recomputed
#' internally). Each output pixel integrates a \code{K × K} neighbourhood
#' of 10 m pixels → effective window of \code{K * 10} metres.
#' \emph{K must be an odd integer \eqn{\geq 3}}. \code{K = 1} is rejected by
#' \code{terra::focal} as a non-meaningful window; in any case, it would be
#' semantically equivalent to \code{block_10_m} (a 1-pixel focal SD on the
#' 10 m grid is the per-pixel value, which \code{block_10_m} already provides).
#'
#' The \code{project} step uses \code{method = "bilinear"} (terra default,
#' made explicit for reproducibility).
#'
#' @param source_raster_path Character. Path to the 1 m source raster
#'   (DSM: \code{rasterize_canopy.vrt}, CHM: \code{chm.tif}).
#' @param mask_raster_path   Character. Path to the 10 m Deciduous_Only mask
#'   (\code{artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi}).
#'   Used as spatial template for \code{project} and as a binary mask.
#' @param fact               Integer. Aggregation factor for \code{method =
#'   "block"}. Default \code{10L} → 10 m output from 1 m source, matching the
#'   legacy. Ignored when \code{method = "focal"}.
#' @param method             Character. Aggregation strategy: \code{"block"}
#'   (default, non-overlapping aggregate) or \code{"focal"} (sliding window).
#' @param K                  Integer or \code{NULL}. Window size in pixels for
#'   \code{method = "focal"}. Must be an odd integer \eqn{\geq 3}. Effective
#'   spatial scale: \code{K * 10} metres. \code{NULL} for \code{method =
#'   "block"}.
#' @param base_raster        \code{SpatRaster} or \code{NULL}. Pre-computed
#'   block-10 m raster to use as input for \code{method = "focal"}, avoiding
#'   redundant computation when the orchestration script already holds it.
#'   If \code{NULL}, the block-10 m base is recomputed internally from
#'   \code{source_raster_path}.
#'
#' @return A masked \code{SpatRaster} at 10 m resolution, aligned to
#'   \code{mask_raster_path}.
#'
#' @references
#' Block legacy: 02_CODES/LiDAR/2.calculate_25m_metrics.R lines 901-916.
#'
#' @examples
#' \dontrun{
#' # Block at 10 m (legacy-equivalent)
#' het10 <- compute_heterogeneity(
#'   source_raster_path = here::here("03_RESULTS","Blois","LiDAR","dsm",
#'                                   "res_1_m","rasterize_canopy.vrt"),
#'   mask_raster_path   = here::here("03_RESULTS","Blois","LiDAR",
#'                                   "Heterogeneity_Masks",
#'    "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi")
#' )
#'
#' # Focal 30 m sliding window using cached block_10m
#' het_f30 <- compute_heterogeneity(
#'   source_raster_path = here::here("03_RESULTS","Blois","LiDAR","dsm",
#'                                   "res_1_m","rasterize_canopy.vrt"),
#'   mask_raster_path   = here::here("03_RESULTS","Blois","LiDAR",
#'                                   "Heterogeneity_Masks",
#'    "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"),
#'   method      = "focal",
#'   K           = 3L,
#'   base_raster = het10
#' )
#' }
#'
#' @export
compute_heterogeneity <- function(source_raster_path, mask_raster_path,
                                  fact        = 10L,
                                  method      = "block",
                                  K           = NULL,
                                  base_raster = NULL) {
  msk <- terra::rast(mask_raster_path)

  if (method == "block") {
    src <- terra::rast(source_raster_path)
    # Block aggregate: sd over fact × fact pixels — mirrors legacy line 901
    agg <- terra::aggregate(src, fact = fact, fun = "sd", na.rm = TRUE)
    # Project onto mask grid — mirrors legacy line 906; method made explicit
    agg <- terra::project(agg, msk, method = "bilinear")
    # Apply mask — mirrors legacy line 907
    agg <- terra::mask(agg, msk)
    return(agg)
  }

  if (method == "focal") {
    if (is.null(K)) {
      stop("compute_heterogeneity: 'K' must be specified for method = 'focal'.")
    }
    if (K < 3L || K %% 2L == 0L) {
      stop(
        "compute_heterogeneity: K = ", K, " is invalid. ",
        "K must be an odd integer >= 3. ",
        "K = 1 is rejected by terra::focal as a non-meaningful window; ",
        "use block_10_m as its proxy (Option C)."
      )
    }

    # Base: block_10_m — use precomputed raster if provided to avoid rework
    if (!is.null(base_raster)) {
      base <- base_raster
    } else {
      src  <- terra::rast(source_raster_path)
      base <- terra::aggregate(src, fact = 10L, fun = "sd", na.rm = TRUE)
      base <- terra::project(base, msk, method = "bilinear")
      base <- terra::mask(base, msk)
    }

    # Sliding-window SD: K × K pixel window on the 10 m block result
    # Effective spatial scale: K * 10 metres
    agg <- terra::focal(base, w = K, fun = "sd", na.rm = TRUE)
    agg <- terra::mask(agg, msk)
    return(agg)
  }

  stop(
    "compute_heterogeneity: unknown method '", method,
    "'. Use 'block' or 'focal'."
  )
}


# ── save_heterogeneity_raster ──────────────────────────────────────────────────

#' @title Write a heterogeneity SpatRaster to disk
#'
#' @description
#' Writes the output of \code{compute_heterogeneity()} to a GeoTIFF file
#' at \code{file.path(out_dir, site, out_filename)}. The output directory
#' \code{file.path(out_dir, site)} is created if absent.
#'
#' The \code{out_filename} argument decouples naming from this function,
#' allowing both the new SM6 passe-2 convention
#' (\code{build_output_filename()}) and the legacy name
#' (\code{"{metric}_sd_res_10_m.tif"}) to be passed without branching.
#'
#' @param het_raster   A \code{SpatRaster}. The heterogeneity raster produced
#'   by \code{compute_heterogeneity()}.
#' @param site         Character. Study site name (e.g. \code{"Blois"}).
#'   Used as a subdirectory under \code{out_dir}.
#' @param out_filename Character. Full filename including extension, e.g.
#'   \code{"dsm_sd_block_10_m.tif"} or \code{"dsm_sd_res_10_m.tif"}.
#' @param out_dir      Character. Root output directory. The subdirectory
#'   \code{file.path(out_dir, site)} is created if absent.
#'
#' @return Character. The full path of the written file (invisibly).
#'
#' @examples
#' \dontrun{
#' save_heterogeneity_raster(
#'   het_raster   = het,
#'   site         = "Blois",
#'   out_filename = build_output_filename("DSM", "block", 10),
#'   out_dir      = here::here("output", "intermediate", "sm6")
#' )
#' }
#'
#' @export
save_heterogeneity_raster <- function(het_raster, site, out_filename, out_dir) {
  dir_path <- file.path(out_dir, site)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  out_path <- file.path(dir_path, out_filename)
  terra::writeRaster(het_raster, out_path, overwrite = TRUE)

  invisible(out_path)
}


# ── window_metrics ────────────────────────────────────────────────────────────

#' @title Sliding-window r/RMSE/Bias/Slope along a covariate axis
#'
#' @description
#' Sorts paired \code{(x_var, y_var)} observations by an ordering covariate
#' (typically \code{dsm_sd}) and computes the four agreement metrics inside
#' overlapping windows of width \code{window_n}, advancing by
#' \code{window_step} in rank. Metrics use the \code{lm(y ~ x)} convention
#' consistent with \code{compute_metrics_by_class()}:
#' \describe{
#'   \item{R}{Pearson \code{cor(x_var, y_var)}}
#'   \item{RMSE}{\code{sqrt(mean((y_var - x_var)^2))}}
#'   \item{Bias}{\code{mean(y_var - x_var)}}
#'   \item{Slope}{\code{coef(lm(y_var ~ x_var))[[2]]}}
#' }
#'
#' Windows containing fewer than \code{min_n_per_window} valid pairs are
#' dropped. The window centre is the median of the ordering covariate inside
#' the window (so the abscissa is on the covariate's own scale).
#'
#' WARNING: adjacent windows overlap by \code{(window_n - window_step) /
#' window_n}, so the resulting series is not independent across windows;
#' do NOT compute naive confidence intervals or significance tests on it.
#'
#' @param x_var            Numeric vector. Predictor (e.g. LAI_ALS).
#' @param y_var            Numeric vector. Response (e.g. LAI_S2).
#' @param dsm_sd           Numeric vector. Ordering covariate.
#' @param window_n         Integer. Window width in points.
#' @param window_step      Integer. Step between window starts in rank.
#' @param min_n_per_window Integer. Skip windows with fewer valid pairs.
#'
#' @return A \code{data.table} with columns
#'   \code{dsm_sd_center, n_window, R, RMSE, Bias, Slope}. Empty
#'   \code{data.table} if fewer than \code{min_n_per_window} valid pairs.
#'
#' @export
window_metrics <- function(x_var, y_var, dsm_sd, window_n, window_step,
                           min_n_per_window) {
  ok  <- !is.na(x_var) & !is.na(y_var) & !is.na(dsm_sd)
  x_var  <- x_var[ok]; y_var <- y_var[ok]; dsm_sd <- dsm_sd[ok]
  o      <- order(dsm_sd)
  x_var  <- x_var[o];  y_var <- y_var[o]; dsm_sd <- dsm_sd[o]
  n_tot  <- length(dsm_sd)

  if (n_tot < min_n_per_window) return(data.table::data.table())

  win_starts <- seq.int(1L, max(1L, n_tot - window_n + 1L), by = window_step)
  rows <- vector("list", length(win_starts))
  for (i in seq_along(win_starts)) {
    s <- win_starts[i]; e <- min(s + window_n - 1L, n_tot)
    n_win <- e - s + 1L
    if (n_win < min_n_per_window) next
    xv <- x_var[s:e]; yv <- y_var[s:e]
    dsm_c <- stats::median(dsm_sd[s:e])

    diff <- yv - xv
    r_val    <- stats::cor(xv, yv)
    rmse_val <- sqrt(mean(diff^2))
    bias_val <- mean(diff)
    fit       <- stats::lm(yv ~ xv)
    slope_val <- stats::coef(fit)[[2L]]

    rows[[i]] <- data.table::data.table(
      dsm_sd_center = dsm_c,
      n_window      = n_win,
      R             = r_val,
      RMSE          = rmse_val,
      Bias          = bias_val,
      Slope         = slope_val
    )
  }
  data.table::rbindlist(rows)
}
