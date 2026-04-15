# ---
# title:  sm6_heterogeneity.R
# desc:   Computation of canopy horizontal heterogeneity metrics for SM6.
#         Exposes compute_heterogeneity() and save_heterogeneity_raster().
#
#         Faithful refactor of 02_CODES/LiDAR/2.calculate_25m_metrics.R
#         (active block lines 901-916) with one extension: the DSM-based
#         metric is now parameterised so the same pipeline applies to either
#         the DSM or the CHM 1 m raster, addressing reviewer comment R2.8
#         (CHM_std vs DSM_std comparison).
#
#         Method: terra::aggregate(fact = 10, fun = "sd", na.rm = TRUE) on
#         a 1 m raster, followed by terra::project onto the 10 m mask grid
#         and terra::mask. Block non-sliding aggregate, no focal window.
# ---

# ── compute_heterogeneity ──────────────────────────────────────────────────────

#' @title Compute block standard deviation of canopy height at 10 m resolution
#'
#' @description
#' Aggregates a 1 m canopy height raster (DSM or CHM) to 10 m by computing
#' the standard deviation over non-overlapping 10 × 10 pixel blocks, then
#' projects the result onto the reference 10 m mask grid and applies the mask.
#'
#' The three-step pipeline — \code{aggregate} → \code{project} → \code{mask}
#' — is a literal refactor of the active DSM-SD computation in
#' \code{02_CODES/LiDAR/2.calculate_25m_metrics.R} lines 901-916.
#' No focal (sliding) window is used: each output pixel represents exactly
#' one non-overlapping 10 × 10 block of the 1 m source raster.
#'
#' The \code{project} step uses \code{method = "bilinear"}, which is terra's
#' default but is made explicit here for reproducibility. The legacy code at
#' line 906 omits \code{method} (same effective behaviour).
#'
#' @param source_raster_path Character. Path to the 1 m source raster
#'   (DSM: \code{rasterize_canopy.vrt}, or CHM: \code{chm.tif}).
#' @param mask_raster_path   Character. Path to the 10 m Deciduous_Only mask
#'   raster
#'   (\code{artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi}).
#'   Combines the BDForêt deciduous filter, low-vegetation exclusion (fCover
#'   majority 90 %), and artifact removal. Used both as a spatial template
#'   for \code{project} and as a binary mask for \code{mask}. Grids of this
#'   raster and \code{mnc_mean_heights_res_10_m.envi} are identical per-site
#'   (verified for Blois: same ext, CRS EPSG:32631, res 10 m).
#' @param fact               Integer. Aggregation factor (number of source
#'   pixels per output pixel along each dimension). Default \code{10L}
#'   reproduces the legacy \code{fact = 10} on 1 m data → 10 m output.
#'
#' @return A \code{SpatRaster} at 10 m resolution with the same extent and
#'   CRS as \code{mask_raster_path}, masked to the footprint of the mask.
#'
#' @references
#' Legacy: 02_CODES/LiDAR/2.calculate_25m_metrics.R lines 901-916
#' (active block aggregate sd on DSM 1 m, project, mask).
#'
#' @examples
#' \dontrun{
#' het <- compute_heterogeneity(
#'   source_raster_path = here::here(
#'     "03_RESULTS", "Blois", "LiDAR", "dsm", "res_1_m",
#'     "rasterize_canopy.vrt"
#'   ),
#'   mask_raster_path = here::here(
#'     "03_RESULTS", "Blois", "LiDAR", "Heterogeneity_Masks",
#'     "artifacts_deciduous_only_low_vegetation_majority_90_p_res_10_m.envi"
#'   )
#' )
#' }
#'
#' @export
compute_heterogeneity <- function(source_raster_path, mask_raster_path,
                                  fact = 10L) {
  src <- terra::rast(source_raster_path)
  msk <- terra::rast(mask_raster_path)

  # Block aggregate: sd over fact × fact pixels — mirrors legacy line 901
  agg <- terra::aggregate(src, fact = fact, fun = "sd", na.rm = TRUE)

  # Project onto mask grid — mirrors legacy line 906; method made explicit
  agg <- terra::project(agg, msk, method = "bilinear")

  # Apply mask — mirrors legacy line 907
  agg <- terra::mask(agg, msk)

  agg
}

# ── save_heterogeneity_raster ──────────────────────────────────────────────────

#' @title Write a heterogeneity SpatRaster to disk
#'
#' @description
#' Writes the output of \code{compute_heterogeneity()} to a GeoTIFF file
#' under \code{out_dir/{site}/{metric}_sd_res_10_m.tif}. The output filename
#' convention mirrors the legacy \code{dsm_sd_res_10_m.tif} for
#' \code{metric = "DSM"} and extends it consistently for \code{metric = "CHM"}.
#'
#' @param het_raster A \code{SpatRaster}. The heterogeneity raster produced
#'   by \code{compute_heterogeneity()}.
#' @param site       Character. Study site name (e.g. \code{"Blois"}).
#'   Used as a subdirectory under \code{out_dir}.
#' @param metric     Character. One of \code{"DSM"} or \code{"CHM"}.
#'   Determines the output filename prefix (lowercased).
#' @param out_dir    Character. Root output directory. The subdirectory
#'   \code{file.path(out_dir, site)} is created if absent.
#'
#' @return Character. The full path of the written file (invisibly).
#'
#' @examples
#' \dontrun{
#' save_heterogeneity_raster(
#'   het_raster = het,
#'   site       = "Blois",
#'   metric     = "DSM",
#'   out_dir    = here::here("revision", "output", "intermediate", "sm6")
#' )
#' }
#'
#' @export
save_heterogeneity_raster <- function(het_raster, site, metric, out_dir) {
  dir_path <- file.path(out_dir, site)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  out_path <- file.path(
    dir_path,
    paste0(tolower(metric), "_sd_res_10_m.tif")
  )

  terra::writeRaster(het_raster, out_path, overwrite = TRUE)

  invisible(out_path)
}
