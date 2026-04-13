# ---
# title:  lidar_normalize.R
# author: Nathan Corroyez, UMR TETIS
# desc:   LiDAR point cloud normalisation wrappers (DSM- and DTM-based).
#         Replaces end_dsm_norm / end_dtm_norm from
#         02_CODES/libraries/functions_lidar.R.
#
#         normalize_dsm() : CHM-based "flattening" — aligns every pixel's
#           canopy top to h_target m, making LAD profiles comparable across
#           cells regardless of absolute canopy height. Matches end_dsm_norm
#           semantics exactly (class-4 shift, 2 m vegetation floor).
#         normalize_dtm() : standard DTM-based normalisation — heights are
#           referenced to the ground. h_min and z_max are arguments to
#           support the reviewer sensitivity analysis.
# ---

library("lidR")
library("cli")

# ── normalize_dsm ──────────────────────────────────────────────────────────────

#' @title Normalise LiDAR point cloud with CHM-based canopy flattening
#'
#' @description
#' Produces a "flattened" point cloud in which every pixel's canopy top is
#' aligned to `h_target` metres, regardless of absolute canopy height. On data
#' pre-classified into two ASPRS classes (2 = ground, 4 = vegetation), this is
#' mathematically equivalent to Z_flat = Z − CHM(x,y) + h_target applied
#' point-by-point (cf. §2.5.3 of manuscript RSE-D-25-04417).
#'
#' This function reproduces the semantics of `end_dsm_norm()` in
#' `02_CODES/libraries/functions_lidar.R` with parameterised CHM resolution,
#' algorithm, vegetation class, and canopy-top anchor height.
#'
#' @details
#' **Role of `h_target`** : `h_target` (default 40 m) is chosen to match the
#' ceiling of the LAD binning grid used downstream by `compute_lad_profile()`
#' (bins from 2.5 m to 39.5 m). Setting the canopy top to exactly 40 m ensures
#' that all canopy-top returns land in the topmost LAD bin.
#'
#' **Physical justification** : Normalising depth from the canopy top
#' rather than from the ground makes within-pixel LAD profiles comparable
#' across cells with different absolute tree heights — the relevant dimension
#' for gap-fraction-based LAI estimation via Beer-Lambert.
#'
#' **Classification assumption** : Input data are assumed to carry exactly two
#' ASPRS classes: class 2 (ground / non-vegetation) and `veg_class`
#' (vegetation, default 4). Points of any other class are left unchanged.
#' A warning is emitted if no points of `veg_class` are found in the tile
#' (degenerate tile — no shift is applied).
#'
#' **Vegetation floor** : After the canopy-top alignment, points with
#' Z < 2 m are reclassified to class 2 and set to Z = 0 (reproduces lines
#' 576-577 of `end_dsm_norm` in `02_CODES/libraries/functions_lidar.R`).
#' This ensures that near-ground sub-canopy returns are treated as gap by
#' `compute_pai()` (refactored `myPAI`) in the downstream LAI estimation.
#' The 2 m floor follows Bouvier et al. (2015) and is currently hardcoded to
#' match the published analysis. If a sensitivity test on the floor threshold
#' is needed, add an `h_min` argument in a future revision.
#'
#' @param las      A `LAS` object (lidR). Must contain `Z` and `Classification`
#'   columns. Returned unchanged (NULL) if empty.
#' @param h_target Numeric. Target height in metres for the canopy top after
#'   flattening. Default 40 m (matches LAD grid ceiling).
#' @param chm_res  Numeric. Resolution in metres for the CHM raster computed
#'   internally. Default 1 m.
#' @param chm_algorithm lidR algorithm object. Algorithm passed to
#'   `lidR::rasterize_canopy()` for CHM computation. Default
#'   `lidR::pitfree()` (Khosravipour et al. 2014).
#' @param veg_class Integer. ASPRS classification code for vegetation points.
#'   Default 4L. The canopy-top shift is applied only to these points.
#'
#' @return A `LAS` object with modified `Z` values, or `NULL` if `las` is
#'   empty. The returned object has the same number of points as the input
#'   after the NA-height filter; no points are removed beyond that step.
#'
#' @references
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Generalizing
#' predictive models of forest inventory attributes using an area-based approach
#' with airborne LiDAR data. Remote Sensing of Environment, 156, 322–334.
#' \doi{10.1016/j.rse.2014.10.004}
#'
#' @examples
#' \dontrun{
#' las <- lidR::readLAS("path/to/tile.las")
#' nlas <- normalize_dsm(las)
#' # With alternative CHM algorithm:
#' nlas <- normalize_dsm(las, chm_algorithm = lidR::p2r(0.2))
#' }
#'
#' @export
normalize_dsm <- function(las,
                          h_target      = 40,
                          chm_res       = 1,
                          chm_algorithm = lidR::pitfree(),
                          veg_class     = 4L) {
  if (is.null(las) || lidR::npoints(las) == 0L) return(NULL)

  # A. Compute CHM and normalise heights per-point (Z <- Z - CHM(x,y))
  chm  <- lidR::rasterize_canopy(las, res = chm_res, pkg = "terra",
                                 algorithm = chm_algorithm)
  nlas <- lidR::normalize_height(las, algo = chm)
  rm(chm)

  # Remove points whose height could not be normalised (NA DSM cells)
  nlas <- lidR::filter_poi(nlas, !is.na(Z))

  # B. Identify vegetation points and apply canopy-top shift
  #    Z[Z > 0] <- 0 : clamp positive outliers (points above local DSM)
  #    Z <- Z + h_target - max(Z) : shift top of canopy to h_target
  cls_veg <- nlas$Classification == veg_class
  if (!any(cls_veg)) {
    cli::cli_warn(
      "normalize_dsm(): no points with Classification == {veg_class} found \\
       in tile. Canopy-top shift skipped — returning post-normalize_height LAS."
    )
    return(nlas)
  }
  nlas$Z[cls_veg & nlas$Z > 0] <- 0
  nlas$Z[cls_veg] <- nlas$Z[cls_veg] + h_target - max(nlas$Z[cls_veg])

  # C. Apply 2 m vegetation floor (Bouvier et al. 2015):
  #    points with Z < 2 → reclassify as class 2, set Z = 0
  nlas <- lidR::classify_poi(nlas, class = 2L, poi = ~Z < 2)
  nlas$Z[nlas$Classification == 2L] <- 0

  gc()
  return(nlas)
}


# ── normalize_dtm ──────────────────────────────────────────────────────────────

#' @title Normalise LiDAR point cloud using DTM (ground-relative normalisation)
#'
#' @description
#' Normalises point heights relative to a Digital Terrain Model (DTM) computed
#' from the point cloud using the TIN algorithm. Heights are expressed as metres
#' above ground level (Z = 0 at ground, positive upward). This is a standard
#' lidR normalisation with no canopy-top shift.
#'
#' Reproduces the semantics of `end_dtm_norm()` in
#' `02_CODES/libraries/functions_lidar.R`, with `h_min` and `z_max` exposed
#' as arguments to support the reviewer sensitivity analysis on the vegetation
#' height threshold.
#'
#' @details
#' **Difference from `normalize_dsm()`** : no canopy-top shift is applied here.
#' Heights are ground-referenced; canopy tops at different absolute heights
#' remain at different Z values. `compute_lad_profile()` with `norm = "dtm"`
#' applies the canopy-top shift internally before computing LAD layers.
#'
#' **Height ceiling** : Before the floor reclassification, points with
#' Z > `z_max` are removed via `filter_poi(Z <= z_max)`. This discards
#' spurious high returns that would otherwise inflate upper-bin counts in
#' `compute_lad_profile()`.
#'
#' **Vegetation floor** : Points with Z < `h_min` are then reclassified to
#' class 2 and set to Z = 0, ensuring they are treated as gap by
#' `compute_pai()` (refactored `myPAI`) downstream. Default `h_min = 2`
#' matches the published analysis (Bouvier et al. 2015). Set
#' `h_min ∈ {2, 3, 5}` for the reviewer sensitivity analysis.
#'
#' @param las    A `LAS` object (lidR).
#' @param h_min  Numeric. Vegetation height floor in metres. Points below this
#'   value are reclassified to class 2 and set to Z = 0. Default 2 m (Bouvier
#'   et al. 2015). For reviewer sensitivity analysis use {2, 3, 5}.
#' @param z_max  Numeric. Maximum height in metres. Points above this value are
#'   removed. Default 40 m.
#'
#' @return A `LAS` object with ground-referenced heights, or `NULL` if `las`
#'   is empty.
#'
#' @references
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Generalizing
#' predictive models of forest inventory attributes using an area-based approach
#' with airborne LiDAR data. Remote Sensing of Environment, 156, 322–334.
#' \doi{10.1016/j.rse.2014.10.004}
#'
#' @examples
#' \dontrun{
#' las <- lidR::readLAS("path/to/tile.las")
#' nlas_dtm <- normalize_dtm(las)
#' # h_min sensitivity:
#' nlas_dtm3 <- normalize_dtm(las, h_min = 3)
#' }
#'
#' @export
normalize_dtm <- function(las, h_min = 2, z_max = 40) {
  if (is.null(las) || lidR::npoints(las) == 0L) return(NULL)

  # A. Compute DTM via TIN interpolation and normalise heights
  dtm  <- lidR::rasterize_terrain(las, res = 1, pkg = "terra",
                                  algorithm = lidR::tin())
  nlas <- lidR::normalize_height(las, dtm)
  rm(dtm)

  # B. Remove points above the ceiling
  nlas <- lidR::filter_poi(nlas, Z <= z_max)

  # C. Apply vegetation floor: points below h_min → class 2, Z = 0
  nlas <- lidR::classify_poi(nlas, class = 2L, poi = ~Z < h_min)
  nlas$Z[nlas$Classification == 2L] <- 0

  gc()
  return(nlas)
}
