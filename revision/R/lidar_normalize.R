# ---
# title:  lidar_normalize.R
# author: Nathan Corroyez, UMR TETIS
# desc:   Refactored LiDAR point cloud normalisation wrappers (DSM- and
#         DTM-based). Replaces end_dsm_norm / end_dtm_norm from
#         02_CODES/libraries/functions_lidar.R.
#         Key change: h_min is now an argument in both functions so the 2 m
#         vegetation floor can be varied for the reviewer sensitivity analysis.
# ---

library("lidR")

# ── normalize_dsm ──────────────────────────────────────────────────────────────

#' @title Normalise LiDAR point cloud using DSM (CHM-based normalisation)
#'
#' @description
#' Normalises point heights relative to the local canopy surface model (DSM),
#' producing the "top-of-canopy" reference frame used in the submitted paper.
#' After normalisation:
#' \enumerate{
#'   \item Ground-class (class 4) points are shifted so that their maximum
#'     height becomes `z_max` m: Z' = Z + z_max - max(Z[class == 4]).
#'   \item Points with normalised height < `h_min` are reclassified to class 2
#'     and their Z set to 0 (sub-canopy / ground floor convention, Bouvier 2015).
#'   \item Points with Z > `z_max` or NA are removed.
#' }
#'
#' Corresponds to `end_dsm_norm()` in `02_CODES/libraries/functions_lidar.R`,
#' with the 2 m floor and 40 m ceiling made into arguments.
#'
#' @param las A `LAS` object (lidR). Must contain a `Classification` column
#'   with class 4 for ground points.
#' @param h_min Numeric. Minimum vegetation height threshold in metres. Points
#'   below this value are set to Z = 0. Default 2 m (Bouvier 2015). For the
#'   reviewer sensitivity analysis use h_min ∈ {2, 3, 5}.
#' @param z_max Numeric. Maximum canopy height anchor in metres. Ground points
#'   are shifted so their top reaches `z_max`. Default 40 m.
#'
#' @return A normalised `LAS` object, or `NULL` if the input is empty.
#'
#' @references
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Remote Sensing of
#' Environment, 156, 322–334. \doi{10.1016/j.rse.2014.10.004}
#'
#' @export
normalize_dsm <- function(las, h_min = 2, z_max = 40) {
  if (is.null(las) || lidR::npoints(las) == 0L) return(NULL)

  # 1. Compute DSM from pit-free algorithm and normalise heights
  dsm  <- lidR::rasterize_canopy(las, res = 1, pkg = "terra",
                                 algorithm = lidR::pitfree())
  nlas <- lidR::normalize_height(las, algo = dsm)
  rm(dsm)

  # 2. Remove NA heights
  nlas <- lidR::filter_poi(nlas, !is.na(Z))

  # 3. Shift class-4 points to top-of-canopy frame (mirrors end_dsm_norm).
  #    Class-4 points that are already > 0 in DSM-normalised frame are noise:
  #    clamp to 0 first, then shift so their maximum reaches z_max.
  cls4 <- nlas$Classification == 4L
  if (any(cls4)) {
    nlas$Z[cls4 & nlas$Z > 0] <- 0
    max_cls4 <- max(nlas$Z[cls4])
    nlas$Z[cls4] <- nlas$Z[cls4] + z_max - max_cls4
  }

  # 4. Apply vegetation floor: points below h_min → class 2, Z = 0
  nlas <- lidR::classify_poi(nlas, class = 2L, poi = ~Z < h_min)
  nlas$Z[nlas$Classification == 2L] <- 0

  gc()
  return(nlas)
}


# ── normalize_dtm ──────────────────────────────────────────────────────────────

#' @title Normalise LiDAR point cloud using DTM (ground-relative normalisation)
#'
#' @description
#' Normalises point heights relative to a Digital Terrain Model (DTM) computed
#' with the TIN algorithm (kriging is available as an alternative but slower).
#' After normalisation:
#' \enumerate{
#'   \item Heights above `z_max` are removed.
#'   \item Points with normalised height < `h_min` are reclassified to class 2
#'     and their Z set to 0.
#' }
#'
#' Corresponds to `end_dtm_norm()` in `02_CODES/libraries/functions_lidar.R`,
#' with the 2 m floor and 40 m ceiling made into arguments.
#'
#' @param las A `LAS` object (lidR).
#' @param h_min Numeric. Minimum vegetation height threshold in metres. Points
#'   below this value are set to Z = 0. Default 2 m. For the reviewer
#'   sensitivity analysis use h_min ∈ {2, 3, 5}.
#' @param z_max Numeric. Maximum height (metres). Points above this are removed.
#'   Default 40 m.
#'
#' @return A normalised `LAS` object, or `NULL` if the input is empty.
#'
#' @references
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Remote Sensing of
#' Environment, 156, 322–334. \doi{10.1016/j.rse.2014.10.004}
#'
#' @export
normalize_dtm <- function(las, h_min = 2, z_max = 40) {
  if (is.null(las) || lidR::npoints(las) == 0L) return(NULL)

  # 1. Compute DTM with TIN interpolation and normalise heights
  dtm  <- lidR::rasterize_terrain(las, res = 1, pkg = "terra",
                                  algorithm = lidR::tin())
  nlas <- lidR::normalize_height(las, dtm)
  rm(dtm)

  # 2. Remove points above z_max
  nlas <- lidR::filter_poi(nlas, Z <= z_max)

  # 3. Apply vegetation floor: points below h_min → class 2, Z = 0
  nlas <- lidR::classify_poi(nlas, class = 2L, poi = ~Z < h_min)
  nlas$Z[nlas$Classification == 2L] <- 0

  gc()
  return(nlas)
}
