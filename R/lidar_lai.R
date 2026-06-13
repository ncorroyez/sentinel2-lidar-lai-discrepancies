# ---
# title:  lidar_lai.R
# author: Nathan Corroyez, UMR TETIS
# desc:   Refactored LAI/LAD computation functions.
#         Replaces myPAI / myPAD / myPAD_dtm from 02_CODES/libraries/functions_lidar.R.
#         Key changes:
#           - compute_pai(): h_min and k are arguments (sensitivity analysis)
#           - compute_lad_profile(): k and z_max are arguments; z0 (binning
#             origin) is fixed at 2 m to preserve comparability with the
#             published profiles. h_min sensitivity is handled upstream in
#             normalize_dsm() / normalize_dtm().
# ---

library("lidR")

# ── compute_pai ────────────────────────────────────────────────────────────────

#' @title Compute Plant Area Index (PAI) from LiDAR gap fraction
#'
#' @description
#' Estimates the Plant Area Index (PAI) of a forest stand using the integrated
#' gap fraction approach and Beer-Lambert law. Points with normalised height
#' below `h_min` are treated as gap returns (ground-level or sub-canopy);
#' the gap fraction P = n_gap / n_total gives PAI = -log(P) / k.
#'
#' The underlying theory (Beer-Lambert applied to ALS vertical gap fraction)
#' originates from MacArthur & Horn (1969) and has been widely applied to
#' airborne LiDAR (Morsdorf et al. 2006; Solberg et al. 2009; Korhonen et al.
#' 2011). The specific implementation here (h_min = 2 m vegetation floor,
#' k = 0.5) follows Bouvier et al. (2015).
#'
#' This function is designed to be called inside `lidR::pixel_metrics()`. It
#' returns a single numeric value per pixel.
#'
#' @param z Numeric vector. Normalised point heights in metres (Z coordinate
#'   from a LiDAR point cloud after CHM- or DTM-based normalisation).
#' @param h_min Numeric. Minimum vegetation height threshold in metres. Points
#'   below this value are treated as non-canopy returns. Default 2 m (Bouvier
#'   et al. 2015 convention). For the reviewer sensitivity analysis use
#'   h_min ∈ {2, 3, 5}.
#' @param k Numeric. Beer-Lambert extinction coefficient (dimensionless).
#'   Default 0.5. For the reviewer sensitivity analysis use k ∈ {0.4, 0.5, 0.6}.
#'
#' @return Numeric scalar. PAI estimate (m² m⁻²). Returns NA if `z` is empty
#'   or contains only NA values. Returns Inf if no point falls below `h_min`
#'   (perfectly closed canopy); downstream scripts should replace Inf with NA.
#'
#' @references
#' MacArthur R.H., Horn H.S. (1969). Foliage profile by vertical measurements.
#' Ecology, 50(5), 802–804. \doi{10.2307/1933693}
#'
#' Morsdorf F., Kötz B., Meier E., Itten K.I., Allgöwer B. (2006). Estimation
#' of LAI and fractional cover from small footprint airborne laser scanning data
#' based on gap fraction. Remote Sensing of Environment, 104, 50–61.
#' \doi{10.1016/j.rse.2006.04.019}
#'
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Generalizing
#' predictive models of forest inventory attributes using an area-based approach
#' with airborne LiDAR data. Remote Sensing of Environment, 156, 322–334.
#' \doi{10.1016/j.rse.2014.10.004}
#'
#' @examples
#' z <- c(0.5, 1.2, 3.4, 8.1, 12.3, 18.5, 22.0)
#' compute_pai(z)                        # h_min = 2, k = 0.5
#' compute_pai(z, h_min = 3, k = 0.4)   # sensitivity analysis
#'
#' @export
compute_pai <- function(z, h_min = 2, k = 0.5) {
  z <- z[!is.na(z)]
  if (length(z) == 0L) return(NA_real_)
  n_gap   <- sum(z < h_min)
  n_total <- length(z)
  if (n_gap == 0L) return(Inf)   # Inf downstream → replace with NA
  -log(n_gap / n_total) / k
}


# ── compute_lad_profile ────────────────────────────────────────────────────────

#' @title Compute Leaf Area Density (LAD) vertical profile from LiDAR
#'
#' @description
#' Computes the layer-wise Leaf Area Density (LAD, m² m⁻² m⁻¹) profile and
#' the cumulative Plant Area Index (PAI, m² m⁻²) from each layer upward,
#' using the gap fraction method in `lidR::LAD()` (Beer-Lambert law,
#' Bouvier et al. 2015).
#'
#' Two normalisation conventions are supported, matching the preprocessing
#' functions in `R/lidar_normalize.R`:
#' \describe{
#'   \item{"dsm"}{CHM-based: point heights are referenced to the local canopy
#'     top (output of `normalize_dsm()`). This is the convention used in the
#'     submitted paper.}
#'   \item{"dtm"}{DTM-based: point heights are referenced to ground level
#'     (output of `normalize_dtm()`). Before computing LAD, heights are shifted
#'     into the top-of-canopy frame via z' = z + z_max - max(z).}
#' }
#'
#' The vertical binning grid is always anchored at z0 = 2 m (bin centres at
#' 2.5, 3.5, ..., z_max - 0.5 m), preserving comparability with the published
#' profiles. Sensitivity to the h_min floor is handled upstream by
#' `normalize_dsm()` / `normalize_dtm()`, which set sub-threshold points to
#' Z = 0 before this function is called.
#'
#' The function is designed to be called inside `lidR::pixel_metrics()`.
#'
#' @param z Numeric vector. Normalised point heights in metres.
#' @param norm Character. Normalisation convention: `"dsm"` (CHM-based, default)
#'   or `"dtm"` (DTM-based).
#' @param z_max Numeric. Maximum canopy height in metres used as the top-of-canopy
#'   anchor. Default 40 m.
#' @param k Numeric. Beer-Lambert extinction coefficient. Default 0.5.
#'   Use {0.4, 0.5, 0.6} for the reviewer k-sensitivity analysis.
#'
#' @return Named list of numeric values intended for use with
#'   `lidR::pixel_metrics()`. Two groups:
#'   \describe{
#'     \item{LAD_Layer_{h}}{LAD value (m² m⁻² m⁻¹) at bin centre h (2.5 to
#'       z_max - 0.5 m, 1 m steps). NA where no points exist in that layer.}
#'     \item{PAI_{h}}{Cumulative PAI from layer h to `z_max` (m² m⁻²). NA
#'       where the corresponding LAD layer is NA.}
#'   }
#'
#' @references
#' Bouvier M., Durrieu S., Fournier R.A., Renaud J.P. (2015). Generalizing
#' predictive models of forest inventory attributes using an area-based approach
#' with airborne LiDAR data. Remote Sensing of Environment, 156, 322–334.
#' \doi{10.1016/j.rse.2014.10.004}
#'
#' @examples
#' set.seed(42)
#' z <- runif(500, 0, 35)
#' prof <- compute_lad_profile(z, norm = "dsm", z_max = 40, k = 0.5)
#' prof[["LAD_Layer_2.5"]]
#' # k sensitivity (grid stays the same, LAD values change):
#' prof_k04 <- compute_lad_profile(z, norm = "dsm", z_max = 40, k = 0.4)
#'
#' @export
compute_lad_profile <- function(z,
                                norm  = c("dsm", "dtm"),
                                z_max = 40,
                                k     = 0.5) {
  norm <- match.arg(norm)

  # DTM convention: shift heights to top-of-canopy frame (same as myPAD_dtm)
  if (norm == "dtm") {
    max_z <- max(z, na.rm = TRUE)
    if (is.infinite(max_z) || is.na(max_z)) {
      return(.empty_lad_profile(z_max))
    }
    z <- z + z_max - max_z
  }

  # Compute LAD using lidR::LAD.
  # z0 = 2 is fixed (matches original myPAD / myPAD_dtm) so the binning grid
  # is always 2.5, 3.5, ..., z_max - 0.5 m regardless of h_min sensitivity.
  # k is passed through for the reviewer k-sensitivity analysis.
  lad_raw <- tryCatch(
    lidR::LAD(z, z0 = 2, dz = 1, k = k),
    error = function(e) NULL
  )
  if (is.null(lad_raw) || nrow(lad_raw) == 0L) {
    return(.empty_lad_profile(z_max))
  }

  # Snap bin centres to 2.5, 3.5, …
  lad_raw$z <- floor(lad_raw$z) + 0.5

  # Full target grid: one bin per metre from 2 m to z_max
  bin_centers <- seq(2.5, z_max - 0.5, by = 1)

  # Align lad_raw to the full grid (vectorised, no for-loop)
  matched  <- match(bin_centers, lad_raw$z)
  lad_vals <- ifelse(is.na(matched), NA_real_, lad_raw$lad[matched])

  # Cumulative PAI from each layer upward (NA propagates)
  pai_vals <- vapply(seq_along(bin_centers), function(i) {
    if (is.na(lad_vals[i])) NA_real_
    else sum(lad_vals[i:length(lad_vals)], na.rm = TRUE)
  }, numeric(1L))

  lad_list <- setNames(as.list(lad_vals), paste0("LAD_Layer_", bin_centers))
  pai_list <- setNames(as.list(pai_vals), paste0("PAI_",       bin_centers))

  c(pai_list, lad_list)
}


# ── internal helper ────────────────────────────────────────────────────────────

#' Return an all-NA LAD profile (used when input is degenerate)
#' @noRd
.empty_lad_profile <- function(z_max) {
  bin_centers <- seq(2.5, z_max - 0.5, by = 1)
  na_list  <- setNames(rep(list(NA_real_), length(bin_centers)), bin_centers)
  lad_list <- setNames(na_list, paste0("LAD_Layer_", bin_centers))
  pai_list <- setNames(na_list, paste0("PAI_",       bin_centers))
  c(pai_list, lad_list)
}
