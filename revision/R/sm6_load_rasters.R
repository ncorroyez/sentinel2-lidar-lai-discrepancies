# ---
# title:  sm6_load_rasters.R
# desc:   Load per-site rasters for SM6b analysis.
#         Exposes load_site_rasters() and the internal pad_filename() helper.
#
#         Faithful refactor of the raster-loading block in:
#         02_CODES/three_factors_analysis_final.R lines 64-86.
# ---

# ── pad_filename ───────────────────────────────────────────────────────────────

#' @title Build PAD raster filename from d_opt value
#'
#' @description
#' Constructs the filename for a PAD (Plant Area Density) raster slice
#' following the naming convention \code{PAD_X_Y.tif}, where \code{Y} is the
#' fixed canopy maximum height and \code{X} is the centre of the 1 m-thick
#' LAD layer at depth \code{dopt_value} measured from the ground:
#'
#' \deqn{X = \text{canopy\_max\_m} - \text{dopt\_value} + 0.5}
#'
#' The \code{+0.5} offset places \code{X} at the centre of the 1 m-thick LAD
#' pixel at that depth (each layer represents a 1 m slice of the canopy;
#' the centre of the layer is 0.5 m from either edge). For example,
#' \code{dopt_value = 4} m from the ground in a 40 m canopy corresponds to
#' the layer spanning [3, 4] m above ground, centred at 3.5 m
#' (\eqn{X = 40 - 4 + 0.5 = 36.5}), giving \code{"PAD_36.5_40.tif"}.
#'
#' \code{format(x_value, nsmall = 1, decimal.mark = ".")} ensures a point
#' decimal separator on all systems (avoids locale-dependent comma).
#'
#' @param dopt_value   Numeric. Optical depth in metres from the ground.
#'   For the submitted manuscript (combined configuration): \code{4},
#'   giving \code{"PAD_36.5_40.tif"}.
#' @param canopy_max_m Numeric. Fixed canopy maximum height in metres.
#'   Default \code{40}.
#'
#' @return Character. Filename string, e.g. \code{"PAD_36.5_40.tif"}.
#'
#' @examples
#' pad_filename(4)    # "PAD_36.5_40.tif"
#' pad_filename(3.5)  # "PAD_37.0_40.tif"
#'
#' @keywords internal
pad_filename <- function(dopt_value, canopy_max_m = 40) {
  x_value <- canopy_max_m - dopt_value + 0.5
  paste0("PAD_",
         format(x_value, nsmall = 1, decimal.mark = "."),
         "_", canopy_max_m, ".tif")
}

# ── load_site_rasters ──────────────────────────────────────────────────────────

#' @title Load all rasters required for SM6b analysis for one site
#'
#' @description
#' Loads seven 10 m SpatRasters needed for the SM6b heterogeneity analysis:
#' \itemize{
#'   \item \strong{lai_als}: Integrated LAI from the ALS vertical LAD profile
#'     (\code{ladstack_classic.tif}). This raster is multi-layer, representing
#'     the vertical LAD profile (one layer per 1 m height bin). Summing across
#'     layers via \code{sum(rast(...), na.rm = TRUE)} produces a single-band
#'     raster where each pixel is the integrated LAI over the full canopy
#'     height (Bouvier et al. 2015 formulation). A \code{stop()} is raised if
#'     \code{terra::nlyr()} returns 1, indicating an unexpected single-layer
#'     file that would silently produce incorrect values.
#'   \item \strong{lai_als_dopt}: ALS LAI at optical depth \code{dopt_value}.
#'     Loaded from \code{PAD_Profiles_dsm_keepTrees/PAD_X_Y.tif} where
#'     \code{X} and \code{Y} are derived via \code{pad_filename(dopt_value)}.
#'     Single-layer.
#'   \item \strong{lai_s2_atbd}: Sentinel-2 LAI from the ATBD
#'     parametrisation (\code{s2lai_summer_atbd_res_10_m.tif}). Single-layer.
#'   \item \strong{lai_s2_opt}: Sentinel-2 LAI from the site-optimal
#'     per-site parametrisation
#'     (\code{s2lai_summer_best_indiv_res_10_m.tif}). Single-layer.
#'   \item \strong{dsm_sd}: DSM-based canopy height standard deviation at
#'     10 m, produced by SM6a (\code{dsm_sd_res_10_m.tif}).
#'   \item \strong{chm_sd}: CHM-based equivalent, produced by SM6a
#'     (\code{chm_sd_res_10_m.tif}).
#'   \item \strong{max_height}: Maximum canopy height at 10 m
#'     (\code{max_res_10_m.tif}).
#' }
#'
#' \strong{Note — lai_s2_common excluded}:
#' \code{s2lai_summer_depth_study_common_res_10_m.tif} is intentionally
#' excluded from this module. It is only used in the 'Sites combined'
#' analysis (legacy lines 316-328 of three_factors_analysis_final.R), which
#' is out of scope for SM6b passe 1. A future SM6c module may reintroduce it.
#'
#' @param site       Character. Site name, one of
#'   \code{c("Aigoual", "Blois", "Mormal")}.
#' @param dopt_value Numeric. Optical depth in metres, used to derive the PAD
#'   filename via \code{pad_filename()}. Default \code{4} (value reported in
#'   the submitted manuscript for the combined configuration).
#' @param sm6a_dir   Character. Root directory of SM6a outputs (the directory
#'   that contains \code{{site}/dsm_sd_res_10_m.tif} and its CHM counterpart).
#'   Typically \code{here::here("revision", "output", "intermediate", "sm6")}.
#' @param ext_dir    Character. Root of the external results directory.
#'   Typically \code{here::here("03_RESULTS")}. The function looks for rasters
#'   under \code{ext_dir/{site}/Metrics/Deciduous_Only/}.
#'
#' @return A named list of seven \code{SpatRaster} objects:
#'   \code{lai_als}, \code{lai_als_dopt}, \code{lai_s2_atbd},
#'   \code{lai_s2_opt}, \code{dsm_sd}, \code{chm_sd}, \code{max_height}.
#'
#' @references
#' Bouvier, M., Durrieu, S., Fournier, R. A., & Renaud, J. P. (2015).
#' Generalizing relationships between canopy structure and light transmission:
#' the case of European temperate forests for a radiative transfer model
#' parametrization. Remote Sensing of Environment, 171, 158-171.
#'
#' Legacy: 02_CODES/three_factors_analysis_final.R lines 64-86.
#'
#' @examples
#' \dontrun{
#' r <- load_site_rasters(
#'   site       = "Blois",
#'   dopt_value = 4,
#'   sm6a_dir   = here::here("revision", "output", "intermediate", "sm6"),
#'   ext_dir    = here::here("03_RESULTS")
#' )
#' terra::nlyr(r$lai_als)  # > 1
#' }
#'
#' @export
load_site_rasters <- function(site, dopt_value = 4, sm6a_dir, ext_dir) {
  dec_only <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")

  # ── ladstack: multi-layer → sum to get LAI_ALS ──────────────────────────────
  ladstack_path <- file.path(dec_only, "ladstack_classic.tif")
  ladstack      <- terra::rast(ladstack_path)
  if (terra::nlyr(ladstack) == 1L) {
    stop(
      "load_site_rasters: ladstack_classic.tif for site '", site, "' has ",
      "only 1 layer. Expected a multi-layer LAD profile stack. ",
      "Check that the correct file is referenced."
    )
  }
  lai_als <- sum(ladstack, na.rm = TRUE)

  # ── PAD at d_opt (single-layer) ──────────────────────────────────────────────
  pad_fn   <- pad_filename(dopt_value)
  pad_path <- file.path(dec_only, "PAD_Profiles_dsm_keepTrees", pad_fn)
  lai_als_dopt <- terra::rast(pad_path)

  # ── S2 rasters (single-layer each) ───────────────────────────────────────────
  lai_s2_atbd <- terra::rast(
    file.path(dec_only, "s2lai_summer_atbd_res_10_m.tif")
  )
  lai_s2_opt <- terra::rast(
    file.path(dec_only, "s2lai_summer_best_indiv_res_10_m.tif")
  )

  # ── SM6a heterogeneity rasters ───────────────────────────────────────────────
  sm6a_site <- file.path(sm6a_dir, site)
  dsm_sd    <- terra::rast(file.path(sm6a_site, "dsm_sd_res_10_m.tif"))
  chm_sd    <- terra::rast(file.path(sm6a_site, "chm_sd_res_10_m.tif"))

  # ── Max height (single-layer) ─────────────────────────────────────────────────
  max_height <- terra::rast(file.path(dec_only, "max_res_10_m.tif"))

  list(
    lai_als      = lai_als,
    lai_als_dopt = lai_als_dopt,
    lai_s2_atbd  = lai_s2_atbd,
    lai_s2_opt   = lai_s2_opt,
    dsm_sd       = dsm_sd,
    chm_sd       = chm_sd,
    max_height   = max_height
  )
}
