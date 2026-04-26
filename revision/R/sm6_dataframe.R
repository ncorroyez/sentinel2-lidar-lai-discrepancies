# ---
# title:  sm6_dataframe.R
# desc:   Build the multi-site pixel-level data.table from per-site rasters.
#         Exposes build_pixel_dt() and build_multisite_dt().
#
#         Values are extracted via terra::values(), replicating the legacy
#         approach in 02_CODES/three_factors_analysis_final.R lines 116-122.
#         NA values are retained globally; filters are applied locally within
#         compute_metrics_by_class().
# ---

#' @title Extract pixel values for one site into a data.table
#'
#' @description
#' Converts a named list of SpatRasters (output of \code{load_site_rasters()})
#' into a \code{data.table} with one row per pixel. Values are extracted via
#' \code{terra::values()}, replicating the legacy approach
#' (\code{as.numeric(values(raster))}) in
#' \code{02_CODES/three_factors_analysis_final.R} lines 116-122.
#'
#' NA values are retained in all columns (no global \code{na.omit()}).
#' Per-combination pair-wise NA filtering is performed later inside
#' \code{compute_metrics_by_class()}.
#'
#' @param rasters   Named list of SpatRasters as returned by
#'   \code{load_site_rasters()}.
#' @param site_name Character. Site label to populate the \code{site} column.
#'
#' @return A \code{data.table} with 8 columns:
#' \describe{
#'   \item{site}{Character. Site name.}
#'   \item{lai_als}{Numeric. Integrated LAI from ALS (sum of LAD layers).
#'     Units: m\eqn{^2}/m\eqn{^2}.}
#'   \item{lai_als_dopt}{Numeric. ALS LAI at optical depth d_opt.}
#'   \item{lai_s2_atbd}{Numeric. Sentinel-2 LAI, ATBD parametrisation.}
#'   \item{lai_s2_opt}{Numeric. Sentinel-2 LAI, site-optimal parametrisation.}
#'   \item{dsm_sd}{Numeric. DSM-based canopy height SD at 10 m (from SM6a).
#'     Units: m.}
#'   \item{chm_sd}{Numeric. CHM-based canopy height SD at 10 m (from SM6a).
#'     Units: m.}
#'   \item{max_height}{Numeric. Maximum canopy height at 10 m. Units: m.}
#' }
#'
#' @examples
#' \dontrun{
#' r  <- load_site_rasters(
#'   "Blois", dopt_value = 4,
#'   sm6a_dir = here::here("revision", "output", "intermediate", "sm6"),
#'   ext_dir  = here::here("03_RESULTS")
#' )
#' dt <- build_pixel_dt(r, "Blois")
#' dt[, lapply(.SD, function(x) sum(!is.na(x)))]  # count non-NA per column
#' }
#'
#' @export
build_pixel_dt <- function(rasters, site_name) {
  data.table::data.table(
    site         = site_name,
    lai_als      = as.numeric(terra::values(rasters$lai_als)),
    lai_als_dopt = as.numeric(terra::values(rasters$lai_als_dopt)),
    lai_s2_atbd  = as.numeric(terra::values(rasters$lai_s2_atbd)),
    lai_s2_opt   = as.numeric(terra::values(rasters$lai_s2_opt)),
    dsm_sd       = as.numeric(terra::values(rasters$dsm_sd)),
    chm_sd       = as.numeric(terra::values(rasters$chm_sd)),
    max_height   = as.numeric(terra::values(rasters$max_height))
  )
}


#' @title Build multi-site pixel-level data.table
#'
#' @description
#' Iterates over a vector of site names, loads per-site rasters via
#' \code{load_site_rasters()}, extracts pixel values via
#' \code{build_pixel_dt()}, and row-binds the results into a single
#' \code{data.table}.
#'
#' @param sites        Character vector. Site names
#'   (e.g. \code{c("Aigoual", "Blois", "Mormal")}).
#' @param dopt_by_site Named numeric vector or single numeric. Optical depth
#'   in metres, per site. If a named vector is supplied (names = site names),
#'   each site uses its own d_opt. If a scalar is supplied, all sites use the
#'   same value. Default \code{4} (backward-compatible).
#' @param sm6a_dir     Character. Root of SM6a output directory.
#'   Typically \code{here::here("revision", "output", "intermediate", "sm6")}.
#' @param ext_dir      Character. Root of external results directory.
#'   Typically \code{here::here("03_RESULTS")}.
#'
#' @return A \code{data.table} with 8 columns (see \code{build_pixel_dt()})
#'   and \code{sum(n_pixels_per_site)} rows.
#'
#' @examples
#' \dontrun{
#' dt <- build_multisite_dt(
#'   sites        = c("Aigoual", "Blois", "Mormal"),
#'   dopt_by_site = c(Aigoual = 4, Blois = 27, Mormal = 13),
#'   sm6a_dir     = here::here("revision", "output", "intermediate", "sm6"),
#'   ext_dir      = here::here("03_RESULTS")
#' )
#' dt[, .N, by = site]
#' }
#'
#' @export
build_multisite_dt <- function(sites, dopt_by_site = 4, sm6a_dir, ext_dir,
                               lai_s2_opt_paths = NULL) {
  dt_list <- vector("list", length(sites))
  names(dt_list) <- sites

  for (site in sites) {
    d_val    <- if (length(dopt_by_site) == 1L) dopt_by_site else dopt_by_site[[site]]
    opt_path <- if (!is.null(lai_s2_opt_paths)) lai_s2_opt_paths[[site]] else NULL
    cli::cli_alert_info("Loading rasters for site: {site} (d_opt = {d_val} m)")
    r               <- load_site_rasters(site, d_val, sm6a_dir, ext_dir,
                                         lai_s2_opt_path = opt_path)
    dt_list[[site]] <- build_pixel_dt(r, site)
    cli::cli_alert_success("  {site}: {nrow(dt_list[[site]])} pixels loaded")
  }

  data.table::rbindlist(dt_list)
}
