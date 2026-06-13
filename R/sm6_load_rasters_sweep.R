# ---
# title:  sm6_load_rasters_sweep.R
# desc:   Load per-site rasters for SM6b passe 2 (multi-scale sweep).
#         Exposes load_site_rasters_sweep(), build_pixel_dt_sweep(), and
#         build_multisite_dt_sweep().
#
#         Extends the SM6b passe 1 loading logic (sm6_load_rasters.R) to
#         consume the 6 heterogeneity variants per metric produced by SM6a
#         passe 2 (07_sm6_compute_heterogeneity.R):
#           block_5, block_10, block_20, block_50, focal_30, focal_50
#         for each metric ∈ c("DSM", "CHM") → 12 SD rasters per site.
#
#         Intentionally self-contained: load_site_rasters() from passe 1 is
#         NOT called or sourced here. Code is partially duplicated to ensure
#         that changes to the passe 1 loader do not silently affect the sweep.
#
#         SM6a passe 2 filename convention (from build_output_filename()):
#           {metric_lower}_sd_{method}_{size_m}_m.tif   (on disk)
#         Column name convention (in the data.table):
#           {metric_lower}_sd_{method}_{size_m}         (no trailing _m)
# ---

# ── pad_filename_sweep (private helper) ────────────────────────────────────────
# Replicates the logic of pad_filename() from sm6_load_rasters.R without
# sourcing it, keeping this module self-contained. Named differently to avoid
# a name clash if both modules are sourced in the same session.
pad_filename_sweep <- function(dopt_value, canopy_max_m = 40) {
  x_value <- canopy_max_m - dopt_value + 0.5
  paste0("PAD_",
         format(x_value, nsmall = 1, decimal.mark = "."),
         "_", canopy_max_m, ".tif")
}


# ── load_site_rasters_sweep ────────────────────────────────────────────────────

#' @title Load rasters for SM6b passe 2 sweep analysis for one site
#'
#' @description
#' Loads all rasters required for the SM6b multi-scale heterogeneity sweep
#' analysis for a single site. Returns a flat named list of
#' \code{SpatRaster} objects.
#'
#' \strong{Fixed rasters} (5): replicate the LAI / max-height rasters from
#' \code{load_site_rasters()} in SM6b passe 1. See that function for detailed
#' documentation of each raster.
#' \itemize{
#'   \item \code{lai_als}: Integrated LAI from the ALS LAD profile stack
#'     (\code{ladstack_classic.tif}). Multi-layer sum via
#'     \code{sum(rast(...), na.rm = TRUE)}. A \code{stop()} is raised if only
#'     one layer is found.
#'   \item \code{lai_als_dopt}: ALS LAI at optical depth \code{dopt_value}
#'     (\code{PAD_Profiles_dsm_keepTrees/PAD_X_Y.tif}).
#'   \item \code{lai_s2_atbd}: Sentinel-2 LAI, ATBD parametrisation.
#'   \item \code{lai_s2_opt}: Sentinel-2 LAI, site-optimal parametrisation.
#'   \item \code{max_height}: Maximum canopy height at 10 m.
#' }
#'
#' \strong{Sweep heterogeneity rasters} (12 per site):
#' For each combination of \code{metric} ∈ \code{c("DSM","CHM")} and
#' \code{variant} ∈ \code{sweep_variants}, loads the SD raster produced by
#' SM6a passe 2 (\code{07_sm6_compute_heterogeneity.R}). File naming on disk
#' follows \code{{metric_lower}_sd_{method}_{size_m}_m.tif}; the key in the
#' returned list drops the trailing \code{_m} to avoid ambiguity with the unit
#' suffix: \code{{metric_lower}_sd_{method}_{size_m}}.
#'
#' All expected raster paths are verified before loading. A single
#' \code{stop()} lists all missing files if any are absent.
#'
#' @param site           Character. Site name, one of
#'   \code{c("Aigoual", "Blois", "Mormal")}.
#' @param dopt_value     Numeric. Optical depth in metres used to derive the
#'   PAD filename. Default \code{4}.
#' @param sweep_variants List of named lists, each with elements
#'   \code{method} (character, \code{"block"} or \code{"focal"}) and
#'   \code{size_m} (integer, effective window size in metres).
#' @param sm6a_dir       Character. Root directory of SM6a passe 2 outputs
#'   (contains \code{{site}/{metric_lower}_sd_{method}_{size_m}_m.tif}).
#' @param ext_dir        Character. Root of external results directory
#'   (\code{03_RESULTS}). Fixed rasters are under
#'   \code{ext_dir/{site}/Metrics/Deciduous_Only/}.
#'
#' @return A named list of \code{SpatRaster} objects with 5 fixed entries and
#'   \code{length(sweep_variants) * 2} sweep entries:
#'   \code{$lai_als}, \code{$lai_als_dopt}, \code{$lai_s2_atbd},
#'   \code{$lai_s2_opt}, \code{$max_height},
#'   \code{$dsm_sd_block_5}, \ldots, \code{$chm_sd_focal_50}.
#'
#' @export
load_site_rasters_sweep <- function(site, dopt_value, sweep_variants,
                                    sm6a_dir, ext_dir,
                                    lai_als_dopt_path = NULL) {
  dec_only  <- file.path(ext_dir, site, "Metrics", "Deciduous_Only")
  sm6a_site <- file.path(sm6a_dir, site)
  metrics   <- c("DSM", "CHM")

  # ── Build expected het raster paths (col_key → disk path) ─────────────────
  het_paths <- list()
  for (metric in metrics) {
    for (variant in sweep_variants) {
      col_key <- paste0(
        tolower(metric), "_sd_", variant$method, "_", variant$size_m
      )
      het_paths[[col_key]] <- file.path(
        sm6a_site,
        paste0(tolower(metric), "_sd_",
               variant$method, "_", variant$size_m, "_m.tif")
      )
    }
  }

  # ── Verify all het rasters exist before loading anything ──────────────────
  missing <- character(0)
  for (key in names(het_paths)) {
    if (!file.exists(het_paths[[key]])) {
      missing <- c(missing, het_paths[[key]])
    }
  }
  if (length(missing) > 0L) {
    stop(
      "load_site_rasters_sweep: ", length(missing),
      " het raster(s) missing for site '", site, "':\n",
      paste0("  ", missing, collapse = "\n"),
      "\nRun 07_sm6_compute_heterogeneity.R (SM6a passe 2) first."
    )
  }

  # ── Fixed LAI rasters (replicates load_site_rasters() passe 1 logic) ─────
  ladstack_path <- file.path(dec_only, "ladstack_classic.tif")
  ladstack      <- terra::rast(ladstack_path)
  if (terra::nlyr(ladstack) == 1L) {
    stop(
      "load_site_rasters_sweep: ladstack_classic.tif for site '", site,
      "' has only 1 layer. Expected a multi-layer LAD profile stack."
    )
  }
  lai_als <- sum(ladstack, na.rm = TRUE)

  pad_fn   <- pad_filename_sweep(dopt_value)
  pad_path <- file.path(dec_only, "PAD_Profiles_dsm_keepTrees", pad_fn)

  if (!is.null(lai_als_dopt_path) && file.exists(lai_als_dopt_path)) {
    lai_als_dopt <- terra::rast(lai_als_dopt_path)
  } else if (file.exists(pad_path)) {
    lai_als_dopt <- terra::rast(pad_path)
  } else {
    stop("load_site_rasters_sweep: LAI_ALS_dopt not found for site '", site,
         "' (d_opt = ", dopt_value, "). Tried:\n  ",
         lai_als_dopt_path, "\n  ", pad_path,
         "\nRun scripts/steps/07_compute_lai_als_dopt.R first.")
  }

  rasters <- list(
    lai_als      = lai_als,
    lai_als_dopt = lai_als_dopt,
    lai_s2_atbd  = terra::rast(
      file.path(sm6a_dir, site, "s2lai_summer_atbd_T_res_10_m.tif")
    ),
    lai_s2_opt   = terra::rast(
      file.path(sm6a_dir, site, "s2lai_summer_opt_common_res_10_m.tif")
    ),
    max_height   = terra::rast(
      file.path(dec_only, "max_res_10_m.tif")
    )
  )

  # ── Sweep het rasters ─────────────────────────────────────────────────────
  for (key in names(het_paths)) {
    rasters[[key]] <- terra::rast(het_paths[[key]])
  }

  rasters
}


# ── build_pixel_dt_sweep ───────────────────────────────────────────────────────

#' @title Extract pixel values for one site into a wide data.table (sweep)
#'
#' @description
#' Converts the named list of \code{SpatRaster} objects returned by
#' \code{load_site_rasters_sweep()} into a \code{data.table} with one row
#' per pixel. The 5 fixed LAI / max-height columns are constructed first;
#' the 12 heterogeneity SD columns (one per metric × variant) are appended
#' by reference via \code{data.table::set()}, following the naming convention
#' \code{{metric_lower}_sd_{method}_{size_m}}.
#'
#' NA values are retained (no global \code{na.omit()}). Per-combination
#' pair-wise NA filtering is handled later in
#' \code{compute_metrics_by_class()}.
#'
#' @param rasters        Named list of \code{SpatRaster} objects as returned
#'   by \code{load_site_rasters_sweep()}.
#' @param site_name      Character. Site label for the \code{site} column.
#' @param sweep_variants List of named lists with \code{method} and
#'   \code{size_m} elements. Must match what was passed to
#'   \code{load_site_rasters_sweep()}.
#' @param metrics        Character vector. Heterogeneity metric sources.
#'   Default \code{c("DSM", "CHM")}.
#'
#' @return A \code{data.table} with 18 columns:
#'   \code{site}, \code{lai_als}, \code{lai_als_dopt}, \code{lai_s2_atbd},
#'   \code{lai_s2_opt}, \code{max_height} (6 fixed), plus 12 SD columns in
#'   metric × variant order.
#'
#' @export
build_pixel_dt_sweep <- function(rasters, site_name, sweep_variants,
                                 metrics = c("DSM", "CHM")) {
  # Fixed columns
  dt <- data.table::data.table(
    site         = site_name,
    lai_als      = as.numeric(terra::values(rasters$lai_als)),
    lai_als_dopt = as.numeric(terra::values(rasters$lai_als_dopt)),
    lai_s2_atbd  = as.numeric(terra::values(rasters$lai_s2_atbd)),
    lai_s2_opt   = as.numeric(terra::values(rasters$lai_s2_opt)),
    max_height   = as.numeric(terra::values(rasters$max_height))
  )

  # Sweep het columns — added by reference to avoid a full copy
  for (metric in metrics) {
    for (variant in sweep_variants) {
      col_key <- paste0(
        tolower(metric), "_sd_", variant$method, "_", variant$size_m
      )
      data.table::set(dt, j = col_key,
                      value = as.numeric(terra::values(rasters[[col_key]])))
    }
  }

  dt
}


# ── build_multisite_dt_sweep ───────────────────────────────────────────────────

#' @title Build multi-site pixel-level data.table for the sweep analysis
#'
#' @description
#' Iterates over a vector of site names, loads per-site rasters via
#' \code{load_site_rasters_sweep()}, extracts pixel values via
#' \code{build_pixel_dt_sweep()}, and row-binds the results into a single
#' wide \code{data.table}.
#'
#' @param sites          Character vector. Site names.
#' @param dopt_value     Numeric. Optical depth in metres. Default \code{4}.
#' @param sweep_variants List of named lists with \code{method} and
#'   \code{size_m} elements.
#' @param sm6a_dir       Character. Root of SM6a passe 2 output directory.
#' @param ext_dir        Character. Root of external results directory.
#'
#' @return A \code{data.table} with 18 columns and
#'   \code{sum(n_pixels_per_site)} rows.
#'
#' @export
build_multisite_dt_sweep <- function(sites, dopt_value, sweep_variants,
                                     sm6a_dir, ext_dir,
                                     lai_als_dopt_paths = NULL) {
  dt_list        <- vector("list", length(sites))
  names(dt_list) <- sites

  for (site in sites) {
    cli::cli_alert_info("Loading rasters for site: {site}")
    als_path <- if (!is.null(lai_als_dopt_paths))
      lai_als_dopt_paths[[site]] else NULL
    r <- load_site_rasters_sweep(
      site, dopt_value, sweep_variants, sm6a_dir, ext_dir,
      lai_als_dopt_path = als_path
    )
    dt_list[[site]] <- build_pixel_dt_sweep(r, site, sweep_variants)
    cli::cli_alert_success(
      "  {site}: {nrow(dt_list[[site]])} pixels loaded"
    )
  }

  data.table::rbindlist(dt_list)
}
