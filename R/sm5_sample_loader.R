# ---
# title:  sm5_sample_loader.R
# desc:   Loads the per-site SM5 stratified-uniform sample (5,000 points per
#         site, h_min=10, LAI-uniform, fCover > 0.9) joined with DSM_sd at
#         each sample position and the four LAI estimators used by the
#         heterogeneity analysis. Used by scripts/steps/20d_* (continuous)
#         and 20e_* (classes on SM5 points).
# ---

# ── build_sm5_site_dt ─────────────────────────────────────────────────────────

#' @title Load the per-site SM5 stratified-uniform sample with DSM_sd + LAI
#'
#' @description
#' For one site, reads (i) the SM5 sampling GPKG with point geometries,
#' (ii) the SM6a DSM_sd raster extracted at sample positions, (iii) the four
#' LAI estimators (LAI_ALS full, LAI_ALS_dopt rescaled to \code{k_select},
#' LAI_S2_ATBD, LAI_S2_opt at \code{column_opt}), and returns a long
#' data.table aligned row-wise with the sample.
#'
#' @param site            Character. Site name ("Aigoual", "Blois", "Mormal").
#' @param d_opt_per_site  Integer. Site-specific optical depth (m).
#' @param column_opt      Character. PROSAIL column name selected at d_opt
#'   (e.g. \code{"LIDFa=2_lai=4_LMA=2_BROWN=1"}).
#' @param k_select        Numeric. Extinction coefficient used in this study
#'   (typically 0.65).
#' @param k_ref           Numeric. Reference k used at PAD precomputation
#'   (typically 0.5). \code{k_scale = k_ref / k_select} rescales PAD LAI.
#' @param h_min           Integer. h_min used at SM5 sampling (typically 10).
#' @param nb_samples      Integer. Number of points per site (typically 5000).
#' @param paths_obj       The \code{paths} list (output/data roots etc.).
#' @param ext_results     Character. Path to external results root (where
#'   \code{ladstack_classic.tif} lives, typically \code{paths$ext_results}).
#' @param dsm_sd_filename Character. DSM_sd raster filename inside SM6a.
#'   Default \code{"dsm_sd_res_10_m.tif"}.
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{site}{character — site name}
#'     \item{dsm_sd}{numeric — DSM_sd extracted at sample position (m)}
#'     \item{lai_als}{numeric — full LAI_ALS rescaled to k_select (m²/m²)}
#'     \item{lai_als_dopt}{numeric — LAI_ALS integrated to d_opt, rescaled
#'                          to k_select (m²/m²)}
#'     \item{lai_s2_atbd}{numeric — S2 ATBD inversion (m²/m²)}
#'     \item{lai_s2_opt}{numeric — S2 OPT inversion at column_opt (m²/m²)}
#'   }
#'
#' @export
build_sm5_site_dt <- function(site, d_opt_per_site, column_opt,
                              k_select, k_ref,
                              h_min            = 10L,
                              nb_samples       = 5000L,
                              paths_obj        = NULL,
                              ext_results      = NULL,
                              dsm_sd_filename  = "dsm_sd_res_10_m.tif") {
  if (is.null(paths_obj))   paths_obj   <- get("paths", envir = .GlobalEnv)
  if (is.null(ext_results)) ext_results <- paths_obj$ext_results

  k_scale   <- k_ref / k_select
  fn_suffix <- sprintf("stratified_uniform_hmin%d_nbSamples_%d", h_min, nb_samples)

  # 1. Sample positions
  gpkg <- file.path(paths_obj$output, "intermediate", "sampling", site,
                    paste0("Sampling_", fn_suffix, ".GPKG"))
  pts <- terra::vect(gpkg)

  # 2. DSM_sd at sample positions
  dsm <- terra::rast(file.path(paths_obj$output, "intermediate", "sm6",
                                site, dsm_sd_filename))
  pts_dsm  <- terra::project(pts, terra::crs(dsm))
  dsm_sd_v <- terra::extract(dsm, pts_dsm)[[2L]]

  # 3. LAI_ATBD (single-column TSV)
  atbd_csv <- file.path(paths_obj$output, "intermediate", "PROSAIL_Models",
                        site, "LIDFa_lai_LMA_BROWN", "atbd",
                        paste0("LAI_estimated_atbd_", fn_suffix, ".csv"))
  atbd_dt     <- data.table::fread(atbd_csv, sep = "\t")
  lai_s2_atbd <- atbd_dt[[1L]]

  # 4. LAI_S2_opt at column_opt (TSV, many columns, pick one)
  common_csv <- file.path(paths_obj$output, "intermediate", "PROSAIL_Models",
                          site, "LIDFa_lai_LMA_BROWN", "common",
                          paste0("LAI_estimated_common_", fn_suffix, ".csv"))
  common_dt  <- data.table::fread(common_csv, sep = "\t")
  if (!(column_opt %in% names(common_dt)))
    stop("Column_opt '", column_opt, "' not in ", common_csv,
         "\nAvailable: ", paste(head(names(common_dt)), collapse = ", "), " …")
  lai_s2_opt <- common_dt[[column_opt]]

  # 5a. LAI_ALS at d_opt (PAD CSV at depth d_opt, rescaled to k_select)
  pad_dopt_csv <- file.path(paths_obj$output, "intermediate", "sampling", site,
                            sprintf("PAD_DSM_Depth_%d_Samples_%s.csv",
                                    d_opt_per_site, fn_suffix))
  pad_dopt_dt  <- data.table::fread(pad_dopt_csv, sep = "\t")
  lai_als_dopt <- pad_dopt_dt$lidar_values * k_scale

  # 5b. Full LAI_ALS — extracted from sum(ladstack_classic.tif). The
  # PAD_Depth_38 CSV is unusable here: only pixels with canopy > 38 m carry
  # a value, so we go back to the raster.
  ladstack_path <- file.path(ext_results, site, "Metrics",
                              "Deciduous_Only", "ladstack_classic.tif")
  ladstack      <- terra::rast(ladstack_path)
  lai_als_rast  <- sum(ladstack, na.rm = TRUE)
  pts_als       <- terra::project(pts, terra::crs(lai_als_rast))
  lai_als       <- terra::extract(lai_als_rast, pts_als)[[2L]] * k_scale

  data.table::data.table(
    site         = site,
    dsm_sd       = dsm_sd_v,
    lai_als      = lai_als,
    lai_als_dopt = lai_als_dopt,
    lai_s2_atbd  = lai_s2_atbd,
    lai_s2_opt   = lai_s2_opt
  )
}
