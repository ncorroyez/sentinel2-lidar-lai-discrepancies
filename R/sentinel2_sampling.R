# ---
# title:  sentinel2_sampling.R
# author: Nathan Corroyez, UMR TETIS
# desc:   Sentinel-2 pixel sampling and LiDAR PAD extraction at sample
#         locations. Replaces get_s2_samples.R (libraries/) and the
#         orchestration scripts Main_s2_02A_extract.R /
#         Main_s2_02B_sample_lidar.R from
#         PROSAIL-Optimization/02_CODES/ for the RSE revision.
#
#         Five public functions, in dependency order:
#           align_and_remove_na_for_aoi()      — clip AOI to valid S2+LiDAR pixels
#           stratified_sampling_uniform()       — LAI-stratified uniform sampling
#           stratified_sampling()               — quantile-stratified sampling (dead code)
#           get_s2_samples()                    — dispatcher: sample S2 reflectances
#           extract_pad_at_samples()            — extract PAD values at sample positions
# ---


# ── align_and_remove_na_for_aoi ───────────────────────────────────────────────

#' @title Clip an AOI vector to pixels valid in both S2 and LiDAR rasters
#'
#' @description
#' Rasterises the study-area AOI polygon onto the Sentinel-2 grid, masks it
#' against the LiDAR LAI raster (removing cells with NA LiDAR coverage), and
#' converts the result back to a polygon. The clipped AOI is optionally written
#' to disk and its path returned, or the object itself is returned if
#' \code{save = FALSE}.
#'
#' This function is a faithful copy of \code{align_and_remove_na_for_aoi()} in
#' \code{get_s2_samples.R}
#' (\code{PROSAIL-Optimization/02_CODES/libraries/}).
#'
#' @details
#' **Path redirection**: the original function received a hard-coded
#' \code{sampling_path} that pointed to
#' \code{01_DATA/{site}/Geo_Files/aoi_sampling_prep.gpkg}, writing into
#' \code{01_DATA/} — a violation of the CLAUDE.md read-only rule for that
#' folder. This refactor replaces the fixed path with an explicit
#' \code{output_path} argument. The caller (orchestration script) supplies a
#' path inside \code{03_RESULTS/} or \code{output/}. When
#' \code{save = FALSE}, nothing is written and \code{output_path} is ignored.
#'
#' **CRS assumption**: \code{terra::project(lidar_r, s2_r)} reprojects the
#' LiDAR raster to the S2 CRS before masking. A hard stop is raised if CRS
#' diverges after reprojection or if extents do not overlap.
#'
#' @param aoi_path        Character. Path to the AOI vector file
#'   (\code{utm_init.shp} or equivalent).
#' @param S2_path         Character. Path to the Sentinel-2 reflectance stack
#'   (ENVI format, no extension; produced by \code{write_reflectance()}).
#' @param lidar_lai_path  Character. Path to the LiDAR LAI raster used to mask
#'   NA pixels (e.g. \code{ladstack.tif}).
#' @param output_path     Character or \code{NULL}. Full path (including
#'   filename and \code{.gpkg} extension) where the clipped AOI polygon will be
#'   written when \code{save = TRUE}. Ignored when \code{save = FALSE}.
#' @param save            Logical. If \code{TRUE} (default), writes the clipped
#'   AOI as a GeoPackage to \code{output_path} and returns \code{output_path}.
#'   If \code{FALSE}, returns the \code{SpatVector} object directly.
#'
#' @return If \code{save = TRUE}: character, the path of the written GeoPackage
#'   (\code{output_path}). If \code{save = FALSE}: a \code{SpatVector} polygon.
#'
#' @examples
#' \dontrun{
#' aoi_sampling_path <- align_and_remove_na_for_aoi(
#'   aoi_path       = here::here("01_DATA","Blois","Geo_Files","utm_init.shp"),
#'   S2_path        = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                               "L2A_T31TCN_A031222_20210614T105443",
#'                               "Reflectance","L2A_T31TCN_A031222_20210614T105443_Refl"),
#'   lidar_lai_path = here::here("01_DATA","Blois","LiDAR",
#'                               "PAD_Profiles_Classic","ladstack.tif"),
#'   output_path    = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                               "sampling","aoi_sampling_prep.gpkg"),
#'   save           = TRUE
#' )
#' }
#'
#' @export
align_and_remove_na_for_aoi <- function(aoi_path, S2_path, lidar_lai_path,
                                         output_path = NULL, save = TRUE) {
  # Read rasters
  s2_r    <- terra::rast(S2_path)
  lidar_r <- terra::rast(lidar_lai_path)
  lidar_r <- terra::project(lidar_r, s2_r)   # ensure CRS alignment
  aoi     <- terra::vect(aoi_path)
  aoi_r   <- terra::rasterize(aoi, s2_r)

  # Spatial extent and CRS assessment
  if (!terra::crs(s2_r) == terra::crs(lidar_r)) {
    stop("CRS mismatch detected among inputs.")
  }
  if (is.null(terra::intersect(terra::ext(s2_r),
                                terra::ext(lidar_r)))) {
    stop("Spatial extents do not overlap among inputs.")
  }

  # Mask and filter NA values
  aoi_masked <- terra::mask(aoi_r, lidar_r)
  aoi_masked <- terra::as.polygons(aoi_masked, na.rm = TRUE)

  # Write the result to file
  if (save) {
    terra::writeVector(x = aoi_masked, filename = output_path, overwrite = TRUE)
  }
  return(output_path)
}


# ── stratified_sampling_uniform ───────────────────────────────────────────────

#' @title Uniform LAI-stratified sampling of Sentinel-2 pixels
#'
#' @description
#' Samples \code{nbSamples} pixel locations from the study area using a
#' uniform-stratified design based on LiDAR-derived LAI. Pixels are first
#' filtered to keep only those with canopy height in (10, 40] m (using
#' \code{max_res_10_m.tif}), then binned along the LAI range
#' [2, 98th-percentile] and drawn proportionally across bins. The approach
#' ensures that the sample covers the full range of observed LAI values rather
#' than over-representing the mode.
#'
#' This function is a copy of \code{stratified_sampling_uniform()} in
#' \code{get_s2_samples.R}
#' (\code{PROSAIL-Optimization/02_CODES/libraries/}) with two explicit
#' corrections.
#'
#' @details
#' **Correction 1 — \code{dplyr::} prefixes**: the original function called
#' \code{filter()}, \code{mutate()}, \code{group_by()}, \code{sample_n()},
#' \code{ungroup()}, \code{slice_sample()}, and \code{cur_group_id()} without
#' namespace prefix, relying on \code{library(dplyr)} being in scope. This
#' refactor adds explicit \code{dplyr::} prefixes throughout.
#'
#' **Correction 2 — \code{nbSamples} argument restored**: the original
#' hardcoded \code{5000} at two locations (computing \code{base_samples} and
#' \code{extra_samples}), silently ignoring the \code{nbSamples} argument.
#' This refactor replaces both occurrences with \code{nbSamples}. For
#' \code{nbSamples = 5000} (the production value used in
#' Main_s2_02A_extract.R), the output is numerically identical to the original.
#'
#' **Height filter**: pixels with \code{max_height <= 10} or
#' \code{max_height > 40} are excluded before sampling. \code{lai_5th} is
#' fixed at 2 (matching the published analysis); \code{lai_95th} is computed
#' as the 98th percentile of filtered LAI values.
#'
#' @param site           Character. Site name (\code{"Aigoual"}, \code{"Blois"},
#'   or \code{"Mormal"}). Passed through to internal logic; currently used only
#'   in commented-out site-specific breaks.
#' @param lidar_lai_path Character. Path to the LiDAR LAI stack raster
#'   (\code{ladstack.tif}). All layers are summed (\code{na.rm = TRUE}) to
#'   obtain total LAI per pixel.
#' @param max_path       Character. Path to the maximum canopy height raster
#'   (\code{max_res_10_m.tif}). Used to filter pixels to the
#'   \code{(10, 40]} m height range.
#' @param nbSamples      Integer. Number of sample points to draw. Default 5000
#'   (the production value). Must be ≤ total number of eligible pixels.
#'
#' @return A \code{SpatVector} of point locations (pixel centres) in the CRS
#'   of \code{lidar_lai_path}, with columns \code{x}, \code{y},
#'   \code{lidar_lai}, \code{max_height}, and \code{lai_bin}.
#'
#' @examples
#' \dontrun{
#' samples <- stratified_sampling_uniform(
#'   site           = "Blois",
#'   lidar_lai_path = here::here("01_DATA","Blois","LiDAR",
#'                               "PAD_Profiles_Classic","ladstack.tif"),
#'   max_path       = here::here("01_DATA","Blois","LiDAR","max_res_10_m.tif"),
#'   nbSamples      = 5000
#' )
#' }
#'
#' @export
stratified_sampling_uniform <- function(site,
                                         lidar_lai_path,
                                         max_path,
                                         fcover_path,
                                         nbSamples,
                                         h_min = 10L) {
  # Open raster
  mraster_r    <- terra::rast(lidar_lai_path)
  mraster      <- sum(mraster_r, na.rm = TRUE)
  max_rast     <- terra::rast(max_path)
  fcover_rast  <- terra::rast(fcover_path)

  # Stack aligned rasters for filtering
  combined_stack        <- c(mraster, max_rast, fcover_rast)
  names(combined_stack) <- c("lidar_lai", "max_height", "fCover")

  # Filter: max_height in [h_min, 40] m AND fCover > 0.9
  data_df     <- as.data.frame(combined_stack, xy = TRUE, na.rm = TRUE)
  filtered_df <- dplyr::filter(data_df,
                                max_height >= h_min & max_height <= 40 &
                                fCover > 0.9)

  # Define breaks within the 5th-95th percentile range
  lai_5th  <- 2
  lai_95th <- floor(quantile(filtered_df$lidar_lai, probs = 0.98, na.rm = TRUE))
  breaks   <- seq(lai_5th, lai_95th, length.out = lai_95th - lai_5th + 1)

  # Compute base samples per bin
  # NOTE: original hardcoded 5000 here (twice) — see @details for correction.
  num_bins     <- length(breaks) - 1
  base_samples <- floor(nbSamples / num_bins)
  extra_samples <- nbSamples %% num_bins

  # Sample
  sampled_points <- dplyr::filter(filtered_df,
                                   lidar_lai >= lai_5th & lidar_lai <= lai_95th) |>
    dplyr::mutate(lai_bin = cut(lidar_lai, breaks = breaks,
                                 include.lowest = TRUE)) |>
    dplyr::group_by(lai_bin) |>
    dplyr::sample_n(size = base_samples + (dplyr::cur_group_id() <= extra_samples),
                    replace = FALSE) |>
    dplyr::ungroup() |>
    dplyr::slice_sample(n = nbSamples)

  samples <- terra::vect(sampled_points, geom = c("x", "y"),
                          crs = terra::crs(mraster_r))
  print(hist(samples$lidar_lai))
  return(samples)
}


# ── stratified_sampling ───────────────────────────────────────────────────────

#' @title Quantile-stratified sampling using mean and L-skewness LiDAR metrics
#'
#' @description
#' Samples \code{nbSamples} pixel locations using a stratified design based on
#' two LiDAR structural metrics: canopy height mean and L-skewness. Strata are
#' defined by quantile-based breaks (Q25, Q50, Q75) applied to field-point
#' distributions of each metric, then crossed into a composite strata raster.
#' Points are drawn proportionally from each stratum via
#' \code{sgsR::sample_strat()}.
#'
#' This function is a faithful copy of \code{stratified_sampling()} in
#' \code{get_s2_samples.R}
#' (\code{PROSAIL-Optimization/02_CODES/libraries/}).
#'
#' @details
#' **DEAD CODE in production**: this sampling method is never invoked by the
#' production pipeline, which uses exclusively
#' \code{stratified_sampling_uniform()} (\code{method = 'stratified_uniform'}).
#' Preserved here for potential future sensitivity analyses (e.g. reviewer
#' requests to test sampling method effects on inversion results). Not validated
#' by the SM3 end-to-end equivalence check.
#'
#' @param mean_path         Character. Path to the raster of mean canopy height
#'   (\code{mean_res_10_m.tif}).
#' @param lskew_path        Character. Path to the raster of canopy height
#'   L-skewness (\code{lskew_res_10_m.tif}).
#' @param field_points_path Character. Path to the field point vector
#'   (\code{data_utm31n.geojson}) used to compute metric quantiles.
#' @param nbSamples         Integer. Number of sample points to draw.
#'
#' @return A \code{SpatVector} of sampled point locations.
#'
#' @examples
#' \dontrun{
#' samples <- stratified_sampling(
#'   mean_path         = here::here("01_DATA","Blois","LiDAR","mean_res_10_m.tif"),
#'   lskew_path        = here::here("01_DATA","Blois","LiDAR","lskew_res_10_m.tif"),
#'   field_points_path = here::here("01_DATA","Blois","Geo_Files","data_utm31n.geojson"),
#'   nbSamples         = 5000
#' )
#' }
#'
#' @export
stratified_sampling <- function(mean_path, lskew_path, field_points_path,
                                 nbSamples) {
  # Open rasters and vector
  files   <- c(mean_path, lskew_path)
  mraster <- terra::rast(files)
  names(mraster) <- c("mean", "lskew")
  field_points  <- terra::vect(field_points_path)
  field_values  <- terra::extract(mraster, field_points, bind = TRUE)

  # Spatial extent and CRS assessment
  if (!terra::crs(mraster[[1]]) == terra::crs(mraster[[2]]) ||
      !terra::crs(mraster[[1]]) == terra::crs(field_points)) {
    stop("CRS mismatch detected among inputs.")
  }
  if (is.null(terra::intersect(terra::ext(mraster[[1]]),
                                terra::ext(mraster[[2]]))) ||
      is.null(terra::intersect(terra::ext(mraster[[1]]),
                                terra::ext(field_points)))) {
    stop("Spatial extents do not overlap among inputs.")
  }

  # Calculate quantiles for mean and lskew
  br.mean  <- quantile(field_values$mean,  probs = c(0.25, 0.50, 0.75),
                        na.rm = TRUE)
  br.lskew <- quantile(field_values$lskew, probs = c(0.25, 0.50, 0.75),
                        na.rm = TRUE)

  # Create stratified raster using quantile-based breaks
  sraster <- sgsR::strat_breaks(mraster,
                                 breaks  = list(br.mean, br.lskew),
                                 plot    = FALSE,
                                 details = FALSE)

  # Sample points stratified by the composite strata
  composite_strata <- sraster[[1]] * 10 + sraster[[2]]
  sampled <- sgsR::sample_strat(sraster  = composite_strata,
                                 nSamp    = nbSamples,
                                 mindist  = NULL,
                                 force    = TRUE,
                                 plot     = FALSE)
  samples <- terra::vect(sampled$geometry)
  return(samples)
}


# ── get_s2_samples ────────────────────────────────────────────────────────────

#' @title Sample Sentinel-2 reflectances at stratified pixel locations
#'
#' @description
#' Dispatcher function that applies one of three spatial sampling strategies
#' (\code{'random'}, \code{'stratified'}, or \code{'stratified_uniform'}) to
#' select \code{nbSamples} pixel locations within the study-area AOI, then
#' extracts Sentinel-2 bidirectional reflectance values at those locations.
#' Returns both the sample positions (as a \code{SpatVector}) and the extracted
#' reflectances (as a \code{data.frame}).
#'
#' This function is a faithful copy of \code{get_s2_samples()} in
#' \code{get_s2_samples.R}
#' (\code{PROSAIL-Optimization/02_CODES/libraries/}) with one modification.
#'
#' @details
#' **Removed \code{above20} argument**: the original signature included an
#' \code{above20} argument (declared but never read in the function body). It
#' has been removed in this refactor. Callers that were passing \code{above20}
#' explicitly will need to update their calls, but no production script passes
#' this argument (verified by grep across the legacy pipeline).
#'
#' **Production method**: only \code{method = 'stratified_uniform'} is used in
#' the production pipeline (Main_s2_02A_extract.R, line 26). The \code{'random'}
#' and \code{'stratified'} branches are preserved for potential future use.
#'
#' **CRS handling**: if the S2 raster CRS differs from the sample vector CRS,
#' the raster is reprojected before extraction (matching original behaviour).
#'
#' @param aoi_path          Character. Path to the (pre-clipped) AOI vector
#'   file produced by \code{align_and_remove_na_for_aoi()}
#'   (\code{aoi_sampling_prep.gpkg}).
#' @param S2_path           Character. Path to the Sentinel-2 reflectance
#'   stack (ENVI format, no extension).
#' @param mean_path         Character. Path to the mean canopy height raster.
#'   Required for \code{method = 'stratified'}; ignored otherwise.
#' @param lskew_path        Character. Path to the L-skewness canopy height
#'   raster. Required for \code{method = 'stratified'}; ignored otherwise.
#' @param lidar_lai_path    Character. Path to the LiDAR LAI stack raster
#'   (\code{ladstack.tif}). Required for \code{method = 'stratified_uniform'}.
#' @param max_path          Character. Path to the maximum canopy height raster.
#'   Required for \code{method = 'stratified_uniform'}.
#' @param field_points_path Character. Path to the field point vector. Required
#'   for \code{method = 'stratified'}; ignored otherwise.
#' @param nbSamples         Integer. Number of pixel locations to sample.
#' @param site              Character. Site name (\code{"Aigoual"},
#'   \code{"Blois"}, or \code{"Mormal"}). Passed to
#'   \code{stratified_sampling_uniform()}.
#' @param method            Character. Sampling method: \code{'random'},
#'   \code{'stratified'}, or \code{'stratified_uniform'} (default
#'   \code{'random'}).
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{sample_location}}{A \code{SpatVector} of sampled point
#'       locations.}
#'     \item{\code{S2_refl}}{A \code{data.frame} of extracted Sentinel-2
#'       reflectances (one row per sample, one column per band).}
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' result <- get_s2_samples(
#'   aoi_path          = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                                  "sampling","aoi_sampling_prep.gpkg"),
#'   S2_path           = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                                  "L2A_T31TCN_A031222_20210614T105443",
#'                                  "Reflectance",
#'                                  "L2A_T31TCN_A031222_20210614T105443_Refl"),
#'   mean_path         = here::here("01_DATA","Blois","LiDAR","mean_res_10_m.tif"),
#'   lskew_path        = here::here("01_DATA","Blois","LiDAR","lskew_res_10_m.tif"),
#'   lidar_lai_path    = here::here("01_DATA","Blois","LiDAR",
#'                                  "PAD_Profiles_Classic","ladstack.tif"),
#'   max_path          = here::here("01_DATA","Blois","LiDAR","max_res_10_m.tif"),
#'   field_points_path = here::here("01_DATA","Blois","Geo_Files",
#'                                  "data_utm31n.geojson"),
#'   nbSamples         = 5000,
#'   site              = "Blois",
#'   method            = "stratified_uniform"
#' )
#' }
#'
#' @export
get_s2_samples <- function(aoi_path, S2_path, mean_path,
                            lskew_path, lidar_lai_path, max_path,
                            fcover_path,
                            field_points_path,
                            nbSamples, site, method = 'random',
                            h_min = 10L) {

  # read vector and raster
  aoi_init <- terra::vect(aoi_path)
  S2_rast  <- terra::rast(S2_path)

  # check if CRS is the same otherwise reproject
  if (terra::crs(aoi_init) == terra::crs(S2_rast)) {
    aoi <- aoi_init
  } else if (!terra::crs(aoi_init) == terra::crs(S2_rast)) {
    aoi <- terra::project(x = aoi_init, y = terra::crs(S2_rast))
  }

  if (method == 'random') {
    samples <- terra::spatSample(x = aoi, size = nbSamples, method = method)
    # Snap samples to pixel centers
    coords     <- terra::geom(samples)[, c("x", "y")]
    cells      <- terra::cellFromXY(S2_rast, coords)
    centers    <- terra::xyFromCell(S2_rast, cells)
    centers_df <- data.frame(x = centers[, 1], y = centers[, 2], layer = 1)

    # Create a SpatVector from the centers_df, specifying the CRS
    samples <- terra::vect(centers_df, geom = c("x", "y"),
                            crs = terra::crs(S2_rast))
  } else if (method == 'stratified') {
    # Automatically at the pixel centers
    samples <- stratified_sampling(mean_path         = mean_path,
                                    lskew_path        = lskew_path,
                                    field_points_path = field_points_path,
                                    nbSamples         = nbSamples)
  } else if (method == 'stratified_uniform') {
    samples <- stratified_sampling_uniform(site           = site,
                                            lidar_lai_path = lidar_lai_path,
                                            max_path       = max_path,
                                            fcover_path    = fcover_path,
                                            nbSamples      = nbSamples,
                                            h_min          = h_min)
  } else {
    stop("Error: method should be either 'random', 'stratified' or 'stratified_uniform.")
  }

  # Extract reflectance from raster
  crs_raster <- terra::crs(S2_rast)
  crs_vector <- terra::crs(samples)
  if (crs_raster != crs_vector) {
    S2_rast <- terra::project(S2_rast, crs_vector)
  }
  S2_refl <- terra::extract(x = S2_rast, y = samples)

  return(list('sample_location' = samples, 'S2_refl' = S2_refl))
}


# ── extract_pad_at_samples ────────────────────────────────────────────────────

#' @title Extract PAD profile values at pre-defined sample locations
#'
#' @description
#' Reads a single PAD raster file and a GeoPackage of sample locations
#' (produced by \code{get_s2_samples()}), extracts the raster value at each
#' point, and writes the result as a TSV CSV file. This is an atomic operation:
#' one PAD raster → one CSV. The orchestration script
#' (\code{03b_extract_lidar_at_samples.R}) handles the outer loop over sites,
#' normalisation types, and depth layers.
#'
#' This function is a refactor of the per-file extraction logic in
#' \code{Main_s2_02B_sample_lidar.R}
#' (\code{PROSAIL-Optimization/02_CODES/Sentinel2/}).
#'
#' @details
#' **Sample ID preservation**: the CSV written contains two columns:
#' \code{lidar_values} (extracted PAD values) and \code{samples_id} (integer
#' index of each sample point, from \code{sample_locations$sample_id}).
#' This matches the column structure of the CSVs written by
#' \code{Main_s2_02B_sample_lidar.R} and expected by the d_opt analysis
#' (\code{Main_s2_04C_analysis_depth.R}).
#'
#' @param sample_gpkg_path Character. Full path to the GeoPackage of sample
#'   locations (output of \code{03a_sample_s2_pixels.R}). Must contain a
#'   \code{sample_id} field.
#' @param pad_rast_path    Character. Full path to the PAD raster file for a
#'   single depth layer (e.g.
#'   \code{01_DATA/{site}/LiDAR/testPADs/PAD_Profiles_DSM_keepTrees/PAD_2.5_40.tif}).
#' @param output_csv_path  Character. Full path (including filename and
#'   \code{.csv} extension) of the output TSV file.
#'
#' @return Invisibly \code{NULL}. Side effect: a TSV CSV file written to
#'   \code{output_csv_path} with columns \code{lidar_values} and
#'   \code{samples_id}.
#'
#' @examples
#' \dontrun{
#' extract_pad_at_samples(
#'   sample_gpkg_path = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                                  "sampling",
#'                                  "Sampling_stratified_uniform_nbSamples_5000.GPKG"),
#'   pad_rast_path    = here::here("01_DATA","Blois","LiDAR","testPADs",
#'                                  "PAD_Profiles_DSM_keepTrees","PAD_2.5_40.tif"),
#'   output_csv_path  = here::here("03_RESULTS","Blois","PROSAIL_Optimization",
#'                                  "sampling",
#'                                  "PAD_DSM_keepTrees_Depth_1_Samples_stratified_uniform_nbSamples_5000.csv")
#' )
#' }
#'
#' @export
extract_pad_at_samples <- function(sample_gpkg_path, pad_rast_path,
                                    output_csv_path) {
  sample_locations <- terra::vect(sample_gpkg_path)
  samples_id       <- sample_locations$sample_id
  lidar_values     <- unlist(
    terra::extract(terra::rast(pad_rast_path), sample_locations)[-1]
  )
  readr::write_delim(
    x     = as.data.frame(cbind(lidar_values, samples_id)),
    file  = output_csv_path,
    delim = '\t'
  )
  invisible(NULL)
}
