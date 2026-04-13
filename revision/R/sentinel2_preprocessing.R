# ---
# title:  sentinel2_preprocessing.R
# author: Nathan Corroyez, UMR TETIS
# desc:   Sentinel-2 preprocessing wrappers — download via CDSE and extraction
#         from local SAFE archives.
#         Replaces Main_s2_00_download.R and Main_s2_01_download.R from
#         PROSAIL-Optimization/02_CODES/Sentinel2/ for the RSE revision.
# ---

library("preprocS2")
library("cli")

# ── download_s2_cdse ──────────────────────────────────────────────────────────

#' @title Download Sentinel-2 L2A imagery from Copernicus Data Space (CDSE)
#'
#' @description
#' Downloads a Sentinel-2 Level-2A product for one site and one acquisition
#' date via the Copernicus Data Space Ecosystem (CDSE) API, using
#' \code{preprocS2::get_s2_raster()}. The geometry of acquisition (solar and
#' view angles) is also extracted and saved (\code{geomAcq = TRUE}).
#' All outputs land in \code{output_dir}.
#'
#' This function reproduces exactly the loop body of
#' \code{Main_s2_01_download.R} (lines 31-47) with \code{here::here()} paths
#' and parameterised credentials.
#'
#' @details
#' **Credentials** : CDSE credentials must be stored in \code{~/.Renviron} as
#' \code{CDSE_ID}, \code{CDSE_SECRET}, and \code{CURL_SSL_BACKEND="openssl"}.
#' They are read at runtime via \code{Sys.getenv()} and passed as a list to
#' \code{get_s2_raster(authentication = ...)}. Never read or display the file
#' \code{01_DATA/copernicus-credentials*.yml}.
#'
#' **CRS** : \code{keepCRS = TRUE} preserves the native UTM projection of the
#' SAFE product (EPSG:32631), consistent with the rest of the pipeline.
#'
#' @param site            Character. Site name (e.g. \code{"Aigoual"},
#'   \code{"Blois"}, \code{"Mormal"}).
#' @param date            Character. Acquisition date in \code{"YYYY-MM-DD"}
#'   format.
#' @param aoi_path        Character. Path to the AOI vector file
#'   (\code{utm_init.shp}). Used to clip the download to the study extent.
#' @param output_dir      Character. Output directory for this site (created
#'   if absent). Matches \code{03_RESULTS/{site}/PROSAIL_Optimization} in the
#'   original pipeline.
#' @param credentials_env Character vector of length 2. Names of the
#'   \code{.Renviron} variables holding the CDSE id and secret. Default
#'   \code{c("CDSE_ID", "CDSE_SECRET")}.
#'
#' @return Invisibly \code{NULL}. Side effects: SAFE product, reflectance
#'   stack, cloud mask, and geometry-of-acquisition rasters written to
#'   \code{output_dir} by \code{preprocS2::get_s2_raster()}.
#'
#' @examples
#' \dontrun{
#' download_s2_cdse(
#'   site       = "Blois",
#'   date       = "2021-06-14",
#'   aoi_path   = here::here("01_DATA", "Blois", "Geo_Files", "utm_init.shp"),
#'   output_dir = here::here("03_RESULTS", "Blois", "PROSAIL_Optimization")
#' )
#' }
#'
#' @export
download_s2_cdse <- function(site,
                             date,
                             aoi_path,
                             output_dir,
                             credentials_env = c("CDSE_ID", "CDSE_SECRET")) {
  id     <- Sys.getenv(credentials_env[1])
  secret <- Sys.getenv(credentials_env[2])
  if (nchar(id) == 0L || nchar(secret) == 0L) {
    stop(
      "download_s2_cdse(): CDSE credentials not found in environment. ",
      "Set ", credentials_env[1], " and ", credentials_env[2],
      " in ~/.Renviron and restart R."
    )
  }
  authentication <- list(id = id, pwd = secret)

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  preprocS2::get_s2_raster(
    aoi_path          = aoi_path,
    datetime          = date,
    output_dir        = output_dir,
    path_S2tilinggrid = NULL,
    overwrite         = TRUE,
    geomAcq           = TRUE,
    authentication    = authentication,
    siteName          = site,
    keepCRS           = TRUE
  )
  invisible(NULL)
}


# ── load_s2_from_safe ─────────────────────────────────────────────────────────

#' @title Extract Sentinel-2 reflectances from a local SAFE archive
#'
#' @description
#' Extracts Sentinel-2 reflectances and associated metadata from a Level-2A
#' SAFE directory already present on disk, using
#' \code{preprocS2::extract_from_S2_L2A()}. Returns the full \code{S2obj}
#' list — containing \code{S2_Stack} (stars object) and \code{S2_Bands}
#' (metadata list) — for downstream cloud-mask generation and reflectance
#' writing.
#'
#' This function reproduces exactly the \code{extract_from_S2_L2A()} call in
#' \code{Main_s2_00_download.R} (lines 30-33) for a single site.
#'
#' @details
#' All bands available at the requested \code{resolution} are extracted;
#' 20 m bands are resampled to 10 m internally by
#' \code{extract_from_S2_L2A()}. Only bands B3, B4, B8 are used for PROSAIL
#' inversion downstream, but the full stack is preserved here for
#' completeness.
#'
#' @param safe_path   Character. Full path to the \code{.SAFE} directory
#'   (e.g. \code{S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_*.SAFE}).
#' @param path_vector Character. Path to the AOI vector file
#'   (\code{utm_init.shp}) used to clip the extracted raster to the study
#'   area.
#' @param resolution  Numeric. Target spatial resolution in metres. Default
#'   \code{10}.
#' @param s2source    Character. Source format identifier passed to
#'   \code{extract_from_S2_L2A()}. Default \code{"SAFE"}.
#'
#' @return A named list with at minimum:
#'   \describe{
#'     \item{\code{S2_Stack}}{A \code{stars} object with all extracted bands.}
#'     \item{\code{S2_Bands}}{A list with elements \code{GRANULE},
#'       \code{metadata}, and \code{metadata_MSI} required by
#'       \code{write_reflectance()}.}
#'   }
#'
#' @examples
#' \dontrun{
#' safe_path <- here::here(
#'   "01_DATA", "Blois", "Sentinel-2", "2021-06-14",
#'   "S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_20210614T140120.SAFE"
#' )
#' aoi_path <- here::here("01_DATA", "Blois", "Geo_Files", "utm_init.shp")
#' s2obj <- load_s2_from_safe(safe_path, aoi_path)
#' }
#'
#' @export
load_s2_from_safe <- function(safe_path,
                               path_vector,
                               resolution = 10,
                               s2source   = "SAFE") {
  s2obj <- preprocS2::extract_from_S2_L2A(
    Path_dir_S2 = safe_path,
    path_vector = path_vector,
    S2source    = s2source,
    resolution  = resolution
  )
  return(s2obj)
}


# ── generate_cloud_mask ───────────────────────────────────────────────────────

#' @title Generate and save a binary cloud mask from a Sentinel-2 stack
#'
#' @description
#' Generates a binary cloud mask from the Sentinel-2 stars stack using
#' \code{preprocS2::save_cloud_s2()}, which derives the mask from the Scene
#' Classification Layer (SCL) embedded in the L2A product. Writes the mask
#' files to \code{Cloud_path}.
#'
#' @details
#' Adapted from \code{functions_sentinel_2.R::generate_cloud_mask} (line 116)
#' for the RSE revision refactor. Behaviour is identical to the original;
#' function body is an exact copy.
#'
#' @param S2_Stack  A \code{stars} object. The Sentinel-2 stack returned as
#'   \code{S2obj$S2_Stack} by \code{load_s2_from_safe()}.
#' @param Cloud_path Character. Directory where cloud mask files will be
#'   written. Must exist or be created before calling this function.
#' @param S2source   Character. Source format identifier (\code{"SAFE"} or
#'   \code{"L2A"}). Default \code{"SAFE"}.
#' @param saveRaw    Logical. Whether to also save the raw SCL raster alongside
#'   the binary mask. Default \code{TRUE}.
#'
#' @return The object returned by \code{preprocS2::save_cloud_s2()}: a list
#'   with paths to the generated mask files (used downstream if the mask path
#'   is needed for later steps).
#'
#' @examples
#' \dontrun{
#' s2obj      <- load_s2_from_safe(safe_path, aoi_path)
#' cloud_path <- here::here("03_RESULTS", "Blois", "PROSAIL_Optimization",
#'                          "granule_name", "CloudMask")
#' dir.create(cloud_path, recursive = TRUE, showWarnings = FALSE)
#' cloudmasks <- generate_cloud_mask(s2obj$S2_Stack, cloud_path)
#' }
#'
#' @export
generate_cloud_mask <- function(S2_Stack, Cloud_path, S2source = "SAFE",
                                saveRaw = TRUE) {
  # Write cloud mask and get filenames
  cloudmasks <- preprocS2::save_cloud_s2(
    S2_stars   = S2_Stack,
    Cloud_path = Cloud_path,
    S2source   = S2source,
    SaveRaw    = saveRaw
  )
  cat("Cloud mask has been successfully created and saved at :",
      Cloud_path, "\n")

  return(cloudmasks)
}


# ── write_reflectance ─────────────────────────────────────────────────────────

#' @title Write Sentinel-2 reflectances to ENVI format
#'
#' @description
#' Writes the Sentinel-2 reflectance stack to disk in ENVI BIL Int16 format
#' using \code{preprocS2::save_reflectance_s2()}. The companion metadata
#' files produced alongside the ENVI image encode the offset and gain applied
#' to the S2 L2A DN values; they are mandatory for correct interpretation
#' by the downstream PROSAIL inversion code.
#'
#' @details
#' Adapted from \code{functions_sentinel_2.R::write_reflectance} (line 140)
#' for the RSE revision refactor. Behaviour is identical to the original;
#' function body is an exact copy.
#'
#' @param S2_Stack  A \code{stars} object. The Sentinel-2 stack
#'   (\code{S2obj$S2_Stack}).
#' @param S2_Bands  A list. The \code{S2_Bands} element of the \code{S2obj}
#'   returned by \code{load_s2_from_safe()}. Must contain \code{GRANULE},
#'   \code{metadata}, and \code{metadata_MSI}.
#' @param refl_path Character. Full output path without extension where the
#'   ENVI image and companion files (\code{.hdr}, \code{.bil}, metadata) will
#'   be written.
#'
#' @return Invisibly \code{NULL}. Side effect: ENVI image and companion
#'   metadata files written to \code{refl_path}.
#'
#' @examples
#' \dontrun{
#' s2obj     <- load_s2_from_safe(safe_path, aoi_path)
#' refl_path <- here::here(
#'   "03_RESULTS", "Blois", "PROSAIL_Optimization",
#'   "granule_name", "Reflectance", "granule_name_Refl"
#' )
#' write_reflectance(s2obj$S2_Stack, s2obj$S2_Bands, refl_path)
#' }
#'
#' @export
write_reflectance <- function(S2_Stack, S2_Bands, refl_path) {
  # Save Reflectance file as ENVI image with BIL interleaves
  # metadata files are important to account for offset applied on S2 L2A products
  tile_S2    <- preprocS2::get_tile(S2_Bands$GRANULE)
  dateAcq_S2 <- preprocS2::get_dateAcq(
    S2product = dirname(dirname(S2_Bands$GRANULE))
  )

  preprocS2::save_reflectance_s2(
    S2_stars   = S2_Stack,
    Refl_path  = refl_path,
    S2Sat      = NULL,
    tile_S2    = tile_S2,
    dateAcq_S2 = dateAcq_S2,
    Format     = "ENVI",
    datatype   = "Int16",
    MTD        = S2_Bands$metadata,
    MTD_MSI    = S2_Bands$metadata_MSI
  )
  cat("Reflectance has been successfully created and saved at :",
      refl_path, "\n")
}
