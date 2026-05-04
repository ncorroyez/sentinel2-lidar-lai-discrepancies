# ---
# title:  prosail_inversion.R
# desc:   Helpers for applying trained SVR ensembles to Sentinel-2 reflectances
#         as part of the PROSAIL hybrid inversion pipeline (SM4).
#
#         Two thin wrappers that isolate the liquidSVM deserialisation and the
#         prosail::prosail_hybrid_apply call from the orchestration script.
#         They reproduce verbatim the corresponding lines of
#         Main_s2_03_PROSAIL_Inversion_bis.R (lines 174-179).
#         Updated for prosail 3.0.0 (2026-04): PROSAIL_Hybrid_Apply →
#         prosail_hybrid_apply; return fields MeanEstimate/StdEstimate →
#         mean_estimate/sd_estimate.
# ---

# ── load_svr_ensemble ──────────────────────────────────────────────────────────

#' @title Load and deserialise a serialised SVR ensemble
#'
#' @description
#' Reads an RDS file produced by \code{revision/R/prosail_lut.R}
#' \code{train_all_simulations()} and deserialises each liquidSVM model it
#' contains, returning a list ready for \code{prosail::prosail_hybrid_apply}.
#'
#' @details
#' This function reproduces verbatim
#' \code{Main_s2_03_PROSAIL_Inversion_bis.R} line 174:
#' \preformatted{
#'   modelSVR <- lapply(readRDS(file = filename),
#'                      liquidSVM::unserialize.liquidSVM)
#' }
#' The serialisation round-trip (via \code{liquidSVM::serialize.liquidSVM}
#' during training and \code{liquidSVM::unserialize.liquidSVM} here) is
#' required because liquidSVM objects cannot be written to disk via base R
#' \code{saveRDS()} without explicit serialisation.
#'
#' @param rds_path Character scalar. Full path to a \code{.rds} file produced
#'   by \code{train_all_simulations()} containing a serialised liquidSVM
#'   ensemble.
#'
#' @return A list of deserialised liquidSVM model objects, suitable as the
#'   \code{RegressionModels} argument of
#'   \code{prosail::prosail_hybrid_apply()}.
#'
#' @examples
#' \dontrun{
#'   svr <- load_svr_ensemble(
#'     here::here("revision", "output", "intermediate",
#'                "PROSAIL_Models", "Blois",
#'                "LIDFa_lai_LMA_BROWN",
#'                "LIDFa=40_lai=2_LMA=0.01_BROWN=0.rds")
#'   )
#' }
load_svr_ensemble <- function(rds_path) {
  lapply(readRDS(file = rds_path), liquidSVM::unserialize.liquidSVM)
}

# ── apply_svr_to_reflectance ───────────────────────────────────────────────────

#' @title Apply an SVR ensemble to a matrix of Sentinel-2 reflectances
#'
#' @description
#' Thin wrapper around \code{prosail::prosail_hybrid_apply}. Takes a
#' deserialised SVR ensemble (produced by \code{load_svr_ensemble()}) and a
#' reflectance matrix sub-selected to the training bands, and returns the raw
#' inversion object containing both the mean LAI estimate and its standard
#' deviation across ensemble members.
#'
#' @details
#' This function reproduces verbatim
#' \code{Main_s2_03_PROSAIL_Inversion_bis.R} line 176:
#' \preformatted{
#'   LAIest <- prosail_hybrid_apply(regression_models = modelSVR,
#'                                  Refl = Refl[, selbands])
#' }
#' The caller is responsible for sub-selecting the reflectance columns to the
#' training bands (via \code{selbands}) before passing
#' \code{reflectance_matrix}. The function does not perform any column
#' filtering internally.
#'
#' The raw return value is intentionally preserved (not split into two
#' vectors) so that the caller can access \code{$mean_estimate} and
#' \code{$sd_estimate} explicitly, mirroring the original script.
#'
#' @param svr_ensemble List. A deserialised SVR ensemble as returned by
#'   \code{load_svr_ensemble()}.
#' @param reflectance_matrix A matrix or data.frame of Sentinel-2 reflectance
#'   values (pixels × bands), already scaled to [0, 1] (i.e. divided by
#'   10 000) and sub-selected to the bands used during training. Typically
#'   \code{Refl[, selbands]} where \code{Refl} is the full reflectance table
#'   read from a \code{S2_reflectance_*.csv} file.
#'
#' @return The raw object returned by \code{prosail::prosail_hybrid_apply},
#'   a list with at least two components:
#'   \describe{
#'     \item{\code{mean_estimate}}{Numeric vector of length \code{nrow(Refl)}.
#'       Mean LAI estimate across ensemble members (m2/m2).}
#'     \item{\code{sd_estimate}}{Numeric vector of length \code{nrow(Refl)}.
#'       Standard deviation of LAI estimates across ensemble members.}
#'   }
#'
#' @examples
#' \dontrun{
#'   svr    <- load_svr_ensemble("path/to/combo.rds")
#'   LAIest <- apply_svr_to_reflectance(
#'     svr_ensemble       = svr,
#'     reflectance_matrix = Refl[, selbands]
#'   )
#'   lai_mean <- LAIest$mean_estimate
#'   lai_sd   <- LAIest$sd_estimate
#' }
apply_svr_to_reflectance <- function(svr_ensemble, reflectance_matrix) {
  prosail::prosail_hybrid_apply(
    regression_models = svr_ensemble,
    refl             = reflectance_matrix
  )
}

# ── predict_lai_full_raster ────────────────────────────────────────────────────

#' @title Apply an SVR ensemble to a full Sentinel-2 reflectance raster
#'
#' @description
#' Loads band-selected layers from a multi-band S2 SpatRaster, scales DN values
#' to [0, 1], applies \code{prosail_hybrid_apply()} to every non-NA pixel, and
#' returns a single-layer SpatRaster of mean LAI estimates. Optionally masks the
#' output with a binary mask raster (1 = keep, NA = exclude).
#'
#' Bands are selected by matching layer names against \code{band_prefixes} using
#' \code{startsWith()}. For the standard S2 preprocS2 output (10 bands at 10 m)
#' the default \code{c("B03", "B04", "B08")} selects layers 2, 3 and 7.
#'
#' @param svr_ensemble     List. Deserialised SVR ensemble from
#'   \code{load_svr_ensemble()}.
#' @param s2_raster        SpatRaster. Multi-band S2 reflectance raster in DN
#'   (values approximately in [0, 10000]).
#' @param band_prefixes    Character vector. Layer name prefixes to select.
#'   Default \code{c("B03", "B04", "B08")}.
#' @param mask_rast        SpatRaster or \code{NULL}. Binary mask (1 = keep).
#'   Applied after prediction via \code{terra::mask()}. Default \code{NULL}.
#' @param dn_scale         Numeric. Divisor to convert DN → reflectance [0, 1].
#'   Default \code{10000}.
#'
#' @return A single-layer SpatRaster with mean LAI estimates (m²/m²). NA pixels
#'   correspond to original NA pixels or masked-out pixels.
#'
#' @export
predict_lai_full_raster <- function(svr_ensemble, s2_raster,
                                    band_prefixes = c("B03", "B04", "B08"),
                                    mask_rast  = NULL,
                                    dn_scale   = 10000,
                                    chunk_size = 100000L) {
  # Select bands by prefix
  band_idx <- which(Reduce(`|`, lapply(band_prefixes, startsWith, x = names(s2_raster))))
  if (length(band_idx) != length(band_prefixes)) {
    stop(
      "predict_lai_full_raster: expected ", length(band_prefixes),
      " bands matching prefixes c(", paste(band_prefixes, collapse = ", "),
      "), found ", length(band_idx), ".\n",
      "Available layer names: ", paste(names(s2_raster), collapse = ", ")
    )
  }
  s2_sel <- s2_raster[[band_idx]]

  # Extract all pixel values: matrix (npix × nbands)
  vals  <- terra::values(s2_sel) / dn_scale
  valid <- stats::complete.cases(vals)

  preds <- rep(NA_real_, nrow(vals))
  if (any(valid)) {
    valid_idx <- which(valid)
    n_valid   <- length(valid_idx)
    n_chunks  <- ceiling(n_valid / chunk_size)
    for (chunk in seq_len(n_chunks)) {
      start <- (chunk - 1L) * chunk_size + 1L
      end   <- min(chunk * chunk_size, n_valid)
      idx   <- valid_idx[start:end]
      result <- prosail::prosail_hybrid_apply(
        regression_models = svr_ensemble,
        refl              = as.matrix(vals[idx, , drop = FALSE])
      )
      preds[idx] <- result$mean_estimate
    }
  }

  # Write predictions back into a single-layer raster template
  out_rast        <- s2_raster[[1L]]
  names(out_rast) <- "lai_s2_opt"
  terra::values(out_rast) <- preds

  if (!is.null(mask_rast)) {
    out_rast <- terra::mask(out_rast, mask_rast)
  }
  out_rast
}
