# ---
# title:  prosail_inversion.R
# desc:   Helpers for applying trained SVR ensembles to Sentinel-2 reflectances
#         as part of the PROSAIL hybrid inversion pipeline (SM4).
#
#         Two thin wrappers that isolate the liquidSVM deserialisation and the
#         prosail::PROSAIL_Hybrid_Apply call from the orchestration script.
#         They reproduce verbatim the corresponding lines of
#         Main_s2_03_PROSAIL_Inversion_bis.R (lines 174-179).
# ---

# ── load_svr_ensemble ──────────────────────────────────────────────────────────

#' @title Load and deserialise a serialised SVR ensemble
#'
#' @description
#' Reads an RDS file produced by \code{revision/R/prosail_lut.R}
#' \code{train_all_simulations()} and deserialises each liquidSVM model it
#' contains, returning a list ready for \code{prosail::PROSAIL_Hybrid_Apply}.
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
#'   \code{prosail::PROSAIL_Hybrid_Apply()}.
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
#' Thin wrapper around \code{prosail::PROSAIL_Hybrid_Apply}. Takes a
#' deserialised SVR ensemble (produced by \code{load_svr_ensemble()}) and a
#' reflectance matrix sub-selected to the training bands, and returns the raw
#' inversion object containing both the mean LAI estimate and its standard
#' deviation across ensemble members.
#'
#' @details
#' This function reproduces verbatim
#' \code{Main_s2_03_PROSAIL_Inversion_bis.R} line 176:
#' \preformatted{
#'   LAIest <- PROSAIL_Hybrid_Apply(RegressionModels = modelSVR,
#'                                  Refl = Refl[, selbands])
#' }
#' The caller is responsible for sub-selecting the reflectance columns to the
#' training bands (via \code{selbands}) before passing
#' \code{reflectance_matrix}. The function does not perform any column
#' filtering internally.
#'
#' The raw return value is intentionally preserved (not split into two
#' vectors) so that the caller can access \code{$MeanEstimate} and
#' \code{$StdEstimate} explicitly, mirroring the original script.
#'
#' @param svr_ensemble List. A deserialised SVR ensemble as returned by
#'   \code{load_svr_ensemble()}.
#' @param reflectance_matrix A matrix or data.frame of Sentinel-2 reflectance
#'   values (pixels × bands), already scaled to [0, 1] (i.e. divided by
#'   10 000) and sub-selected to the bands used during training. Typically
#'   \code{Refl[, selbands]} where \code{Refl} is the full reflectance table
#'   read from a \code{S2_reflectance_*.csv} file.
#'
#' @return The raw object returned by \code{prosail::PROSAIL_Hybrid_Apply},
#'   a list with at least two components:
#'   \describe{
#'     \item{\code{MeanEstimate}}{Numeric vector of length \code{nrow(Refl)}.
#'       Mean LAI estimate across ensemble members (m2/m2).}
#'     \item{\code{StdEstimate}}{Numeric vector of length \code{nrow(Refl)}.
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
#'   lai_mean <- LAIest$MeanEstimate
#'   lai_sd   <- LAIest$StdEstimate
#' }
apply_svr_to_reflectance <- function(svr_ensemble, reflectance_matrix) {
  prosail::PROSAIL_Hybrid_Apply(
    RegressionModels = svr_ensemble,
    Refl             = reflectance_matrix
  )
}
