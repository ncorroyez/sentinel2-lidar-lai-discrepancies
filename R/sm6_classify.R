# ---
# title:  sm6_classify.R
# desc:   Classify pixels into heterogeneity classes and define the three
#         LAI combination configurations used throughout SM6b.
#         Exposes classify_heterogeneity() and get_lai_combinations().
# ---

#' @title Add heterogeneity class column to a pixel data.table
#'
#' @description
#' Classifies each pixel into Low, Medium, or High heterogeneity using
#' calibrated thresholds from \code{calibrate_thresholds()}. The new column
#' \code{het_class} is added to the data.table \emph{by reference}
#' (\code{data.table::set}), modifying the object in place without a copy.
#'
#' Classification uses \code{base::cut()} with
#' \code{breaks = c(-Inf, low_threshold, high_threshold, Inf)},
#' \code{right = FALSE}, \code{include.lowest = TRUE}, matching the legacy
#' interval convention (lower-closed, upper-open bins) from
#' \code{02_CODES/three_factors_analysis_final.R} lines 216-222. Pixels
#' where \code{metric_col} is \code{NA} receive \code{NA} in
#' \code{het_class}.
#'
#' @param dt             A \code{data.table} as produced by
#'   \code{build_multisite_dt()}. Modified by reference.
#' @param metric_col     Character. Name of the heterogeneity metric column
#'   (\code{"dsm_sd"} or \code{"chm_sd"}).
#' @param low_threshold  Numeric. Threshold between Low and Medium classes.
#' @param high_threshold Numeric. Threshold between Medium and High classes.
#'
#' @return The input \code{data.table}, invisibly, with an added
#'   \code{het_class} column (character: \code{"Low"}, \code{"Medium"},
#'   \code{"High"}, or \code{NA} where \code{metric_col} is \code{NA}).
#'
#' @examples
#' \dontrun{
#' thr <- calibrate_thresholds(dt, "dsm_sd",
#'                             sites = c("Aigoual", "Blois", "Mormal"))
#' classify_heterogeneity(dt, "dsm_sd",
#'                        low_threshold  = thr$low,
#'                        high_threshold = thr$high)
#' dt[, .N, by = .(site, het_class)]
#' }
#'
#' @export
classify_heterogeneity <- function(dt, metric_col, low_threshold,
                                   high_threshold) {
  cls <- cut(
    dt[[metric_col]],
    breaks         = c(-Inf, low_threshold, high_threshold, Inf),
    labels         = c("Low", "Medium", "High"),
    right          = FALSE,
    include.lowest = TRUE
  )
  data.table::set(dt, j = "het_class", value = as.character(cls))
  invisible(dt)
}


#' @title Return the three LAI combination definitions for SM6b
#'
#' @description
#' Returns a list of three named lists, each describing one LAI comparison
#' used in the SM6b analysis. The combinations replicate the three
#' configurations from the legacy
#' \code{02_CODES/three_factors_analysis_final.R} lines 58-61:
#' \enumerate{
#'   \item \code{("LAI_S2_ATBD", "LAI_ALS")}: ATBD S2 vs. full ALS LAI.
#'   \item \code{("LAI_S2_ATBD", "LAI_ALS_dopt")}: ATBD S2 vs. ALS at d_opt.
#'   \item \code{("LAI_S2_opt", "LAI_ALS_dopt")}: optimal S2 vs. ALS at d_opt.
#' }
#'
#' Each element has three fields:
#' \describe{
#'   \item{name}{Short identifier used in the output \code{combination} column.}
#'   \item{lidar_col}{Name of the LiDAR LAI column in the data.table. This is
#'     the \emph{predictor} variable in \code{lm(s2 ~ lidar)} — see the
#'     regression convention documented in \code{compute_metrics_by_class()}.}
#'   \item{s2_col}{Name of the Sentinel-2 LAI column in the data.table.
#'     This is the \emph{response} in \code{lm(s2 ~ lidar)}.}
#' }
#'
#' @return A list of three named lists, each with elements
#'   \code{name}, \code{lidar_col}, \code{s2_col}.
#'
#' @examples
#' combos <- get_lai_combinations()
#' combos[[1]]$name      # "ATBD_vs_ALS"
#' combos[[1]]$lidar_col # "lai_als"
#' combos[[1]]$s2_col    # "lai_s2_atbd"
#'
#' @export
get_lai_combinations <- function() {
  list(
    list(
      name      = "ATBD_vs_ALS",
      lidar_col = "lai_als",
      s2_col    = "lai_s2_atbd"
    ),
    list(
      name      = "ATBD_vs_ALS_dopt",
      lidar_col = "lai_als_dopt",
      s2_col    = "lai_s2_atbd"
    ),
    list(
      name      = "opt_vs_ALS_dopt",
      lidar_col = "lai_als_dopt",
      s2_col    = "lai_s2_opt"
    )
  )
}
