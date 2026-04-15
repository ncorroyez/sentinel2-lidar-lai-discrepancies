# ---
# title:  sm6_metrics.R
# desc:   Compute per-class LAI comparison metrics for SM6b.
#         Exposes compute_metrics_by_class().
#
#         Refactor of compute_metrics v2 in:
#         02_CODES/three_factors_analysis_final.R lines 210-296.
#
#         Regression convention: lm(lidar ~ s2) — lidar is the response,
#         s2 is the predictor. Matches legacy line 270. This is the INVERSE
#         of the SM5 convention (lm(s2 ~ lidar)) and is intentionally
#         preserved here to faithfully reproduce the SM6 legacy behaviour.
#         Do not harmonise with SM5 without explicit instruction.
#
#         Extension relative to legacy: R2 = r^2 is added (not in legacy).
# ---

#' @title Compute LAI comparison metrics by heterogeneity class
#'
#' @description
#' For each combination of (site, LAI combination, heterogeneity class),
#' computes seven statistical metrics comparing LiDAR-derived LAI (response)
#' against Sentinel-2 LAI (predictor):
#'
#' \enumerate{
#'   \item \strong{R}: Pearson correlation coefficient,
#'     \code{cor(lidar, s2, use = "complete.obs")}.
#'   \item \strong{R2}: Coefficient of determination, \eqn{R^2 = r^2}.
#'     Extension relative to the legacy (not present in
#'     \code{three_factors_analysis_final.R}).
#'   \item \strong{RMSE}: Root mean squared error,
#'     \eqn{\sqrt{\overline{(lidar - s2)^2}}}.
#'   \item \strong{Bias}: Mean signed difference,
#'     \eqn{\overline{lidar - s2}}.
#'   \item \strong{Bias_pvalue}: Two-sided p-value from a one-sample t-test
#'     on the pairwise differences, \code{t.test(lidar - s2)$p.value}.
#'     Tests \eqn{H_0: \mu_{lidar - s2} = 0}.
#'   \item \strong{Slope}: Slope \eqn{\hat\beta_1} of \code{lm(lidar ~ s2)}.
#'   \item \strong{Slope_pvalue}: Two-sided p-value for
#'     \eqn{H_0: \beta_1 = 1} (not \eqn{\beta_1 = 0}), computed as
#'     \eqn{2 P\!\left(|t_{n-2}| > |(\hat\beta_1 - 1) / SE(\hat\beta_1)|\right)}.
#' }
#'
#' \strong{Regression convention}: \code{lm(lidar ~ s2)}, where \code{lidar}
#' is the \emph{response} and \code{s2} the \emph{predictor}. This matches
#' \code{02_CODES/three_factors_analysis_final.R} \strong{line 270} and is
#' the \emph{inverse} of the SM5 convention (\code{lm(s2 ~ lidar)}).
#' The asymmetry is inherited from the legacy and is intentionally not
#' corrected here (passe 1). Document the discrepancy when reporting results.
#'
#' Groups with fewer than 10 complete (non-NA) pairs are retained in the
#' output with \code{NA} for all metric columns and \code{n} set to the
#' actual (insufficient) count.
#'
#' @param dt            A \code{data.table} with column \code{"het_class"}
#'   (added by \code{classify_heterogeneity()}) and the LAI columns
#'   referenced by \code{combinations}.
#' @param combinations  List of combination definitions as returned by
#'   \code{get_lai_combinations()}.
#' @param metric_source Character. Label for the heterogeneity metric source,
#'   \code{"DSM"} or \code{"CHM"}. Stored in the output
#'   \code{metric_source} column.
#' @param downsample_n  Integer or \code{NULL}. Maximum number of complete
#'   pairs to use per (site × combination × class) group before computing
#'   metrics. \code{NULL} (default) uses all available pixels — this is the
#'   recommended default for the refactored pipeline. Pass
#'   \code{downsample_n = 5000L} to replicate the legacy sampling strategy
#'   from \code{three_factors_analysis_final.R} lines 229-234
#'   (\code{slice_sample(n = 5000)}) for direct comparability with results
#'   in the submitted manuscript. The \code{n >= 10} filter is applied
#'   \emph{after} downsampling.
#' @param seed          Integer. Random seed set via \code{set.seed()} before
#'   any downsampling. Ignored when \code{downsample_n} is \code{NULL}.
#'   Default \code{42L} (matches \code{set.seed(42)} in the legacy).
#'
#' @return A \code{data.table} with one row per
#'   (site × metric_source × combination × het_class):
#' \describe{
#'   \item{site}{Character.}
#'   \item{metric_source}{Character. \code{"DSM"} or \code{"CHM"}.}
#'   \item{combination}{Character. One of \code{"ATBD_vs_ALS"},
#'     \code{"ATBD_vs_ALS_dopt"}, \code{"opt_vs_ALS_dopt"}.}
#'   \item{het_class}{Character. \code{"Low"}, \code{"Medium"}, or
#'     \code{"High"}.}
#'   \item{n}{Integer. Number of complete pairs used.}
#'   \item{R}{Numeric. Pearson r.}
#'   \item{R2}{Numeric. \eqn{r^2}.}
#'   \item{RMSE}{Numeric.}
#'   \item{Bias}{Numeric.}
#'   \item{Bias_pvalue}{Numeric.}
#'   \item{Slope}{Numeric.}
#'   \item{Slope_pvalue}{Numeric.}
#' }
#' Total rows: \code{|sites| × |combinations| × 3 classes = 27} (before
#' filtering rows where \code{n < 10}).
#'
#' @references
#' Legacy: 02_CODES/three_factors_analysis_final.R lines 210-296
#' (\code{compute_metrics} v2 with p-values).
#'
#' @examples
#' \dontrun{
#' combos <- get_lai_combinations()
#' result <- compute_metrics_by_class(dt, combos, metric_source = "DSM")
#' result[het_class == "Low"]
#' }
#'
#' @export
compute_metrics_by_class <- function(dt, combinations, metric_source,
                                     downsample_n = NULL, seed = 42L) {
  sites       <- unique(dt[["site"]])
  het_classes <- c("Low", "Medium", "High")

  if (!is.null(downsample_n)) {
    set.seed(seed)
  }

  n_rows <- length(sites) * length(combinations) * length(het_classes)
  rows   <- vector("list", n_rows)
  idx    <- 0L

  for (site in sites) {
    dt_site <- dt[dt[["site"]] == site, ]

    for (combo in combinations) {
      lidar_col  <- combo$lidar_col
      s2_col     <- combo$s2_col
      combo_name <- combo$name

      for (cls in het_classes) {
        idx <- idx + 1L

        # ── Pair-wise complete cases for this class ───────────────────────────
        in_class <- dt_site[["het_class"]] == cls
        in_class[is.na(in_class)] <- FALSE

        lidar_raw <- dt_site[[lidar_col]][in_class]
        s2_raw    <- dt_site[[s2_col]][in_class]
        complete  <- !is.na(lidar_raw) & !is.na(s2_raw)

        lidar_vec <- lidar_raw[complete]
        s2_vec    <- s2_raw[complete]

        # ── Optional downsampling ─────────────────────────────────────────────
        # downsample_n = NULL  → use all pixels (default refactor behaviour)
        # downsample_n = 5000  → replicates legacy slice_sample(n = 5000)
        if (!is.null(downsample_n) && length(lidar_vec) > downsample_n) {
          idx_sample <- sample(length(lidar_vec), downsample_n)
          lidar_vec  <- lidar_vec[idx_sample]
          s2_vec     <- s2_vec[idx_sample]
        }

        n_obs     <- length(lidar_vec)

        if (n_obs < 10L) {
          rows[[idx]] <- data.table::data.table(
            site          = site,
            metric_source = metric_source,
            combination   = combo_name,
            het_class     = cls,
            n             = n_obs,
            R             = NA_real_,
            R2            = NA_real_,
            RMSE          = NA_real_,
            Bias          = NA_real_,
            Bias_pvalue   = NA_real_,
            Slope         = NA_real_,
            Slope_pvalue  = NA_real_
          )
          next
        }

        # ── R and R2 ──────────────────────────────────────────────────────────
        r_val  <- stats::cor(lidar_vec, s2_vec, use = "complete.obs")
        r2_val <- r_val ^ 2

        # ── RMSE ─────────────────────────────────────────────────────────────
        rmse_val <- sqrt(mean((lidar_vec - s2_vec) ^ 2, na.rm = TRUE))

        # ── Bias and Bias_pvalue ──────────────────────────────────────────────
        diff_vec <- lidar_vec - s2_vec
        bias_val <- mean(diff_vec, na.rm = TRUE)
        bias_p   <- if (sum(!is.na(diff_vec)) > 1L) {
          stats::t.test(diff_vec)$p.value
        } else {
          NA_real_
        }

        # ── Slope lm(lidar ~ s2) and Slope_pvalue (H0: slope = 1) ───────────
        # Convention: lm(lidar ~ s2) — lidar is response, s2 is predictor.
        # Matches legacy line 270 of three_factors_analysis_final.R.
        # INVERSE of SM5 convention lm(s2 ~ lidar) — do not harmonise.
        lm_fit    <- stats::lm(lidar_vec ~ s2_vec)
        slope_val <- stats::coef(lm_fit)[2L]
        slope_p   <- if (n_obs > 2L) {
          se_slope <- summary(lm_fit)$coefficients[2L, 2L]
          t_stat   <- (slope_val - 1) / se_slope
          2 * stats::pt(abs(t_stat),
                        df         = stats::df.residual(lm_fit),
                        lower.tail = FALSE)
        } else {
          NA_real_
        }

        rows[[idx]] <- data.table::data.table(
          site          = site,
          metric_source = metric_source,
          combination   = combo_name,
          het_class     = cls,
          n             = n_obs,
          R             = r_val,
          R2            = r2_val,
          RMSE          = rmse_val,
          Bias          = bias_val,
          Bias_pvalue   = bias_p,
          Slope         = slope_val,
          Slope_pvalue  = slope_p
        )
      }
    }
  }

  data.table::rbindlist(rows)
}
