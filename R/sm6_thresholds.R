# ---
# title:  sm6_thresholds.R
# desc:   Dynamic threshold calibration for canopy heterogeneity classes.
#         Exposes calibrate_thresholds().
#
#         Replaces the fixed thresholds (c(-Inf, 2.5, 5, Inf)) used in the
#         legacy 02_CODES/three_factors_analysis_final.R lines 217-222.
# ---

#' @title Calibrate Low/Medium/High heterogeneity thresholds dynamically
#'
#' @description
#' Searches a grid of quantile-derived threshold pairs \code{(low, high)} and
#' selects the pair that (1) ensures a minimum pixel count in every
#' (site × class) cell, and (2) maximises class balance on the pooled
#' 3-site distribution.
#'
#' \strong{Algorithm}:
#' \enumerate{
#'   \item Build a candidate grid from quantiles of the pooled (3-site)
#'     non-NA distribution:
#'     \code{low_quantiles  = seq(0.05, 0.50, by = 0.01)},
#'     \code{high_quantiles = seq(0.50, 0.95, by = 0.01)}.
#'     Convert to absolute values via \code{stats::quantile()}. Keep only
#'     pairs with \code{low_val < high_val} (strict).
#'   \item For each candidate pair, classify the full data.table into
#'     Low / Medium / High via \code{base::cut()} and compute a
#'     \eqn{|sites| \times 3} count matrix.
#'   \item Retain pairs where \code{min(counts) >= target_n_per_class} across
#'     all (site × class) cells. If no pair qualifies, raise an informative
#'     error asking to relax the constraint.
#'   \item Among qualifying pairs, select the one minimising the balance
#'     score: \code{max(abs(pool_props - 1/3))}, where \code{pool_props} is
#'     the 3-class proportional distribution on the pooled sample.
#'     A score of 0 means perfectly equal thirds; lower is better.
#' }
#'
#' @param df               A \code{data.table} with at least columns
#'   \code{"site"} and the column named by \code{metric_col}.
#' @param metric_col       Character. Name of the heterogeneity metric column
#'   (e.g. \code{"dsm_sd"} or \code{"chm_sd"}).
#' @param sites            Character vector. Site names expected in
#'   \code{df[["site"]]}. Determines the rows of the count matrix.
#' @param target_n_per_class Integer. Minimum pixel count required in every
#'   (site × class) cell. Default \code{5000L}.
#'
#' @return A named list:
#' \describe{
#'   \item{low}{Numeric. Lower threshold (Low | Medium boundary).}
#'   \item{high}{Numeric. Upper threshold (Medium | High boundary).}
#'   \item{balance_score}{Numeric. \code{max(abs(pool_props - 1/3))} for the
#'     selected pair. Lower is better.}
#'   \item{counts}{Integer matrix (\code{|sites| × 3}). Row names are site
#'     names; column names are \code{c("Low", "Medium", "High")}.}
#' }
#'
#' @examples
#' \dontrun{
#' thr <- calibrate_thresholds(
#'   dt, "dsm_sd",
#'   sites = c("Aigoual", "Blois", "Mormal")
#' )
#' cat("Low:", thr$low, "  High:", thr$high,
#'     "  balance:", round(thr$balance_score, 3), "\n")
#' print(thr$counts)
#' }
#'
#' @export
calibrate_thresholds <- function(df, metric_col, sites,
                                 target_n_per_class = 5000L) {
  vals        <- df[[metric_col]]
  vals_non_na <- vals[!is.na(vals)]

  if (length(vals_non_na) == 0L) {
    stop("calibrate_thresholds: '", metric_col,
         "' contains no non-NA values.")
  }

  # ── 1. Candidate grid ────────────────────────────────────────────────────────
  low_q  <- seq(0.05, 0.50, by = 0.01)
  high_q <- seq(0.50, 0.95, by = 0.01)

  low_vals  <- stats::quantile(vals_non_na, probs = low_q,  names = FALSE)
  high_vals <- stats::quantile(vals_non_na, probs = high_q, names = FALSE)

  grid <- expand.grid(
    low_val  = low_vals,
    high_val = high_vals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$low_val < grid$high_val, , drop = FALSE]

  n_cand <- nrow(grid)
  if (n_cand == 0L) {
    stop("calibrate_thresholds: no valid (low < high) threshold pairs found.")
  }

  # ── 2. Evaluate every candidate pair ─────────────────────────────────────────
  balance_scores  <- numeric(n_cand)
  min_counts_vec  <- integer(n_cand)
  count_matrices  <- vector("list", n_cand)

  site_col   <- df[["site"]]
  metric_vec <- df[[metric_col]]
  cls_labels <- c("Low", "Medium", "High")

  for (i in seq_len(n_cand)) {
    lv <- grid$low_val[i]
    hv <- grid$high_val[i]

    cls <- cut(
      metric_vec,
      breaks         = c(-Inf, lv, hv, Inf),
      labels         = cls_labels,
      right          = FALSE,
      include.lowest = TRUE
    )

    # Per-site × per-class counts
    counts_mat <- matrix(
      0L,
      nrow     = length(sites),
      ncol     = 3L,
      dimnames = list(sites, cls_labels)
    )
    for (s in sites) {
      idx <- site_col == s & !is.na(cls)
      tbl <- table(cls[idx])
      for (cl in cls_labels) {
        counts_mat[s, cl] <- if (cl %in% names(tbl)) as.integer(tbl[[cl]]) else 0L
      }
    }

    # Pooled balance score
    pool_tbl <- table(cls[!is.na(cls)])
    pool_vec <- as.integer(pool_tbl)
    # Pad to length 3 if a class is absent from pool
    if (length(pool_vec) < 3L) {
      pool_vec <- c(pool_vec, rep(0L, 3L - length(pool_vec)))
    }
    pool_total <- sum(pool_vec)
    pool_props <- if (pool_total > 0L) pool_vec / pool_total else rep(0, 3L)

    balance_scores[i]  <- max(abs(pool_props - 1 / 3))
    min_counts_vec[i]  <- min(counts_mat)
    count_matrices[[i]] <- counts_mat
  }

  # ── 3. Filter by target_n_per_class ─────────────────────────────────────────
  valid <- which(min_counts_vec >= target_n_per_class)
  if (length(valid) == 0L) {
    stop(
      "calibrate_thresholds: no threshold pair satisfies ",
      "target_n_per_class = ", target_n_per_class,
      " across all ", length(sites), " sites \u00d7 3 classes. ",
      "Consider relaxing the constraint (lower target_n_per_class)."
    )
  }

  # ── 4. Select best balance ────────────────────────────────────────────────────
  best_idx <- valid[which.min(balance_scores[valid])]

  list(
    low           = grid$low_val[best_idx],
    high          = grid$high_val[best_idx],
    balance_score = balance_scores[best_idx],
    counts        = count_matrices[[best_idx]]
  )
}
