# ---
# title:  sm6_metrics.R
# desc:   Compute per-class LAI comparison metrics for SM6b.
#         Exposes compute_metrics_by_class() and sample_indices_uniform_lai().
#
#         Within each (site × heterogeneity class), pixels are sampled
#         uniformly on LAI_ALS in [2, class-specific max LAI_ALS]; the same
#         indices are then used to evaluate all three LAI combinations.
#         This mirrors the stratified-uniform LAI_ALS sampling described in
#         the manuscript SM6 methodology.
#
#         Regression convention: lm(s2 ~ lidar) — s2 is the response,
#         lidar is the predictor. Harmonized with SM5 convention as of
#         2026-04-15.
#
#         Extension relative to legacy: R2 = r^2 is added (not in legacy).
# ---

# ── sample_indices_uniform_lai ────────────────────────────────────────────────

#' @title Uniform-LAI sample indices
#'
#' @description Build row indices for a sample with approximately uniform
#'   LAI_ALS distribution between \code{lai_min} and \code{floor(quantile(
#'   lai_values, lai_max_quantile))}. Bins are \code{bin_width}-wide
#'   (1 unit by default, matching the manuscript). Values outside
#'   [\code{lai_min}, upper bound] or NA are excluded.
#'
#' @param lai_values        Numeric vector of LAI_ALS values.
#' @param n_target          Target sample size.
#' @param lai_min           Lower bound (default 2 m²/m²).
#' @param lai_max_quantile  Quantile used to derive the upper bound
#'                          (default \code{0.98}). The upper bound is
#'                          \code{floor(quantile)}. Set to \code{1} or
#'                          \code{NULL} to use the observed max.
#' @param bin_width         LAI bin width (default 1.0; produces integer
#'                          breaks from \code{lai_min} to \code{floor(p)}).
#' @param seed              RNG seed.
#'
#' @return Integer vector of row indices into \code{lai_values}.
#' @export
sample_indices_uniform_lai <- function(lai_values, n_target,
                                       lai_min = 2.0,
                                       lai_max_quantile = 0.98,
                                       bin_width = 1.0,
                                       seed = 42L) {
  if (length(lai_values) == 0L) return(integer(0))
  set.seed(seed)

  valid_lo <- !is.na(lai_values) & lai_values >= lai_min
  if (sum(valid_lo) == 0L) return(integer(0))

  if (is.null(lai_max_quantile) || lai_max_quantile >= 1) {
    lai_max <- floor(max(lai_values[valid_lo]))
  } else {
    q <- as.numeric(stats::quantile(lai_values[valid_lo],
                                     probs = lai_max_quantile,
                                     na.rm = TRUE))
    lai_max <- floor(q)
  }

  if (lai_max <= lai_min) return(which(valid_lo))

  bin_edges <- seq(lai_min, lai_max, by = bin_width)
  n_bins    <- length(bin_edges) - 1L
  if (n_bins < 1L) return(which(valid_lo))

  valid     <- valid_lo & lai_values <= lai_max
  bin_idx   <- cut(lai_values, breaks = bin_edges,
                   include.lowest = TRUE, labels = FALSE)
  n_per_bin <- ceiling(n_target / n_bins)

  out <- integer(0)
  for (b in seq_len(n_bins)) {
    pool <- which(valid & bin_idx == b)
    if (length(pool) == 0L) next
    take <- min(length(pool), n_per_bin)
    out  <- c(out, if (length(pool) == take) pool else sample(pool, take))
  }
  out
}


# ── sample_uniform_lai_dt ─────────────────────────────────────────────────────

#' @title Per-site uniform-LAI sample of a multi-site data.table
#'
#' @description Returns a sub-set of \code{dt} where, for each site,
#'   pixels are sampled to approximate a uniform LAI_ALS distribution in
#'   [\code{lai_min}, \code{floor(quantile(lai_als, lai_max_quantile))}].
#'   This is the manuscript's site-level uniform-LAI sampling used as the
#'   input to the SM6 heterogeneity classification.
#'
#' @param dt                 \code{data.table} with columns \code{site}
#'                           and \code{lai_als}.
#' @param n_target_per_site  Target sample size per site.
#' @param lai_min            Lower bound (default 2).
#' @param lai_max_quantile   Quantile defining the upper bound
#'                           (default 0.98 → \code{floor(p98)}).
#' @param bin_width          Bin width in LAI units (default 1).
#' @param seed               Base RNG seed (offset per site).
#'
#' @return A \code{data.table} with the same columns as \code{dt}.
#' @export
sample_uniform_lai_dt <- function(dt, n_target_per_site,
                                  lai_min = 2.0,
                                  lai_max_quantile = 0.98,
                                  bin_width = 1.0,
                                  seed = 42L) {
  sites <- unique(dt[["site"]])
  out   <- vector("list", length(sites))
  for (i in seq_along(sites)) {
    sub <- dt[site == sites[i]]
    idx <- sample_indices_uniform_lai(
      lai_values       = sub$lai_als,
      n_target         = n_target_per_site,
      lai_min          = lai_min,
      lai_max_quantile = lai_max_quantile,
      bin_width        = bin_width,
      seed             = seed + i - 1L
    )
    out[[i]] <- sub[idx]
  }
  data.table::rbindlist(out)
}

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
#'     \eqn{\overline{s2 - lidar}}.
#'     Positive Bias = S2 overestimates LiDAR.
#'   \item \strong{Bias_pvalue}: Two-sided p-value from a one-sample t-test
#'     on the pairwise differences, \code{t.test(s2 - lidar)$p.value}.
#'     Tests \eqn{H_0: \mu_{s2 - lidar} = 0}.
#'   \item \strong{Slope}: Slope \eqn{\hat\beta_1} of \code{lm(s2 ~ lidar)}.
#'   \item \strong{Slope_pvalue}: Two-sided p-value for
#'     \eqn{H_0: \beta_1 = 1} (not \eqn{\beta_1 = 0}), computed as
#'     \eqn{2 P\!\left(|t_{n-2}| > |(\hat\beta_1 - 1) / SE(\hat\beta_1)|\right)}.
#' }
#'
#' \strong{Regression convention}: \code{lm(s2 ~ lidar)}, where \code{s2}
#' is the \emph{response} and \code{lidar} the \emph{predictor}. Harmonized
#' with the SM5 convention as of 2026-04-15.
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
#' @param uniform_sample_n  Integer or \code{NULL}. If \code{NULL}
#'   (recommended for SM6: sample once site-level via
#'   \code{sample_uniform_lai_dt()} BEFORE this call, then pass the
#'   sampled \code{dt} here), no further sampling is done. If an integer,
#'   additionally sample uniformly on \code{lai_als} within each class
#'   using \code{sample_indices_uniform_lai()}.
#' @param lai_min          Lower bound on LAI_ALS for the per-class
#'   sample (default \code{2.0}). Ignored when
#'   \code{uniform_sample_n} is \code{NULL}.
#' @param lai_max_quantile Quantile defining the upper bound for the
#'   per-class sample (default \code{0.98}).
#' @param bin_width        Per-class bin width (default \code{1}).
#' @param seed             Integer. RNG seed (default \code{42L}).
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
                                     uniform_sample_n = NULL,
                                     lai_min = 2.0,
                                     lai_max_quantile = 0.98,
                                     bin_width = 1.0,
                                     seed    = 42L,
                                     include_total = TRUE) {
  sites       <- unique(dt[["site"]])
  het_classes <- c("Low", "Medium", "High")

  n_rows <- length(sites) * length(combinations) * length(het_classes)
  rows   <- vector("list", n_rows)
  idx    <- 0L

  # Accumulator: for each site, pooled sampled indices across the 3 classes.
  # Used by the Total block so that Total = union(Low + Medium + High) on the
  # SAMPLED pixels, not the full unsampled raster.
  pooled_sampled <- vector("list", length(sites))
  names(pooled_sampled) <- sites

  for (site_name in sites) {
    dt_site <- dt[site == site_name]
    pooled_sampled[[site_name]] <- integer(0)

    for (cls in het_classes) {

      # ── Build a single uniform-LAI sample per (site × class) ────────────────
      in_class_idx <- which(!is.na(dt_site[["het_class"]]) &
                              dt_site[["het_class"]] == cls)

      if (length(in_class_idx) == 0L) {
        for (combo in combinations) {
          idx <- idx + 1L
          rows[[idx]] <- data.table::data.table(
            site = site_name, metric_source = metric_source,
            combination = combo$name, het_class = cls,
            n = 0L, R = NA_real_, R2 = NA_real_, RMSE = NA_real_,
            Bias = NA_real_, Bias_pvalue = NA_real_,
            Slope = NA_real_, Slope_pvalue = NA_real_
          )
        }
        next
      }

      lai_class <- dt_site[["lai_als"]][in_class_idx]

      if (is.null(uniform_sample_n)) {
        sampled_idx <- in_class_idx
      } else {
        sub <- sample_indices_uniform_lai(
          lai_values       = lai_class,
          n_target         = uniform_sample_n,
          lai_min          = lai_min,
          lai_max_quantile = lai_max_quantile,
          bin_width        = bin_width,
          seed             = seed + match(cls, het_classes) - 1L
        )
        sampled_idx <- in_class_idx[sub]
      }

      # Accumulate for the Total row
      pooled_sampled[[site_name]] <- c(pooled_sampled[[site_name]], sampled_idx)

      # ── Evaluate each combo on the shared sample ────────────────────────────
      for (combo in combinations) {
        idx <- idx + 1L
        lidar_col  <- combo$lidar_col
        s2_col     <- combo$s2_col
        combo_name <- combo$name

        lidar_raw <- dt_site[[lidar_col]][sampled_idx]
        s2_raw    <- dt_site[[s2_col]][sampled_idx]
        complete  <- !is.na(lidar_raw) & !is.na(s2_raw)

        lidar_vec <- lidar_raw[complete]
        s2_vec    <- s2_raw[complete]
        n_obs     <- length(lidar_vec)

        if (n_obs < 10L) {
          rows[[idx]] <- data.table::data.table(
            site          = site_name,
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

        # ── Bias and Bias_pvalue — Harmonized with lm(s2 ~ lidar) convention as of 2026-04-16.
        # Bias = mean(s2 - lidar): positive = S2 overestimates LiDAR.
        diff_vec <- s2_vec - lidar_vec
        bias_val <- mean(diff_vec, na.rm = TRUE)
        bias_p   <- if (sum(!is.na(diff_vec)) > 1L) {
          stats::t.test(diff_vec)$p.value
        } else {
          NA_real_
        }

        # ── Slope lm(s2 ~ lidar) and Slope_pvalue (H0: slope = 1) ──────────
        # Convention: lm(s2 ~ lidar) — s2 is response, lidar is predictor.
        # Harmonized with SM5 convention as of 2026-04-15.
        lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
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
          site          = site_name,
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

  # ── Total rows: pooled sampled pixels across (Low + Medium + High) ───────
  # Same sampling discipline as the per-class rows: Total = union of the
  # three class samples. This keeps Total comparable with the class rows
  # (n ≈ sum of class n) rather than the full unsampled raster.
  if (include_total) {
    total_rows <- vector("list", length(sites) * length(combinations))
    tidx <- 0L
    for (site_name in sites) {
      dt_site     <- dt[site == site_name]
      sampled_idx <- pooled_sampled[[site_name]]
      for (combo in combinations) {
        tidx       <- tidx + 1L
        lidar_col  <- combo$lidar_col
        s2_col     <- combo$s2_col
        combo_name <- combo$name

        lidar_raw <- dt_site[[lidar_col]][sampled_idx]
        s2_raw    <- dt_site[[s2_col]][sampled_idx]
        complete  <- !is.na(lidar_raw) & !is.na(s2_raw)
        lidar_vec <- lidar_raw[complete]
        s2_vec    <- s2_raw[complete]
        n_obs     <- length(lidar_vec)

        if (n_obs < 10L) {
          total_rows[[tidx]] <- data.table::data.table(
            site = site_name, metric_source = metric_source,
            combination = combo_name, het_class = "Total",
            n = n_obs, R = NA_real_, R2 = NA_real_, RMSE = NA_real_,
            Bias = NA_real_, Bias_pvalue = NA_real_,
            Slope = NA_real_, Slope_pvalue = NA_real_
          )
          next
        }

        r_val     <- stats::cor(lidar_vec, s2_vec, use = "complete.obs")
        rmse_val  <- sqrt(mean((lidar_vec - s2_vec)^2, na.rm = TRUE))
        diff_vec  <- s2_vec - lidar_vec
        bias_val  <- mean(diff_vec, na.rm = TRUE)
        bias_p    <- stats::t.test(diff_vec)$p.value
        lm_fit    <- stats::lm(s2_vec ~ lidar_vec)
        slope_val <- stats::coef(lm_fit)[2L]
        se_slope  <- summary(lm_fit)$coefficients[2L, 2L]
        t_stat    <- (slope_val - 1) / se_slope
        slope_p   <- 2 * stats::pt(abs(t_stat),
                                   df         = stats::df.residual(lm_fit),
                                   lower.tail = FALSE)

        total_rows[[tidx]] <- data.table::data.table(
          site          = site_name,
          metric_source = metric_source,
          combination   = combo_name,
          het_class     = "Total",
          n             = n_obs,
          R             = r_val,
          R2            = r_val ^ 2,
          RMSE          = rmse_val,
          Bias          = bias_val,
          Bias_pvalue   = bias_p,
          Slope         = slope_val,
          Slope_pvalue  = slope_p
        )
      }
    }
    rows <- c(rows, total_rows)
  }

  data.table::rbindlist(rows)
}
