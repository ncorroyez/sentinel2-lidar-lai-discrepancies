# ---
# title:  sm5_dopt.R
# desc:   Selection of the optimal canopy integration depth (d_opt) for the
#         d_opt analysis (SM5). Implements four selection methods: pearson
#         (which.max R), rmse (which.min RMSE), bias (which.min |Bias|), and
#         pareto (multi-criteria score: mean of normalised criteria).
#
#         Addresses reviewer comments R3.minor.13 and R4.spec.2, which
#         requested testing selection criteria alternative to Pearson
#         correlation alone.
# ---

# ── compute_pareto_front ───────────────────────────────────────────────────────

#' @title Compute the Pareto-optimal canopy depth using four agreement metrics
#'
#' @description
#' Given the metrics table for a single (Site × Norm × Column) triplet across
#' all available depths (typically 1–38), identifies the Pareto front (non-
#' dominated depths) and selects within that front the depth that minimises
#' the MEAN of normalised minimisable criteria (\code{Avg_score}).
#' The available criteria are:
#' \describe{
#'   \item{m1}{\code{1 - R} (high Pearson r is desirable)}
#'   \item{m2}{\code{RMSE}}
#'   \item{m3}{\code{abs(Bias)}}
#'   \item{m4}{\code{abs(Slope - 1)}}
#' }
#' By default all four are used; passing \code{criteria = c("R","RMSE","Bias")}
#' drops \code{|Slope - 1|}.
#'
#' Three-stage procedure:
#' \enumerate{
#'   \item Each criterion is min-max normalised over all valid candidates
#'         (\code{(m - min) / (max - min)}); constant criteria are set to
#'         0. \code{Avg_score = 1 - mean(m_norm)} ("higher = better",
#'         1 = best, 0 = worst over the valid range).
#'   \item Dominance identifies Front-1 over the same candidate set.
#'         Dominance is invariant to a positive-monotone rescaling, so it
#'         is tested on the raw minimisable criteria.
#'   \item The selected depth is \code{argmax(Avg_score)}, which by
#'         construction lies in Front-1 (a dominated candidate is strictly
#'         worse on at least one criterion, hence has strictly smaller
#'         Avg_score). Ties: \code{which.max()} returns the smallest depth.
#' }
#' A depth d1 dominates d2 if and only if all criterion values at d1 are
#' \eqn{\le} those at d2, with at least one strict inequality.
#'
#' @param metrics_subset A \code{data.table} with columns \code{Depth}
#'   (integer), \code{R2} (numeric), \code{RMSE} (numeric), \code{Bias}
#'   (numeric), \code{Slope} (numeric). Pre-filtered to one
#'   (Site × Norm × Column) triplet. Typically 38 rows (one per depth).
#'
#' @return A named list:
#' \describe{
#'   \item{d_opt}{Integer. Front-1 depth with maximum \code{Avg_score}.
#'                \code{NA_integer_} if fewer than 2 non-NA depths available.}
#'   \item{Avg_score_best}{Numeric. \code{Avg_score} at \code{d_opt}
#'                        (\code{1 - mean(m_norm)} over criteria normalised
#'                        on the full valid candidate set; 1 = best,
#'                        0 = worst). \code{NA_real_} if too few candidates.}
#'   \item{pareto_front_depths}{Integer vector of non-dominated depths.
#'                              \code{NA_integer_} if too few non-NA depths.}
#'   \item{all_distances}{Named numeric vector of \code{Avg_score} for every
#'                        candidate. Names are depth indices as character.
#'                        Higher = better. \code{NULL} if too few candidates.}
#' }
#'
#' @references
#' R3.minor.13, R4.spec.2: reviewers requested multi-criteria d_opt selection
#' and reporting of the full Pareto front.
#'
#' @examples
#' sub <- data.table::data.table(
#'   Depth = 1:5,
#'   R2    = c(0.60, 0.70, 0.80, 0.75, 0.65),
#'   RMSE  = c(1.20, 1.00, 0.80, 0.90, 1.10),
#'   Bias  = c(0.30, 0.10, -0.05, 0.20, 0.40),
#'   Slope = c(0.90, 0.95,  1.00, 0.97, 0.85)
#' )
#' compute_pareto_front(sub)
#'
#' @param criteria Character vector of criteria to include in the Pareto
#'   dominance test AND the Pareto score. Subset of
#'   \code{c("R", "RMSE", "Bias", "Slope")}. Default uses all four. Pass
#'   \code{c("R", "RMSE", "Bias")} to drop \code{|Slope - 1|}.
#'
#' @export
compute_pareto_front <- function(metrics_subset,
                                 criteria    = c("R", "RMSE", "Bias", "Slope"),
                                 norm_method = c("max", "minmax"),
                                 norm_scope  = c("front1", "all")) {
  criteria    <- match.arg(criteria, c("R", "RMSE", "Bias", "Slope"),
                            several.ok = TRUE)
  norm_method <- match.arg(norm_method)
  norm_scope  <- match.arg(norm_scope)

  na_result <- list(
    d_opt               = NA_integer_,
    Avg_score_best      = NA_real_,
    pareto_front_depths = NA_integer_,
    all_distances       = NULL
  )

  depths <- metrics_subset[["Depth"]]

  # Round input metrics to 2 dp BEFORE building the criteria matrix. This
  # makes Pareto dominance and normalised score invariant to sub-0.01
  # floating-point noise: two configs differing only at the 3rd decimal
  # are treated as ties on that criterion.
  R_r     <- round(metrics_subset[["R"]],     2)
  RMSE_r  <- round(metrics_subset[["RMSE"]],  2)
  Bias_r  <- round(metrics_subset[["Bias"]],  2)
  Slope_r <- round(metrics_subset[["Slope"]], 2)

  # Build the metric matrix from the requested criteria (minimisable form)
  m_cols <- list(
    R     = 1 - R_r,
    RMSE  = RMSE_r,
    Bias  = abs(Bias_r),
    Slope = abs(Slope_r - 1)
  )
  m_mat_full <- do.call(cbind, m_cols[criteria])

  valid <- stats::complete.cases(m_mat_full)
  if (sum(valid) < 2L) return(na_result)

  depths_v <- depths[valid]
  m_mat    <- m_mat_full[valid, , drop = FALSE]

  # Pareto front: non-dominated depths over ALL valid depths.
  # d1 dominates d2 iff all m_mat[d2, ] <= m_mat[d1, ] with >= 1 strict ineq.
  # Dominance is invariant to positive-monotone rescaling, so it can be
  # tested directly on the raw minimisable criteria (no need to normalise).
  n         <- nrow(m_mat)
  dominated <- logical(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      if (all(m_mat[j, ] <= m_mat[i, ]) &&
          any(m_mat[j, ] <  m_mat[i, ])) {
        dominated[i] <- TRUE
        break
      }
    }
  }
  front_idx    <- which(!dominated)
  front_depths <- depths_v[front_idx]

  # Normalise criteria min-max:
  #   norm_scope = "all"    → over every valid candidate (use case: d_opt
  #                            selection across depths in [1, dthr]).
  #   norm_scope = "front1" → over Front-1 candidates only (use case:
  #                            PROSAIL OPT selection, where dozens of
  #                            absurd configurations on the ~225 grid
  #                            would otherwise stretch min/max bounds
  #                            and compress the realistic configs).
  #   minmax → (m - min) / (max - min)  ∈ [0, 1]
  #   max    → m / max                  ∈ [min/max, 1]
  # A constant criterion is set to 0.
  norm_set_idx <- if (norm_scope == "all") seq_len(nrow(m_mat)) else front_idx
  m_ref <- m_mat[norm_set_idx, , drop = FALSE]

  proj_one <- function(x, xref) {
    hi <- max(xref)
    if (norm_method == "max") {
      if (hi == 0) return(rep(0, length(x)))
      return(x / hi)
    }
    lo <- min(xref)
    if (hi == lo) return(rep(0, length(x)))
    (x - lo) / (hi - lo)
  }
  m_norm <- vapply(seq_len(ncol(m_mat)), function(j) {
    proj_one(m_mat[, j], m_ref[, j])
  }, numeric(nrow(m_mat)))
  if (!is.matrix(m_norm)) m_norm <- matrix(m_norm, ncol = ncol(m_mat))

  # Avg_score: HIGHER = better. 1 - mean(minmax-normalised), rounded 2 dp.
  # Ties broken by which.max() → smallest depth.
  avg_v <- round(1 - rowMeans(m_norm), 2)

  if (norm_scope == "front1") {
    # Restrict the argmax search to Front-1 (dominated candidates may have
    # Avg_score outside [0, 1] under this projection — not meaningful).
    candidate_idx  <- front_idx
    sub_avg        <- avg_v[candidate_idx]
    idx_in_front   <- which.max(sub_avg)
    d_opt          <- depths_v[candidate_idx[idx_in_front]]
    Avg_score_best <- sub_avg[[idx_in_front]]
    # all_distances: keep only Front-1 values; NA elsewhere.
    all_distances        <- rep(NA_real_, length(depths_v))
    all_distances[front_idx] <- avg_v[front_idx]
  } else {
    idx_max        <- which.max(avg_v)
    d_opt          <- depths_v[idx_max]
    Avg_score_best <- avg_v[[idx_max]]
    all_distances  <- avg_v
  }
  names(all_distances) <- as.character(depths_v)

  list(
    d_opt               = d_opt,
    Avg_score_best      = Avg_score_best,
    pareto_front_depths = front_depths,
    all_distances       = all_distances
  )
}

# ── select_dopt_single (internal) ─────────────────────────────────────────────
#
# Applies one selection method to a (Site × Norm × Column) subset (all depths).
# Returns a one-row data.table with columns:
#   method_dopt, d_opt, R, R2, RMSE, Bias, Slope, Avg_score, pareto_front.
# Avg_score and pareto_front are NA for methods other than "pareto".
# pareto_front is encoded as "3,4,5,7" (comma-separated sorted depths).
#
select_dopt_single <- function(metrics_subset, method,
                                pareto_criteria   = c("R","RMSE","Bias","Slope"),
                                pareto_norm_method = c("max", "minmax"),
                                pareto_norm_scope  = c("front1", "all")) {
  pareto_norm_method <- match.arg(pareto_norm_method)
  pareto_norm_scope  <- match.arg(pareto_norm_scope)
  depths <- metrics_subset[["Depth"]]

  d_opt        <- NA_integer_
  Avg_score    <- NA_real_
  pareto_front <- NA_character_

  if (method == "pearson") {
    idx <- which.max(metrics_subset[["R"]])
    if (length(idx) > 0L) d_opt <- depths[idx]

  } else if (method == "rmse") {
    idx <- which.min(metrics_subset[["RMSE"]])
    if (length(idx) > 0L) d_opt <- depths[idx]

  } else if (method == "bias") {
    idx <- which.min(abs(metrics_subset[["Bias"]]))
    if (length(idx) > 0L) d_opt <- depths[idx]

  } else if (method == "slope") {
    idx <- which.min(abs(metrics_subset[["Slope"]] - 1))
    if (length(idx) > 0L) d_opt <- depths[idx]

  } else if (method == "pareto") {
    pf    <- compute_pareto_front(metrics_subset,
                                   criteria    = pareto_criteria,
                                   norm_method = pareto_norm_method,
                                   norm_scope  = pareto_norm_scope)
    d_opt <- pf$d_opt

    Avg_score <- if (is.na(pf$d_opt)) NA_real_ else pf$Avg_score_best

    pareto_front <- if (anyNA(pf$pareto_front_depths)) {
      NA_character_
    } else {
      paste(sort(pf$pareto_front_depths), collapse = ",")
    }
  }

  # Extract metrics at d_opt; all NA if d_opt itself is NA.
  # Round to 2 dp for consistency with the Pareto computation and CSVs.
  if (is.na(d_opt)) {
    r_val     <- NA_real_
    r2_val    <- NA_real_
    rmse_val  <- NA_real_
    bias_val  <- NA_real_
    slope_val <- NA_real_
  } else {
    row_d     <- metrics_subset[Depth == d_opt]
    r_val     <- round(row_d[["R"]][[1L]],     2)
    r2_val    <- round(row_d[["R2"]][[1L]],    2)
    rmse_val  <- round(row_d[["RMSE"]][[1L]],  2)
    bias_val  <- round(row_d[["Bias"]][[1L]],  2)
    slope_val <- round(row_d[["Slope"]][[1L]], 2)
  }

  data.table::data.table(
    method_dopt  = method,
    d_opt        = d_opt,
    R            = r_val,
    R2           = r2_val,
    RMSE         = rmse_val,
    Bias         = bias_val,
    Slope        = slope_val,
    Avg_score    = Avg_score,
    pareto_front = pareto_front
  )
}

# ── select_prosail_opt ────────────────────────────────────────────────────────

#' @title Select optimal PROSAIL configuration at the reference d_opt
#'
#' @description
#' Two-stage selection for SM5 passe 2, run under two d_opt modes:
#' \describe{
#'   \item{\code{"per_site"}}{Each individual site uses its own d_opt from
#'     \code{dopt_dt}. The Pareto selection over PROSAIL columns is performed
#'     on the metrics of that site at its own optimal depth.}
#'   \item{\code{"all_sites"}}{The d_opt from the \code{"All_sites"} row of
#'     \code{dopt_dt} is fixed and applied identically to each individual site.
#'     The Pareto selection is then performed on per-site metrics at that
#'     common depth, not on the aggregated All_sites data.}
#' }
#' Both modes are always computed and returned in a single table distinguished
#' by the \code{d_opt_source} column. \code{compute_pareto_front()} is reused
#' with PROSAIL columns playing the role of "depths" in criterion space.
#'
#' @param metrics_dt    \code{data.table} from \code{build_metrics_table()}
#'   (SM5 passe 1). Must contain Site, Norm, Depth, Column, R2, RMSE, Bias,
#'   Slope. The "All_sites" rows are used only to supply pooled metrics; the
#'   Pareto selection always runs on individual-site rows.
#' @param dopt_dt       \code{data.table} from \code{dopt_reference.csv}
#'   (ATBD only). Must contain Site, Norm, method_dopt, d_opt. Must include
#'   both individual-site rows and a row with \code{Site == "All_sites"} for
#'   the \code{"all_sites"} mode.
#' @param dopt_method   Character scalar. Which \code{method_dopt} to use for
#'   fixing depth. Default \code{"pareto"}.
#' @param norm_filter   Character scalar passed to \code{grepl()} to subset
#'   Norm values. Default \code{"keepTrees"}. Set to \code{NULL} for all norms.
#' @param individual_sites Character vector. Sites to process individually.
#'   Default \code{c("Aigoual", "Blois", "Mormal")}.
#' @param fixed_d_opt   Integer or \code{NULL}. If provided, adds a third mode
#'   \code{"fixed"} where all sites use this depth (e.g. \code{4L} for the
#'   submitted-paper reference). Pareto selection runs on each site's individual
#'   metrics at that fixed depth. Default \code{NULL} (mode omitted).
#'
#' @return A \code{data.table} with one row per (Site × Norm × d_opt_source):
#' \describe{
#'   \item{d_opt_source}{character — \code{"per_site"} or \code{"all_sites"}}
#'   \item{Site, Norm}{character}
#'   \item{d_opt}{integer — depth used (site-specific or All_sites value)}
#'   \item{method_dopt}{character}
#'   \item{Column_opt}{character — best PROSAIL column at d_opt}
#'   \item{ATBD}{logical}
#'   \item{R, R2, RMSE, Bias, Slope}{numeric — metrics at d_opt for Column_opt}
#'   \item{Avg_score}{numeric — mean of Front-1-normalised criteria for
#'                    Column_opt (1 = best within front, 0 = worst)}
#'   \item{n_pareto_front}{integer — number of non-dominated columns}
#' }
#'
#' @references
#' R3.minor.13, R4.spec.2: multi-criteria PROSAIL configuration selection.
#'
#' @export
select_prosail_opt <- function(metrics_dt,
                                dopt_dt,
                                dopt_method        = "pareto",
                                norm_filter        = "keepTrees",
                                individual_sites   = c("Aigoual", "Blois", "Mormal"),
                                fixed_d_opt        = NULL,
                                pareto_criteria    = c("R","RMSE","Bias","Slope"),
                                pareto_norm_method = c("max", "minmax"),
                                pareto_norm_scope  = c("front1", "all")) {
  pareto_norm_method <- match.arg(pareto_norm_method)
  pareto_norm_scope  <- match.arg(pareto_norm_scope)   # default front1
                                                       # for PROSAIL OPT

  # ── Apply norm filter ────────────────────────────────────────────────────────
  if (!is.null(norm_filter) && nzchar(norm_filter)) {
    metrics_dt <- metrics_dt[grepl(norm_filter, Norm)]
    dopt_dt    <- dopt_dt[grepl(norm_filter, Norm)]
  }

  dopt_ref <- dopt_dt[method_dopt == dopt_method, .(Site, Norm, d_opt)]

  # Normalise combined-site labels produced by script 06 ("Sites averaged",
  # "Sites combined", etc.) to the canonical "All_sites" key expected below.
  combined_labels <- c("Sites averaged", "Sites_averaged",
                       "Sites combined",  "Sites_combined")
  dopt_ref[Site %in% combined_labels, Site := "All_sites"]

  if (nrow(dopt_ref) == 0L) {
    cli::cli_abort(
      "select_prosail_opt: no rows in dopt_dt with method_dopt == {.val {dopt_method}}."
    )
  }

  # ── Internal helper: Pareto over PROSAIL columns for one (site, norm, depth) ─
  pareto_over_columns <- function(site, norm, d_val, source_label) {
    sub <- metrics_dt[Site == site & Norm == norm & Depth == d_val]
    if (nrow(sub) == 0L) {
      cli::cli_warn(
        "select_prosail_opt: no rows for Site={site} Norm={norm} Depth={d_val} ",
        "[{source_label}] — skipped."
      )
      return(NULL)
    }

    sub_pf <- data.table::copy(sub)
    sub_pf[, tmp_depth := .I]

    pf <- compute_pareto_front(
      data.table::data.table(
        Depth = sub_pf[["tmp_depth"]],
        R     = sub_pf[["R"]],
        RMSE  = sub_pf[["RMSE"]],
        Bias  = sub_pf[["Bias"]],
        Slope = sub_pf[["Slope"]]
      ),
      criteria    = pareto_criteria,
      norm_method = pareto_norm_method,
      norm_scope  = pareto_norm_scope
    )

    if (is.na(pf$d_opt)) {
      cli::cli_warn(
        "select_prosail_opt: Pareto returned NA for Site={site} Norm={norm} ",
        "[{source_label}] — skipped."
      )
      return(NULL)
    }

    best_row <- sub_pf[tmp_depth == pf$d_opt]
    n_pf     <- if (anyNA(pf$pareto_front_depths)) NA_integer_ else
                  length(pf$pareto_front_depths)

    data.table::data.table(
      d_opt_source  = source_label,
      Site          = site,
      Norm          = norm,
      d_opt         = d_val,
      method_dopt   = dopt_method,
      Column_opt    = best_row[["Column"]],
      ATBD          = best_row[["ATBD"]],
      R              = round(best_row[["R"]],     2),
      R2             = round(best_row[["R2"]],    2),
      RMSE           = round(best_row[["RMSE"]],  2),
      Bias           = round(best_row[["Bias"]],  2),
      Slope          = round(best_row[["Slope"]], 2),
      Avg_score      = pf$Avg_score_best,
      n_pareto_front = n_pf
    )
  }

  results <- list()

  # ── Mode 1: per_site — each site uses its own d_opt ─────────────────────────
  cli::cli_alert_info("select_prosail_opt: mode = per_site")
  site_dopt <- dopt_ref[Site %in% individual_sites]
  for (i in seq_len(nrow(site_dopt))) {
    row  <- site_dopt[i]
    res  <- pareto_over_columns(row$Site, row$Norm, row$d_opt, "per_site")
    if (!is.null(res)) results <- c(results, list(res))
  }

  # ── Mode 2: all_sites — All_sites d_opt applied to each individual site ─────
  cli::cli_alert_info("select_prosail_opt: mode = all_sites")
  norms_available <- unique(dopt_ref[Site == "All_sites", Norm])
  for (norm in norms_available) {
    d_all <- dopt_ref[Site == "All_sites" & Norm == norm, d_opt]
    if (length(d_all) == 0L || is.na(d_all)) next
    for (site in individual_sites) {
      res <- pareto_over_columns(site, norm, d_all, "all_sites")
      if (!is.null(res)) results <- c(results, list(res))
    }
  }

  # ── Mode 3: fixed — user-supplied depth applied to each individual site ────────
  if (!is.null(fixed_d_opt)) {
    cli::cli_alert_info(
      "select_prosail_opt: mode = fixed (d_opt = {fixed_d_opt})"
    )
    norms_available <- unique(metrics_dt[grepl(
      if (!is.null(norm_filter) && nzchar(norm_filter)) norm_filter else ".",
      Norm), Norm])
    for (norm in norms_available) {
      for (site in individual_sites) {
        res <- pareto_over_columns(site, norm, fixed_d_opt,
                                   paste0("fixed_", fixed_d_opt))
        if (!is.null(res)) results <- c(results, list(res))
      }
    }
  }

  if (length(results) == 0L) {
    cli::cli_warn("select_prosail_opt: no results — returning empty table.")
    return(data.table::data.table(
      d_opt_source = character(), Site = character(), Norm = character(),
      d_opt = integer(), method_dopt = character(), Column_opt = character(),
      ATBD = logical(), R = numeric(), R2 = numeric(), RMSE = numeric(),
      Bias = numeric(), Slope = numeric(),
      Avg_score = numeric(), n_pareto_front = integer()
    ))
  }

  data.table::rbindlist(results, use.names = TRUE)
}

# ── select_dopt ───────────────────────────────────────────────────────────────

#' @title Select optimal canopy integration depth for all (Site × Norm × Column)
#'   triplets
#'
#' @description
#' Applies four d_opt selection methods to all triplets present in the SM5
#' metrics table. For each (Site × Norm × Column) triplet, the optimal depth
#' is computed independently per method:
#' \describe{
#'   \item{pearson}{\code{d_opt = Depth[which.max(R)]}}
#'   \item{rmse}{\code{d_opt = Depth[which.min(RMSE)]}}
#'   \item{bias}{\code{d_opt = Depth[which.min(|Bias|)]}}
#'   \item{pareto}{Multi-criteria. See \code{\link{compute_pareto_front}}.}
#' }
#' The "All_sites" rows produced by SM5 PART 2 are processed alongside
#' per-site rows. The \code{Method} column (sampling strategy) is treated
#' as a pass-through field; grouping is strictly over (Site, Norm, Column).
#'
#' @param metrics_dt A \code{data.table} with the schema produced by
#'   \code{build_metrics_table()} (SM5 passe 1):
#'   Site (character), Norm (character), Depth (integer), Method (character),
#'   Column (character), R (numeric), R2 (numeric), RMSE (numeric),
#'   Bias (numeric), Slope (numeric), ATBD (logical).
#' @param methods Character vector. Selection methods to apply.
#'   Default \code{c("pearson", "rmse", "bias", "pareto")}.
#' @param max_depth Integer or \code{NULL}. Upper bound on depth considered for
#'   d_opt selection. Depths > \code{max_depth} are excluded before running any
#'   selection method. When \code{NULL} (default) all depths are used. Set to
#'   \code{h_min} to enforce the redundancy-zone constraint: pixels filtered at
#'   h_min have heights > h_min, so depths > h_min repeat the same canopy crown
#'   material and provide no additional information.
#' @param prosail_filter Character scalar or vector controlling which PROSAIL
#'   columns are included before selection:
#'   \describe{
#'     \item{"ATBD"}{(default) Keep only \code{ATBD == TRUE} rows.
#'                   Corresponds to the canonical PROSAIL baseline of the
#'                   submitted paper.}
#'     \item{"all"}{No filter. All columns present in the CSV.}
#'     \item{character vector}{Keep rows where
#'                             \code{Column \%in\% prosail_filter}.}
#'   }
#'
#' @return A \code{data.table} with columns (in order):
#' \describe{
#'   \item{Site}{character}
#'   \item{Norm}{character}
#'   \item{Column}{character}
#'   \item{ATBD}{logical}
#'   \item{method_dopt}{character — "pearson", "rmse", "bias", or "pareto"}
#'   \item{d_opt}{integer}
#'   \item{R}{numeric — value at d_opt}
#'   \item{R2}{numeric — value at d_opt}
#'   \item{RMSE}{numeric — value at d_opt}
#'   \item{Bias}{numeric — value at d_opt}
#'   \item{Slope}{numeric — value at d_opt}
#'   \item{Avg_score}{numeric — mean of Front-1-normalised criteria at d_opt
#'                    (1 = best, 0 = worst); NA unless method_dopt == "pareto"}
#'   \item{pareto_front}{character — NA unless method_dopt == "pareto";
#'                       sorted depths as "3,4,5,7"}
#' }
#'
#' @references
#' R3.minor.13, R4.spec.2: reviewers requested testing criteria alternative to
#' Pearson correlation alone for d_opt selection.
#'
#' @examples
#' mt <- data.table::data.table(
#'   Site   = "TestSite",
#'   Norm   = "DSM",
#'   Depth  = 1:3,
#'   Method = "stratified_uniform",
#'   Column = "LIDFa=1_lai=1_LMA=1_BROWN=1",
#'   R      = c(0.80, 0.90, 0.85),
#'   R2     = c(0.64, 0.81, 0.72),
#'   RMSE   = c(1.20, 0.90, 1.05),
#'   Bias   = c(0.30, 0.05, 0.15),
#'   Slope  = c(0.85, 0.98, 0.92),
#'   ATBD   = TRUE
#' )
#' select_dopt(mt, methods = c("pearson", "pareto"), prosail_filter = "ATBD")
#'
#' @export
select_dopt <- function(metrics_dt,
                        methods         = c("pearson", "rmse", "bias", "pareto"),
                        max_depth       = NULL,
                        prosail_filter  = "ATBD",
                        pareto_criteria = c("R","RMSE","Bias","Slope"),
                        pareto_norm_method = c("max", "minmax"),
                        pareto_norm_scope  = c("front1", "all")) {
  pareto_norm_method <- match.arg(pareto_norm_method)
  pareto_norm_scope  <- match.arg(pareto_norm_scope)

  # ── Filter PROSAIL columns ───────────────────────────────────────────────────
  if (identical(prosail_filter, "ATBD")) {
    metrics_f <- metrics_dt[ATBD == TRUE]
  } else if (identical(prosail_filter, "all")) {
    metrics_f <- metrics_dt
  } else {
    metrics_f <- metrics_dt[Column %in% prosail_filter]
  }

  # ── Enforce redundancy-zone constraint ──────────────────────────────────────
  # Pixels filtered at h_min have heights > h_min; depths > h_min integrate
  # the same crown material and add no new signal → exclude before selection.
  if (!is.null(max_depth)) {
    metrics_f <- metrics_f[Depth <= max_depth]
  }

  # ── Unique (Site, Norm, Column) triplets ─────────────────────────────────────
  triplets   <- unique(metrics_f[, .(Site, Norm, Column, ATBD)])
  n_triplets <- nrow(triplets)

  if (n_triplets == 0L) {
    cli::cli_warn(
      "select_dopt: no triplets after filtering — returning empty table."
    )
    return(data.table::data.table(
      Site = character(), Norm = character(), Column = character(),
      ATBD = logical(), method_dopt = character(), d_opt = integer(),
      R = numeric(), R2 = numeric(), RMSE = numeric(),
      Bias = numeric(), Slope = numeric(),
      Avg_score = numeric(), pareto_front = character()
    ))
  }

  # Pre-allocate list for all (triplet × method) combinations
  all_rows <- vector("list", n_triplets * length(methods))
  idx      <- 0L

  for (i in seq_len(n_triplets)) {
    site <- triplets[i, Site]
    norm <- triplets[i, Norm]
    col  <- triplets[i, Column]
    atbd <- triplets[i, ATBD]

    sub <- metrics_f[Site == site & Norm == norm & Column == col]

    for (method in methods) {
      idx          <- idx + 1L
      row          <- select_dopt_single(sub, method,
                                          pareto_criteria    = pareto_criteria,
                                          pareto_norm_method = pareto_norm_method,
                                          pareto_norm_scope  = pareto_norm_scope)
      row[, Site   := site]
      row[, Norm   := norm]
      row[, Column := col]
      row[, ATBD   := atbd]
      all_rows[[idx]] <- row
    }
  }

  result <- data.table::rbindlist(all_rows[seq_len(idx)], use.names = TRUE)
  data.table::setcolorder(result, c(
    "Site", "Norm", "Column", "ATBD", "method_dopt",
    "d_opt", "R", "R2", "RMSE", "Bias", "Slope",
    "Avg_score", "pareto_front"
  ))
  result
}
