# ---
# title:  sm5_dopt.R
# desc:   Selection of the optimal canopy integration depth (d_opt) for the
#         d_opt analysis (SM5). Implements four selection methods: pearson
#         (which.max R), rmse (which.min RMSE), bias (which.min |Bias|), and
#         pareto (multi-criteria distance to utopian point).
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
#' all available depths (typically 1–38), selects the depth that minimises the
#' Euclidean distance to the utopian point in a four-dimensional normalised
#' criterion space. The four criteria to minimise are:
#' \describe{
#'   \item{m1}{\code{1 - R2} (high R² is desirable)}
#'   \item{m2}{\code{RMSE}}
#'   \item{m3}{\code{abs(Bias)}}
#'   \item{m4}{\code{abs(Slope - 1)}}
#' }
#' Each criterion is normalised to \[0, 1\] across the valid depths for the
#' triplet: \code{m_norm = (m - min(m)) / (max(m) - min(m))}. A constant
#' criterion (max == min) is set to 0 for all depths. Depths with NA on any
#' of the four criteria are excluded before normalisation. The utopian point
#' is \code{(0, 0, 0, 0)}.
#'
#' In addition to the selected depth, the function returns the complete set of
#' Pareto-non-dominated depths for use in sensitivity figures (R3.minor.13,
#' R4.spec.2). A depth d1 dominates d2 if and only if
#' \code{m_norm[d1, i] <= m_norm[d2, i]} for all i, with at least one strict
#' inequality.
#'
#' @param metrics_subset A \code{data.table} with columns \code{Depth}
#'   (integer), \code{R2} (numeric), \code{RMSE} (numeric), \code{Bias}
#'   (numeric), \code{Slope} (numeric). Pre-filtered to one
#'   (Site × Norm × Column) triplet. Typically 38 rows (one per depth).
#'
#' @return A named list:
#' \describe{
#'   \item{d_opt}{Integer. Depth with minimum distance to utopian point.
#'                \code{NA_integer_} if fewer than 2 non-NA depths available.}
#'   \item{dist_min}{Numeric. Distance to utopian point at \code{d_opt}.
#'                   \code{NA_real_} if fewer than 2 non-NA depths.}
#'   \item{pareto_front_depths}{Integer vector of non-dominated depths.
#'                              \code{NA_integer_} if fewer than 2 non-NA
#'                              depths.}
#'   \item{all_distances}{Named numeric vector of utopian distances for all
#'                        non-NA depths, names are depth indices as character.
#'                        \code{NULL} if fewer than 2 non-NA depths.}
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
#' @export
compute_pareto_front <- function(metrics_subset) {
  na_result <- list(
    d_opt               = NA_integer_,
    dist_min            = NA_real_,
    pareto_front_depths = NA_integer_,
    all_distances       = NULL
  )

  depths <- metrics_subset[["Depth"]]

  # Four criteria to minimise — literal specification
  m1 <- 1 - metrics_subset[["R2"]]           # m1 = 1 - R2
  m2 <- metrics_subset[["RMSE"]]             # m2 = RMSE
  m3 <- abs(metrics_subset[["Bias"]])        # m3 = abs(Bias)
  m4 <- abs(metrics_subset[["Slope"]] - 1)  # m4 = abs(Slope - 1)

  # Exclude depths with NA on any criterion
  valid <- !is.na(m1) & !is.na(m2) & !is.na(m3) & !is.na(m4)
  if (sum(valid) < 2L) return(na_result)

  depths_v <- depths[valid]
  m_mat    <- cbind(m1[valid], m2[valid], m3[valid], m4[valid])

  # Normalise each criterion to [0, 1] locally over valid depths.
  # Constant criterion (max == min) → 0 for all depths (no discrimination).
  m_norm <- apply(m_mat, 2L, function(x) {
    lo <- min(x); hi <- max(x)
    if (hi == lo) return(rep(0, length(x)))
    (x - lo) / (hi - lo)
  })

  # Euclidean distance to utopian point (0, 0, 0, 0)
  dist_v         <- sqrt(rowSums(m_norm^2))
  names(dist_v)  <- as.character(depths_v)

  # d_opt: depth with minimum utopian distance
  idx_min  <- which.min(dist_v)
  d_opt    <- depths_v[idx_min]
  dist_min <- dist_v[[idx_min]]

  # Pareto front: non-dominated depths.
  # d1 dominates d2 iff all m_norm[d1,] <= m_norm[d2,] with >= 1 strict ineq.
  n          <- nrow(m_norm)
  dominated  <- logical(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      if (all(m_norm[j, ] <= m_norm[i, ]) &&
          any(m_norm[j, ] <  m_norm[i, ])) {
        dominated[i] <- TRUE
        break
      }
    }
  }
  front_depths <- depths_v[!dominated]

  list(
    d_opt               = d_opt,
    dist_min            = dist_min,
    pareto_front_depths = front_depths,
    all_distances       = dist_v
  )
}

# ── select_dopt_single (internal) ─────────────────────────────────────────────
#
# Applies one selection method to a (Site × Norm × Column) subset (all depths).
# Returns a one-row data.table with columns:
#   method_dopt, d_opt, R, R2, RMSE, Bias, Slope, dist_utopia, pareto_front.
# dist_utopia and pareto_front are NA for methods other than "pareto".
# pareto_front is encoded as "3,4,5,7" (comma-separated sorted depths).
#
select_dopt_single <- function(metrics_subset, method) {
  depths <- metrics_subset[["Depth"]]

  d_opt        <- NA_integer_
  dist_utopia  <- NA_real_
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
    pf    <- compute_pareto_front(metrics_subset)
    d_opt <- pf$d_opt

    dist_utopia <- if (is.na(pf$d_opt)) NA_real_ else pf$dist_min

    pareto_front <- if (anyNA(pf$pareto_front_depths)) {
      NA_character_
    } else {
      paste(sort(pf$pareto_front_depths), collapse = ",")
    }
  }

  # Extract metrics at d_opt; all NA if d_opt itself is NA
  if (is.na(d_opt)) {
    r_val     <- NA_real_
    r2_val    <- NA_real_
    rmse_val  <- NA_real_
    bias_val  <- NA_real_
    slope_val <- NA_real_
  } else {
    row_d     <- metrics_subset[Depth == d_opt]
    r_val     <- row_d[["R"]][[1L]]
    r2_val    <- row_d[["R2"]][[1L]]
    rmse_val  <- row_d[["RMSE"]][[1L]]
    bias_val  <- row_d[["Bias"]][[1L]]
    slope_val <- row_d[["Slope"]][[1L]]
  }

  data.table::data.table(
    method_dopt  = method,
    d_opt        = d_opt,
    R            = r_val,
    R2           = r2_val,
    RMSE         = rmse_val,
    Bias         = bias_val,
    Slope        = slope_val,
    dist_utopia  = dist_utopia,
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
#'   \item{dist_utopia}{numeric — Pareto distance for Column_opt}
#'   \item{n_pareto_front}{integer — number of non-dominated columns}
#' }
#'
#' @references
#' R3.minor.13, R4.spec.2: multi-criteria PROSAIL configuration selection.
#'
#' @export
select_prosail_opt <- function(metrics_dt,
                                dopt_dt,
                                dopt_method      = "pareto",
                                norm_filter      = "keepTrees",
                                individual_sites = c("Aigoual", "Blois", "Mormal"),
                                fixed_d_opt      = NULL) {

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

    pf <- compute_pareto_front(data.table::data.table(
      Depth = sub_pf[["tmp_depth"]],
      R2    = sub_pf[["R2"]],
      RMSE  = sub_pf[["RMSE"]],
      Bias  = sub_pf[["Bias"]],
      Slope = sub_pf[["Slope"]]
    ))

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
      R             = best_row[["R"]],
      R2            = best_row[["R2"]],
      RMSE          = best_row[["RMSE"]],
      Bias          = best_row[["Bias"]],
      Slope         = best_row[["Slope"]],
      dist_utopia   = pf$dist_min,
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
      dist_utopia = numeric(), n_pareto_front = integer()
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
#'   \item{dist_utopia}{numeric — NA unless method_dopt == "pareto"}
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
                        methods        = c("pearson", "rmse", "bias", "pareto"),
                        max_depth      = NULL,
                        prosail_filter = "ATBD") {

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
      dist_utopia = numeric(), pareto_front = character()
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
      row          <- select_dopt_single(sub, method)
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
    "dist_utopia", "pareto_front"
  ))
  result
}
