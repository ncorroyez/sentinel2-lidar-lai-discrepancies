# ---
# title:  sm5_aggregate.R
# desc:   Orchestration of the metrics table build for the d_opt analysis (SM5).
#         Exposes build_metrics_table() and the loader helper
#         read_sampling_and_estimated().
#
#         Faithful refactor of Main_s2_04C_analysis_depth.R (December 2025,
#         PROSAIL-Optimization/02_CODES/Sentinel2/) with the following
#         intentional changes vs legacy:
#           - R2 column added (not in legacy).
#           - name_strategy cleaned ("LIDFa_lai_LMA_BROWN", no _Agg_10m suffix).
#           - Four norm_methods active instead of only DSM_keepTrees.
#           - S2 LAI CSV read from SM4 output (revision/output/intermediate/...)
#             not from 03_RESULTS/...PROSAIL_Models/ (legacy path).
#           - Single fwrite() at script level; no double-write with append.
#           - List accumulation instead of repeated rbindlist() for performance.
#           - Known incohérence samples_id PART 1 vs PART 2 reproduced faithfully.
# ---

# ── detect_atbd (internal) ────────────────────────────────────────────────────
#
# Counts how many parameters in a column name are at their first-level index
# (the ATBD baseline), and returns TRUE if all parms2test are at level 1.
#
# Legacy 04C (line 101) used stringr::str_count(colname, "\\.1(\\D|$)") against
# the old dot-separated format "LIDFa.1_lai.1_LMA.1_BROWN.1".
# SM4 refactored column names use "=" separator: "LIDFa=1_lai=1_LMA=1_BROWN=1".
# Separator is hardcoded to "=" (SM4 output format).
#
detect_atbd <- function(colname, n_params) {
  pattern  <- "=1(\\D|$)"
  m        <- gregexpr(pattern, colname, perl = TRUE)[[1L]]
  n_matches <- if (m[[1L]] == -1L) 0L else length(m)
  isTRUE(n_matches == n_params)
}

# ── read_sampling_and_estimated ────────────────────────────────────────────────

#' @title Load LiDAR PAD CSV and S2 LAI estimated CSV for one combination
#'
#' @description
#' Loads the two input files required to compute metrics for a single
#' (site × norm × depth × sampling_method) combination.
#'
#' \strong{LiDAR PAD CSV} (SM3 output):
#' \code{03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#' PAD_{norm}_Depth_{depth}_Samples_{method}_nbSamples_{n}.csv}
#' Columns: \code{lidar_values} (numeric PAD sum), \code{samples_id} (integer).
#'
#' \strong{S2 LAI estimated CSV} (SM4 output):
#' \code{revision/output/intermediate/PROSAIL_Models/{site}/{strategy}/
#' LAI_estimated_{method}_nbSamples_{n}.csv}
#' Columns: one per PROSAIL combination, named \code{LIDFa=X_lai=Y_LMA=Z_BROWN=W}.
#'
#' Note: the S2 LAI path differs from the legacy 04C path
#' (\code{03_RESULTS/{site}/PROSAIL_Optimization/PROSAIL_Models/}).
#'
#' @param site            Character. Study site name.
#' @param norm            Character. Normalisation type (e.g. "DSM_keepTrees").
#' @param depth           Integer. Canopy depth index (1–38).
#' @param sampling_method Character. Sampling method (e.g. "stratified_uniform").
#' @param name_strategy   Character. LUT strategy name.
#' @param nb_samples      Integer. Number of samples (e.g. 5000).
#'
#' @return A named list \code{list(lidar_dt, s2_dt)} of \code{data.table}s,
#'   or \code{NULL} if either file is absent.
#'
#' @export
read_sampling_and_estimated <- function(site, norm, depth,
                                         sampling_method, name_strategy,
                                         nb_samples) {
  lidar_csv <- here::here(
    "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
    paste0("PAD_", norm, "_Depth_", depth,
           "_Samples_", sampling_method,
           "_nbSamples_", nb_samples, ".csv")
  )
  s2_csv <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy,
    paste0("LAI_estimated_", sampling_method,
           "_nbSamples_", nb_samples, ".csv")
  )

  if (!file.exists(lidar_csv) || !file.exists(s2_csv)) return(NULL)

  list(
    lidar_dt = data.table::fread(lidar_csv, header = TRUE, sep = "\t"),
    s2_dt    = data.table::fread(s2_csv,    header = TRUE, sep = "\t")
  )
}

# ── build_metrics_table ────────────────────────────────────────────────────────

#' @title Build full metrics table for the d_opt analysis (SM5)
#'
#' @description
#' Iterates over all site × norm × depth × sampling_method × PROSAIL column
#' combinations and computes R, R2, RMSE, Bias, Slope between the S2 LAI
#' estimate and the LiDAR PAD-based LAI. Returns a single \code{data.table}
#' covering both per-site rows and "All_sites" aggregated rows.
#'
#' \strong{PART 1 — per-site} (mirrors 04C lines 63–125):
#' For each (site × norm × depth × method), loads the two CSVs and computes
#' metrics for every PROSAIL column. No samples_id alignment.
#'
#' \strong{PART 2 — All_sites} (mirrors 04C lines 140–211):
#' Aggregates the three sites' data for each (norm × depth × method), with
#' a samples_id alignment that generates a spurious \code{Column = "samples_id"}
#' row per (norm × depth × method). Behaviour reproduced faithfully in passe 1;
#' correction planned for passe 2.
#'
#' @param sites           Character vector. Sites to process.
#'   Default \code{c("Aigoual", "Blois", "Mormal")}.
#' @param norm_methods    Character vector. Normalisation types.
#'   Default \code{c("DSM", "DTM", "DSM_keepTrees", "DTM_keepTrees")}.
#'   Legacy 04C had only \code{"DSM_keepTrees"} active; all four activated here.
#' @param depths          Integer vector. Depth indices. Default \code{1:38}.
#' @param sampling_methods Character vector. Default \code{"stratified_uniform"}.
#' @param name_strategy   Character. LUT strategy. Default \code{"LIDFa_lai_LMA_BROWN"}.
#' @param parms2test      Character vector. Parameter names for ATBD detection.
#'   Default \code{c("LIDFa", "lai", "LMA", "BROWN")}.
#' @param nb_samples      Integer. Number of samples. Default \code{5000L}.
#'
#' @return A \code{data.table} with columns (in order):
#'   Site, Norm, Depth, Method, Column, R, R2, RMSE, Bias, Slope, ATBD.
#'   Types: character / character / integer / character / character /
#'          numeric × 5 / logical.
#'
#' @references
#' Legacy: PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_04C_analysis_depth.R
#'
#' @export
build_metrics_table <- function(
    sites            = c("Aigoual", "Blois", "Mormal"),
    norm_methods     = c("DSM", "DTM", "DSM_keepTrees", "DTM_keepTrees"),
    depths           = 1:38,
    sampling_methods = "stratified_uniform",
    name_strategy    = "LIDFa_lai_LMA_BROWN",
    parms2test       = c("LIDFa", "lai", "LMA", "BROWN"),
    nb_samples       = 5000L) {

  n_params  <- length(parms2test)
  all_rows  <- list()   # accumulate dt chunks; single rbindlist at the end

  # ── PART 1: per-site (mirrors 04C lines 63–125) ─────────────────────────────
  cli::cli_h2("SM5 PART 1 — per-site metrics")

  for (site in sites) {
    cli::cli_alert_info("Site: {site}")
    for (norm in norm_methods) {
      for (depth in depths) {
        for (sampling_method in sampling_methods) {

          dat <- read_sampling_and_estimated(
            site, norm, depth, sampling_method, name_strategy, nb_samples
          )
          if (is.null(dat)) next

          lidar_dt <- dat$lidar_dt
          s2_dt    <- dat$s2_dt

          # Row-count guard (mirrors 04C line 91-93)
          if (nrow(lidar_dt) != nrow(s2_dt)) {
            stop("SM5: row mismatch — site=", site, " norm=", norm,
                 " depth=", depth)
          }

          # NOTE: no samples_id alignment here.
          # Mirrors 04C PART 1 where alignment is commented out (lines 87-89).
          # Known incohérence vs PART 2 — correction planned for passe 1.5.

          chunk <- data.table::rbindlist(lapply(colnames(s2_dt), function(col) {
            m <- compute_metrics_s2_lidar(
              s2_vals    = s2_dt[[col]],
              lidar_vals = lidar_dt[["lidar_values"]]
            )
            data.table::data.table(
              Site   = site,
              Norm   = norm,
              Depth  = depth,
              Method = sampling_method,
              Column = col,
              R      = round(m$R,     3L),
              R2     = round(m$R2,    3L),
              RMSE   = round(m$RMSE,  3L),
              Bias   = round(m$Bias,  3L),
              Slope  = round(m$Slope, 3L),
              ATBD   = detect_atbd(col, n_params)
            )
          }))
          all_rows <- c(all_rows, list(chunk))
        }
      }
    }
  }

  # ── PART 2: All_sites aggregated (mirrors 04C lines 140–211) ─────────────────
  cli::cli_h2("SM5 PART 2 — All_sites aggregated metrics")

  for (norm in norm_methods) {
    for (depth in depths) {
      for (sampling_method in sampling_methods) {

        agg_lidar <- data.table::data.table()
        agg_s2    <- data.table::data.table()

        for (site in sites) {
          dat <- read_sampling_and_estimated(
            site, norm, depth, sampling_method, name_strategy, nb_samples
          )
          if (is.null(dat)) next

          lidar_dt <- dat$lidar_dt
          s2_dt    <- dat$s2_dt

          # Align by samples_id — mirrors 04C PART 2 lines 167-168.
          # Creates a synthetic samples_id column (row index) in s2_dt that
          # persists into agg_s2, generating a spurious "samples_id" row in
          # the metrics output. Reproduced faithfully.
          s2_dt[, samples_id := .I]
          s2_dt <- s2_dt[samples_id %in% lidar_dt[["samples_id"]]]

          agg_lidar <- data.table::rbindlist(
            list(agg_lidar, lidar_dt)
          )
          agg_s2 <- data.table::rbindlist(
            list(agg_s2, s2_dt)
          )
        }

        if (nrow(agg_lidar) == 0L) next

        chunk <- data.table::rbindlist(lapply(colnames(agg_s2), function(col) {
          vals_s2 <- agg_s2[[col]]

          # Pre-filter complete cases (mirrors 04C PART 2 lines 183-184)
          valid <- stats::complete.cases(agg_lidar[["lidar_values"]], vals_s2)

          m <- compute_metrics_s2_lidar(
            s2_vals    = vals_s2[valid],
            lidar_vals = agg_lidar[["lidar_values"]][valid]
          )
          data.table::data.table(
            Site   = "All_sites",
            Norm   = norm,
            Depth  = depth,
            Method = sampling_method,
            Column = col,
            R      = round(m$R,     3L),
            R2     = round(m$R2,    3L),
            RMSE   = round(m$RMSE,  3L),
            Bias   = round(m$Bias,  3L),
            Slope  = round(m$Slope, 3L),
            ATBD   = detect_atbd(col, n_params)
          )
        }))
        all_rows <- c(all_rows, list(chunk))
      }
    }
  }

  # ── Combine all chunks ────────────────────────────────────────────────────────
  if (length(all_rows) == 0L) {
    cli::cli_warn("SM5: no data — all input files missing. Returning empty table.")
    return(data.table::data.table(
      Site = character(), Norm = character(), Depth = integer(),
      Method = character(), Column = character(),
      R = numeric(), R2 = numeric(), RMSE = numeric(),
      Bias = numeric(), Slope = numeric(), ATBD = logical()
    ))
  }

  data.table::rbindlist(all_rows, use.names = TRUE)
}
