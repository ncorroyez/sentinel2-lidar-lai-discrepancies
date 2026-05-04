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
#' @param h_min           Integer. Minimum canopy height filter in metres
#'   (e.g. 10L, 15L, 20L). Used to build h_min-specific filenames.
#'
#' @return A named list \code{list(lidar_dt, s2_dt)} of \code{data.table}s,
#'   or \code{NULL} if either file is absent.
#'
#' @export
read_sampling_and_estimated <- function(site, norm, depth,
                                         sampling_method, name_strategy,
                                         nb_samples, h_min = 10L,
                                         lai_scenario = "per_site") {
  lidar_csv <- file.path(
    paths$ext_results, site, "PROSAIL_Optimization", "sampling",
    paste0("PAD_", norm, "_Depth_", depth,
           "_Samples_", sampling_method,
           "_hmin", h_min,
           "_nbSamples_", nb_samples, ".csv")
  )
  s2_csv <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy, lai_scenario,
    paste0("LAI_estimated_", lai_scenario, "_", sampling_method,
           "_hmin", h_min,
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
#' @param h_min           Integer. Minimum canopy height filter in metres.
#'   Default \code{10L}. Added as a column \code{h_min} in the output.
#' @param lai_scenario    Character. LAI source scenario for SVR training
#'   (\code{"per_site"} or \code{"common"}). Default \code{"per_site"}.
#'   Added as a column \code{lai_scenario} in the output.
#' @param lidar_scale     Numeric. Multiplicative scale applied to all LiDAR
#'   PAD values before computing metrics. Use \code{(k_ref / k_select) *
#'   cos(theta_select * pi / 180)} when k or scan angle differ from the
#'   reference embedded in the PAD CSVs (k_ref = 0.5, theta = 0°). Default
#'   \code{1} (no rescaling).
#'
#' @return A \code{data.table} with columns (in order):
#'   Site, Norm, Depth, Method, h_min, lai_scenario, Column, R, R2, RMSE,
#'   Bias, Slope, ATBD.
#'   Types: character / character / integer / character / integer / character /
#'          character / numeric × 5 / logical.
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
    nb_samples       = 5000L,
    h_min            = 10L,
    lai_scenario     = "per_site",
    lidar_scale      = 1) {

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
            site, norm, depth, sampling_method, name_strategy, nb_samples,
            h_min = h_min, lai_scenario = lai_scenario
          )
          if (is.null(dat)) next

          lidar_dt <- dat$lidar_dt
          s2_dt    <- dat$s2_dt

          # Align by samples_id: S2 rows are in GPKG feature order (row index =
          # implicit samples_id 1..n). Reorder lidar_values to match S2 row order
          # via match(), making alignment explicit and robust to any CSV reordering.
          s2_ids         <- seq_len(nrow(s2_dt))
          lidar_aligned  <- lidar_dt[["lidar_values"]][
            match(s2_ids, lidar_dt[["samples_id"]])
          ] * lidar_scale
          if (anyNA(lidar_aligned)) {
            warning("SM5: unmatched samples_id — site=", site, " norm=", norm,
                    " depth=", depth, ". Rows with NA lidar dropped.")
          }

          chunk <- data.table::rbindlist(lapply(colnames(s2_dt), function(col) {
            valid <- !is.na(lidar_aligned) & !is.na(s2_dt[[col]])
            m <- compute_metrics_s2_lidar(
              s2_vals    = s2_dt[[col]][valid],
              lidar_vals = lidar_aligned[valid]
            )
            data.table::data.table(
              Site         = site,
              Norm         = norm,
              Depth        = depth,
              Method       = sampling_method,
              h_min        = h_min,
              lai_scenario = lai_scenario,
              Column       = col,
              R            = round(m$R,     3L),
              R2           = round(m$R2,    3L),
              RMSE         = round(m$RMSE,  3L),
              Bias         = round(m$Bias,  3L),
              Slope        = round(m$Slope, 3L),
              ATBD         = detect_atbd(col, n_params)
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
            site, norm, depth, sampling_method, name_strategy, nb_samples,
            h_min = h_min, lai_scenario = lai_scenario
          )
          if (is.null(dat)) next

          lidar_dt <- dat$lidar_dt
          s2_dt    <- dat$s2_dt

          # Align by samples_id: assign row-index as samples_id to S2, then
          # keep only rows present in LiDAR (inner join on samples_id).
          s2_dt[, samples_id := .I]
          lidar_dt <- lidar_dt[samples_id %in% s2_dt[["samples_id"]]]
          s2_dt    <- s2_dt[samples_id %in% lidar_dt[["samples_id"]]]
          # Reorder both tables to identical samples_id order before stacking.
          data.table::setkey(lidar_dt, samples_id)
          data.table::setkey(s2_dt,    samples_id)

          agg_lidar <- data.table::rbindlist(
            list(agg_lidar, lidar_dt)
          )
          agg_s2 <- data.table::rbindlist(
            list(agg_s2, s2_dt)
          )
        }

        if (nrow(agg_lidar) == 0L) next

        s2_cols <- setdiff(colnames(agg_s2), "samples_id")

        chunk <- data.table::rbindlist(lapply(s2_cols, function(col) {
          vals_s2 <- agg_s2[[col]]

          # Pre-filter complete cases (mirrors 04C PART 2 lines 183-184)
          valid <- stats::complete.cases(agg_lidar[["lidar_values"]], vals_s2)

          m <- compute_metrics_s2_lidar(
            s2_vals    = vals_s2[valid],
            lidar_vals = agg_lidar[["lidar_values"]][valid] * lidar_scale
          )
          data.table::data.table(
            Site         = "All_sites",
            Norm         = norm,
            Depth        = depth,
            Method       = sampling_method,
            h_min        = h_min,
            lai_scenario = lai_scenario,
            Column       = col,
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
      Method = character(), h_min = integer(), lai_scenario = character(),
      Column = character(), R = numeric(), R2 = numeric(), RMSE = numeric(),
      Bias = numeric(), Slope = numeric(), ATBD = logical()
    ))
  }

  data.table::rbindlist(all_rows, use.names = TRUE)
}
