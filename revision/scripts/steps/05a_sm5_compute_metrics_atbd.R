# ---
# title:  05a_sm5_compute_metrics_atbd.R
# desc:   Orchestration — computes the metrics table (R², RMSE, Bias, Slope)
#         for the ATBD-only inversion across all depths, sites, norms and h_min
#         values. Output feeds directly into 06 (d_opt selection).
#
#         Reuses build_metrics_table() from revision/R/sm5_aggregate.R, called
#         with the ATBD-specific CSV filenames produced by 04a.
#
# Prerequisites:
#   03b — PAD CSVs at 03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   04a — LAI_estimated_atbd_*.csv at
#         revision/output/intermediate/PROSAIL_Models/{site}/
#         LIDFa_lai_LMA_BROWN/atbd/
#
# Output:
#   revision/output/intermediate/sm5/
#   all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05a_sm5_compute_metrics_atbd.R")
# ---

library("here")
library("data.table")
library("cli")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "sm5_metrics.R"))
source(here::here("revision", "R", "sm5_aggregate.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites            <- c("Aigoual", "Blois", "Mormal")
norm_methods     <- c(
                      # "DSM", "DTM", 
                      "DSM_keepTrees", "DTM_keepTrees"
                      )
depths           <- 1:38
sampling_methods <- "stratified_uniform"
name_strategy    <- "LIDFa_lai_LMA_BROWN"
parms2test       <- c("LIDFa", "lai", "LMA", "BROWN")
nb_samples       <- 5000L
h_min_values     <- c(10L, 15L, 20L)

# ── Output ─────────────────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "intermediate", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_dir,
                     paste0("all_results_atbd_", name_strategy, ".csv"))

# ── Build table (ATBD LAI estimated CSVs, atbd/ sub-directory) ────────────────
# read_sampling_and_estimated() is called with a custom s2_csv path that
# points to the atbd/ sub-directory produced by 04a.

build_metrics_table_atbd <- function(sites, norm_methods, depths,
                                      sampling_methods, name_strategy,
                                      parms2test, nb_samples, h_min_values) {

  n_params <- length(parms2test)
  all_rows <- list()

  cli::cli_h2("SM5a PART 1 — per-site (ATBD)")

  for (site in sites) {
    cli::cli_alert_info("Site: {site}")
    for (norm in norm_methods) {
      for (depth in depths) {
        for (sampling_method in sampling_methods) {
          for (h_min in h_min_values) {

            lidar_csv <- file.path(
              paths$ext_results, site, "PROSAIL_Optimization", "sampling",
              paste0("PAD_", norm, "_Depth_", depth,
                     "_Samples_", sampling_method,
                     "_hmin", h_min,
                     "_nbSamples_", nb_samples, ".csv")
            )
            s2_csv <- here::here(
              "revision", "output", "intermediate", "PROSAIL_Models",
              site, name_strategy, "atbd",
              paste0("LAI_estimated_atbd_", sampling_method,
                     "_hmin", h_min,
                     "_nbSamples_", nb_samples, ".csv")
            )

            if (!file.exists(lidar_csv) || !file.exists(s2_csv)) next

            lidar_dt <- data.table::fread(lidar_csv, header = TRUE, sep = "\t")
            s2_dt    <- data.table::fread(s2_csv,    header = TRUE, sep = "\t")

            s2_ids        <- seq_len(nrow(s2_dt))
            lidar_aligned <- lidar_dt[["lidar_values"]][
              match(s2_ids, lidar_dt[["samples_id"]])
            ]

            col <- "LAI_atbd"
            valid <- !is.na(lidar_aligned) & !is.na(s2_dt[[col]])
            m <- compute_metrics_s2_lidar(
              s2_vals    = s2_dt[[col]][valid],
              lidar_vals = lidar_aligned[valid]
            )
            all_rows <- c(all_rows, list(data.table::data.table(
              Site   = site,
              Norm   = norm,
              Depth  = depth,
              Method = sampling_method,
              h_min  = h_min,
              Column = col,
              R      = round(m$R,     2L),
              R2     = round(m$R2,    2L),
              RMSE   = round(m$RMSE,  2L),
              Bias   = round(m$Bias,  2L),
              Slope  = round(m$Slope, 2L),
              ATBD   = TRUE
            )))
          }
        }
      }
    }
  }

  cli::cli_h2("SM5a PART 2 — All_sites aggregated (ATBD)")

  for (norm in norm_methods) {
    for (depth in depths) {
      for (sampling_method in sampling_methods) {
        for (h_min in h_min_values) {

          agg_lidar <- data.table::data.table()
          agg_s2    <- data.table::data.table()

          for (site in sites) {
            lidar_csv <- file.path(
              paths$ext_results, site, "PROSAIL_Optimization", "sampling",
              paste0("PAD_", norm, "_Depth_", depth,
                     "_Samples_", sampling_method,
                     "_hmin", h_min,
                     "_nbSamples_", nb_samples, ".csv")
            )
            s2_csv <- here::here(
              "revision", "output", "intermediate", "PROSAIL_Models",
              site, name_strategy, "atbd",
              paste0("LAI_estimated_atbd_", sampling_method,
                     "_hmin", h_min,
                     "_nbSamples_", nb_samples, ".csv")
            )
            if (!file.exists(lidar_csv) || !file.exists(s2_csv)) next

            lidar_dt <- data.table::fread(lidar_csv, header = TRUE, sep = "\t")
            s2_dt    <- data.table::fread(s2_csv,    header = TRUE, sep = "\t")

            s2_dt[, samples_id := .I]
            lidar_dt <- lidar_dt[samples_id %in% s2_dt[["samples_id"]]]
            s2_dt    <- s2_dt[samples_id %in% lidar_dt[["samples_id"]]]
            data.table::setkey(lidar_dt, samples_id)
            data.table::setkey(s2_dt,    samples_id)

            agg_lidar <- data.table::rbindlist(list(agg_lidar, lidar_dt))
            agg_s2    <- data.table::rbindlist(list(agg_s2,    s2_dt))
          }

          if (nrow(agg_lidar) == 0L) next

          col   <- "LAI_atbd"
          valid <- stats::complete.cases(agg_lidar[["lidar_values"]],
                                         agg_s2[[col]])
          m <- compute_metrics_s2_lidar(
            s2_vals    = agg_s2[[col]][valid],
            lidar_vals = agg_lidar[["lidar_values"]][valid]
          )
          all_rows <- c(all_rows, list(data.table::data.table(
            Site   = "All_sites",
            Norm   = norm,
            Depth  = depth,
            Method = sampling_method,
            h_min  = h_min,
            Column = col,
            R      = round(m$R,     2L),
            R2     = round(m$R2,    2L),
            RMSE   = round(m$RMSE,  2L),
            Bias   = round(m$Bias,  2L),
            Slope  = round(m$Slope, 2L),
            ATBD   = TRUE
          )))
        }
      }
    }
  }

  if (length(all_rows) == 0L) {
    cli::cli_warn("SM5a: no data found — check 04a outputs.")
    return(data.table::data.table())
  }
  data.table::rbindlist(all_rows, use.names = TRUE)
}

# ── Run ────────────────────────────────────────────────────────────────────────

t0 <- proc.time()

combined <- build_metrics_table_atbd(
  sites            = sites,
  norm_methods     = norm_methods,
  depths           = depths,
  sampling_methods = sampling_methods,
  name_strategy    = name_strategy,
  parms2test       = parms2test,
  nb_samples       = nb_samples,
  h_min_values     = h_min_values
)

elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

data.table::fwrite(combined, out_csv)

cli::cli_h1("SM5a — résumé")
cli::cli_bullets(c(
  "v" = "Fichier écrit : {out_csv}",
  "i" = "Lignes totales : {nrow(combined)}",
  "i" = "Durée : {elapsed} min"
))
