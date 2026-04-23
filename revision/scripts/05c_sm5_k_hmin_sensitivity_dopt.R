# ---
# title:  05c_sm5_k_hmin_sensitivity_dopt.R
# desc:   Sensitivity of Pareto d_opt to the extinction coefficient k
#         (0.4, 0.5, 0.6) and to the minimum height threshold h_min
#         (10, 15, 20 m).
#
#         k sensitivity: Beer-Lambert gives LAI_ALS proportional to 1/k,
#         so LAI_ALS_k = LAI_ALS_{k_ref} x (k_ref / k). Metrics are
#         recomputed for each k by scaling the existing PAD sampling CSVs;
#         S2 estimates are unchanged. Each PAD file is read once and metrics
#         are computed for all k values in the same pass.
#
#         h_min sensitivity: the existing all_results_atbd CSV already covers
#         h_min in {10, 15, 20}. select_dopt is applied per h_min with
#         max_depth = h_min (redundancy-zone constraint).
#
# Prerequisites:
#   03b — PAD CSVs at 03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   04a — ATBD LAI estimated CSVs
#   05a — all_results_atbd_LIDFa_lai_LMA_BROWN.csv (h_min sensitivity only)
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05c_sm5_k_hmin_sensitivity_dopt.R")
# ---

library("here")
library("data.table")
library("cli")

source(here::here("revision", "R", "sm5_metrics.R"))
source(here::here("revision", "R", "sm5_dopt.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites           <- c("Aigoual", "Blois", "Mormal")
norms_select    <- c("DSM_keepTrees", "DTM_keepTrees")
depths          <- 1:38
sampling_method <- "stratified_uniform"
name_strategy   <- "LIDFa_lai_LMA_BROWN"
nb_samples      <- 5000L
h_min_values    <- c(10L, 15L, 20L)
k_values        <- c(0.4, 0.5, 0.6)
k_ref           <- 0.5          # k used when computing the existing PAD CSVs

# ── Part 1: k sensitivity — recompute metrics for all k in one file pass ───────

cli::cli_h1("Part 1: k sensitivity")

all_rows <- list()

for (hm in h_min_values) {
  cli::cli_alert_info("h_min = {hm} m")
  for (site in sites) {
    for (norm in norms_select) {

      s2_csv <- here::here(
        "revision", "output", "intermediate", "PROSAIL_Models",
        site, name_strategy, "atbd",
        paste0("LAI_estimated_atbd_", sampling_method,
               "_hmin", hm, "_nbSamples_", nb_samples, ".csv")
      )
      if (!file.exists(s2_csv)) next
      s2_dt <- data.table::fread(s2_csv, header = TRUE, sep = "\t")

      for (depth in depths) {
        lidar_csv <- here::here(
          "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
          paste0("PAD_", norm, "_Depth_", depth,
                 "_Samples_", sampling_method,
                 "_hmin", hm, "_nbSamples_", nb_samples, ".csv")
        )
        if (!file.exists(lidar_csv)) next

        lidar_dt <- data.table::fread(lidar_csv, header = TRUE, sep = "\t")

        # Align by sample ID (same logic as 05a)
        s2_ids   <- seq_len(nrow(s2_dt))
        lidar_ref <- lidar_dt[["lidar_values"]][
          match(s2_ids, lidar_dt[["samples_id"]])
        ]

        # Compute metrics for all k values using the same aligned pair
        for (k_val in k_values) {
          lidar_k <- lidar_ref * (k_ref / k_val)
          valid   <- !is.na(lidar_k) & !is.na(s2_dt[["LAI_atbd"]])
          if (sum(valid) < 5L) next
          m <- compute_metrics_s2_lidar(
            s2_vals    = s2_dt[["LAI_atbd"]][valid],
            lidar_vals = lidar_k[valid]
          )
          all_rows <- c(all_rows, list(data.table::data.table(
            k      = k_val,
            h_min  = hm,
            Site   = site,
            Norm   = norm,
            Depth  = as.integer(depth),
            R      = round(m$R,     3L),
            R2     = round(m$R2,    3L),
            RMSE   = round(m$RMSE,  3L),
            Bias   = round(m$Bias,  3L),
            Slope  = round(m$Slope, 3L),
            ATBD   = TRUE,
            Column = "LAI_atbd"
          )))
        }
      }
    }
  }
}

k_dt <- data.table::rbindlist(all_rows, use.names = TRUE)
cli::cli_alert_success("Metrics computed: {nrow(k_dt)} rows ({nrow(k_dt) / length(k_values) / length(depths)} site x norm x h_min combos)")

# ── d_opt (Pareto) for each (k x h_min) ───────────────────────────────────────

dopt_k_rows <- list()
for (k_val in k_values) {
  for (hm in h_min_values) {
    sub <- k_dt[k == k_val & h_min == hm & Norm %in% norms_select]
    if (nrow(sub) == 0L) next
    d <- select_dopt(sub, methods = "pareto", max_depth = hm,
                     prosail_filter = "ATBD")
    d[, k     := k_val]
    d[, h_min := hm]
    dopt_k_rows <- c(dopt_k_rows, list(d[method_dopt == "pareto"]))
  }
}
dopt_k_dt <- data.table::rbindlist(dopt_k_rows, use.names = TRUE, fill = TRUE)

# Wide: one column per k value, rows = Site x Norm x h_min
wide_k <- data.table::dcast(
  dopt_k_dt,
  Site + Norm + h_min ~ paste0("k=", k),
  value.var = "d_opt"
)
data.table::setorder(wide_k, h_min, Norm, Site)

cli::cli_h2("d_opt (Pareto) by k | max_depth = h_min")
print(wide_k)

# ── Part 2: h_min sensitivity — from existing metrics CSV ─────────────────────

cli::cli_h1("Part 2: h_min sensitivity")

atbd_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
)
if (!file.exists(atbd_csv))
  stop("Run 05a first:\n  ", atbd_csv)

metrics_dt <- data.table::fread(atbd_csv)
metrics_dt <- metrics_dt[ATBD == TRUE & Norm %in% norms_select &
                           h_min %in% h_min_values]

dopt_hmin_rows <- list()
for (hm in h_min_values) {
  sub <- metrics_dt[h_min == hm]
  d   <- select_dopt(sub, methods = "pareto", max_depth = hm,
                     prosail_filter = "ATBD")
  d[, h_min := hm]
  dopt_hmin_rows <- c(dopt_hmin_rows, list(d[method_dopt == "pareto"]))
}
dopt_hmin_dt <- data.table::rbindlist(dopt_hmin_rows, use.names = TRUE, fill = TRUE)

wide_hmin <- data.table::dcast(
  dopt_hmin_dt,
  Site + Norm ~ paste0("hmin=", h_min),
  value.var = "d_opt"
)
data.table::setorder(wide_hmin, Norm, Site)

cli::cli_h2("d_opt (Pareto) by h_min | max_depth = h_min, k = 0.5")
print(wide_hmin)

# ── Joint summary: k x h_min ──────────────────────────────────────────────────

cli::cli_h1("Joint summary: d_opt (Pareto) by k x h_min")

wide_joint <- data.table::dcast(
  dopt_k_dt,
  Site + Norm ~ paste0("k=", k, "_hmin=", h_min),
  value.var = "d_opt"
)
data.table::setorder(wide_joint, Norm, Site)
print(wide_joint)

cat("Done.\n")
