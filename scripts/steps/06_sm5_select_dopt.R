# ---
# title:  06_sm5_select_dopt.R
# desc:   Orchestration — selects the optimal canopy integration depth (d_opt)
#         by four methods (pearson, rmse, bias, pareto) and writes two CSV
#         files to output/intermediate/sm5/:
#           - dopt_reference.csv  (ATBD == TRUE columns only)
#           - dopt_all.csv        (all PROSAIL columns)
#
#         Calls select_dopt() from R/sm5_dopt.R, which addresses
#         reviewer comments R3.minor.13 and R4.spec.2.
#
#         Prerequisites:
#           SM5 passe 1 — CSV must exist at:
#             output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#           Run scripts/05_sm5_compute_metrics.R first.
#
# Run from the project root (NC_Full/):
#   source("scripts/06_sm5_select_dopt.R")
# ---

# ── Source ─────────────────────────────────────────────────────────────────────

library(here)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_dopt.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites_individual <- c("Aigoual", "Blois", "Mormal")

# h_min to analyse. Ignored when k_select/theta_select are non-reference
# (05c CSV only has h_min = 10).
h_min_select  <- 10L

# k_select = 0.5 and theta_select = 0 → use 05a's CSV
# any other (k, theta)                → use 05c's kt CSV (h_min = 10 only)
k_ref         <- 0.5
k_select      <- 0.65
theta_select  <- 0

# "pool"    — use pooled-pixel combined row ("Sites combined")
# "average" — compute site-averaged metrics post-hoc ("Sites averaged")
combined_mode <- "average"

# Pareto criteria — subset of c("R", "RMSE", "Bias", "Slope").
# Default uses all 4. Drop "Slope" if you want |Slope-1| to NOT influence
# Pareto selection (avoids picking compromise configs that have small
# |Slope-1| but large Bias).
pareto_criteria <- c("R", "RMSE", "Bias", "Slope")
# pareto_criteria <- c("R", "RMSE", "Bias")

# Pareto normalisation method — "max" (default) or "minmax".
#   "max"    : m / max(m)              → preserves absolute ratios; a narrow-
#                                         range criterion stays narrow.
#   "minmax" : (m - min) / (max - min) → strict [0,1] but AMPLIFIES narrow
#                                         ranges (a 0.03 gap on R becomes as
#                                         wide as a 1.0 gap on RMSE).
pareto_norm_method <- "minmax"
# pareto_norm_method <- "max"

# ── Load data ──────────────────────────────────────────────────────────────────

use_reference <- isTRUE(k_select == k_ref && theta_select == 0)

if (use_reference) {
  sm5_csv <- here::here(
    "output", "intermediate", "sm5",
    "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(sm5_csv))
    stop("SM5a ATBD CSV not found:\n  ", sm5_csv,
         "\nRun scripts/05a_sm5_compute_metrics_atbd.R first.")
  metrics_dt <- data.table::fread(sm5_csv)
  metrics_dt <- metrics_dt[h_min == h_min_select & grepl("keepTrees", Norm)]
} else {
  kt_csv <- here::here(
    "output", "intermediate", "sm5",
    "kt_sensitivity_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(kt_csv))
    stop("kt CSV not found:\n  ", kt_csv,
         "\nRun scripts/05c_sm5_k_zmin_sensitivity_dopt.R first.")
  dt_kt <- data.table::fread(kt_csv)
  metrics_dt <- dt_kt[k == k_select & theta == theta_select &
                        grepl("keepTrees", Norm)]
  if (nrow(metrics_dt) == 0L)
    stop("No rows found for k=", k_select, ", theta=", theta_select)
  metrics_dt[, h_min := h_min_select]
  cli::cli_alert_info(
    "Using 05c CSV — k={k_select}, theta={theta_select}°, h_min={h_min_select}"
  )
}

# ── Apply combined_mode ────────────────────────────────────────────────────────

combined_row_names <- c("All_sites", "Sites_combined", "Sites_averaged",
                        "Sites combined", "Sites averaged")
dt_individual <- metrics_dt[!Site %in% combined_row_names]

if (combined_mode == "average") {
  combined_label <- "Sites averaged"
  avg_dt <- dt_individual[Site %in% sites_individual,
    .(R     = mean(R,     na.rm = TRUE),
      R2    = mean(R2,    na.rm = TRUE),
      RMSE  = mean(RMSE,  na.rm = TRUE),
      Bias  = mean(Bias,  na.rm = TRUE),
      Slope = mean(Slope, na.rm = TRUE),
      ATBD  = TRUE, Column = "LAI_atbd",
      h_min = h_min_select),
    by = .(Norm, Depth)
  ]
  avg_dt[, Site := combined_label]
  metrics_dt <- data.table::rbindlist(list(dt_individual, avg_dt),
                                       use.names = TRUE, fill = TRUE)
} else {
  combined_label <- "Sites combined"
  pool_dt <- metrics_dt[Site %in% combined_row_names][, Site := combined_label]
  metrics_dt <- data.table::rbindlist(list(dt_individual, pool_dt),
                                       use.names = TRUE, fill = TRUE)
}

# ── Select d_opt ───────────────────────────────────────────────────────────────

t0 <- proc.time()

dopt_reference <- select_dopt(
  metrics_dt,
  methods         = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth       = h_min_select,
  prosail_filter  = "ATBD",
  pareto_criteria    = pareto_criteria,
  pareto_norm_method = pareto_norm_method
)

dopt_all <- select_dopt(
  metrics_dt,
  methods         = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth       = h_min_select,
  prosail_filter  = "all",
  pareto_criteria    = pareto_criteria,
  pareto_norm_method = pareto_norm_method
)

# No max_depth constraint — d_opt free to go up to 38 m
dopt_reference_free <- select_dopt(
  metrics_dt,
  methods         = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth       = NULL,
  prosail_filter  = "ATBD",
  pareto_criteria    = pareto_criteria,
  pareto_norm_method = pareto_norm_method
)

dopt_all_free <- select_dopt(
  metrics_dt,
  methods         = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth       = NULL,
  prosail_filter  = "all",
  pareto_criteria    = pareto_criteria,
  pareto_norm_method = pareto_norm_method
)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Write CSV ──────────────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "intermediate", "sm5")
data.table::fwrite(dopt_reference, file.path(out_dir, "dopt_reference.csv"))
data.table::fwrite(dopt_all,       file.path(out_dir, "dopt_all.csv"))

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("SM5 d_opt — résumé")
cli::cli_bullets(c(
  "v" = "dopt_reference.csv — {nrow(dopt_reference)} lignes (ATBD)",
  "v" = "dopt_all.csv       — {nrow(dopt_all)} lignes (toutes configs)",
  "i" = "Durée : {elapsed} s"
))

make_wide <- function(dopt_dt) {
  data.table::dcast(
    dopt_dt,
    Site + Norm ~ method_dopt,
    value.var     = "d_opt",
    fun.aggregate = function(x) x[[1L]]
  )
}

w_constrained <- make_wide(dopt_reference)
w_free        <- make_wide(dopt_reference_free)

cli::cli_h2("1. d_opt par Norm — max_depth = {h_min_select} m (ATBD)")
data.table::setorder(w_constrained, Norm, Site)
print(w_constrained)

cli::cli_h2("2. d_opt par Site — max_depth = {h_min_select} m (ATBD)")
data.table::setorder(w_constrained, Site, Norm)
print(w_constrained)

cli::cli_h2("3. d_opt par Norm — sans contrainte de profondeur (ATBD)")
data.table::setorder(w_free, Norm, Site)
print(w_free)

cli::cli_h2("4. d_opt par Site — sans contrainte de profondeur (ATBD)")
data.table::setorder(w_free, Site, Norm)
print(w_free)

# ── Métriques au d_opt Pareto (constrained) — par site et norm ────────────────

cli::cli_h2("5. Métriques au d_opt Pareto (max_depth = {h_min_select} m, ATBD)")

pareto_metrics <- dopt_reference[
  method_dopt == "pareto",
  .(Site, Norm, d_opt,
    R    = round(R,     2),
    R2   = round(R2,    2),
    RMSE = round(RMSE,  2),
    Bias = round(Bias,  2),
    Slope = round(Slope, 2))
]

norm_short <- c(DSM_keepTrees = "CHM", DTM_keepTrees = "DTM")
pareto_metrics[, Norm_short := norm_short[Norm]]

sites_order <- c(sites_individual, combined_label)
pareto_metrics[, Site := factor(Site, levels = sites_order)]
data.table::setorder(pareto_metrics, Site, Norm)

for (s in sites_order) {
  sub <- pareto_metrics[Site == s]
  if (nrow(sub) == 0L) next
  cli::cli_h3("{s}")
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    cli::cli_text(
      "  {r$Norm_short} (d_opt = {r$d_opt} m) — ",
      "r = {r$R}  R2 = {r$R2}  RMSE = {r$RMSE}  Bias = {r$Bias}  Slope = {r$Slope}"
    )
  }
}
