# ---
# title:  06_sm5_select_dopt.R
# desc:   Orchestration — selects the optimal canopy integration depth (d_opt)
#         by four methods (pearson, rmse, bias, pareto) and writes two CSV
#         files to revision/output/intermediate/sm5/:
#           - dopt_reference.csv  (ATBD == TRUE columns only)
#           - dopt_all.csv        (all PROSAIL columns)
#
#         Calls select_dopt() from revision/R/sm5_dopt.R, which addresses
#         reviewer comments R3.minor.13 and R4.spec.2.
#
#         Prerequisites:
#           SM5 passe 1 — CSV must exist at:
#             revision/output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#           Run revision/scripts/05_sm5_compute_metrics.R first.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/06_sm5_select_dopt.R")
# ---

# ── Prerequisites ──────────────────────────────────────────────────────────────

sm5_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
)

if (!file.exists(sm5_csv)) {
  stop(
    "SM5a ATBD CSV not found:\n  ", sm5_csv,
    "\nRun revision/scripts/05a_sm5_compute_metrics_atbd.R first."
  )
}

# ── Source and load ────────────────────────────────────────────────────────────

source(here::here("revision", "R", "sm5_dopt.R"))

# h_min to analyse (must match a value produced by 05_sm5_compute_metrics.R).
# Change to 15L or 20L to run the sensitivity analysis for other thresholds.
h_min_select <- 10L

metrics_dt <- data.table::fread(sm5_csv)

# Filter to the selected h_min threshold
metrics_dt <- metrics_dt[h_min == h_min_select]

# Keep only keepTrees normalisation methods (DSM_keepTrees, DTM_keepTrees)
metrics_dt <- metrics_dt[grepl("keepTrees", Norm)]

# ── Select d_opt ───────────────────────────────────────────────────────────────

t0 <- proc.time()

dopt_reference <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth      = h_min_select,
  prosail_filter = "ATBD"
)

dopt_all <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth      = h_min_select,
  prosail_filter = "all"
)

# No max_depth constraint — d_opt free to go up to 38 m
dopt_reference_free <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth      = NULL,
  prosail_filter = "ATBD"
)

dopt_all_free <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "slope", "pareto"),
  max_depth      = NULL,
  prosail_filter = "all"
)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Write CSV ──────────────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "intermediate", "sm5")
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
