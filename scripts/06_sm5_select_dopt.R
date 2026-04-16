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
  "all_results_combined_LIDFa_lai_LMA_BROWN.csv"
)

if (!file.exists(sm5_csv)) {
  stop(
    "SM5 passe 1 CSV not found:\n  ", sm5_csv,
    "\nRun revision/scripts/05_sm5_compute_metrics.R first."
  )
}

# ── Source and load ────────────────────────────────────────────────────────────

source(here::here("revision", "R", "sm5_dopt.R"))

metrics_dt <- data.table::fread(sm5_csv)

# ── Select d_opt ───────────────────────────────────────────────────────────────

t0 <- proc.time()

dopt_reference <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "pareto"),
  prosail_filter = "ATBD"
)

dopt_all <- select_dopt(
  metrics_dt,
  methods        = c("pearson", "rmse", "bias", "pareto"),
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

# Wide table: d_opt by (Site × Norm) × method_dopt for the ATBD reference.
# With the standard "LIDFa_lai_LMA_BROWN" LUT there is exactly one ATBD column,
# so aggregation below is stable (fun.aggregate takes first value if >1).
wide <- data.table::dcast(
  dopt_reference,
  Site + Norm ~ method_dopt,
  value.var    = "d_opt",
  fun.aggregate = function(x) x[[1L]]
)

cli::cli_h2("d_opt par Site × Norm — méthodes comparées (ATBD)")
print(wide)
