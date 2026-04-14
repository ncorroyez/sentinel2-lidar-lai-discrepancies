# ---
# title:  05_sm5_compute_metrics.R
# desc:   Orchestration — computes the full metrics table for the d_opt
#         analysis (SM5) and writes it to revision/output/intermediate/sm5/.
#
#         Calls build_metrics_table() from revision/R/sm5_aggregate.R,
#         which is a faithful refactor of:
#           PROSAIL-Optimization/02_CODES/Sentinel2/Main_s2_04C_analysis_depth.R
#
#         Prerequisites:
#           SM3 — PAD CSVs must exist at:
#             03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#             PAD_{norm}_Depth_{d}_Samples_stratified_uniform_nbSamples_5000.csv
#           SM4 — LAI estimated CSVs must exist at:
#             revision/output/intermediate/PROSAIL_Models/{site}/
#             LIDFa_lai_LMA_BROWN/
#             LAI_estimated_stratified_uniform_nbSamples_5000.csv
#
#         If either set of prerequisites is absent, build_metrics_table()
#         skips the missing combinations and warns; the output CSV is written
#         but will be empty or partial.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05_sm5_compute_metrics.R")
# ---

source(here::here("revision", "R", "sm5_metrics.R"))
source(here::here("revision", "R", "sm5_aggregate.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites            <- c("Aigoual", "Blois", "Mormal")

# All four normalisation types active (legacy 04C had only DSM_keepTrees active;
# all four enabled here per Nathan's decision).
norm_methods     <- c("DSM", "DTM", "DSM_keepTrees", "DTM_keepTrees")

depths           <- 1:38          # same as legacy 04C line 30
sampling_methods <- "stratified_uniform"
name_strategy    <- "LIDFa_lai_LMA_BROWN"  # no _Agg_10m suffix (§4)
parms2test       <- c("LIDFa", "lai", "LMA", "BROWN")
nb_samples       <- 5000L

# ── Output path ────────────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "intermediate", "sm5")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
out_csv <- file.path(out_dir,
                     paste0("all_results_combined_", name_strategy, ".csv"))

# ── Build table ────────────────────────────────────────────────────────────────

t0 <- proc.time()

combined <- build_metrics_table(
  sites            = sites,
  norm_methods     = norm_methods,
  depths           = depths,
  sampling_methods = sampling_methods,
  name_strategy    = name_strategy,
  parms2test       = parms2test,
  nb_samples       = nb_samples
)

elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

# ── Write CSV (single fwrite — no append, no double-write) ────────────────────

data.table::fwrite(combined, out_csv)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("SM5 — résumé")
cli::cli_bullets(c(
  "v" = "Fichier écrit : {out_csv}",
  "i" = "Lignes totales : {nrow(combined)}",
  "i" = paste0(
    "Sites couverts : ",
    paste(sort(unique(combined$Site)), collapse = ", ")
  ),
  "i" = paste0(
    "Plage Depth : ",
    min(combined$Depth), " – ", max(combined$Depth)
  ),
  "i" = paste0(
    "Norms couvertes : ",
    paste(sort(unique(combined$Norm)), collapse = ", ")
  ),
  "i" = "Durée : {elapsed} min"
))
