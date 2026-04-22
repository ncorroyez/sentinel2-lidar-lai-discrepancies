# ---
# title:  10_sm5_compute_metrics_full.R
# desc:   Orchestration — computes the full metrics table for all 270 PROSAIL
#         configurations × depths × h_min × LAI scenarios, and writes a single
#         combined CSV to revision/output/intermediate/sm5/.
#
#         Calls build_metrics_table() from revision/R/sm5_aggregate.R for each
#         (h_min × lai_scenario) combination and stacks the results.
#
# Prerequisites:
#   03b — PAD CSVs at 03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   09  — LAI estimated CSVs at
#         revision/output/intermediate/PROSAIL_Models/{site}/
#         LIDFa_lai_LMA_BROWN/{scenario}/
#         LAI_estimated_{scenario}_stratified_uniform_hmin{h_min}_nbSamples_5000.csv
#
# Output:
#   revision/output/intermediate/sm5/
#   all_results_combined_LIDFa_lai_LMA_BROWN.csv
#   (columns: Site, Norm, Depth, Method, h_min, lai_scenario, Column,
#             R, R2, RMSE, Bias, Slope, ATBD)
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/10_sm5_compute_metrics_full.R")
# ---

library("here")
library("data.table")
library("cli")

source(here::here("revision", "R", "sm5_metrics.R"))
source(here::here("revision", "R", "sm5_aggregate.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites            <- c("Aigoual", "Blois", "Mormal")
norm_methods     <- c("DSM", "DTM", "DSM_keepTrees", "DTM_keepTrees")
depths           <- 1:38
sampling_methods <- "stratified_uniform"
name_strategy    <- "LIDFa_lai_LMA_BROWN"
parms2test       <- c("LIDFa", "lai", "LMA", "BROWN")
nb_samples       <- 5000L
h_min_values     <- c(10L, 15L, 20L)
lai_scenarios    <- c("per_site", "common", "fixed_4")

# ── Output ─────────────────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "intermediate", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(out_dir,
                     paste0("all_results_combined_", name_strategy, ".csv"))

# ── Build table ────────────────────────────────────────────────────────────────

t0 <- proc.time()

chunks <- list()
for (lai_scenario in lai_scenarios) {
  for (h_min in h_min_values) {
    cli::cli_alert_info("h_min = {h_min} m | scenario = {lai_scenario}")
    chunks <- c(chunks, list(
      build_metrics_table(
        sites            = sites,
        norm_methods     = norm_methods,
        depths           = depths,
        sampling_methods = sampling_methods,
        name_strategy    = name_strategy,
        parms2test       = parms2test,
        nb_samples       = nb_samples,
        h_min            = h_min,
        lai_scenario     = lai_scenario
      )
    ))
  }
}
combined <- data.table::rbindlist(chunks, use.names = TRUE)

elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

data.table::fwrite(combined, out_csv)

cli::cli_h1("SM5 full — résumé")
cli::cli_bullets(c(
  "v" = "Fichier écrit : {out_csv}",
  "i" = "Lignes totales : {nrow(combined)}",
  "i" = paste0("Sites : ",       paste(sort(unique(combined$Site)),         collapse = ", ")),
  "i" = paste0("Norms : ",       paste(sort(unique(combined$Norm)),         collapse = ", ")),
  "i" = paste0("h_min : ",       paste(sort(unique(combined$h_min)),        collapse = ", "), " m"),
  "i" = paste0("Scénarios LAI : ", paste(sort(unique(combined$lai_scenario)), collapse = ", ")),
  "i" = paste0("Profondeurs : ", min(combined$Depth), " – ", max(combined$Depth)),
  "i" = "Durée : {elapsed} min"
))
