# ---
# title:  09_apply_prosail_full.R
# desc:   Orchestration — applies all 270 trained SVR ensembles (from 08) to
#         Sentinel-2 reflectances sampled by 03a, for all sites × h_min values
#         × LAI scenarios (per_site / common / fixed_4).
#
# Prerequisite (08): SVR RDS files at
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/{scenario}/*.rds
#
# Prerequisite (03a): S2 reflectance CSVs at
#   03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   S2_reflectance_stratified_uniform_hmin{h_min}_nbSamples_5000.csv
#
# Output: two CSVs per site × scenario × h_min at:
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/{scenario}/
#     LAI_estimated_{scenario}_stratified_uniform_hmin{h_min}_nbSamples_5000.csv
#     LAI_estimated_{scenario}_stratified_uniform_hmin{h_min}_nbSamples_5000_SD.csv
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/09_apply_prosail_full.R")
# ---

library("here")
library("readr")
library("cli")

source(here::here("revision", "R", "prosail_inversion.R"))

# ── Parameters ─────────────────────────────────────────────----------------------------------------------------------------

sites              <- c("Aigoual", "Blois", "Mormal")
Samplingmethods    <- "stratified_uniform"
name_strategy      <- "LIDFa_lai_LMA_BROWN"
nbSamples_sampling <- 5000
h_min_values       <- c(10L, 15L, 20L)

# Scenarios produced by 07 then trained by 08:
#   per_site — each site's own Pareto d_opt
#   common   — combined-sites d_opt (Sites averaged or Sites combined) applied
#              to all 3 sites
# k_select, theta_select, combined_mode must match 06 and 07.
lai_scenarios <- c("per_site", "common")

S2BandSelect_2 <- c("B03", "B04", "B08")

# ── Simulation strategy — combination labels ──────────────────────────────────

strategy_path <- here::here(
  "revision", "output", "intermediate", "PROSAIL_Models",
  "Simulations_Strategy", name_strategy, "simulation_strategy.rds"
)
simulations <- readRDS(strategy_path)

combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

cat("Combinations:", length(combination_labels), "\n")

# ── Inversion loop ────────────────────────────────────────────────────────────

t_global <- Sys.time()

for (lai_scenario in lai_scenarios) {

  cat("\n══ Scenario:", lai_scenario, "══\n")

  for (site in sites) {

    cat("\n── Site:", site, "──\n")
    t_site <- Sys.time()

    models_dir   <- here::here(
      "revision", "output", "intermediate", "PROSAIL_Models",
      site, name_strategy, lai_scenario
    )
    filename_SVR <- file.path(models_dir, paste0(combination_labels, ".rds"))

    missing <- !file.exists(filename_SVR)
    if (any(missing)) {
      warning(sum(missing), " SVR files missing for site=", site,
              " scenario=", lai_scenario, " — run 08 first. Skipping.")
      next
    }

    for (h_min in h_min_values) {

      cat("  h_min:", h_min, "m\n")

      for (Samplingmethod in Samplingmethods) {

        Sampling_path <- here::here(
          "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
          paste0("S2_reflectance_", Samplingmethod,
                 "_hmin", h_min,
                 "_nbSamples_", nbSamples_sampling, ".csv")
        )

        if (!file.exists(Sampling_path)) {
          warning("S2 CSV not found for site=", site, " h_min=", h_min,
                  " — run 03a first. Skipping.")
          next
        }

        Refl <- readr::read_delim(Sampling_path, delim = "\t",
                                   show_col_types = FALSE)

        selbands <- match(S2BandSelect_2, sub(" .*", "", names(Refl)))
        Refl     <- Refl / 10000

        LAI_Mean <- data.frame(matrix(NA, nrow = nrow(Refl),
                                       ncol = length(combination_labels)))
        LAI_SD   <- data.frame(matrix(NA, nrow = nrow(Refl),
                                       ncol = length(combination_labels)))

        cli::cli_progress_bar("Inversion SVR", total = length(combination_labels))

        for (i in seq_along(combination_labels)) {
          modelSVR      <- load_svr_ensemble(filename_SVR[i])
          LAIest        <- apply_svr_to_reflectance(modelSVR, Refl[, selbands])
          LAI_Mean[, i] <- LAIest$mean_estimate
          LAI_SD[, i]   <- LAIest$sd_estimate
          cli::cli_progress_update()
        }
        cli::cli_progress_done()

        names(LAI_Mean) <- combination_labels
        names(LAI_SD)   <- combination_labels

        mean_path <- file.path(
          models_dir,
          paste0("LAI_estimated_", lai_scenario, "_", Samplingmethod,
                 "_hmin", h_min, "_nbSamples_", nbSamples_sampling, ".csv")
        )
        sd_path <- file.path(
          models_dir,
          paste0("LAI_estimated_", lai_scenario, "_", Samplingmethod,
                 "_hmin", h_min, "_nbSamples_", nbSamples_sampling, "_SD.csv")
        )

        readr::write_delim(round(LAI_Mean, 4), mean_path, delim = "\t")
        readr::write_delim(round(LAI_SD,   4), sd_path,   delim = "\t")
        cat("  Written:", basename(mean_path), "\n")

        rm(LAI_Mean, LAI_SD)
        gc()
      }
    }

    cat("Site", site, "done —",
        round(difftime(Sys.time(), t_site, units = "mins"), 1), "min\n")
  }
}

cat("\n── 09 done — total:",
    round(difftime(Sys.time(), t_global, units = "mins"), 1), "min ──\n")
