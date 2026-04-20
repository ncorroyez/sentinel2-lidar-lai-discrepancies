# ---
# title:  04_apply_prosail_inversion.R
# desc:   Orchestration — applies trained SVR ensembles (SM2) to Sentinel-2
#         reflectances sampled by SM3, producing per-combination LAI estimates
#         for all three study sites.
#
#         Reproduces faithfully the inversion section of
#         Main_s2_03_PROSAIL_Inversion_bis.R (lines 151-201), refactored with
#         here::here() paths and revision/R/prosail_inversion.R.
#
# Prerequisite (SM2): SVR model .rds files must be present at
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/*.rds
#   (produced by revision/scripts/02_train_prosail_lut.R)
#
# Prerequisite (SM3): S2 reflectance CSVs must be present at
#   03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   S2_reflectance_stratified_uniform_nbSamples_5000.csv
#   (produced by revision/scripts/03a_sample_s2_pixels.R)
#
# Prerequisite (SM2): simulation strategy RDS must be present at
#   revision/output/intermediate/PROSAIL_Models/Simulations_Strategy/
#   LIDFa_lai_LMA_BROWN/simulation_strategy.rds
#   (produced by revision/scripts/02_train_prosail_lut.R via
#    get_parameter_combinations())
#
# Output: two CSVs per site per sampling method, at:
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/
#     LAI_estimated_stratified_uniform_nbSamples_5000.csv
#     LAI_estimated_stratified_uniform_nbSamples_5000_SD.csv
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/04_apply_prosail_inversion.R")
# ---

source(here::here("revision", "R", "prosail_inversion.R"))

# ── Global parameters — recopie fidèle de _03_bis_ ────────────────────────────

sites             <- c("Aigoual", "Blois", "Mormal")
Samplingmethods   <- "stratified_uniform"   # one method in production (_03_bis_ l.34)
name_strategy     <- "LIDFa_lai_LMA_BROWN"  # without _Agg_10m suffix; cohérent with SM2
nbSamples_sampling <- 5000                  # number of S2 reflectance samples (_03_bis_ nbSamplesS2Refl)

# S2BandSelect = prosail-internal band names used by SM2 training.
# S2BandSelect_2 = CSV column names written by preprocS2 / SM3, used here to
#   match selbands from the reflectance table. Both are required and distinct.
S2BandSelect   <- c("B3",  "B4",  "B8")   # training bands (prosail internal notation)
S2BandSelect_2 <- c("B03", "B04", "B08")  # CSV column prefix (preprocS2 notation, _03_bis_ l.30)

# ── Simulation strategy — combination labels ───────────────────────────────────

strategy_path <- here::here(
  "revision", "output", "intermediate", "PROSAIL_Models",
  "Simulations_Strategy", name_strategy, "simulation_strategy.rds"
)
simulations <- readRDS(strategy_path)

combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

cat("Combination labels:", length(combination_labels), "combinations\n")

# ── SVR model directories (SM2 output) ────────────────────────────────────────
# Points to revision/output/intermediate/ (SM2 output), NOT 03_RESULTS/
# (historical path of _03_bis_). SM4 reads the .rds produced by SM2 refactorisé.

prosail_models_dir <- setNames(
  vapply(sites, function(s) {
    here::here("revision", "output", "intermediate",
               "PROSAIL_Models", s, name_strategy)
  }, character(1)),
  sites
)

# ── Global chronometer ────────────────────────────────────────────────────────

t_global <- Sys.time()

# ── Inversion loop — mirrors _03_bis_ lines 126-203 ──────────────────────────

for (site in sites) {

  cat("\n── Inversion PROSAIL for site:", site, "──\n")
  t_site <- Sys.time()

  # Paths to SVR .rds files for this site (mirrors _03_bis_ l.129)
  filename_SVR <- file.path(
    prosail_models_dir[site],
    paste0(combination_labels, ".rds")
  )

  for (Samplingmethod in Samplingmethods) {

    # i. Read S2 reflectance CSV (mirrors _03_bis_ l.153-157)
    Sampling_path <- here::here(
      "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
      paste0("S2_reflectance_", Samplingmethod,
             "_nbSamples_", nbSamples_sampling, ".csv")
    )
    Refl <- readr::read_delim(
      file           = Sampling_path,
      delim          = "\t",
      show_col_types = FALSE
    )

    # ii. Band index for prosail_hybrid_apply (mirrors _03_bis_ l.158-161)
    #     Matches S2BandSelect_2 against the column-name prefix
    #     (preprocS2 writes "B03 ..." with extra metadata after the space).
    selbands <- match(
      x     = S2BandSelect_2,
      table = sub(" .*", "", names(Refl))
    )

    # iii. Scale reflectances to [0, 1] (mirrors _03_bis_ l.162)
    Refl <- Refl / 10000

    # iv. Initialise output data.frames 5000 × 270 (mirrors _03_bis_ l.168-170)
    LAI_Mean <- data.frame(matrix(
      data = NA,
      nrow = nrow(Refl),
      ncol = length(combination_labels)
    ))
    LAI_SD <- data.frame(matrix(
      data = NA,
      nrow = nrow(Refl),
      ncol = length(combination_labels)
    ))

    # v. Progress bar (cli replaces progress::progress_bar — not in whitelist)
    cli::cli_progress_bar(
      "Inversion SVR",
      total = length(combination_labels)
    )

    # vi. Inner loop over 270 combinations (mirrors _03_bis_ l.171-180)
    for (i in seq_along(combination_labels)) {

      modelSVR <- load_svr_ensemble(filename_SVR[i])
      LAIest   <- apply_svr_to_reflectance(
        svr_ensemble       = modelSVR,
        reflectance_matrix = Refl[, selbands]
      )
      LAI_Mean[, i] <- LAIest$mean_estimate
      LAI_SD[, i]   <- LAIest$sd_estimate

      cli::cli_progress_update()
    }

    cli::cli_progress_done()

    # vii. Rename columns (mirrors _03_bis_ l.182-183)
    names(LAI_Mean) <- combination_labels
    names(LAI_SD)   <- combination_labels

    # viii. Write output CSVs (mirrors _03_bis_ l.185-197, path adapted to SM2 output)
    output_dir <- prosail_models_dir[site]
    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    mean_path <- file.path(
      output_dir,
      paste0("LAI_estimated_", Samplingmethod,
             "_nbSamples_", nbSamples_sampling, ".csv")
    )
    sd_path <- file.path(
      output_dir,
      paste0("LAI_estimated_", Samplingmethod,
             "_nbSamples_", nbSamples_sampling, "_SD.csv")
    )

    readr::write_delim(
      x     = round(x = LAI_Mean, digits = 4),
      file  = mean_path,
      delim = "\t"
    )
    readr::write_delim(
      x     = round(x = LAI_SD, digits = 4),
      file  = sd_path,
      delim = "\t"
    )

    cat("Written:", basename(mean_path), "\n")
    cat("Written:", basename(sd_path),   "\n")

    # ix. Defensive cleanup (mirrors _03_bis_ l.198-200)
    rm(LAI_Mean)
    rm(LAI_SD)
    gc()
  }

  cat("Site", site, "done —",
      round(difftime(Sys.time(), t_site, units = "mins"), 1), "min\n")
}

cat("\n── All sites done — total:",
    round(difftime(Sys.time(), t_global, units = "mins"), 1), "min ──\n")
