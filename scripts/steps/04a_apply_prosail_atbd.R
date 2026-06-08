# ---
# title:  04a_apply_prosail_atbd.R
# desc:   Orchestration — applies the ATBD SVR ensemble (from 02a) to the
#         Sentinel-2 reflectances sampled by 03a, for all sites × h_min values.
#         Produces LAI_S2_atbd estimates used by 05a to determine d_opt.
#
# Prerequisite (02a): ATBD SVR RDS must be present at
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/atbd/LIDFa=1_lai=1_LMA=1_BROWN=1.rds
#
# Prerequisite (03a): S2 reflectance CSVs must be present at
#   03_RESULTS/{site}/PROSAIL_Optimization/sampling/
#   S2_reflectance_stratified_uniform_hmin{h_min}_nbSamples_5000.csv
#
# Output: one CSV per site × h_min at:
#   revision/output/intermediate/PROSAIL_Models/{site}/
#   LIDFa_lai_LMA_BROWN/atbd/
#   LAI_estimated_atbd_stratified_uniform_hmin{h_min}_nbSamples_5000.csv
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/04a_apply_prosail_atbd.R")
# ---

library("here")
library("readr")
library("cli")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "prosail_inversion.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites              <- c("Aigoual", "Blois", "Mormal")
Samplingmethod     <- "stratified_uniform"
name_strategy      <- "LIDFa_lai_LMA_BROWN"
nbSamples_sampling <- 5000
h_min_values       <- c(10L, 15L, 20L)
S2BandSelect_2     <- c("B03", "B04", "B08")

atbd_label <- "LIDFa=1_lai=1_LMA=1_BROWN=1"

# ── Inversion loop ────────────────────────────────────────────────────────────

t_global <- Sys.time()

for (site in sites) {

  cat("\n── ATBD inversion for site:", site, "──\n")

  models_dir   <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy, "atbd"
  )
  filename_svr <- file.path(models_dir, paste0(atbd_label, ".rds"))

  if (!file.exists(filename_svr)) {
    warning("ATBD SVR not found for site ", site, " — run 02a first. Skipping.")
    next
  }

  modelSVR <- load_svr_ensemble(filename_svr)

  for (h_min in h_min_values) {

    Sampling_path <- file.path(
      paths$ext_results, site, "PROSAIL_Optimization", "sampling",
      paste0("S2_reflectance_", Samplingmethod,
             "_hmin", h_min,
             "_nbSamples_", nbSamples_sampling, ".csv")
    )

    if (!file.exists(Sampling_path)) {
      warning("S2 reflectance CSV not found for site ", site,
              " h_min=", h_min, " — run 03a first. Skipping.")
      next
    }

    Refl <- readr::read_delim(Sampling_path, delim = "\t", show_col_types = FALSE)

    selbands <- match(S2BandSelect_2, sub(" .*", "", names(Refl)))
    Refl     <- Refl / 10000

    LAIest <- apply_svr_to_reflectance(
      svr_ensemble       = modelSVR,
      reflectance_matrix = Refl[, selbands]
    )

    out_csv <- file.path(
      models_dir,
      paste0("LAI_estimated_atbd_", Samplingmethod,
             "_hmin", h_min,
             "_nbSamples_", nbSamples_sampling, ".csv")
    )

    readr::write_delim(
      x     = round(data.frame(LAI_atbd = LAIest$mean_estimate), 4),
      file  = out_csv,
      delim = "\t"
    )

    cat("  h_min", h_min, "— written:", basename(out_csv), "\n")
    gc()
  }
}

cat("\n── 04a done — total:",
    round(difftime(Sys.time(), t_global, units = "mins"), 1), "min ──\n")
