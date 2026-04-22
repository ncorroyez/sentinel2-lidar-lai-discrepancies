# ---
# title:  prosail_lut.R
# author: Nathan Corroyez, UMR TETIS
# desc:   PROSAIL LUT generation and SVR ensemble training wrappers.
#         Replaces get_combination.R, define_parm_combinations.R,
#         PROSAIL_train_sensitivity.R, and PROSAIL_train_sensitivity_subset.R
#         from PROSAIL-Optimization/02_CODES/libraries/ for the RSE revision.
#
#         Five public functions, in dependency order:
#           get_parameter_combinations() — build simulation grid (270 combos)
#           build_prosail_lut()          — generate LUT + apply SRF + add noise
#           train_svr_ensemble()         — fit liquidSVM ensemble on one LUT
#           train_one_simulation()       — per-chunk worker (maps over grid rows)
#           train_all_simulations()      — orchestrator: sequential or parallel
# ---


# ── get_parameter_combinations ────────────────────────────────────────────────

#' @title Define PROSAIL parameter sampling strategies and simulation grid
#'
#' @description
#' Builds the cross-product simulation grid used to train the PROSAIL SVR
#' ensemble. For each parameter in \code{parms2test}, a list of sampling
#' variants (ATBD default or alternative Gaussian/uniform distributions) is
#' looked up from the hard-coded table — ported verbatim from
#' \code{get_combination.R} in
#' \code{PROSAIL-Optimization/02_CODES/libraries/}. An \code{expand.grid()}
#' over the variant indices produces the full simulation grid. The resulting
#' \code{simulations} list is saved to \code{output_dir} as an RDS and
#' returned.
#'
#' @details
#' **Simulation grid size**: with \code{parms2test = c("LIDFa","lai","LMA","BROWN")}
#' and their respective 5, 6, 3, 3 variants, the grid has 5 × 6 × 3 × 3 = 270
#' rows. Each row is passed to \code{train_one_simulation()} as one SVR model
#' to fit.
#'
#' **Embedded distributions**: the hard-coded parameter table exactly matches
#' \code{get_combination.R} (active, non-commented section). Any change to the
#' distributions must be made here and documented in the commit.
#'
#' **Output path**: the RDS lands at
#' \code{output_dir/Simulations_Strategy/<name_strategy>/simulation_strategy.rds},
#' mirroring the original pipeline's convention.
#'
#' @param parms2test   Character vector. Names of PROSAIL parameters to vary.
#'   Supported: \code{"LIDFa"}, \code{"lai"}, \code{"LMA"}, \code{"BROWN"}.
#'   Default strategy: \code{c("LIDFa","lai","LMA","BROWN")}.
#' @param output_dir   Character. Root output directory (e.g.
#'   \code{revision/output/intermediate/PROSAIL_Models/<site>}). Sub-directory
#'   \code{Simulations_Strategy/<name_strategy>/} is created inside.
#' @param name_strategy Character. Label for this parameter strategy (e.g.
#'   \code{"LIDFa_lai_LMA_BROWN"}). Used as sub-directory name and in output
#'   file names.
#' @param overwrite    Logical. If \code{FALSE} (default), stops with an error
#'   if the RDS already exists. Set \code{TRUE} to overwrite.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{grid_simu}}{A \code{data.frame} with one row per simulation
#'       and one column per parameter; values are 1-based indices into each
#'       parameter's variant list.}
#'     \item{\code{combination}}{A named list, one element per parameter in
#'       \code{parms2test}. Each element is itself a list of variant
#'       specifications (either the string \code{"ATBD"} or a named numeric
#'       vector with \code{distribution}, \code{min}, \code{max}, etc.).}
#'   }
#'
#' @examples
#' \dontrun{
#' out_dir <- here::here("revision", "output", "intermediate",
#'                       "PROSAIL_Models", "Blois")
#' simulations <- get_parameter_combinations(
#'   parms2test   = c("LIDFa", "lai", "LMA", "BROWN"),
#'   output_dir   = out_dir,
#'   name_strategy = "LIDFa_lai_LMA_BROWN"
#' )
#' nrow(simulations$grid_simu)  # 270
#' }
#'
# ── sample_lai_stratified_uniform ─────────────────────────────────────────────

#' @title Stratified uniform sampling of LiDAR LAI values for SVR training
#'
#' @description
#' Draws \code{n_samples} values from a vector of LiDAR-derived LAI values
#' using a stratified uniform strategy: the range [\code{lai_min},
#' \code{quantile(x, q_max)}] is divided into 1 m²/m² bins, and an equal
#' number of observations is drawn from each bin (with replacement when a bin
#' has fewer observations than required). This avoids over-representing the
#' modal LAI class in the SVR training set and ensures the model is exposed to
#' the full range of observed LAI values.
#'
#' @param x Numeric vector. LAI values (already masked: all values >= lai_min).
#' @param n_samples Integer. Total number of samples to draw. Default 5000.
#' @param lai_min Numeric. Lower bound of the sampling range (default 2).
#' @param q_max Numeric in (0,1). Quantile used as upper bound (default 0.95).
#'
#' @return Numeric vector of length \code{n_samples}.
#'
#' @examples
#' \dontrun{
#'   lai_vals <- terra::values(lidar_lai_rast, na.rm = TRUE)
#'   pool <- sample_lai_stratified_uniform(lai_vals, n_samples = 5000)
#' }
#' @export
sample_lai_stratified_uniform <- function(x, n_samples = 5000,
                                          lai_min = 2, q_max = 0.98) {
  lai_max  <- quantile(x, q_max, na.rm = TRUE)
  breaks   <- seq(lai_min, ceiling(lai_max), by = 1)
  n_bins   <- length(breaks) - 1
  per_bin  <- floor(n_samples / n_bins)
  extra    <- n_samples - per_bin * n_bins   # distribute remainder to first bins

  out <- vector("numeric", n_samples)
  idx <- 1L
  for (b in seq_len(n_bins)) {
    lo   <- breaks[b]
    hi   <- breaks[b + 1]
    pool <- x[x >= lo & x < hi]
    k    <- per_bin + (b <= extra)
    if (length(pool) == 0L) pool <- x  # fallback: draw from full range
    out[idx:(idx + k - 1L)] <- sample(pool, size = k, replace = length(pool) < k)
    idx  <- idx + k
  }
  out
}


#' @export
get_parameter_combinations <- function(parms2test,
                                        output_dir,
                                        name_strategy,
                                        overwrite = FALSE) {
  # Internal helper — verbatim port of get_combination.R (active section only).
  # Returns the list of sampling variants for a single PROSAIL parameter.
  get_combination <- function(parm) {
    combo <- list()
    combo$LIDFa <- list(
      "ATBD",
      c(min = 20, max = 50, distribution = "gauss", mean = 35, sd = 20),
      c(min = 20, max = 60, distribution = "gauss", mean = 40, sd = 20),
      c(min = 30, max = 50, distribution = "gauss", mean = 40, sd = 20),
      c(min = 30, max = 60, distribution = "gauss", mean = 45, sd = 20)
    )
    combo$lai <- list(
      "ATBD",
      c(min = 0,  max = 9,  distribution = "gauss", mean = 4, sd = 3),
      c(min = 0,  max = 11, distribution = "gauss", mean = 5, sd = 3),
      c(min = 0,  max = 13, distribution = "gauss", mean = 6, sd = 3),
      c(distribution = "LiDAR_LAI"),
      c(distribution = "LiDAR_LAI_Best_Site_Depth")
    )
    combo$LMA <- list(
      "ATBD",
      c(min = 0.003, max = 0.02, distribution = "gauss", mean = 0.010, sd = 0.005),
      c(min = 0.003, max = 0.03, distribution = "gauss", mean = 0.015, sd = 0.005)
    )
    combo$BROWN <- list(
      "ATBD",
      c(min = 0, max = 1, distribution = "gauss",   mean = 0, sd = 0.1),
      c(min = 0, max = 0, distribution = "uniform")
    )
    combo[[parm]]
  }

  combination <- parm_vals <- list()
  for (parm in parms2test) {
    combination[[parm]] <- get_combination(parm)
    parm_vals[[parm]]   <- seq_len(length(combination[[parm]]))
  }
  grid_simu   <- expand.grid(parm_vals)
  simulations <- list(grid_simu = grid_simu, combination = combination)

  output_subdir <- file.path(output_dir, "Simulations_Strategy", name_strategy)
  filename      <- file.path(output_subdir, "simulation_strategy.rds")
  if (file.exists(filename) && !overwrite) {
    stop(
      "get_parameter_combinations(): simulation strategy already exists:\n  ",
      filename, "\nSet overwrite = TRUE to replace it."
    )
  }
  dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(simulations, file = filename)
  simulations
}


# ── build_prosail_lut ─────────────────────────────────────────────────────────

#' @title Generate a PROSAIL bidirectional reflectance LUT and apply sensor SRF
#'
#' @description
#' Generates a Look-Up Table (LUT) of Sentinel-2 bidirectional reflectances by
#' running the 4SAIL canopy RTM coupled with PROSPECT-5B via
#' \code{prosail::Generate_LUT_4SAIL()}, convolving the 1 nm outputs to
#' Sentinel-2 band centres with \code{prosail::applySensorCharacteristics()},
#' and adding ATBD-calibrated noise with \code{prosail::apply_noise_atbd()}.
#' Returns only the rows corresponding to \code{bands_to_select}.
#'
#' @details
#' This function consolidates the three-step LUT pipeline that appears inline
#' in \code{PROSAIL_train_sensitivity_subset.R}:
#' (1) \code{generate_lut_4sail()} (prosail ≥ 3.0, replaces \code{Generate_LUT_4SAIL}),
#' (2) \code{apply_sensor_characteristics()},
#' (3) \code{apply_noise_atbd()}.
#'
#' Updated for prosail 3.0.0 (2026-04): uses \code{generate_lut_4sail()} with
#' \code{spec_soil_ossl} (48-spectrum OSSL soil database, replaces 3-spectrum
#' \code{SpecSOIL}) and \code{spec_atm} (snake_case column names, replaces
#' \code{SpecATM}). The \code{surface_refl} component of the returned list
#' is the full coupled BRF that \code{apply_sensor_characteristics()} expects.
#' Passing \code{SpecATM} (camelCase columns) instead of \code{spec_atm} silently
#' produces an empty \code{surface_refl} — hence the explicit snake_case argument.
#'
#' The function does not write any files; it is called from within
#' \code{train_one_simulation()}.
#'
#' @param input_prosail A \code{data.frame} of PROSAIL input parameters as
#'   returned by \code{prosail::get_atbd_v3_lut_input()}, one row per LUT sample.
#' @param srf           A sensor spectral response function object as returned
#'   by \code{prosail::get_srf_sensor("Sentinel_2")}.
#' @param bands_to_select Integer vector. Row indices (into the BRF matrix)
#'   corresponding to the bands selected for inversion
#'   (e.g. \code{match(c("B3","B4","B8"), srf$Spectral_Bands)}).
#'
#' @return A numeric matrix of dimensions \code{length(bands_to_select) ×
#'   nrow(input_prosail)}: noisy simulated reflectances for the selected bands,
#'   one column per LUT sample.
#'
#' @examples
#' \dontrun{
#' srf  <- prosail::get_srf_sensor("Sentinel_2")
#' bsel <- match(c("B3","B4","B8"), srf$Spectral_Bands)
#' ip   <- prosail::get_atbd_v3_lut_input(nb_samples = 500)
#' lut  <- build_prosail_lut(ip, srf, bsel)
#' dim(lut)  # 3 x 500
#' }
#'
#' @export
build_prosail_lut <- function(input_prosail, srf, bands_to_select) {
  res <- prosail::generate_lut_4sail(
    input_prosail = input_prosail,
    spec_prospect = prospect::spec_prospect_full_range,
    spec_soil     = prosail::spec_soil_ossl,
    spec_atm      = prosail::spec_atm
  )
  brf_lut_1nm <- res$surface_refl
  brf_lut <- prosail::apply_sensor_characteristics(
    wvl  = prospect::spec_prospect_full_range$lambda,
    refl = brf_lut_1nm,
    srf  = srf
  )
  brf_lut_noise <- prosail::apply_noise_atbd(brf_lut)
  brf_lut_noise[bands_to_select, ]
}


# ── train_svr_ensemble ────────────────────────────────────────────────────────

#' @title Train a liquidSVM ensemble on a PROSAIL LUT
#'
#' @description
#' Fits a Support Vector Regression (SVR) ensemble of
#' \code{round(n_samples / 100)} models to predict LAI from simulated
#' Sentinel-2 reflectances, using \code{prosail::PROSAIL_Hybrid_Train()}.
#' Each model in the ensemble is immediately serialised via
#' \code{liquidSVM::serialize.liquidSVM()} for safe RDS storage and
#' cross-session restoration.
#'
#' @details
#' The 1/100 scaling of ensemble size relative to sample count reproduces
#' the heuristic from \code{PROSAIL_train_sensitivity_subset.R}
#' (\code{nbmodels <- round(nbSamples/100)}). With \code{n_samples = 1000}
#' this yields 10 ensemble members.
#'
#' The returned list is ready to be passed to \code{saveRDS()} and later
#' restored via \code{lapply(readRDS(file), liquidSVM::unserialize.liquidSVM)}
#' for application in \code{prosail::prosail_hybrid_apply()}.
#'
#' Updated for prosail 3.0.0 (2026-04): uses \code{prosail_hybrid_train()}
#' (replaces \code{PROSAIL_Hybrid_Train}) with renamed arguments
#' \code{refl_lut}/\code{input_variables}/\code{nb_bagg}.
#'
#' @param brf_lut_noise Numeric matrix. Noisy LUT reflectances as returned by
#'   \code{build_prosail_lut()}: dimensions \code{n_bands × n_samples}.
#' @param input_lai     Numeric vector of length \code{n_samples}. LAI values
#'   from \code{input_prosail$lai}: the regression target variable.
#' @param n_samples     Integer. Number of LUT samples (= \code{nrow(input_prosail)}).
#'   Controls ensemble size via \code{round(n_samples / 100)}.
#'
#' @return A list of \code{round(n_samples / 100)} serialised liquidSVM model
#'   objects (output of \code{liquidSVM::serialize.liquidSVM()}).
#'
#' @examples
#' \dontrun{
#' srf   <- prosail::get_srf_sensor("Sentinel_2")
#' bsel  <- match(c("B3","B4","B8"), srf$Spectral_Bands)
#' ip    <- prosail::get_atbd_v3_lut_input(nb_samples = 1000)
#' lut   <- build_prosail_lut(ip, srf, bsel)
#' svr   <- train_svr_ensemble(lut, ip$lai, n_samples = 1000)
#' length(svr)  # 10
#' }
#'
#' @export
train_svr_ensemble <- function(brf_lut_noise, input_lai, n_samples) {
  nb_models <- round(n_samples / 100)
  model_svr <- prosail::prosail_hybrid_train(
    refl_lut        = brf_lut_noise,
    input_variables = input_lai,
    nb_bagg         = nb_models
  )
  lapply(model_svr, liquidSVM::serialize.liquidSVM)
}


# ── train_one_simulation ──────────────────────────────────────────────────────

#' @title Train SVR models for one chunk of the PROSAIL simulation grid
#'
#' @description
#' Worker function called by \code{train_all_simulations()} — either via
#' \code{mapply()} (sequential) or \code{future.apply::future_mapply()}
#' (parallel). For each row of \code{simuset} it overrides the ATBD baseline
#' PROSAIL parameter distributions according to the combination index, generates
#' a LUT with \code{build_prosail_lut()}, trains an SVR ensemble with
#' \code{train_svr_ensemble()}, and writes the serialised model to the
#' corresponding element of \code{filename}.
#'
#' @details
#' **ATBD v3 baseline**: \code{prosail::get_atbd_v3_lut_input()} (prosail 3.0.0,
#' replaces \code{get_InputPROSAIL(atbd=TRUE)}) is called once per
#' \code{train_one_simulation()} call. In sequential mode each chunk is a single
#' row of the grid, so the baseline is always fresh. In parallel mode a chunk may
#' contain multiple rows; the overrides accumulate within a chunk (matching the
#' original \code{PROSAIL_train_sensitivity_subset} behaviour).
#'
#' **Parameter name mapping**: \code{simulations$grid_simu} column names follow
#' legacy camelCase (\code{"LIDFa"}, \code{"LMA"}, \code{"BROWN"}) while the
#' \code{get_atbd_v3_lut_input()} output uses snake_case (\code{"lidf_a"},
#' \code{"lma"}, \code{"brown"}). An internal \code{parm_col_map} vector
#' translates grid names to data.frame column names before overriding.
#'
#' **LiDAR LAI sampling**: \code{lidar_lai_rast} and \code{lidar_lai_site_rast}
#' must be pre-extracted numeric vectors (use \code{terra::values(rast, na.rm=TRUE)}
#' in the calling script). \code{lidar_lai_common_rast} and
#' \code{lidar_lai_3m_rast} may be passed as \code{SpatRaster} objects;
#' \code{terra::values()} is called inside this function. For the standard
#' 270-combination grid \code{("LIDFa","lai","LMA","BROWN")}, only
#' \code{LiDAR_LAI} and \code{LiDAR_LAI_Best_Site_Depth} are active; the other
#' two variants are dead code retained for forward compatibility.
#'
#' This function is a faithful refactor of
#' \code{PROSAIL_train_sensitivity_subset()} from
#' \code{PROSAIL-Optimization/02_CODES/libraries/}, updated for prosail 3.0.0.
#'
#' @param simuset             Data frame. One or more rows from
#'   \code{simulations$grid_simu}; each column is a parameter and each cell is
#'   the 1-based index of the variant to apply.
#' @param filename            Character vector, same length as \code{nrow(simuset)}.
#'   Full paths (including \code{.rds} extension) where each trained model will
#'   be saved.
#' @param combinations        Named list. \code{simulations$combination}: variant
#'   specifications for each parameter.
#' @param lidar_lai_rast      Numeric vector. Pre-extracted non-NA LAI values
#'   from the full LAD stack raster (\code{terra::values(sum(rast(...)), na.rm=TRUE)}).
#' @param lidar_lai_site_rast Numeric vector. Pre-extracted non-NA LAI values
#'   from the site-best-depth raster.
#' @param lidar_lai_common_rast \code{SpatRaster} or \code{NULL}. LAI raster for
#'   the common-best-depth variant. Values extracted internally. Dead code for
#'   the standard grid.
#' @param lidar_lai_3m_rast   \code{SpatRaster} or \code{NULL}. LAI raster for
#'   the 3 m threshold variant. Values extracted internally. Dead code for the
#'   standard grid.
#' @param geom_s2             Named list with elements \code{MinAngle} and
#'   \code{MaxAngle}, each a named numeric vector with fields \code{vza},
#'   \code{sza}, \code{psi} (Sentinel-2 geometry of acquisition).
#' @param nbSamples           Integer. Number of LUT samples per simulation.
#'   Default 2000.
#' @param p                   A \code{progressr} progressor function, or
#'   \code{NULL} (default) if no progress reporting is needed.
#' @param S2BandSelect        Character vector. Sentinel-2 band names to include
#'   in the LUT for inversion. Default \code{c("B3","B4","B8")}.
#'
#' @return Invisibly \code{NULL}. Side effect: one RDS file per row of
#'   \code{simuset} written to the corresponding element of \code{filename}.
#'
#' @export
train_one_simulation <- function(simuset, filename, combinations,
                                  lidar_lai_rast,
                                  lidar_lai_site_rast,
                                  lidar_lai_common_rast,
                                  lidar_lai_3m_rast,
                                  geom_s2,
                                  nbSamples    = 5000,
                                  p            = NULL,
                                  S2BandSelect = c("B3", "B4", "B8")) {
  # Reshape geometry of acquisition into geom_acq list expected by get_atbd_v3_lut_input
  geom_acq <- list(
    min = data.frame(
      tto = geom_s2$MinAngle["vza"],
      tts = geom_s2$MinAngle["sza"],
      psi = geom_s2$MinAngle["psi"]
    ),
    max = data.frame(
      tto = geom_s2$MaxAngle["vza"],
      tts = geom_s2$MaxAngle["sza"],
      psi = geom_s2$MaxAngle["psi"]
    )
  )

  # ATBD v3 baseline distribution — shared starting point for all rows in simuset.
  # Uses get_atbd_v3_lut_input() (prosail 3.0.0) in place of get_InputPROSAIL(atbd=TRUE).
  # The v3 baseline uses updated parameter priors and spec_soil_ossl (48 OSSL spectra).
  input_prosail <- prosail::get_atbd_v3_lut_input(
    nb_samples         = nbSamples,
    geom_acq           = geom_acq,
    codistribution_lai = FALSE
  )

  # Sensor spectral response and band index — computed once per chunk
  srf             <- prosail::get_srf_sensor("Sentinel_2")
  bands_to_select <- match(S2BandSelect, srf$Spectral_Bands)

  # Map from grid parameter names (legacy camelCase from get_combination) to the
  # snake_case column names used by get_atbd_v3_lut_input() output data.frame.
  parm_col_map <- c(LIDFa = "lidf_a", lai = "lai", LMA = "lma", BROWN = "brown")

  # One SVR model per row of simuset
  for (i in seq_len(nrow(simuset))) {
    simuparm <- simuset[i, , drop = FALSE]

    # Apply per-parameter distribution overrides (non-ATBD only)
    for (parm in names(simuparm)) {
      valparm  <- combinations[[parm]][[as.numeric(simuparm[[parm]])]]
      col_name <- parm_col_map[[parm]]
      if (valparm[1] != "ATBD") {
        if (valparm["distribution"] == "uniform") {
          distparm <- stats::runif(
            n   = nrow(input_prosail),
            min = as.numeric(valparm["min"]),
            max = as.numeric(valparm["max"])
          )
        } else if (valparm["distribution"] == "gauss") {
          distparm <- truncnorm::rtruncnorm(
            n    = nrow(input_prosail),
            a    = as.numeric(valparm["min"]),
            b    = as.numeric(valparm["max"]),
            mean = as.numeric(valparm["mean"]),
            sd   = as.numeric(valparm["sd"])
          )
        } else if (valparm["distribution"] == "LiDAR_LAI") {
          # lidar_lai_rast is a pre-extracted numeric vector (terra::values())
          distparm <- sample(lidar_lai_rast, size = nrow(input_prosail))
        } else if (valparm["distribution"] == "LiDAR_LAI_Best_Site_Depth") {
          # lidar_lai_site_rast is a pre-extracted numeric vector (terra::values())
          distparm <- sample(lidar_lai_site_rast, size = nrow(input_prosail))
        } else if (valparm["distribution"] == "LiDAR_LAI_Best_Common_Depth") {
          distparm <- sample(
            terra::values(lidar_lai_common_rast, na.rm = TRUE),
            size = nrow(input_prosail)
          )
        } else if (valparm["distribution"] == "LiDAR_LAI_3m?") {
          distparm <- sample(
            terra::values(lidar_lai_3m_rast, na.rm = TRUE),
            size = nrow(input_prosail)
          )
        }
        input_prosail[[col_name]] <- distparm
      }
    }

    # Generate LUT → apply SRF → add noise → train SVR → serialise → save
    brf_lut_noise <- build_prosail_lut(input_prosail, srf, bands_to_select)
    model_svr     <- train_svr_ensemble(brf_lut_noise, input_prosail$lai, nbSamples)
    saveRDS(model_svr, file = filename[i])
    rm(model_svr)
    if (!is.null(p)) p()
  }
  gc()
  invisible(NULL)
}


# ── train_all_simulations ─────────────────────────────────────────────────────

#' @title Orchestrate PROSAIL SVR training across the full simulation grid
#'
#' @description
#' Dispatches \code{train_one_simulation()} across all rows of the simulation
#' grid, either sequentially (with \code{progressr} progress bar) or in
#' parallel using \code{future}/\code{future.apply} cluster workers. The
#' 270 trained SVR models are written to \code{filename} (one RDS per row of
#' \code{simulations$grid_simu}).
#'
#' This function is a faithful refactor of \code{PROSAIL_train_sensitivity()}
#' from \code{PROSAIL-Optimization/02_CODES/libraries/}. The refactored version
#' additionally forwards \code{lidar_lai_common_rast} and
#' \code{lidar_lai_3m_rast} to the worker for forward compatibility (they are
#' dead code for the standard 270-combination grid but present in the original
#' function signature).
#'
#' @details
#' **Sequential mode** (\code{nbCPU = 1}): the simulation grid is split into
#' 270 single-row chunks and iterated with \code{mapply()}.
#'
#' **Parallel mode** (\code{nbCPU > 1}): the grid is split into at most
#' \code{nbCPU} chunks, each processed by a \code{future} cluster worker via
#' \code{future.apply::future_mapply()}. Packages \code{terra} and
#' \code{prosail} are exported to each worker automatically.
#'
#' @param simulations         Named list. Output of
#'   \code{get_parameter_combinations()}; must contain \code{$grid_simu} and
#'   \code{$combination}.
#' @param geom_s2             Named list. Sentinel-2 geometry of acquisition for
#'   the site, as returned by \code{get_s2_angles()}.
#' @param lidar_lai_rast      Numeric vector or \code{NULL}. Pre-extracted LAI
#'   values from the full LAD stack (via \code{terra::values(sum(rast(...)),
#'   na.rm=TRUE)} in the calling script). Passed to \code{train_one_simulation()}.
#' @param lidar_lai_site_rast Numeric vector or \code{NULL}. Pre-extracted LAI
#'   values from the site-best-depth raster.
#' @param lidar_lai_common_rast \code{SpatRaster} or \code{NULL}. Dead code for
#'   standard grid; retained for forward compatibility.
#' @param lidar_lai_3m_rast   \code{SpatRaster} or \code{NULL}. Dead code for
#'   standard grid; retained for forward compatibility.
#' @param nbSamples           Integer. Number of LUT samples per model.
#'   Default 2000.
#' @param filename            Character vector, length equal to
#'   \code{nrow(simulations$grid_simu)}. Full paths (with \code{.rds}) where
#'   models will be saved by \code{train_one_simulation()}.
#' @param nbCPU               Integer. Number of CPU cores to use. Default 1
#'   (sequential). Values > 1 enable \code{future} parallel execution.
#' @param S2BandSelect        Character vector. Sentinel-2 band names for
#'   inversion. Default \code{c("B3","B4","B8")}.
#'
#' @return Invisibly \code{NULL}. Side effect: 270 RDS files written to the
#'   paths given in \code{filename}.
#'
#' @examples
#' \dontrun{
#' simulations <- get_parameter_combinations(
#'   parms2test    = c("LIDFa","lai","LMA","BROWN"),
#'   output_dir    = here::here("revision","output","intermediate","PROSAIL_Models","Blois"),
#'   name_strategy = "LIDFa_lai_LMA_BROWN"
#' )
#' train_all_simulations(
#'   simulations = simulations,
#'   geom_s2     = geom_s2[["Blois"]],
#'   lidar_lai_rast      = lai_vals,
#'   lidar_lai_site_rast = lai_site_vals,
#'   nbSamples   = 1000,
#'   filename    = filename_svr,
#'   nbCPU       = 4
#' )
#' }
#'
#' @export
train_all_simulations <- function(simulations,
                                   geom_s2,
                                   lidar_lai_rast        = NULL,
                                   lidar_lai_site_rast   = NULL,
                                   lidar_lai_common_rast = NULL,
                                   lidar_lai_3m_rast     = NULL,
                                   nbSamples             = 5000,
                                   filename,
                                   nbCPU                 = 1,
                                   S2BandSelect          = c("B3", "B4", "B8")) {
  nb_split <- if (nbCPU == 1L) {
    nrow(simulations$grid_simu)
  } else {
    min(nrow(simulations$grid_simu), nbCPU)
  }

  filename_split <- suppressWarnings(split(x = filename,
                                           f = seq(nb_split)))
  simu_split     <- suppressWarnings(split(x = simulations$grid_simu,
                                           f = seq(nb_split)))

  more_args <- list(
    combinations          = simulations$combination,
    lidar_lai_rast        = lidar_lai_rast,
    lidar_lai_site_rast   = lidar_lai_site_rast,
    lidar_lai_common_rast = lidar_lai_common_rast,
    lidar_lai_3m_rast     = lidar_lai_3m_rast,
    geom_s2               = geom_s2,
    S2BandSelect          = S2BandSelect,
    nbSamples             = nbSamples
  )

  if (nbCPU == 1L) {
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(simulations$grid_simu))
      mapply(
        FUN      = train_one_simulation,
        simuset  = simu_split,
        filename = filename_split,
        MoreArgs = c(more_args, list(p = p)),
        SIMPLIFY = FALSE
      )
    })
  } else {
    cl <- parallel::makeCluster(nbCPU)
    future::plan("cluster", workers = cl)
    progressr::handlers("cli")
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(simulations$grid_simu))
      future.apply::future_mapply(
        FUN             = train_one_simulation,
        simuset         = simu_split,
        filename        = filename_split,
        MoreArgs        = c(more_args, list(p = p)),
        future.packages = c("terra", "prosail"),
        future.seed     = TRUE,
        SIMPLIFY        = FALSE
      )
    })
    parallel::stopCluster(cl)
    future::plan("sequential")
  }
  invisible(NULL)
}
