# ---
# title:  pipeline_launcher.R
# desc:   Generates pipeline.RData at the project root.
#         Load it once in any R session to get convenience run_*() functions
#         without having to remember script paths.
#
#         Quick start:
#           load(here::here("pipeline.RData"))   # or double-click in RStudio
#           run_all()                             # full pipeline
#           run_prepare()                         # Phase 1 only
#           run_dopt()                            # Phase 2 only
#           run_phase("03_prosail_compute")       # any individual script
#
#         The .RData contains only function closures — no embedded data.
#         config.yml is re-read at each call, so switching local ↔ smb
#         profile takes effect immediately.
#
#         To regenerate pipeline.RData after editing this file:
#           source("revision/pipeline_launcher.R")
#
# Run from project root (NC_Full/):
#   source("revision/pipeline_launcher.R")
# ---

library(here)
library(cli)

# ── Internal helper ────────────────────────────────────────────────────────────

.run_script <- function(rel_path) {
  path  <- here::here(rel_path)
  label <- sub("\\.R$", "", basename(path))
  cli::cli_rule(left = "{label}")
  t0 <- proc.time()
  source(path, local = FALSE)
  elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
  cli::cli_alert_success("{label} — {elapsed} min")
  invisible(elapsed)
}

# ── Public functions ───────────────────────────────────────────────────────────

#' Run the complete pipeline (all 5 phases)
run_all <- function() {
  .run_script("revision/scripts/00_run_all.R")
}

#' Phase 1 — SVR ATBD training, S2 pixel sampling, PAD extraction
run_prepare <- function() {
  .run_script("revision/scripts/01_prepare.R")
}

#' Phase 2 — d_opt selection + LAI_ALS_dopt rasters
#' @param figures  If TRUE (default), also run figure scripts
run_dopt <- function(figures = TRUE) {
  .run_script("revision/scripts/02_dopt_compute.R")
  if (figures) .run_script("revision/scripts/02_dopt_figures.R")
}

#' Phase 3 — PROSAIL full inversion + LAI_S2_opt rasters
#' @param figures  If TRUE (default), also run figure scripts
run_prosail <- function(figures = TRUE) {
  .run_script("revision/scripts/03_prosail_compute.R")
  if (figures) .run_script("revision/scripts/03_prosail_figures.R")
}

#' Phase 4 — Heterogeneity analysis (SM6)
#' @param figures  If TRUE (default), also run figure scripts
run_het <- function(figures = TRUE) {
  .run_script("revision/scripts/04_het_compute.R")
  if (figures) .run_script("revision/scripts/04_het_figures.R")
}

#' Phase 5 — Sensitivity analyses (k, fCover, h_min, joint)
#' @param figures  If TRUE (default), also run figure scripts
run_sensitivity <- function(figures = TRUE) {
  .run_script("revision/scripts/05_sensitivity.R")
  if (figures) .run_script("revision/scripts/05_sensitivity_figures.R")
}

#' Run an arbitrary pipeline script by number or name
#' @param script  Number (e.g. 3) or script name without .R
#'                (e.g. "03_prosail_compute")
run_phase <- function(script) {
  scripts_dir <- here::here("revision/scripts")
  if (is.numeric(script)) {
    pattern <- sprintf("^%02d_", script)
    candidates <- list.files(scripts_dir, pattern = pattern, full.names = FALSE)
    if (length(candidates) == 0L)
      stop("No script found matching phase number ", script)
    if (length(candidates) > 1L) {
      cli::cli_bullets(c("i" = "Multiple matches: {candidates}"))
      cli::cli_abort("Specify a full name to disambiguate.")
    }
    .run_script(file.path("revision/scripts", candidates[[1L]]))
  } else {
    name <- if (grepl("\\.R$", script)) script else paste0(script, ".R")
    .run_script(file.path("revision/scripts", name))
  }
}

#' List available pipeline scripts with their descriptions
pipeline_help <- function() {
  scripts_dir <- here::here("revision/scripts")
  files <- list.files(scripts_dir, pattern = "^0[0-9].*\\.R$",
                      full.names = TRUE)
  cli::cli_h1("Pipeline scripts")
  for (f in sort(files)) {
    hdr <- readLines(f, n = 5L, warn = FALSE)
    desc_line <- grep("^# desc:", hdr, value = TRUE)
    desc <- if (length(desc_line) > 0L)
      sub("^# desc:\\s*", "", desc_line[[1L]])
    else ""
    cli::cli_bullets(c("*" = "{basename(f)}  {.emph {desc}}"))
  }
  cli::cli_text("")
  cli::cli_text("Usage: run_phase(1) or run_phase(\"02_dopt_compute\")")
}

# ── Save to pipeline.RData ─────────────────────────────────────────────────────

launcher_env <- new.env(parent = emptyenv())
launcher_env$.run_script    <- .run_script
launcher_env$run_all        <- run_all
launcher_env$run_prepare    <- run_prepare
launcher_env$run_dopt       <- run_dopt
launcher_env$run_prosail    <- run_prosail
launcher_env$run_het        <- run_het
launcher_env$run_sensitivity <- run_sensitivity
launcher_env$run_phase      <- run_phase
launcher_env$pipeline_help  <- pipeline_help

out_rdata <- here::here("pipeline.RData")
save(list = ls(launcher_env), envir = launcher_env, file = out_rdata)

cli::cli_alert_success("pipeline.RData written — {out_rdata}")
cli::cli_text("")
cli::cli_text("Load with: {.code load(here::here(\"pipeline.RData\"))}")
cli::cli_text("Then call: {.code run_all()}, {.code run_dopt()}, {.code pipeline_help()}, etc.")
