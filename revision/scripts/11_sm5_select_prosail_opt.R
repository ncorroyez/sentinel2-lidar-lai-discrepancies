# ---
# title:  06b_sm5_select_prosail_opt.R
# desc:   Orchestration — SM5 passe 2.
#         Selects the optimal PROSAIL configuration (LAI_S2_opt) under two
#         d_opt modes, both at the same time:
#
#           per_site  : each site uses its own d_opt (from dopt_reference.csv)
#                       → Pareto selection over PROSAIL columns in that site's
#                         data at its site-specific optimal depth.
#
#           all_sites : the d_opt from the "All_sites" row of dopt_reference.csv
#                       is fixed and applied to each individual site separately
#                       → Pareto selection over PROSAIL columns in each site's
#                         data at the common cross-site depth.
#
#         Both modes are written to the same CSV, distinguished by the
#         d_opt_source column. The selected Column_opt becomes LAI_S2_opt for
#         downstream analyses (SM5 passe 3, SM6b).
#
#         Writes:
#           revision/output/intermediate/sm5/prosail_opt.csv
#             One row per (d_opt_source × Site × Norm).
#
#         Prerequisites:
#           - revision/output/intermediate/sm5/dopt_reference.csv
#             (script 06, ATBD reference d_opt — must include All_sites rows)
#           - revision/output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#             (script 05, full metrics table)
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/06b_sm5_select_prosail_opt.R")
# ---

# ── Configuration ──────────────────────────────────────────────────────────────

dopt_method      <- "pareto"    # one of: "pearson", "rmse", "bias", "slope", "pareto"
norm_filter      <- "keepTrees" # grepl() pattern — set to NULL to keep all norms
individual_sites <- c("Aigoual", "Blois", "Mormal")
fixed_d_opt      <- 4L          # reference depth from submitted paper; NULL to skip
# h_min and lai_scenario to analyse — change for sensitivity analysis.
h_min_select     <- 10L
lai_scenario_select <- "per_site"  # one of: "per_site", "common", "fixed_4"

# ── Prerequisites ──────────────────────────────────────────────────────────────

dopt_ref_csv <- here::here(
  "revision", "output", "intermediate", "sm5", "dopt_reference.csv"
)
sm5_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "all_results_combined_LIDFa_lai_LMA_BROWN.csv"
)

for (f in c(dopt_ref_csv, sm5_csv)) {
  if (!file.exists(f)) {
    stop("Required file not found:\n  ", f,
         "\nRun scripts 06 and 10 first.")
  }
}

# ── Source and load ────────────────────────────────────────────────────────────

source(here::here("revision", "R", "sm5_dopt.R"))

dopt_ref   <- data.table::fread(dopt_ref_csv)
metrics_dt <- data.table::fread(sm5_csv)

# Filter to the selected h_min and LAI scenario
metrics_dt <- metrics_dt[h_min == h_min_select & lai_scenario == lai_scenario_select]

# ── Select optimal PROSAIL configuration (both modes) ─────────────────────────

t0 <- proc.time()

prosail_opt <- select_prosail_opt(
  metrics_dt       = metrics_dt,
  dopt_dt          = dopt_ref,
  dopt_method      = dopt_method,
  norm_filter      = norm_filter,
  individual_sites = individual_sites,
  fixed_d_opt      = fixed_d_opt
)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Write CSV ──────────────────────────────────────────────────────────────────

out_dir  <- here::here("revision", "output", "intermediate", "sm5")
out_path <- file.path(out_dir, "prosail_opt.csv")
data.table::fwrite(prosail_opt, out_path)

# ── Console summary ────────────────────────────────────────────────────────────

cli::cli_h1("SM5 passe 2 — PROSAIL optimal configuration")
cli::cli_bullets(c(
  "v" = "prosail_opt.csv — {nrow(prosail_opt)} lignes",
  "i" = "d_opt method    : {dopt_method}",
  "i" = "Norm filter     : {norm_filter %||% 'none'}",
  "i" = "Durée           : {elapsed} s"
))

display <- data.table::copy(prosail_opt)
data.table::setorder(display, d_opt_source, Norm, Site)

cli::cli_h2("Mode : per_site (d_opt propre à chaque site)")
print(display[d_opt_source == "per_site",
              .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, dist_utopia)])

cli::cli_h2("Mode : all_sites (d_opt All_sites appliqué à chaque site)")
print(display[d_opt_source == "all_sites",
              .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, dist_utopia)])

if (!is.null(fixed_d_opt)) {
  cli::cli_h2("Mode : fixed_{fixed_d_opt} (d_opt={fixed_d_opt} — référence papier soumis)")
  print(display[d_opt_source == paste0("fixed_", fixed_d_opt),
                .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, dist_utopia)])
}
