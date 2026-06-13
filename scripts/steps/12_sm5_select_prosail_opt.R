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
#           output/intermediate/sm5/prosail_opt.csv
#             One row per (d_opt_source × Site × Norm).
#
#         Prerequisites:
#           - output/intermediate/sm5/dopt_reference.csv
#             (script 06, ATBD reference d_opt — must include All_sites rows)
#           - output/intermediate/sm5/
#             all_results_combined_LIDFa_lai_LMA_BROWN.csv
#             (script 05, full metrics table)
#
# Run from the project root (NC_Full/):
#   source("scripts/06b_sm5_select_prosail_opt.R")
# ---

library(here)
library(data.table)
library(cli)

source(here::here("R", "paths.R"))

# ── Configuration ──────────────────────────────────────────────────────────────

dopt_method      <- "pareto"    # one of: "pearson", "rmse", "bias", "slope", "pareto"
norm_filter      <- "keepTrees" # grepl() pattern — set to NULL to keep all norms
individual_sites <- c("Aigoual", "Blois", "Mormal")
fixed_d_opt      <- 4L          # reference depth from submitted paper; NULL to skip
# h_min and lai_scenario to analyse — change for sensitivity analysis.
h_min_select     <- 10L
lai_scenario_select <- "common"    # one of: "per_site", "common", "fixed_4"

# Pareto criteria + normalisation method —
# MUST match scripts/steps/06_sm5_select_dopt.R for consistency.
pareto_criteria <- c("R", "RMSE", "Bias", "Slope")
# pareto_criteria <- c("R", "RMSE", "Bias")
pareto_norm_method <- "minmax"
# pareto_norm_method <- "max"

# ── Prerequisites ──────────────────────────────────────────────────────────────

dopt_ref_csv <- file.path(paths$output, "intermediate", "sm5", "dopt_reference.csv"
)
sm5_csv <- file.path(paths$output, "intermediate", "sm5",
  "all_results_combined_LIDFa_lai_LMA_BROWN.csv"
)

for (f in c(dopt_ref_csv, sm5_csv)) {
  if (!file.exists(f)) {
    stop("Required file not found:\n  ", f,
         "\nRun scripts 06 and 10 first.")
  }
}

# ── Source and load ────────────────────────────────────────────────────────────

source(here::here("R", "sm5_dopt.R"))

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
  fixed_d_opt      = fixed_d_opt,
  pareto_criteria    = pareto_criteria,
  pareto_norm_method = pareto_norm_method
)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)

# ── Write CSV ──────────────────────────────────────────────────────────────────

out_dir  <- file.path(paths$output, "intermediate", "sm5")
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
              .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, Avg_score)])

cli::cli_h2("Mode : all_sites (d_opt All_sites appliqué à chaque site)")
print(display[d_opt_source == "all_sites",
              .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, Avg_score)])

if (!is.null(fixed_d_opt)) {
  cli::cli_h2("Mode : fixed_{fixed_d_opt} (d_opt={fixed_d_opt} — référence papier soumis)")
  print(display[d_opt_source == paste0("fixed_", fixed_d_opt),
                .(Site, Norm, d_opt, Column_opt, R2, RMSE, Bias, Slope, Avg_score)])
}

# ── Markdown summary (persisté pour relecture rapide) ────────────────────────
# Decodes the "LIDFa=X_lai=Y_LMA=Z_BROWN=W" string into human-readable
# parameter priors (matches the level definitions in R/prosail_lut.R).
md_path <- file.path(out_dir, "prosail_opt_summary.md")

LIDFa_lbl <- c(
  "1" = "ATBD",
  "2" = "gauss(35, sd=20, [20,50])",
  "3" = "gauss(40, sd=20, [25,55])",
  "4" = "gauss(40, sd=20, [30,55])",
  "5" = "gauss(45, sd=20, [30,60])"
)
lai_lbl <- c(
  "1" = "ATBD",
  "2" = "gauss(4, sd=3, [0,9])",
  "3" = "gauss(5, sd=3, [0,11])",
  "4" = "LiDAR_LAI",
  "5" = "LiDAR_LAI_Best_Site_Depth"
)
LMA_lbl <- c(
  "1" = "ATBD",
  "2" = "gauss(0.010, sd=0.005, [0.003,0.02])",
  "3" = "gauss(0.015, sd=0.005, [0.003,0.03])"
)
BROWN_lbl <- c(
  "1" = "ATBD",
  "2" = "gauss(0, sd=0.1, [0,1])",
  "3" = "uniform(0,0)"
)

decode_col <- function(col_opt) {
  m <- regmatches(col_opt, regexec(
    "LIDFa=([0-9]+)_lai=([0-9]+)_LMA=([0-9]+)_BROWN=([0-9]+)", col_opt))[[1L]]
  list(LIDFa = m[2], lai = m[3], LMA = m[4], BROWN = m[5])
}

write_mode_table <- function(con, mode_label, mode_key) {
  cat("\n## Mode : ", mode_label, "\n\n", sep = "", file = con, append = TRUE)
  cat("| Site | d_opt (m) | LIDFa | lai | LMA | BROWN | R | RMSE | Bias | Slope | Avg score |\n",
      file = con, append = TRUE)
  cat("|---|---|---|---|---|---|---|---|---|---|---|\n",
      file = con, append = TRUE)
  sub <- display[d_opt_source == mode_key & Norm == norm_filter]
  for (i in seq_len(nrow(sub))) {
    row <- sub[i]; d <- decode_col(row$Column_opt)
    cat(sprintf(
      "| %s | %d | %s (lvl=%s) | %s (lvl=%s) | %s (lvl=%s) | %s (lvl=%s) | %.2f | %.2f | %.2f | %.2f | %.3f |\n",
      row$Site, row$d_opt,
      LIDFa_lbl[d$LIDFa], d$LIDFa,
      lai_lbl[d$lai],     d$lai,
      LMA_lbl[d$LMA],     d$LMA,
      BROWN_lbl[d$BROWN], d$BROWN,
      row$R, row$RMSE, row$Bias, row$Slope, row$Avg_score
    ), file = con, append = TRUE)
  }
}

if (file.exists(md_path)) file.remove(md_path)
cat("# PROSAIL Pareto selection — summary\n\n",
    "Generated by `scripts/steps/12_sm5_select_prosail_opt.R`.  \n",
    "Norm filter: `", norm_filter %||% "none", "`  \n",
    "d_opt method: `", dopt_method, "`\n",
    sep = "", file = md_path, append = TRUE)

write_mode_table(md_path, "per_site (each site uses its own d_opt)", "per_site")
write_mode_table(md_path, "all_sites (d_opt averaged across sites)", "all_sites")
if (!is.null(fixed_d_opt))
  write_mode_table(md_path,
                   paste0("fixed_", fixed_d_opt, " (paper-submitted reference)"),
                   paste0("fixed_", fixed_d_opt))

cat("\n## Parameter-level legend\n\n", file = md_path, append = TRUE)
cat("**LIDFa** (mean leaf angle, deg):\n", file = md_path, append = TRUE)
for (k in names(LIDFa_lbl))
  cat(sprintf("- lvl %s : %s\n", k, LIDFa_lbl[[k]]), file = md_path, append = TRUE)
cat("\n**lai**:\n", file = md_path, append = TRUE)
for (k in names(lai_lbl))
  cat(sprintf("- lvl %s : %s\n", k, lai_lbl[[k]]), file = md_path, append = TRUE)
cat("\n**LMA** (g/cm²):\n", file = md_path, append = TRUE)
for (k in names(LMA_lbl))
  cat(sprintf("- lvl %s : %s\n", k, LMA_lbl[[k]]), file = md_path, append = TRUE)
cat("\n**BROWN**:\n", file = md_path, append = TRUE)
for (k in names(BROWN_lbl))
  cat(sprintf("- lvl %s : %s\n", k, BROWN_lbl[[k]]), file = md_path, append = TRUE)

cli::cli_alert_success("Markdown summary written: {md_path}")
