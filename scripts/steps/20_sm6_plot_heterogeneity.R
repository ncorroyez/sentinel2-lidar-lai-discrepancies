# ---
# title:  20_sm6_plot_heterogeneity.R
# desc:   Visualisation — SM6b heterogeneity class analysis.
#         Style inspired by plot_all_metrics_grid() in
#         02_CODES/three_factors_analysis_final.R.
#
#         For each metric_source (DSM, CHM): a 4-row × 3-column facet grid
#         (rows: r, RMSE, Bias, Slope; columns: Aigoual, Blois, Mormal)
#         with x = heterogeneity class (Low/Medium/High) and three coloured
#         lines for the three LAI comparison combinations:
#
#           ATBD_vs_ALS       LAI[S2_ATBD] vs LAI[ALS]
#           ATBD_vs_ALS_dopt  LAI[S2_ATBD] vs LAI[ALS_dopt]
#           opt_vs_ALS_dopt   LAI[S2_opt]  vs LAI[ALS_dopt]
#
#         Reference hlines: Bias → y = 0, Slope → y = 1.
#
#         Reads:
#           output/intermediate/sm6/heterogeneity_analysis.csv
#
#         Figures written to output/figures/:
#           het_metrics_grid_DSM.pdf
#           het_metrics_grid_CHM.pdf
#           het_metrics_grid_DSM.png
#           het_metrics_grid_CHM.png
#
# Run from the project root (NC_Full/):
#   source("scripts/20_sm6_plot_heterogeneity.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("cli")

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_plot_heterogeneity.R"))

# ── d_opt modes to plot — produces one figure per mode ────────────────────────
dopt_modes_to_plot <- c("per_site", "all_sites")

# Helper: load the SM6a per-mode CSV produced by script 15
load_mode_long <- function(mode) {
  fn <- file.path(paths$output, "intermediate", "sm6",
                   paste0("heterogeneity_analysis_", mode, ".csv"))
  load_heter_csv_long(fn)
}

# ── Output directory ───────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Produce one figure per (mode × metric_source) ─────────────────────────────

for (current_mode in dopt_modes_to_plot) {
  cli::cli_h2("Mode: {current_mode}")
  long_dt <- load_mode_long(current_mode)
  cli::cli_alert_info("Loaded {nrow(long_dt)} rows for mode '{current_mode}'.")

  for (ms in c("DSM", "CHM")) {
    sub_ms <- long_dt[metric_source == ms]
    if (nrow(sub_ms) == 0L) {
      cli::cli_warn("No data for metric_source = {ms} (mode = {current_mode}) — skipping.")
      next
    }
    p_ms <- plot_het_grid(sub_ms, ms_short = ms)

    out_pdf <- file.path(out_dir,
      paste0("het_metrics_grid_", ms, "_", current_mode, ".pdf"))
    out_png <- file.path(out_dir,
      paste0("het_metrics_grid_", ms, "_", current_mode, ".png"))

    ggplot2::ggsave(out_pdf, p_ms, width = 26, height = 22, units = "cm",
                    device = cairo_pdf)
    ggplot2::ggsave(out_png, p_ms, width = 26, height = 22, units = "cm",
                    device = "png", dpi = 600, bg = "white")

    cli::cli_alert_success("Written: {out_pdf}")
    cli::cli_alert_success("Written: {out_png}")
  }
}

cli::cli_h1("Done — SM6 heterogeneity figures")

# ── Console summary: raw metrics per metric_source x combination x site x level ─

fmt_pval <- function(p) {
  if (is.na(p) || length(p) == 0L) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

for (current_mode in dopt_modes_to_plot) {
  fn <- file.path(paths$output, "intermediate", "sm6",
                   paste0("heterogeneity_analysis_", current_mode, ".csv"))
  if (!file.exists(fn)) next
  dt_mode <- data.table::fread(fn)
  for (ms in unique(dt_mode$metric_source)) {
    cli::cli_h1("[{current_mode}] Raw metrics -- metric_source: {ms}")

    for (s in site_levels) {
      cli::cli_h2("{s}")
      dt_s <- dt_mode[metric_source == ms & site == s]

      for (combo in combo_levels) {
        cli::cli_h3("{combo}")
        dt_c <- dt_s[combination == combo]

        for (hc in het_levels) {
          r <- dt_c[het_class == hc]
          if (nrow(r) == 0L) next
          line <- sprintf(
            "  %-7s (n=%s)  r=%.3f  R2=%.3f  RMSE=%.3f  Bias=%.3f (p %s)  Slope=%.3f (p %s)",
            hc,
            formatC(r$n, format = "d", big.mark = ","),
            r$R, r$R2, r$RMSE,
            r$Bias,  fmt_pval(r$Bias_pvalue),
            r$Slope, fmt_pval(r$Slope_pvalue)
          )
          message(line)
        }
      }
    }
  }
}
