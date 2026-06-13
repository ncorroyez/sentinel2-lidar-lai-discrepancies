# ---
# title:  20e_plot_sm5_heter_classes.R
# desc:   Plot the SM5-sample discrete heterogeneity analysis produced by
#         scripts/steps/20e_sm5_heter_classes.R. Uses the shared
#         plot_het_grid() helper so the figure style is identical to
#         scripts/steps/20_sm6_plot_heterogeneity.R (which uses the SM6a
#         pixel sample).
#
#         Reads:
#           output/intermediate/sm6/heterogeneity_analysis_sm5_{mode}.csv
#
#         Figures written to output/figures/:
#           het_metrics_grid_sm5_DSM_per_site.{pdf,png}
#           het_metrics_grid_sm5_DSM_all_sites.{pdf,png}
#
# Run from project root (NC_Full/):
#   source("scripts/steps/20e_plot_sm5_heter_classes.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_plot_heterogeneity.R"))

dopt_modes_to_plot <- c("per_site", "all_sites")

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (current_mode in dopt_modes_to_plot) {
  cli::cli_h2("Mode: {current_mode}")

  csv_path <- file.path(paths$output, "intermediate", "sm6",
                        paste0("heterogeneity_analysis_sm5_",
                               current_mode, ".csv"))
  long_dt  <- load_heter_csv_long(csv_path)
  cli::cli_alert_info("Loaded {nrow(long_dt)} rows for mode '{current_mode}'.")

  for (ms in unique(long_dt$metric_source)) {
    sub_ms <- long_dt[metric_source == ms]
    if (nrow(sub_ms) == 0L) next
    p_ms <- plot_het_grid(sub_ms, ms_short = ms)

    out_pdf <- file.path(out_dir,
      paste0("het_metrics_grid_sm5_", ms, "_", current_mode, ".pdf"))
    out_png <- file.path(out_dir,
      paste0("het_metrics_grid_sm5_", ms, "_", current_mode, ".png"))

    ggplot2::ggsave(out_pdf, p_ms, width = 26, height = 22, units = "cm",
                    device = cairo_pdf)
    ggplot2::ggsave(out_png, p_ms, width = 26, height = 22, units = "cm",
                    device = "png", dpi = 600, bg = "white")

    cli::cli_alert_success("Written: {out_pdf}")
    cli::cli_alert_success("Written: {out_png}")
  }
}

cli::cli_h1("Done -- SM5-sample heterogeneity class figures")
