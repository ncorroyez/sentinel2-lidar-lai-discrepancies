# ---
# title:  18c_fcover_sensitivity_plot.R
# desc:   Supplementary figure — sensitivity of LAI comparison metrics
#         (R, RMSE, Bias, Slope) to fCover threshold (0.80 / 0.90 / 0.95).
#         Reads fcover_sensitivity_atbd.csv produced by 18_fcover_sensitivity.R.
#
#         Output:
#           revision/output/figures/sm_fcover_sensitivity.pdf
#           revision/output/figures/sm_fcover_sensitivity.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/18c_fcover_sensitivity_plot.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

# ── I/O ────────────────────────────────────────────────────────────────────────

in_csv  <- here::here(
  "revision", "output", "intermediate", "reviewers",
  "fcover_sensitivity_atbd.csv"
)
out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dt <- data.table::fread(in_csv)

# ── Prepare ────────────────────────────────────────────────────────────────────

site_order  <- c("Aigoual", "Blois", "Mormal")
site_colors <- c(Aigoual = "#E69F00", Blois = "#56B4E9", Mormal = "#009E73")

dt[, site           := factor(site, levels = site_order)]
dt[, fcover_pct     := fcover_threshold * 100]
dt[, abs_Bias       := abs(Bias)]
dt[, abs_Slope1     := abs(Slope - 1)]

ref_thr <- 90   # reference fCover threshold (%)

# Metrics to display and their y-axis labels
metric_list <- list(
  list(col = "R",          label = expression(italic(r))),
  list(col = "RMSE",       label = expression(RMSE~(m^2~m^{-2}))),
  list(col = "Bias",       label = expression(Bias~(m^2~m^{-2}))),
  list(col = "Slope",      label = expression(Slope))
)

# ── Build one panel per metric ─────────────────────────────────────────────────

make_panel <- function(col, ylab, show_legend = FALSE) {
  p <- ggplot2::ggplot(
    dt,
    ggplot2::aes(
      x     = fcover_pct,
      y     = .data[[col]],
      colour = site,
      group  = site
    )
  ) +
    ggplot2::geom_vline(
      xintercept = ref_thr,
      linetype   = "dashed",
      colour     = "grey60",
      linewidth  = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
    ggplot2::scale_x_continuous(
      breaks = c(80, 90, 95),
      labels = c("80%", "90%", "95%")
    ) +
    ggplot2::scale_colour_manual(
      values = site_colors,
      name   = "Site"
    ) +
    ggplot2::labs(x = "fCover threshold", y = ylab) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor    = ggplot2::element_blank(),
      legend.position     = if (show_legend) "right" else "none",
      axis.title.y        = ggplot2::element_text(size = 9)
    )

  # Add horizontal reference line at y = 0 for Bias, y = 1 for Slope
  if (col == "Bias")  p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3)
  if (col == "Slope") p <- p + ggplot2::geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3)

  p
}

panels <- lapply(seq_along(metric_list), function(i) {
  m <- metric_list[[i]]
  make_panel(m$col, m$label, show_legend = (i == length(metric_list)))
})

# ── Assemble 2×2 ──────────────────────────────────────────────────────────────

fig <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Sensitivity of LAI_ALS vs LAI_S2 metrics to fCover threshold",
    subtitle = "Dashed line: reference threshold (90%). ATBD parametrisation.",
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40")
    )
  ) &
  ggplot2::theme(legend.position = "right")

# ── Save ───────────────────────────────────────────────────────────────────────

w_cm <- 22
h_cm <- 16

ggplot2::ggsave(
  file.path(out_dir, "sm_fcover_sensitivity.pdf"),
  fig, width = w_cm, height = h_cm, units = "cm"
)
ggplot2::ggsave(
  file.path(out_dir, "sm_fcover_sensitivity.png"),
  fig, width = w_cm, height = h_cm, units = "cm", dpi = 300
)

cli::cli_alert_success("Saved sm_fcover_sensitivity.pdf/png ({w_cm}×{h_cm} cm)")
