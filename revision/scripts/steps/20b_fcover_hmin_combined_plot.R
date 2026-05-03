# ---
# title:  20b_fcover_hmin_combined_plot.R
# desc:   Supplementary figure (B) — two independent 1D sensitivity analyses
#         displayed on a single plate:
#           Row 1: sensitivity to fCover threshold (0.80 / 0.90 / 0.95)
#           Row 2: sensitivity to h_min (2 / 3 / 4 / 5 m)
#         Each row shows 4 comparison metrics: r, RMSE, Bias, Slope.
#         Both axes reflect independent parameter sweeps at the other
#         parameter's reference value (fCover = 0.90; h_min = 2 m).
#
#         Reads:
#           revision/output/intermediate/reviewers/fcover_sensitivity_atbd.csv
#           revision/output/intermediate/reviewers/h_min_sensitivity_atbd.csv
#
#         Output:
#           revision/output/figures/sm_fcover_hmin_combined.pdf
#           revision/output/figures/sm_fcover_hmin_combined.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/20b_fcover_hmin_combined_plot.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

# ── I/O ────────────────────────────────────────────────────────────────────────

int_dir <- file.path(paths$output, "intermediate", "reviewers")
out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dt_fc   <- data.table::fread(file.path(int_dir, "fcover_sensitivity_atbd.csv"))
dt_hmin <- data.table::fread(file.path(int_dir, "h_min_sensitivity_atbd.csv"))

# ── Shared aesthetics ──────────────────────────────────────────────────────────

site_order  <- c("Aigoual", "Blois", "Mormal")
site_colors <- c(Aigoual = "#E69F00", Blois = "#56B4E9", Mormal = "#009E73")

metric_list <- list(
  list(col = "R",     label = expression(italic(r))),
  list(col = "RMSE",  label = expression(RMSE~(m^2~m^{-2}))),
  list(col = "Bias",  label = expression(Bias~(m^2~m^{-2}))),
  list(col = "Slope", label = expression(Slope))
)

theme_sens <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor  = ggplot2::element_blank(),
    legend.position   = "none",
    axis.title.y      = ggplot2::element_text(size = 9)
  )

# ── Row 1: fCover sensitivity ──────────────────────────────────────────────────

dt_fc[, site      := factor(site, levels = site_order)]
dt_fc[, fcover_pct := fcover_threshold * 100]

make_fcover_panel <- function(col, ylab, show_legend = FALSE) {
  p <- ggplot2::ggplot(
    dt_fc,
    ggplot2::aes(x = fcover_pct, y = .data[[col]],
                 colour = site, group = site)
  ) +
    ggplot2::geom_vline(xintercept = 90, linetype = "dashed",
                        colour = "grey60", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
    ggplot2::scale_x_continuous(breaks = c(80, 90, 95),
                                labels = c("80%", "90%", "95%")) +
    ggplot2::scale_colour_manual(values = site_colors, name = "Site") +
    ggplot2::labs(x = "fCover threshold", y = ylab) +
    theme_sens +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")

  if (col == "Bias")  p <- p + ggplot2::geom_hline(yintercept = 0,
                                colour = "grey40", linewidth = 0.3)
  if (col == "Slope") p <- p + ggplot2::geom_hline(yintercept = 1,
                                colour = "grey40", linewidth = 0.3)
  p
}

# ── Row 2: h_min sensitivity ───────────────────────────────────────────────────

dt_hmin[, site := factor(site, levels = site_order)]

make_hmin_panel <- function(col, ylab, show_legend = FALSE) {
  p <- ggplot2::ggplot(
    dt_hmin,
    ggplot2::aes(x = h_min_value, y = .data[[col]],
                 colour = site, group = site)
  ) +
    ggplot2::geom_vline(xintercept = 2, linetype = "dashed",
                        colour = "grey60", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
    ggplot2::scale_x_continuous(breaks = c(2, 3, 4, 5)) +
    ggplot2::scale_colour_manual(values = site_colors, name = "Site") +
    ggplot2::labs(x = expression(h[min]~(m)), y = ylab) +
    theme_sens +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")

  if (col == "Bias")  p <- p + ggplot2::geom_hline(yintercept = 0,
                                colour = "grey40", linewidth = 0.3)
  if (col == "Slope") p <- p + ggplot2::geom_hline(yintercept = 1,
                                colour = "grey40", linewidth = 0.3)
  p
}

# ── Assemble 2 × 4 plate ──────────────────────────────────────────────────────

row1 <- lapply(seq_along(metric_list), function(i) {
  m <- metric_list[[i]]
  make_fcover_panel(m$col, m$label, show_legend = (i == length(metric_list)))
})

row2 <- lapply(seq_along(metric_list), function(i) {
  m <- metric_list[[i]]
  make_hmin_panel(m$col, m$label, show_legend = (i == length(metric_list)))
})

fig <- (
  (row1[[1]] | row1[[2]] | row1[[3]] | row1[[4]]) /
  (row2[[1]] | row2[[2]] | row2[[3]] | row2[[4]])
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation() &
  ggplot2::theme(legend.position = "right")

# ── Save ───────────────────────────────────────────────────────────────────────

w_cm <- 28
h_cm <- 14

ggplot2::ggsave(file.path(out_dir, "sm_fcover_hmin_combined.pdf"),
                fig, width = w_cm, height = h_cm, units = "cm", device = cairo_pdf)
ggplot2::ggsave(file.path(out_dir, "sm_fcover_hmin_combined.png"),
                fig, width = w_cm, height = h_cm, units = "cm", dpi = 600, bg = "white")

cli::cli_alert_success(
  "Saved sm_fcover_hmin_combined.pdf/png ({w_cm}×{h_cm} cm)"
)
