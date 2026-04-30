# ---
# title:  20c_fcover_hmin_joint_plot.R
# desc:   Supplementary figure (A) — joint fCover × h_min sensitivity,
#         crossed design. Shows how each comparison metric responds to the
#         combination of both thresholds, with site as facet.
#
#         Layout: 2×2 panels (one per metric: r, RMSE, Bias, Slope).
#         Each panel: x = fCover threshold (3 values), lines coloured by
#         h_min (4 values), facets = site (3 columns).
#         Reference cell (fCover = 0.90, h_min = 2 m) is marked with a
#         filled symbol on each line.
#
#         Reads:
#           revision/output/intermediate/reviewers/
#             joint_fcover_hmin_sensitivity_atbd.csv
#
#         Output:
#           revision/output/figures/sm_fcover_hmin_joint.pdf
#           revision/output/figures/sm_fcover_hmin_joint.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/20c_fcover_hmin_joint_plot.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

# ── I/O ────────────────────────────────────────────────────────────────────────

in_csv  <- here::here(
  "revision", "output", "intermediate", "reviewers",
  "joint_fcover_hmin_sensitivity_atbd.csv"
)
out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dt <- data.table::fread(in_csv)

# ── Prepare ────────────────────────────────────────────────────────────────────

site_order  <- c("Aigoual", "Blois", "Mormal")
# 4 h_min values: colourblind-safe sequential palette (light to dark)
hmin_colors <- c("2" = "#d9f0a3", "3" = "#78c679", "4" = "#238443",
                 "5" = "#004529")
hmin_labels <- c("2" = "h_min = 2 m", "3" = "h_min = 3 m",
                 "4" = "h_min = 4 m", "5" = "h_min = 5 m")

dt[, site        := factor(site, levels = site_order)]
dt[, fcover_pct  := fcover_threshold * 100]
dt[, hmin_fac    := factor(h_min_value)]
# Flag the reference cell for filled symbol
dt[, is_ref      := fcover_threshold == 0.90 & h_min_value == 2]

metric_list <- list(
  list(col = "R",     label = expression(italic(r))),
  list(col = "RMSE",  label = expression(RMSE~(m^2~m^{-2}))),
  list(col = "Bias",  label = expression(Bias~(m^2~m^{-2}))),
  list(col = "Slope", label = expression(Slope))
)

theme_joint <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor    = ggplot2::element_blank(),
    strip.background    = ggplot2::element_rect(fill = "grey92"),
    strip.text          = ggplot2::element_text(size = 8, face = "bold"),
    legend.position     = "none",
    axis.title.y        = ggplot2::element_text(size = 8),
    axis.title.x        = ggplot2::element_text(size = 8)
  )

# ── Panel builder ──────────────────────────────────────────────────────────────

make_joint_panel <- function(col, ylab, show_legend = FALSE) {
  p <- ggplot2::ggplot(
    dt,
    ggplot2::aes(x     = fcover_pct,
                 y     = .data[[col]],
                 colour = hmin_fac,
                 group  = hmin_fac)
  ) +
    ggplot2::facet_wrap(~ site, nrow = 1, scales = "fixed") +
    ggplot2::geom_vline(xintercept = 90, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.75) +
    # Open circles for all points
    ggplot2::geom_point(size = 2.2, shape = 21, fill = "white", stroke = 1.0) +
    # Filled circle at the reference cell (fCover = 0.90, h_min = 2 m)
    ggplot2::geom_point(
      data   = dt[is_ref == TRUE & get(col) == get(col)],
      ggplot2::aes(x = fcover_pct, y = .data[[col]], colour = hmin_fac),
      size   = 3.0, shape = 16
    ) +
    ggplot2::scale_x_continuous(breaks = c(80, 90, 95),
                                labels = c("80%", "90%", "95%")) +
    ggplot2::scale_colour_manual(
      values = hmin_colors,
      labels = hmin_labels,
      name   = expression(h[min]~(m))
    ) +
    ggplot2::labs(x = "fCover threshold", y = ylab) +
    theme_joint +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")

  if (col == "Bias")  p <- p + ggplot2::geom_hline(yintercept = 0,
                                colour = "grey40", linewidth = 0.3)
  if (col == "Slope") p <- p + ggplot2::geom_hline(yintercept = 1,
                                colour = "grey40", linewidth = 0.3)
  p
}

# ── Assemble 2×2 plate ────────────────────────────────────────────────────────

panels <- lapply(seq_along(metric_list), function(i) {
  m <- metric_list[[i]]
  make_joint_panel(m$col, m$label, show_legend = (i == length(metric_list)))
})

fig <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Joint sensitivity of LAI_ALS vs LAI_S2 metrics to fCover × h_min",
    subtitle = paste0(
      "x-axis: fCover threshold; colours: h_min. ",
      "Facets: site. Filled symbol = reference (fCover = 90%, h_min = 2 m). ",
      "ATBD parametrisation."
    ),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40")
    )
  ) &
  ggplot2::theme(legend.position = "right")

# ── Save ───────────────────────────────────────────────────────────────────────

w_cm <- 28
h_cm <- 18

ggplot2::ggsave(file.path(out_dir, "sm_fcover_hmin_joint.pdf"),
                fig, width = w_cm, height = h_cm, units = "cm")
ggplot2::ggsave(file.path(out_dir, "sm_fcover_hmin_joint.png"),
                fig, width = w_cm, height = h_cm, units = "cm", dpi = 300)

cli::cli_alert_success(
  "Saved sm_fcover_hmin_joint.pdf/png ({w_cm}×{h_cm} cm)"
)
