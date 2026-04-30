# ---
# title:  17b_k_sensitivity_plot.R
# desc:   Figure — sensitivity of LAI_ALS vs LAI_S2_ATBD metrics to the
#         extinction coefficient k (0.3 to 0.8).
#
#         Layout: 3 panels (RMSE | Bias | Slope), x = k, one coloured line
#         per site at theta = 0 degrees. Shaded ribbon = range across all
#         theta values (0 to 30 degrees) to show the secondary theta effect.
#         Reference: vertical line at k = 0.5 (Bouvier et al. 2015 default).
#         R is not shown: it is invariant to scale transformations of LAI_ALS.
#
#         Reads:
#           revision/output/intermediate/reviewers/
#             k_theta_sensitivity_atbd.csv
#
#         Outputs:
#           revision/output/figures/reviewers/sm_k_sensitivity.pdf / .png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/17b_k_sensitivity_plot.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(scales)
library(cli)

# ── I/O ────────────────────────────────────────────────────────────────────────

in_csv  <- here::here(
  "revision", "output", "intermediate", "reviewers",
  "k_theta_sensitivity_atbd.csv"
)
out_dir <- file.path(paths$output, "figures", "reviewers")

if (!file.exists(in_csv))
  stop("Missing: ", in_csv, "\nRun revision/scripts/17_k_sensitivity.R first.")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dt <- data.table::fread(in_csv)
cli::cli_alert_success("Loaded {nrow(dt)} rows")

# ── Prepare ────────────────────────────────────────────────────────────────────

site_levels   <- c("Aigoual", "Blois", "Mormal")
site_colours  <- c(Aigoual = "#1b7837", Blois = "#762a83", Mormal = "#d6604d")
k_ref         <- 0.6

metrics_sel <- c("R", "RMSE", "Bias", "Slope")

ref_hlines <- data.table::data.table(
  metric     = factor(c("R", "Bias", "Slope"), levels = metrics_sel),
  yintercept = c(1, 0, 1)
)

dt[, site := factor(site, levels = site_levels)]

# Theta = 0 line + ribbon (min–max across all theta) per (site, k)
dt_theta0 <- dt[theta_deg == 0L, .(site, k_value, R, RMSE, Bias, Slope)]

# Flatten range columns
dt_ribbon <- dt[, .(
  R_lo     = min(R,     na.rm = TRUE),
  R_hi     = max(R,     na.rm = TRUE),
  RMSE_lo  = min(RMSE,  na.rm = TRUE),
  RMSE_hi  = max(RMSE,  na.rm = TRUE),
  Bias_lo  = min(Bias,  na.rm = TRUE),
  Bias_hi  = max(Bias,  na.rm = TRUE),
  Slope_lo = min(Slope, na.rm = TRUE),
  Slope_hi = max(Slope, na.rm = TRUE)
), by = .(site, k_value)]
dt_ribbon[, site := factor(site, levels = site_levels)]

# Long format for line data
long_line <- data.table::melt(
  dt_theta0,
  id.vars       = c("site", "k_value"),
  measure.vars  = metrics_sel,
  variable.name = "metric",
  value.name    = "value"
)
long_line[, metric := factor(metric, levels = metrics_sel)]

# Long format for ribbon (need lo and hi separately, then merge)
long_ribbon <- data.table::melt(
  dt_ribbon,
  id.vars      = c("site", "k_value"),
  measure.vars = patterns("_lo$", "_hi$"),
  value.name   = c("lo", "hi"),
  variable.name = "metric"
)
long_ribbon[, metric := factor(
  gsub("_lo$", "", names(dt_ribbon)[
    grep("_lo$", names(dt_ribbon))
  ])[metric],
  levels = metrics_sel
)]

# ── Build plot ─────────────────────────────────────────────────────────────────

theme_k <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor  = ggplot2::element_blank(),
    strip.background  = ggplot2::element_rect(fill = "grey92"),
    strip.text        = ggplot2::element_text(face = "bold", size = 10),
    legend.position   = "bottom",
    legend.title      = ggplot2::element_blank(),
    legend.key.width  = ggplot2::unit(1.4, "cm")
  )

p <- ggplot2::ggplot() +
  # Theta range ribbon
  ggplot2::geom_ribbon(
    data = long_ribbon,
    ggplot2::aes(x = k_value, ymin = lo, ymax = hi,
                 fill = site, group = site),
    alpha = 0.15, show.legend = FALSE
  ) +
  # Reference lines
  ggplot2::geom_hline(
    data        = ref_hlines,
    ggplot2::aes(yintercept = yintercept),
    linetype = "dashed", colour = "grey50", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  # k = 0.5 reference
  ggplot2::geom_vline(
    xintercept = k_ref,
    linetype   = "dotted", colour = "grey40", linewidth = 0.6
  ) +
  # Theta = 0 lines
  ggplot2::geom_line(
    data = long_line,
    ggplot2::aes(x = k_value, y = value, colour = site, group = site),
    linewidth = 1.1
  ) +
  ggplot2::geom_point(
    data = long_line,
    ggplot2::aes(x = k_value, y = value, colour = site),
    size = 2.5, shape = 16
  ) +
  ggplot2::facet_wrap(
    ~ metric,
    nrow   = 1,
    scales = "free_y",
    labeller = ggplot2::as_labeller(
      c(R     = "italic(r)",
        RMSE  = "RMSE~(m^2~m^{-2})",
        Bias  = "Bias~(m^2~m^{-2})",
        Slope = "Slope"),
      default = ggplot2::label_parsed
    )
  ) +
  ggplot2::scale_colour_manual(values = site_colours) +
  ggplot2::scale_fill_manual(values   = site_colours) +
  ggplot2::scale_x_continuous(
    name   = expression(italic(k) ~ (extinction ~ coefficient)),
    breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  ) +
  ggplot2::labs(
    y       = NULL,
    caption = paste0(
      "Solid line: theta = 0 deg (nadir). ",
      "Shaded band: range across theta in {0, 5, ..., 30} deg. ",
      "Dotted vertical: k = 0.6 (value used in this study). ",
      "Dashed horizontal: Bias = 0 / Slope = 1."
    )
  ) +
  theme_k +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(size = 8, colour = "grey40", hjust = 0)
  )

# ── Save ───────────────────────────────────────────────────────────────────────

w_cm <- 30
h_cm <- 10

out_pdf <- file.path(out_dir, "sm_k_sensitivity.pdf")
out_png <- file.path(out_dir, "sm_k_sensitivity.png")

ggplot2::ggsave(out_pdf, p, width = w_cm, height = h_cm, units = "cm")
ggplot2::ggsave(out_png, p, width = w_cm, height = h_cm, units = "cm",
                device = "png", dpi = 300)

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
