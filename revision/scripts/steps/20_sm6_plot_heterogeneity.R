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
#           revision/output/intermediate/sm6/heterogeneity_analysis.csv
#
#         Figures written to revision/output/figures/sm6/:
#           het_metrics_grid_DSM.pdf
#           het_metrics_grid_CHM.pdf
#           het_metrics_grid_DSM.png
#           het_metrics_grid_CHM.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/20_sm6_plot_heterogeneity.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("cli")

# ── Load data ──────────────────────────────────────────────────────────────────

het_csv <- here::here(
  "revision", "output", "intermediate", "sm6", "heterogeneity_analysis.csv"
)
if (!file.exists(het_csv))
  stop("Missing: ", het_csv, "\nRun scripts 14 and 15 first.")

dt <- data.table::fread(het_csv)
cli::cli_alert_success(
  "Loaded {nrow(dt)} rows from heterogeneity_analysis.csv"
)

# ── Factor levels ──────────────────────────────────────────────────────────────

het_levels   <- c("Low", "Medium", "High", "Total")
combo_levels <- c("ATBD_vs_ALS", "ATBD_vs_ALS_dopt", "opt_vs_ALS_dopt")
site_levels  <- c("Aigoual", "Blois", "Mormal")

# Metric labels used as facet row labels (parsed as plotmath expressions)
metric_parsed_levels <- c("italic(r)", "RMSE~(m^2/m^2)", "Bias~(m^2/m^2)",
                           "Slope")

dt[, het_class   := factor(het_class,   levels = het_levels)]
dt[, combination := factor(combination, levels = combo_levels)]
dt[, site        := factor(site,        levels = site_levels)]

# ── Colour palette (green gradient, matching three_factors_analysis_final.R) ──

combo_colours <- c(
  ATBD_vs_ALS      = "#004d00",
  ATBD_vs_ALS_dopt = "#008000",
  opt_vs_ALS_dopt  = "#4CBB17"
)

combo_labels <- c(
  ATBD_vs_ALS      = expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS]),
  ATBD_vs_ALS_dopt = expression(LAI[S2_ATBD] ~ "vs" ~ LAI["ALS_dopt"]),
  opt_vs_ALS_dopt  = expression(LAI[S2_opt]  ~ "vs" ~ LAI["ALS_dopt"])
)

# ── Pivot to long format ───────────────────────────────────────────────────────

long_dt <- data.table::melt(
  dt,
  id.vars       = c("site", "metric_source", "combination", "het_class", "n"),
  measure.vars  = c("R", "RMSE", "Bias", "Slope"),
  variable.name = "metric",
  value.name    = "value"
)
long_dt[, metric := factor(
  as.character(metric),
  levels = c("R", "RMSE", "Bias", "Slope"),
  labels = metric_parsed_levels
)]

# ── Reference horizontal lines ─────────────────────────────────────────────────

ref_lines <- data.table::data.table(
  metric     = factor(
    c("Bias~(m^2/m^2)", "Slope"),
    levels = metric_parsed_levels
  ),
  yintercept = c(0, 1)
)

# ── Plot function ──────────────────────────────────────────────────────────────

plot_het_grid <- function(data_ms, ms_label, base_size = 14) {
  rmse_zero <- data.frame(
    site   = site_levels,
    metric = factor("RMSE~(m^2/m^2)", levels = metric_parsed_levels),
    value  = 0
  )

  # Split normal classes vs Total
  normal_dt <- data_ms[het_class != "Total"]
  total_dt  <- data_ms[het_class == "Total"]

  # Jitter Total stars per combination around x = 4 (position of "Total")
  x_offsets <- setNames(c(-0.2, 0.0, 0.2), combo_levels)
  total_dt[, x_num := 4 + x_offsets[as.character(combination)]]

  ggplot2::ggplot(
    data_ms,
    ggplot2::aes(
      x      = het_class,
      y      = value,
      colour = combination,
      group  = combination
    )
  ) +
    ggplot2::geom_blank(
      data        = rmse_zero,
      ggplot2::aes(y = value),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      data        = ref_lines,
      ggplot2::aes(yintercept = yintercept),
      linetype    = "dashed",
      colour      = "grey40",
      linewidth   = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(data = normal_dt, linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_point(data = normal_dt, size = 3.0, shape = 16, na.rm = TRUE) +
    ggplot2::geom_point(
      data        = total_dt,
      ggplot2::aes(x = x_num, y = value, colour = combination),
      shape       = 8,
      size        = 3.5,
      stroke      = 1.5,
      inherit.aes = FALSE,
      na.rm       = TRUE
    ) +
    ggplot2::facet_grid(
      metric ~ site,
      scales   = "free_y",
      labeller = ggplot2::labeller(
        metric = ggplot2::label_parsed,
        site   = ggplot2::label_value
      )
    ) +
    ggplot2::scale_x_discrete(limits = het_levels) +
    ggplot2::scale_colour_manual(
      values = combo_colours,
      labels = combo_labels,
      name   = "Relationship",
      guide  = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 3.5, linewidth = 1.1)
      )
    ) +
    ggplot2::labs(
      x = bquote("Horizontal Heterogeneity (" * .(ms_label)[SD] * ")"),
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ── Output directory ───────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures", "sm6")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Produce one figure per metric_source ──────────────────────────────────────

for (ms in c("DSM", "CHM")) {
  sub_ms <- long_dt[metric_source == ms]
  if (nrow(sub_ms) == 0L) {
    cli::cli_warn("No data for metric_source = {ms} — skipping.")
    next
  }
  p_ms <- plot_het_grid(sub_ms, ms_label = ms)

  out_pdf <- file.path(out_dir, paste0("het_metrics_grid_", ms, ".pdf"))
  out_png <- file.path(out_dir, paste0("het_metrics_grid_", ms, ".png"))

  ggplot2::ggsave(out_pdf, p_ms, width = 26, height = 22, units = "cm")
  ggplot2::ggsave(out_png, p_ms, width = 26, height = 22, units = "cm",
                  device = "png", dpi = 300)

  cli::cli_alert_success("Written: {out_pdf}")
  cli::cli_alert_success("Written: {out_png}")
}

cli::cli_h1("Done — SM6 heterogeneity figures")

# ── Console summary: raw metrics per metric_source x combination x site x level ─

fmt_pval <- function(p) {
  if (is.na(p) || length(p) == 0L) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

for (ms in unique(dt$metric_source)) {
  cli::cli_h1("Raw metrics -- metric_source: {ms}")

  for (s in site_levels) {
    cli::cli_h2("{s}")
    dt_s <- dt[metric_source == ms & site == s]

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
