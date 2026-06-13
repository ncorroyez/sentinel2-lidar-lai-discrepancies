# ---
# title:  sm6_plot_heterogeneity.R
# desc:   Shared constants and the plot_het_grid() function used by both
#         scripts/steps/20_sm6_plot_heterogeneity.R (SM6a pixels) and
#         scripts/steps/20e_plot_sm5_heter_classes.R (SM5 sample).
# ---

# ── Factor levels ──────────────────────────────────────────────────────────────

het_levels   <- c("Low", "Medium", "High", "Total")
combo_levels <- c("ATBD_vs_ALS", "ATBD_vs_ALS_dopt", "opt_vs_ALS_dopt")
site_levels  <- c("Aigoual", "Blois", "Mormal")

metric_parsed_levels <- c("italic(r)", "RMSE~(m^2/m^2)", "Bias~(m^2/m^2)",
                           "Slope")

# ── Colour palette ─────────────────────────────────────────────────────────────

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

# ── Reference horizontal lines ─────────────────────────────────────────────────

het_ref_lines <- data.table::data.table(
  metric     = factor(
    c("Bias~(m^2/m^2)", "Slope"),
    levels = metric_parsed_levels
  ),
  yintercept = c(0, 1)
)

# ── load_heter_csv_long ───────────────────────────────────────────────────────

#' @title Load a heterogeneity CSV and pivot to long format for plotting
#'
#' @param csv_path Character. Path to a CSV with the schema produced by
#'   scripts/steps/15_sm6_analysis.R or 20e_sm5_heter_classes.R.
#' @return A long-format \code{data.table} with the factor levels above.
#' @export
load_heter_csv_long <- function(csv_path) {
  if (!file.exists(csv_path))
    stop("Missing: ", csv_path)
  dt <- data.table::fread(csv_path)
  dt[, het_class   := factor(het_class,   levels = het_levels)]
  dt[, combination := factor(combination, levels = combo_levels)]
  dt[, site        := factor(site,        levels = site_levels)]
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
  long_dt
}

# ── plot_het_grid ─────────────────────────────────────────────────────────────

#' @title 4-row × 3-column metric-vs-class facet plot
#'
#' @description Rows = R/RMSE/Bias/Slope, columns = 3 sites, x = heterogeneity
#' class, three coloured lines = LAI comparison combinations. Stars at x =
#' "Total" mark the pooled-sample metric. Reference dashed lines at Bias = 0
#' and Slope = 1.
#'
#' @param data_ms  data.table from \code{load_heter_csv_long()}, already
#'   filtered to one \code{metric_source} (e.g. "DSM" or "CHM").
#' @param ms_short Character. Short tag used in the x-axis label
#'   (\code{ms_short[SD]}, plotmath).
#' @param base_size Numeric. ggplot base font size.
#' @return A \code{ggplot} object.
#' @export
plot_het_grid <- function(data_ms, ms_short, base_size = 14) {
  rmse_zero <- data.frame(
    site   = site_levels,
    metric = factor("RMSE~(m^2/m^2)", levels = metric_parsed_levels),
    value  = 0
  )

  normal_dt <- data_ms[het_class != "Total"]
  total_dt  <- data_ms[het_class == "Total"]

  x_offsets <- setNames(c(-0.2, 0.0, 0.2), combo_levels)
  total_dt[, x_num := 4 + x_offsets[as.character(combination)]]

  ggplot2::ggplot(
    data_ms,
    ggplot2::aes(x = het_class, y = value,
                 colour = combination, group = combination)
  ) +
    ggplot2::geom_blank(
      data        = rmse_zero,
      ggplot2::aes(y = value),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      data        = het_ref_lines,
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
      shape       = 8, size = 3.5, stroke = 1.5,
      inherit.aes = FALSE, na.rm = TRUE
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
      x = bquote("Horizontal Heterogeneity (" * .(ms_short)[SD] * ")"),
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
