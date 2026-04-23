# ---
# title:  05b_sm5_plot_dopt_metrics.R
# desc:   Visualisation — one figure per agreement metric (R, RMSE, |Bias|,
#         |Slope-1|, Pareto distance to utopia), each showing all
#         Site × h_min combinations as facets (facet_grid h_min ~ Site).
#
#         Within each panel:
#           - Lines coloured by Norm (DSM_keepTrees blue solid,
#             DTM_keepTrees red dashed) — matches Main_01 palette.
#           - Gray shaded rectangle for depths > h_min (redundancy zone):
#             canopy heights > h_min mean integrating beyond h_min adds no
#             new information, so d_opt selection ignores these depths.
#           - Vertical lines mark d_opt per selection method and per Norm
#             (line type matches method, colour matches norm).
#
#         Reads:
#           revision/output/intermediate/sm5/all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#           revision/output/intermediate/sm5/dopt_reference.csv  (optional)
#
#         Figures written to revision/output/figures/sm5/.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05b_sm5_plot_dopt_metrics.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal")
norms_select  <- c("DSM_keepTrees", "DTM_keepTrees")
h_min_values  <- c(10L, 15L, 20L)

norm_colours  <- c(DSM_keepTrees = "#2166ac", DTM_keepTrees = "#d6604d")
norm_labels   <- c(DSM_keepTrees = "DSM", DTM_keepTrees = "DTM")
norm_linetypes <- c(DSM_keepTrees = "solid",  DTM_keepTrees = "dashed")

# ── Load data ─────────────────────────────────────────────────────────────────

atbd_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
)
if (!file.exists(atbd_csv))
  stop("Run 05a_sm5_compute_metrics_atbd.R first:\n  ", atbd_csv)

dt <- data.table::fread(atbd_csv)
dt <- dt[ATBD == TRUE & Norm %in% norms_select & h_min %in% h_min_values]
dt[, h_min_label := paste0("h_min = ", h_min, " m")]
dt[, h_min_label := factor(h_min_label,
                            levels = paste0("h_min = ", h_min_values, " m"))]

# ── Compute Pareto distance to utopia ─────────────────────────────────────────

compute_dist_utopia <- function(sub) {
  m1 <- 1 - sub$R2
  m2 <- sub$RMSE
  m3 <- abs(sub$Bias)
  m4 <- abs(sub$Slope - 1)

  norm01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / diff(r)
  }

  d <- sqrt(norm01(m1)^2 + norm01(m2)^2 + norm01(m3)^2 + norm01(m4)^2)
  data.table::data.table(Depth = sub$Depth, dist_utopia = d)
}

dist_dt <- dt[, compute_dist_utopia(.SD), by = .(Site, Norm, h_min)]
dt      <- merge(dt, dist_dt, by = c("Site", "Norm", "Depth", "h_min"))

# ── Plot helper ───────────────────────────────────────────────────────────────

make_metric_figure <- function(plot_dt, yvar, ylab,
                                ylim = NULL, show_legend = FALSE) {

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = Depth, y = .data[[yvar]],
                 colour   = Norm,
                 linetype = Norm)
  ) +
    # Redundancy zone: gray rect for depths > h_min (per facet row)
    # Each h_min_label facet has a different h_min value — use geom_rect
    # with a data.frame that maps h_min_label to its numeric value.
    ggplot2::geom_rect(
      data = data.table::data.table(
        h_min_label = factor(paste0("h_min = ", h_min_values, " m"),
                             levels = paste0("h_min = ", h_min_values, " m")),
        h_min_val   = h_min_values
      ),
      ggplot2::aes(xmin = h_min_val, xmax = Inf, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill  = "grey85",
      alpha = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = norm_colours, labels = norm_labels) +
    ggplot2::scale_linetype_manual(values = norm_linetypes, labels = norm_labels) +
    ggplot2::scale_x_continuous(breaks = c(1, seq(5, 38, by = 5))) +
    ggplot2::labs(x = "Depth from canopy top (m)", y = ylab,
                  colour = NULL, linetype = NULL) +
    ggplot2::facet_grid(h_min_label ~ Site) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position  = if (show_legend) "bottom" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92")
    )

  if (!is.null(ylim))
    p <- p + ggplot2::coord_cartesian(ylim = ylim)

  p
}

# ── Output directory ──────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Build derived columns ─────────────────────────────────────────────────────

dt[, abs_Bias  := abs(Bias)]
dt[, abs_Slope1 := abs(Slope - 1)]

# ── Individual metric figures ─────────────────────────────────────────────────

metrics_spec <- list(
  list(var = "R",           ylab = "Pearson R",           file = "R"),
  list(var = "RMSE",        ylab = "RMSE (m²/m²)",        file = "RMSE"),
  list(var = "abs_Bias",    ylab = "|Bias| (m²/m²)",       file = "absBias"),
  list(var = "abs_Slope1",  ylab = "|Slope - 1|",          file = "absSlope1"),
  list(var = "dist_utopia", ylab = "Pareto dist. to utopia", file = "Pareto")
)

for (spec in metrics_spec) {

  cat("Plotting metric:", spec$var, "\n")

  fig <- make_metric_figure(
    plot_dt     = dt,
    yvar        = spec$var,
    ylab        = spec$ylab,
    show_legend = (spec$var == "dist_utopia")
  ) +
    ggplot2::labs(
      title    = paste0("d_opt metric — ", spec$ylab, " (ATBD)"),
      subtitle = "Gray zone: depth > h_min (redundancy zone, excluded from d_opt selection)"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7,  colour = "grey40"),
      legend.position = "bottom"
    )

  out_file <- file.path(out_dir, paste0("dopt_metric_", spec$file, "_ATBD.pdf"))
  ggplot2::ggsave(out_file, fig, width = 24, height = 18, units = "cm")
  cat("  Written:", out_file, "\n")
}

# ── Combined 5-panel figure ───────────────────────────────────────────────────

cat("Plotting combined 5-panel figure...\n")

p_R      <- make_metric_figure(dt, "R",           "Pearson R")
p_RMSE   <- make_metric_figure(dt, "RMSE",        "RMSE (m²/m²)")
p_Bias   <- make_metric_figure(dt, "abs_Bias",    "|Bias| (m²/m²)")
p_Slope  <- make_metric_figure(dt, "abs_Slope1",  "|Slope - 1|")
p_Pareto <- make_metric_figure(dt, "dist_utopia", "Pareto dist.",
                                show_legend = TRUE)

fig_combined <-
  (p_R | p_RMSE | p_Bias) / (p_Slope | p_Pareto | patchwork::plot_spacer()) +
  patchwork::plot_annotation(
    title    = "d_opt — five agreement metrics vs canopy integration depth (ATBD)",
    subtitle = "Gray: redundancy zone (depth > h_min)",
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8,  colour = "grey40")
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

out_combined <- file.path(out_dir, "dopt_metrics_combined_ATBD.pdf")
ggplot2::ggsave(out_combined, fig_combined, width = 36, height = 30, units = "cm")
cat("  Written:", out_combined, "\n")

cat("Done.\n")
