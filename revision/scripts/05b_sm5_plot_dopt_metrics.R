# ---
# title:  05b_sm5_plot_dopt_metrics.R
# desc:   Visualisation — R (Pearson r) vs canopy integration depth, faceted
#         by h_min × Site, matching the style of Fig. 5 from the submitted
#         paper (Main_01_analysis_optim_depth.R).
#
#         Per panel:
#           - Curves coloured by Norm (DSM_keepTrees blue, DTM_keepTrees red)
#           - Gray shaded rectangle for depths > h_min (redundancy zone)
#           - Downward arrows at Pareto d_opt per norm, coloured by norm
#           - Dashed black horizontal: r_all = R at Depth=38 (total LAI ref)
#           - Dotted black horizontal: r_uniform = R at Depth=h_min (sampled
#             LAI ref — boundary of redundancy zone)
#
#         Reads:
#           revision/output/intermediate/sm5/all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#
#         Figures written to revision/output/figures/sm5/.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05b_sm5_plot_dopt_metrics.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("cli")

source(here::here("revision", "R", "sm5_dopt.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites         <- c("Aigoual", "Blois", "Mormal", "Sites combined")
norms_select  <- c("DSM_keepTrees", "DTM_keepTrees")
h_min_values  <- c(10L, 15L, 20L)

norm_colours   <- c(DSM_keepTrees = "#2166ac", DTM_keepTrees = "#d6604d")
norm_labels    <- c(DSM_keepTrees = "DSM", DTM_keepTrees = "DTM")
norm_linetypes <- c(DSM_keepTrees = "solid",  DTM_keepTrees = "dashed")

hmin_levels <- paste0("h_min = ", h_min_values, " m")

# ── Load and prepare data ──────────────────────────────────────────────────────

atbd_csv <- here::here(
  "revision", "output", "intermediate", "sm5",
  "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
)
if (!file.exists(atbd_csv))
  stop("Run 05a_sm5_compute_metrics_atbd.R first:\n  ", atbd_csv)

dt <- data.table::fread(atbd_csv)
dt <- dt[ATBD == TRUE & Norm %in% norms_select & h_min %in% h_min_values]

# Rename "All_sites" → "Sites combined" to match figure caption
dt[Site == "All_sites", Site := "Sites combined"]

dt[, Site := factor(Site, levels = sites)]
dt[, h_min_label := factor(paste0("h_min = ", h_min, " m"), levels = hmin_levels)]

# ── Derived columns used in multi-metric figures ───────────────────────────────

dt[, abs_Bias   := abs(Bias)]
dt[, abs_Slope1 := abs(Slope - 1)]

# Pareto distance to utopia (one value per Site × Norm × h_min × Depth)
compute_dist_utopia <- function(sub) {
  m1 <- 1 - sub$R2; m2 <- sub$RMSE
  m3 <- abs(sub$Bias); m4 <- abs(sub$Slope - 1)
  norm01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / diff(r)
  }
  data.table::data.table(
    Depth = sub$Depth,
    dist_utopia = sqrt(norm01(m1)^2 + norm01(m2)^2 + norm01(m3)^2 + norm01(m4)^2)
  )
}
dist_dt <- dt[, compute_dist_utopia(.SD), by = .(Site, Norm, h_min)]
dt      <- merge(dt, dist_dt, by = c("Site", "Norm", "Depth", "h_min"))

# ── Pareto d_opt per (Site × Norm × h_min) ────────────────────────────────────
# Computed here so arrows can be drawn at the chosen d_opt.

cli::cli_alert_info("Computing Pareto d_opt for each h_min…")
dopt_by_hmin <- data.table::rbindlist(lapply(h_min_values, function(hm) {
  sub <- dt[h_min == hm]
  d   <- select_dopt(sub, methods = "pareto", max_depth = hm,
                     prosail_filter = "ATBD")
  d[method_dopt == "pareto"][, h_min := hm]
}))

dopt_by_hmin[, h_min_label := factor(
  paste0("h_min = ", h_min, " m"), levels = hmin_levels
)]

# Arrow data: join d_opt with metric values at that depth
arrow_dt <- merge(
  dopt_by_hmin[, .(Site, Norm, h_min, h_min_label, d_opt, R, RMSE, Bias, Slope)],
  dt[, .(Site, Norm, Depth, h_min, dist_utopia)],
  by.x = c("Site", "Norm", "d_opt", "h_min"),
  by.y = c("Site", "Norm", "Depth", "h_min"),
  all.x = TRUE
)
arrow_dt[, abs_Bias   := abs(Bias)]
arrow_dt[, abs_Slope1 := abs(Slope - 1)]
arrow_dt[, Site := factor(Site, levels = sites)]

# ── Reference horizontal lines for the R figure ───────────────────────────────
# r_all     (dashed): mean R at Depth=38 across norms — "total LAI" reference
# r_uniform (dotted): mean R at Depth=h_min across norms — "sampled LAI" ref,
#            i.e. the last valid depth before the redundancy zone.

rref_all <- dt[Depth == 38L & Norm %in% norms_select,
               .(r_all = mean(R, na.rm = TRUE)),
               by = .(Site, h_min, h_min_label)]

rref_uniform <- data.table::rbindlist(lapply(h_min_values, function(hm) {
  sub <- dt[Depth == hm & h_min == hm & Norm %in% norms_select,
            .(r_uniform = mean(R, na.rm = TRUE)),
            by = .(Site, h_min)]
  sub[, h_min_label := factor(paste0("h_min = ", hm, " m"), levels = hmin_levels)]
  sub
}))

hline_dt <- merge(rref_all, rref_uniform, by = c("Site", "h_min", "h_min_label"),
                  all = TRUE)
hline_dt[, Site := factor(Site, levels = sites)]

# ── Plot helper ───────────────────────────────────────────────────────────────

make_metric_figure <- function(plot_dt, yvar, ylab,
                                ylim = NULL, show_legend = FALSE,
                                arrow_dt = NULL) {

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = Depth, y = .data[[yvar]],
                 colour = Norm, linetype = Norm)
  ) +
    ggplot2::geom_rect(
      data = data.table::data.table(
        h_min_label = factor(hmin_levels, levels = hmin_levels),
        h_min_val   = h_min_values
      ),
      ggplot2::aes(xmin = h_min_val, xmax = Inf, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey85", alpha = 0.5
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

  # Downward arrows at d_opt per norm — same colour as the respective curve
  if (!is.null(arrow_dt) && nrow(arrow_dt) > 0 && yvar %in% names(arrow_dt)) {
    y_range <- range(plot_dt[[yvar]], na.rm = TRUE)
    offset  <- diff(y_range) * 0.08

    for (norm_i in norms_select) {
      sub_a <- data.table::copy(arrow_dt[Norm == norm_i])
      if (nrow(sub_a) == 0L) next
      sub_a[, y_end   := get(yvar)]
      sub_a[, y_start := y_end + offset]
      p <- p + ggplot2::geom_segment(
        data = sub_a,
        ggplot2::aes(x = d_opt, xend = d_opt,
                     y = y_start, yend = y_end),
        colour    = norm_colours[[norm_i]],
        linewidth = 0.6,
        arrow     = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"),
                                   type = "closed"),
        inherit.aes = FALSE
      )
    }
  }

  p
}

# ── Output directory ──────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Figure 5 — R only, with reference lines and arrows ───────────────────────

cat("Plotting R figure (Fig. 5 style)…\n")

p_R_fig <- make_metric_figure(
  plot_dt     = dt,
  yvar        = "R",
  ylab        = "Pearson r",
  show_legend = TRUE,
  arrow_dt    = arrow_dt
) +
  ggplot2::geom_hline(
    data        = hline_dt,
    ggplot2::aes(yintercept = r_all),
    linetype    = "dashed",
    colour      = "grey20",
    linewidth   = 0.5
  ) +
  ggplot2::geom_hline(
    data        = hline_dt,
    ggplot2::aes(yintercept = r_uniform),
    linetype    = "dotted",
    colour      = "grey20",
    linewidth   = 0.5
  ) +
  ggplot2::labs(
    title    = "Pearson r — LAI_S2 (ATBD) vs LAI_ALS, by depth (ATBD)",
    subtitle = paste(
      "Gray: redundancy zone | Arrows: Pareto d_opt per norm",
      "| dashed: r_all (full canopy, depth=38)",
      "| dotted: r_uniform (LAI at h_min boundary)"
    )
  ) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(size = 11, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7,  colour = "grey40"),
    legend.position = "bottom"
  )

out_r <- file.path(out_dir, "dopt_R_ATBD.pdf")
ggplot2::ggsave(out_r, p_R_fig, width = 28, height = 18, units = "cm")
cat("  Written:", out_r, "\n")

# ── Multi-metric exploratory figures ─────────────────────────────────────────

metrics_spec <- list(
  list(var = "R",           ylab = "Pearson R",             file = "R"),
  list(var = "RMSE",        ylab = "RMSE (m²/m²)",          file = "RMSE"),
  list(var = "abs_Bias",    ylab = "|Bias| (m²/m²)",        file = "absBias"),
  list(var = "abs_Slope1",  ylab = "|Slope - 1|",           file = "absSlope1"),
  list(var = "dist_utopia", ylab = "Pareto dist. to utopia", file = "Pareto")
)

for (spec in metrics_spec) {
  cat("Plotting metric:", spec$var, "\n")
  fig <- make_metric_figure(
    plot_dt     = dt,
    yvar        = spec$var,
    ylab        = spec$ylab,
    show_legend = (spec$var == "dist_utopia"),
    arrow_dt    = arrow_dt
  ) +
    ggplot2::labs(
      title    = paste0("d_opt metric — ", spec$ylab, " (ATBD)"),
      subtitle = "Gray zone: depth > h_min (redundancy zone) | Arrows: Pareto d_opt"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7,  colour = "grey40"),
      legend.position = "bottom"
    )
  out_file <- file.path(out_dir, paste0("dopt_metric_", spec$file, "_ATBD.pdf"))
  ggplot2::ggsave(out_file, fig, width = 28, height = 18, units = "cm")
  cat("  Written:", out_file, "\n")
}

cat("Done.\n")
