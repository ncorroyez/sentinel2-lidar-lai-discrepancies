# ---
# title:  07_sm5_scatter_lai_atbd.R
# desc:   Scatter plot 1 row x 3 columns (one per site): LAI_ALS (y) vs
#         LAI_S2_ATBD (x). Both norms overlaid with colour coding (CHM blue,
#         DTM red). Uses Pareto d_opt from dopt_reference.csv (h_min = 10 m).
#
#         Top-left annotation per norm (colour coded):
#           CHM: y = a*x + b  r = ...  R2 = ...  RMSE = ...  Bias = ...
#           DTM: y = a*x + b  ...
#           (equation from lm(LAI_ALS ~ LAI_S2); Bias = mean(S2 - LiDAR))
#
#         Bottom-right inset: overlapping density curves for LAI_S2 (black),
#         LAI_ALS CHM (blue), LAI_ALS DTM (red).
#
# Prerequisites:
#   06_sm5_select_dopt.R  -> revision/output/intermediate/sm5/dopt_reference.csv
#   04a                   -> LAI_estimated_atbd_*.csv
#   03b                   -> PAD_*_Depth_*_Samples_*.csv
#
# Run from project root (NC_Full/):
#   source("revision/scripts/07_sm5_scatter_lai_atbd.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("cli")

source(here::here("revision", "R", "sm5_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites           <- c("Aigoual", "Blois", "Mormal")
norms_select    <- c("DSM_keepTrees", "DTM_keepTrees")
sampling_method <- "stratified_uniform"
name_strategy   <- "LIDFa_lai_LMA_BROWN"
nb_samples      <- 5000L
h_min_pixel     <- 10L

norm_colours <- c(DSM_keepTrees = "#2166ac", DTM_keepTrees = "#d6604d")
norm_labels  <- c(DSM_keepTrees = "CHM",     DTM_keepTrees = "DTM")

# ── Load Pareto d_opt per (Site, Norm) ────────────────────────────────────────

dopt_csv <- here::here("revision", "output", "intermediate", "sm5",
                        "dopt_reference.csv")
if (!file.exists(dopt_csv))
  stop("Run 06_sm5_select_dopt.R first:\n  ", dopt_csv)

dopt_pareto <- data.table::fread(dopt_csv)[
  method_dopt == "pareto" & Site %in% sites,
  .(Site, Norm, d_opt)
]

# ── Load paired (LAI_ALS, LAI_S2) for each Site x Norm ───────────────────────

load_paired <- function(site, norm, d_opt_val) {
  s2_path <- here::here(
    "revision", "output", "intermediate", "PROSAIL_Models",
    site, name_strategy, "atbd",
    paste0("LAI_estimated_atbd_", sampling_method,
           "_hmin", h_min_pixel, "_nbSamples_", nb_samples, ".csv")
  )
  lidar_path <- here::here(
    "03_RESULTS", site, "PROSAIL_Optimization", "sampling",
    paste0("PAD_", norm, "_Depth_", d_opt_val,
           "_Samples_", sampling_method,
           "_hmin", h_min_pixel, "_nbSamples_", nb_samples, ".csv")
  )
  if (!file.exists(s2_path) || !file.exists(lidar_path)) {
    cli::cli_warn("Files missing for {site} / {norm}, skipping")
    return(NULL)
  }
  s2_dt    <- data.table::fread(s2_path,    header = TRUE, sep = "\t")
  lidar_dt <- data.table::fread(lidar_path, header = TRUE, sep = "\t")

  s2_ids  <- seq_len(nrow(s2_dt))
  als_v   <- lidar_dt[["lidar_values"]][
    match(s2_ids, lidar_dt[["samples_id"]])
  ]
  s2_v  <- s2_dt[["LAI_atbd"]]
  valid <- !is.na(als_v) & !is.na(s2_v)

  data.table::data.table(
    Site    = site,
    Norm    = norm,
    d_opt   = d_opt_val,
    LAI_S2  = s2_v[valid],
    LAI_ALS = als_v[valid]
  )
}

cli::cli_h1("Loading paired data...")
pairs <- vector("list", nrow(dopt_pareto))
for (i in seq_len(nrow(dopt_pareto))) {
  r <- dopt_pareto[i]
  cli::cli_alert_info("{r$Site} / {r$Norm} at d_opt = {r$d_opt} m")
  pairs[[i]] <- load_paired(r$Site, r$Norm, r$d_opt)
}
dt <- data.table::rbindlist(pairs, use.names = TRUE)

# ── Compute per-(Site x Norm) statistics ──────────────────────────────────────

stats_dt <- dt[, {
  m   <- compute_metrics_s2_lidar(LAI_S2, LAI_ALS)
  fit <- stats::lm(LAI_ALS ~ LAI_S2)      # line plotted on scatter (y = ALS)
  a   <- round(stats::coef(fit)[[2L]], 2)
  b   <- round(stats::coef(fit)[[1L]], 2)
  b_str <- if (b >= 0) paste0("+", b) else as.character(b)
  list(
    a    = a,
    b    = b,
    b_str = b_str,
    r    = round(m$R,    2),
    R2   = round(m$R2,   2),
    RMSE = round(m$RMSE, 2),
    Bias = round(m$Bias, 2)
  )
}, by = .(Site, Norm)]

# Build the multi-line annotation label per (Site, Norm)
stats_dt[, label := paste0(
  norm_labels[Norm], ": y = ", a, "x ", b_str, "\n",
  "r = ", r, "  R\u00b2 = ", R2, "\n",
  "RMSE = ", RMSE, "  Bias = ", Bias
)]

cli::cli_h2("Metrics summary")
print(stats_dt[, .(Site, Norm, a, b, r, R2, RMSE, Bias)])

# ── Global axis limits (square, shared across panels) ─────────────────────────

xy_lim <- c(floor(min(dt$LAI_S2, dt$LAI_ALS, na.rm = TRUE)),
             ceiling(max(dt$LAI_S2, dt$LAI_ALS, na.rm = TRUE)))

# ── Panel builder ─────────────────────────────────────────────────────────────

make_panel <- function(site_name, show_y = TRUE) {

  site_dt    <- dt[Site == site_name]
  site_stats <- stats_dt[Site == site_name]

  # Annotation y positions: CHM block first, DTM block below
  x_ann   <- xy_lim[1] + diff(xy_lim) * 0.03
  y_start <- xy_lim[2] - diff(xy_lim) * 0.02
  block_h <- diff(xy_lim) * 0.24          # vertical space per 3-line block

  # ── Main scatter ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    site_dt,
    ggplot2::aes(x = LAI_S2, y = LAI_ALS, colour = Norm)
  ) +
    ggplot2::geom_point(alpha = 0.20, size = 0.5, shape = 16) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                          se = FALSE, linewidth = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dashed", colour = "grey45",
                          linewidth = 0.5) +
    ggplot2::scale_colour_manual(
      values = norm_colours, labels = norm_labels, name = NULL
    ) +
    ggplot2::coord_equal(xlim = xy_lim, ylim = xy_lim) +
    ggplot2::scale_x_continuous(
      breaks = pretty(xy_lim, n = 5)
    ) +
    ggplot2::scale_y_continuous(
      breaks = pretty(xy_lim, n = 5)
    ) +
    ggplot2::labs(
      title = site_name,
      x     = expression(LAI[S2]~"(ATBD)"),
      y     = if (show_y) expression(LAI[ALS]) else NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title.y     = if (!show_y) ggplot2::element_blank() else
                           ggplot2::element_text()
    )

  # ── Coloured text annotation (top-left), one block per norm ─────────────────
  for (j in seq_along(norms_select)) {
    norm_j <- norms_select[[j]]
    st_j   <- site_stats[Norm == norm_j]
    if (nrow(st_j) == 0L) next
    y_ann <- y_start - block_h * (j - 1L)
    p <- p + ggplot2::annotate(
      "text",
      x       = x_ann,
      y       = y_ann,
      label   = st_j$label,
      colour  = norm_colours[[norm_j]],
      size    = 3.0,
      hjust   = 0,
      vjust   = 1,
      family  = "mono",
      lineheight = 1.05
    )
  }

  # ── Inset: density curves (bottom-right) ─────────────────────────────────────
  # S2 values are the same for both norms; take from DSM to avoid duplicates
  s2_dens  <- data.table::data.table(
    x = site_dt[Norm == norms_select[[1L]], LAI_S2], grp = "S2"
  )
  chm_dens <- data.table::data.table(
    x = site_dt[Norm == "DSM_keepTrees", LAI_ALS], grp = "CHM"
  )
  dtm_dens <- data.table::data.table(
    x = site_dt[Norm == "DTM_keepTrees", LAI_ALS], grp = "DTM"
  )
  hist_dt <- data.table::rbindlist(list(s2_dens, chm_dens, dtm_dens))
  hist_dt[, grp := factor(grp, levels = c("S2", "CHM", "DTM"))]

  dens_clr <- c(S2  = "grey20",
                CHM = norm_colours[["DSM_keepTrees"]],
                DTM = norm_colours[["DTM_keepTrees"]])
  dens_lty <- c(S2 = "solid", CHM = "dashed", DTM = "dotted")

  h <- ggplot2::ggplot(hist_dt,
                        ggplot2::aes(x = x, colour = grp, linetype = grp)) +
    ggplot2::geom_density(linewidth = 0.5, adjust = 1.2) +
    ggplot2::scale_colour_manual(values = dens_clr, name = NULL) +
    ggplot2::scale_linetype_manual(values = dens_lty, name = NULL) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      legend.position      = c(0.97, 0.97),
      legend.justification = c(1, 1),
      legend.key.height    = ggplot2::unit(0.30, "cm"),
      legend.key.width     = ggplot2::unit(0.55, "cm"),
      legend.text          = ggplot2::element_text(size = 6.5),
      legend.background    = ggplot2::element_rect(fill  = "white",
                                                    colour = NA),
      legend.margin        = ggplot2::margin(1, 2, 1, 2),
      panel.grid           = ggplot2::element_blank(),
      axis.text            = ggplot2::element_text(size = 6),
      axis.ticks           = ggplot2::element_line(linewidth = 0.3),
      plot.background      = ggplot2::element_rect(fill   = "white",
                                                    colour = "grey65",
                                                    linewidth = 0.4)
    )

  # Inset at bottom-right: 55-100% x, 0-40% y of panel
  p + patchwork::inset_element(h, left = 0.55, right = 1.0,
                                 bottom = 0.0,  top  = 0.40,
                                 align_to = "panel")
}

# ── Assemble & add shared norm legend ─────────────────────────────────────────

cli::cli_h1("Building figure...")

panel_list <- list(
  make_panel("Aigoual", show_y = TRUE),
  make_panel("Blois",   show_y = FALSE),
  make_panel("Mormal",  show_y = FALSE)
)

# Shared CHM / DTM colour legend via a thin guide panel
guide_dt <- data.table::data.table(
  x    = c(1, 2),
  y    = c(1, 1),
  Norm = norms_select
)
p_guide <- ggplot2::ggplot(guide_dt,
                            ggplot2::aes(x = x, y = y, colour = Norm)) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::scale_colour_manual(values = norm_colours,
                                labels = norm_labels, name = NULL) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      direction = "horizontal",
      override.aes = list(linewidth = 1.5)
    )
  ) +
  ggplot2::theme_void(base_size = 13) +
  ggplot2::theme(
    legend.position = "top",
    legend.key.width = ggplot2::unit(1.2, "cm"),
    legend.text = ggplot2::element_text(size = 11)
  )

fig <- patchwork::wrap_plots(panel_list, nrow = 1) /
  p_guide +
  patchwork::plot_layout(heights = c(20, 1)) +
  patchwork::plot_annotation(
    caption = paste0(
      "Grey dashed: 1:1 line. Solid lines: OLS (y = LAI_ALS ~ LAI_S2). ",
      "Bias = mean(S2 - LiDAR)."
    )
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_dir  <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "scatter_LAI_ALS_vs_S2_ATBD.pdf")
ggplot2::ggsave(out_file, fig, width = 28, height = 12, units = "cm",
                 device = "pdf")
cli::cli_alert_success("Written: {out_file}")
cat("Done.\n")
