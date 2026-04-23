# ---
# title:  07_sm5_scatter_lai_atbd.R
# desc:   2-row x 3-column figure (one column per site).
#         Row 1: 2D count scatter — LAI_ALS DTM (y) vs LAI_S2_ATBD (x),
#                axes fixed [0, 15]. Top-left: y = ax+b, r, R2, RMSE, Bias.
#         Row 2: overlapping histograms of LAI_ALS and LAI_S2, full [0, 15] range.
#
#         Source rasters (all forest pixels, not masked, no keepTrees):
#           LAI_ALS: 03_RESULTS/{site}/Metrics/Not_Masked/
#                    lidarlai_optim_depth_res_10_m.tif  (DTM, d_opt)
#           LAI_S2:  03_RESULTS/{site}/Metrics/Not_Masked/
#                    s2lai_{date}_atbd_res_10_m.tif
#
# Run from project root (NC_Full/):
#   source("revision/scripts/07_sm5_scatter_lai_atbd.R")
# ---

library("here")
library("data.table")
library("terra")
library("ggplot2")
library("patchwork")
library("cli")

source(here::here("revision", "R", "sm5_metrics.R"))

# ── Site metadata ─────────────────────────────────────────────────────────────

site_meta <- list(
  Aigoual = list(s2_date = "2021-07-11"),
  Blois   = list(s2_date = "2021-06-14"),
  Mormal  = list(s2_date = "2021-06-14")
)
sites  <- names(site_meta)
xy_lim <- c(0, 15)

# ── Load all-pixel paired data from rasters ───────────────────────────────────

cli::cli_h1("Loading raster data (all pixels)...")

load_raster_pairs <- function(site) {
  meta <- site_meta[[site]]
  base <- here::here("03_RESULTS", site, "Metrics", "Not_Masked")

  als_path <- file.path(base, "lidarlai_res_10_m.tif")
  s2_path  <- file.path(base,
    paste0("s2lai_", meta$s2_date, "_atbd_res_10_m.tif"))

  if (!file.exists(als_path)) {
    cli::cli_warn("ALS raster missing: {als_path}")
    return(NULL)
  }
  if (!file.exists(s2_path)) {
    cli::cli_warn("S2 raster missing: {s2_path}")
    return(NULL)
  }

  r_als <- terra::rast(als_path)
  r_s2  <- terra::rast(s2_path)

  # Align S2 to ALS grid if needed
  if (!terra::compareGeom(r_als, r_s2, stopOnError = FALSE)) {
    r_s2 <- terra::resample(r_s2, r_als, method = "bilinear")
  }

  als_v <- terra::values(r_als, mat = FALSE)
  s2_v  <- terra::values(r_s2,  mat = FALSE)

  valid <- !is.na(als_v) & !is.na(s2_v) & als_v >= 0 & s2_v >= 0
  cli::cli_alert_success("{site}: {format(sum(valid), big.mark=',')} valid pixels")

  data.table::data.table(
    Site    = site,
    LAI_ALS = als_v[valid],
    LAI_S2  = s2_v[valid]
  )
}

dt <- data.table::rbindlist(lapply(sites, load_raster_pairs), use.names = TRUE)

# ── Statistics ────────────────────────────────────────────────────────────────

stats_dt <- dt[, {
  m   <- compute_metrics_s2_lidar(LAI_S2, LAI_ALS)
  fit <- stats::lm(LAI_ALS ~ LAI_S2)
  a   <- round(stats::coef(fit)[[2L]], 2)
  b   <- round(stats::coef(fit)[[1L]], 2)
  b_str <- if (b >= 0) paste0("+", b) else as.character(b)
  list(
    label = paste0(
      "y = ", a, "x ", b_str, "\n",
      "r = ",   round(m$R,    2), "   R2 = ", round(m$R2,   2), "\n",
      "RMSE = ", round(m$RMSE, 2), "   Bias = ",   round(m$Bias,  2)
    )
  )
}, by = Site]

cli::cli_h2("Metrics")
print(stats_dt[, .(Site, label)])

# ── Count colour scale ─────────────────────────────────────────────────────────

count_breaks <- c(1, 10, 100, 1000, 5000)

# ── Panel builders ─────────────────────────────────────────────────────────────

make_scatter <- function(site_name, show_y = TRUE, show_legend = FALSE) {
  site_dt  <- dt[Site == site_name]
  lbl      <- stats_dt[Site == site_name, label]
  n_obs    <- format(nrow(site_dt), big.mark = ",")

  p <- ggplot2::ggplot(site_dt,
                        ggplot2::aes(x = LAI_S2, y = LAI_ALS)) +
    ggplot2::geom_bin2d(binwidth = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                          se = FALSE, colour = "firebrick",
                          linewidth = 0.85) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dashed", colour = "grey75",
                          linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(
      name   = "Count",
      trans  = "log10",
      breaks = count_breaks,
      labels = scales::comma(count_breaks),
      option = "plasma",
      guide  = if (show_legend)
                 ggplot2::guide_colourbar(
                   barheight = ggplot2::unit(3.5, "cm"),
                   barwidth  = ggplot2::unit(0.4, "cm"),
                   title.position = "top"
                 )
               else "none"
    ) +
    ggplot2::coord_equal(xlim = xy_lim, ylim = xy_lim, expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::annotate("text",
      x = 0.4, y = 14.7,
      label  = lbl, hjust = 0, vjust = 1,
      size   = 3.0, colour = "white",
      family = "mono", lineheight = 1.05
    ) +
    ggplot2::annotate("text",
      x = 14.7, y = 0.3,
      label  = paste0("n = ", n_obs),
      hjust  = 1, vjust = 0,
      size   = 2.8, colour = "white"
    ) +
    ggplot2::labs(
      title = site_name,
      x     = expression(LAI[S2]~"(ATBD)"),
      y     = if (show_y) expression(LAI[ALS]~"(DTM)") else NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title.y     = if (!show_y) ggplot2::element_blank() else
                           ggplot2::element_text()
    )
  p
}

make_hist <- function(site_name, show_y = TRUE) {
  site_dt <- dt[Site == site_name]

  long_dt <- data.table::rbindlist(list(
    data.table::data.table(val = site_dt$LAI_S2,  var = "LAI[S2]"),
    data.table::data.table(val = site_dt$LAI_ALS, var = "LAI[ALS]")
  ))
  long_dt[, var := factor(var, levels = c("LAI[S2]", "LAI[ALS]"))]

  var_colours <- c("LAI[S2]"  = "#444444", "LAI[ALS]" = "#2166ac")

  ggplot2::ggplot(long_dt, ggplot2::aes(x = val, fill = var)) +
    ggplot2::geom_histogram(
      binwidth = 0.3, position = "identity",
      alpha = 0.55, colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = var_colours,
      labels = scales::parse_format(),
      name   = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = xy_lim,
      breaks = c(0, 5, 10, 15),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = expression(LAI~value),
      y = if (show_y) "Count" else NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position    = "bottom",
      legend.key.size    = ggplot2::unit(0.4, "cm"),
      legend.text        = ggplot2::element_text(size = 10),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.y       = if (!show_y) ggplot2::element_blank() else
                             ggplot2::element_text()
    )
}

# ── Assemble ───────────────────────────────────────────────────────────────────

cli::cli_h1("Building figure...")

s1 <- make_scatter("Aigoual", show_y = TRUE,  show_legend = FALSE)
s2 <- make_scatter("Blois",   show_y = FALSE, show_legend = FALSE)
s3 <- make_scatter("Mormal",  show_y = FALSE, show_legend = TRUE)

h1 <- make_hist("Aigoual", show_y = TRUE)
h2 <- make_hist("Blois",   show_y = FALSE)
h3 <- make_hist("Mormal",  show_y = FALSE)

fig <- (s1 | s2 | s3) / (h1 | h2 | h3) +
  patchwork::plot_layout(heights = c(3, 1.6))

# ── Save ───────────────────────────────────────────────────────────────────────

out_dir  <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "scatter_LAI_ALS_vs_S2_ATBD.pdf")
ggplot2::ggsave(out_file, fig, width = 30, height = 20, units = "cm",
                 device = "pdf")
cli::cli_alert_success("Written: {out_file}")
cat("Done.\n")
