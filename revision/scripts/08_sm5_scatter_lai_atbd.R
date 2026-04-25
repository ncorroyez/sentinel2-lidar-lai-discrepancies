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

source(here::here("revision", "R", "paths.R"))
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
  base <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")

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
  list(
    a     = a,
    b     = b,
    b_str = if (b >= 0) paste0(" + ", b) else paste0(" - ", abs(b)),
    r     = round(m$R,    2),
    R2    = round(m$R2,   2),
    RMSE  = round(m$RMSE, 2),
    Bias  = round(m$Bias, 2)
  )
}, by = Site]

cli::cli_h2("Metrics")
print(stats_dt)

# ── Count colour scale ─────────────────────────────────────────────────────────

# ── Panel builders ─────────────────────────────────────────────────────────────

make_hist_inset <- function(site_name) {
  site_dt <- dt[Site == site_name]

  long_dt <- data.table::rbindlist(list(
    data.table::data.table(val = site_dt$LAI_S2,  var = "LAI[S2_ATBD]"),
    data.table::data.table(val = site_dt$LAI_ALS, var = "LAI[ALS]")
  ))
  long_dt[, var := factor(var, levels = c("LAI[S2_ATBD]", "LAI[ALS]"))]

  var_colours <- c("LAI[S2_ATBD]" = "#444444", "LAI[ALS]" = "#2166ac")

  ggplot2::ggplot(long_dt, ggplot2::aes(x = val, fill = var)) +
    ggplot2::geom_histogram(
      binwidth = 0.5, position = "identity",
      alpha = 0.60, colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = var_colours,
      labels = scales::parse_format(),
      name   = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = xy_lim, breaks = c(0, 5, 10, 15), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position      = c(0.97, 0.97),
      legend.justification = c(1, 1),
      legend.key.height    = ggplot2::unit(0.35, "cm"),
      legend.key.width     = ggplot2::unit(0.50, "cm"),
      legend.text          = ggplot2::element_text(size = 9),
      legend.background    = ggplot2::element_rect(fill = "white",
                                                    colour = NA),
      legend.margin        = ggplot2::margin(2, 3, 2, 3),
      panel.grid           = ggplot2::element_blank(),
      axis.text            = ggplot2::element_text(size = 8.5),
      axis.ticks           = ggplot2::element_line(linewidth = 0.3),
      plot.background      = ggplot2::element_rect(fill = "white",
                                                    colour = "grey60",
                                                    linewidth = 0.4)
    )
}

make_scatter <- function(site_name, show_y = TRUE, show_x = TRUE,
                          show_legend = FALSE) {
  site_dt  <- dt[Site == site_name]
  st       <- stats_dt[Site == site_name]
  n_obs    <- format(nrow(site_dt), big.mark = ",")

  # Annotation: top-right, one metric per line
  x_ann <- 14.6
  y_top <- 14.65
  dy    <- 1.00          # inter-line spacing
  sz    <- 4.5           # text size

  p <- ggplot2::ggplot(site_dt,
                        ggplot2::aes(x = LAI_S2, y = LAI_ALS)) +
    ggplot2::geom_bin2d(bins = 400) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                          se = FALSE, colour = "firebrick",
                          linewidth = 0.85) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dashed", colour = "grey55",
                          linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(
      name   = "Count",
      option = "viridis",
      guide  = if (show_legend)
                 ggplot2::guide_colourbar(
                   barheight      = ggplot2::unit(3.5, "cm"),
                   barwidth       = ggplot2::unit(0.4, "cm"),
                   title.position = "top"
                 )
               else "none"
    ) +
    ggplot2::coord_equal(xlim = xy_lim, ylim = xy_lim, expand = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::annotate("text",
      x = x_ann, y = y_top,
      label = paste0("y = ", st$a, "x", st$b_str),
      hjust = 1, vjust = 1, size = sz, colour = "black"
    ) +
    ggplot2::annotate("text",
      x = x_ann, y = y_top - dy,
      label = paste0("RMSE = ", st$RMSE),
      hjust = 1, vjust = 1, size = sz, colour = "black"
    ) +
    ggplot2::annotate("text",
      x = x_ann, y = y_top - 2 * dy,
      label = paste0("Bias = ", st$Bias),
      hjust = 1, vjust = 1, size = sz, colour = "black"
    ) +
    ggplot2::annotate("text",
      x = x_ann, y = y_top - 3 * dy,
      label = paste0("italic(r)~'= ", st$r, "'"),
      hjust = 1, vjust = 1, size = sz, colour = "black",
      parse = TRUE
    ) +
    ggplot2::annotate("text",
      x = x_ann, y = y_top - 4 * dy,
      label = paste0("R^2~'= ", st$R2, "'"),
      hjust = 1, vjust = 1, size = sz, colour = "black",
      parse = TRUE
    ) +
    ggplot2::annotate("text",
      x = 14.6, y = 0.3,
      label = paste0("n = ", n_obs),
      hjust = 1, vjust = 0, size = 3.2, colour = "black"
    ) +
    ggplot2::labs(
      title = site_name,
      x     = if (show_x) expression(LAI[S2_ATBD]) else NULL,
      y     = if (show_y) expression(LAI[ALS]) else NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title.x     = if (!show_x) ggplot2::element_blank() else
                           ggplot2::element_text(),
      axis.title.y     = if (!show_y) ggplot2::element_blank() else
                           ggplot2::element_text()
    )

  # Inset histogram: bottom-right (55-100% x, 0-48% y of panel area)
  h <- make_hist_inset(site_name)
  p + patchwork::inset_element(h, left = 0.487, right = 1.0,
                                 bottom = 0.0,  top  = 0.513,
                                 align_to = "panel")
}

# ── Assemble ───────────────────────────────────────────────────────────────────

cli::cli_h1("Building figure...")

p1 <- make_scatter("Aigoual", show_y = TRUE,  show_x = FALSE, show_legend = FALSE)
p2 <- make_scatter("Blois",   show_y = FALSE, show_x = TRUE,  show_legend = FALSE)
p3 <- make_scatter("Mormal",  show_y = FALSE, show_x = FALSE, show_legend = TRUE)

fig <- (p1 | p2 | p3)

# ── Save ───────────────────────────────────────────────────────────────────────

out_dir  <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "scatter_LAI_ALS_vs_S2_ATBD.png")
ggplot2::ggsave(out_file, fig, width = 30, height = 12, units = "cm",
                 device = "png", dpi = 300)
cli::cli_alert_success("Written: {out_file}")
cat("Done.\n")
