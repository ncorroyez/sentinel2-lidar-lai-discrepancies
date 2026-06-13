# ---
# title:  08b_sm5_scatter_lai_atbd_lowcanopy.R
# desc:   Companion to 08_sm5_scatter_lai_atbd.R restricted to LOW-CANOPY
#         pixels (max_height <= 10 m). Same layout: one column per site,
#         2D bin-count scatter LAI_ALS (x) vs LAI_S2_ATBD (y), regression
#         + 1:1 line + stats annotation at top-left.
#
#         Sample sizes (deciduous-only, h_max <= 10 m, k = 0.65):
#           Aigoual ≈ 1.3k, Blois ≈ 12k, Mormal ≈ 2.5k pixels.
#
#         Source rasters:
#           LAI_ALS:    data/results/{site}/Metrics/Deciduous_Only/
#                       lidarlai_res_10_m.tif  (CHM_keepTrees, rescaled to
#                       k_select = 0.65)
#           LAI_S2:     output/intermediate/sm6/{site}/
#                       s2lai_summer_atbd_T_res_10_m.tif  (prosail 3.0.0 ATBD_T)
#           max height: data/results/{site}/Metrics/Deciduous_Only/
#                       max_res_10_m.tif (used to keep pixels with h <= 10 m)
#
#         Outputs:
#           output/figures/scatter_LAI_ALS_vs_S2_ATBD_h_le_10.{pdf,png}
#
# Run from project root (NC_Full/):
#   source("scripts/steps/08b_sm5_scatter_lai_atbd_lowcanopy.R")
# ---

library("here")
library("data.table")
library("terra")
library("ggplot2")
library("patchwork")
library("cli")

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_metrics.R"))

# ── Site metadata ─────────────────────────────────────────────────────────────

site_meta <- list(
  Aigoual = list(s2_date = "2021-07-11"),
  Blois   = list(s2_date = "2021-06-14"),
  Mormal  = list(s2_date = "2021-06-14")
)
sites    <- names(site_meta)
xy_lim   <- c(0, 15)
h_max_m  <- 10   # keep pixels with max canopy height <= h_max_m

# k rescaling — same as 08
k_ref    <- 0.5
k_select <- 0.65
k_scale  <- k_ref / k_select   # ≈ 0.7692

# ── Load all-pixel paired data + height filter ────────────────────────────────

cli::cli_h1("Loading raster data (max height <= {h_max_m} m)...")

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")

load_raster_pairs <- function(site) {
  base <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")

  als_path <- file.path(base, "lidarlai_res_10_m.tif")
  hmax_path <- file.path(base, "max_res_10_m.tif")
  s2_path  <- file.path(sm6a_dir, site, "s2lai_summer_atbd_T_res_10_m.tif")

  for (p in c(als_path, hmax_path, s2_path)) {
    if (!file.exists(p)) {
      cli::cli_warn("Missing raster for {site}: {p}")
      return(NULL)
    }
  }

  r_als  <- terra::rast(als_path)
  r_hmax <- terra::rast(hmax_path)
  r_s2   <- terra::rast(s2_path)

  if (!terra::compareGeom(r_als, r_s2, stopOnError = FALSE)) {
    r_s2 <- terra::resample(r_s2, r_als, method = "bilinear")
  }
  if (!terra::compareGeom(r_als, r_hmax, stopOnError = FALSE)) {
    r_hmax <- terra::resample(r_hmax, r_als, method = "near")
  }

  als_v  <- terra::values(r_als,  mat = FALSE) * k_scale   # rescale to k_select
  s2_v   <- terra::values(r_s2,   mat = FALSE)
  hmax_v <- terra::values(r_hmax, mat = FALSE)

  valid <- !is.na(als_v) & !is.na(s2_v) & !is.na(hmax_v) &
           als_v >= 0 & s2_v >= 0 & hmax_v <= h_max_m

  cli::cli_alert_success(
    "{site}: {format(sum(valid), big.mark=',')} pixels with h <= {h_max_m} m"
  )

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
  fit <- stats::lm(LAI_S2 ~ LAI_ALS)
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

cli::cli_h2("Metrics (low-canopy, h_max <= {h_max_m} m)")
print(stats_dt)

# ── Panel builder ─────────────────────────────────────────────────────────────

make_scatter <- function(site_name, show_y = TRUE, show_x = TRUE,
                          show_legend = FALSE) {
  site_dt <- dt[Site == site_name]
  st      <- stats_dt[Site == site_name]

  x_ann <- 0.4
  y_top <- 14.65
  dy    <- 1.20
  sz    <- 5.5

  ggplot2::ggplot(site_dt, ggplot2::aes(x = LAI_ALS, y = LAI_S2)) +
    # Coarser binning (50 vs 400) — sample sizes here are ~100x smaller than
    # the all-pixel scatter (08), so bins=400 would leave most cells at
    # count=1 with no contrast.
    ggplot2::geom_bin2d(bins = 50) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                          se = FALSE, colour = "firebrick", linewidth = 0.85) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dashed", colour = "grey55", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(
      name   = "Count",
      option = "viridis",
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
    ggplot2::annotate("text", x = x_ann, y = y_top,
      label = paste0("y = ", st$a, "x", st$b_str),
      hjust = 0, vjust = 1, size = sz) +
    ggplot2::annotate("text", x = x_ann, y = y_top - dy,
      label = paste0("RMSE = ", st$RMSE),
      hjust = 0, vjust = 1, size = sz) +
    ggplot2::annotate("text", x = x_ann, y = y_top - 2 * dy,
      label = paste0("Bias = ", st$Bias),
      hjust = 0, vjust = 1, size = sz) +
    ggplot2::annotate("text", x = x_ann, y = y_top - 3 * dy,
      label = paste0("italic(r)~'= ", st$r, "'"),
      hjust = 0, vjust = 1, size = sz, parse = TRUE) +
    ggplot2::annotate("text", x = x_ann, y = y_top - 4 * dy,
      label = paste0("R^2~'= ", st$R2, "'"),
      hjust = 0, vjust = 1, size = sz, parse = TRUE) +
    ggplot2::labs(
      title = site_name,
      x     = if (show_x) expression(LAI[ALS]) else NULL,
      y     = if (show_y) expression(LAI[S2_ATBD]) else NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title.x     = if (!show_x) ggplot2::element_blank() else ggplot2::element_text(),
      axis.title.y     = if (!show_y) ggplot2::element_blank() else ggplot2::element_text()
    )
}

# ── Assemble & save ───────────────────────────────────────────────────────────

cli::cli_h1("Building figure...")

p1 <- make_scatter("Aigoual", show_y = TRUE,  show_x = FALSE, show_legend = FALSE)
p2 <- make_scatter("Blois",   show_y = FALSE, show_x = TRUE,  show_legend = FALSE)
p3 <- make_scatter("Mormal",  show_y = FALSE, show_x = FALSE, show_legend = TRUE)

fig <- (p1 | p2 | p3)

out_dir  <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_pdf <- file.path(out_dir, "scatter_LAI_ALS_vs_S2_ATBD_h_le_10.pdf")
out_png <- file.path(out_dir, "scatter_LAI_ALS_vs_S2_ATBD_h_le_10.png")

ggplot2::ggsave(out_pdf, fig, width = 30, height = 12, units = "cm",
                device = cairo_pdf)
ggplot2::ggsave(out_png, fig, width = 30, height = 12, units = "cm",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
