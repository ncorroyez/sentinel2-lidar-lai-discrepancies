# ---
# title:  26_plot_als_metrics_scatter.R
# desc:   Figure — density scatterplots of ALS-derived metrics for Aigoual,
#         Blois, and Mormal.
#           (a) Max Height (m) vs DSM_SD
#           (b) Max Height (m) vs LAI_ALS
#           (c) LAI_ALS        vs DSM_SD
#         Layout: 3 rows (a/b/c) × 3 columns (sites, via facet_wrap).
#         Uses geom_hex (log-count viridis fill) + lm smooth.
#
#         Note: the raster file is dsm_sd_res_10_m.tif (DSM-based SD, not CHM).
#         Labels were incorrectly written as CHM_SD in the submitted paper;
#         this script uses DSM_SD throughout.
#
#         Reads (read-only):
#           revision/output/intermediate/sm6/{site}/
#             dsm_sd_res_10_m.tif          (produced by script 14)
#           {ext_results}/{site}/Metrics/Deciduous_Only/
#             lidarlai_res_10_m.tif
#             max_res_10_m.tif
#
#         Outputs:
#           revision/output/figures/reviewers/als_metrics_scatter.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/26_plot_als_metrics_scatter.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)
library(terra)
library(cli)

source(here::here("revision", "R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites <- c("Aigoual", "Blois", "Mormal")

out_dir <- file.path(paths$output, "figures", "reviewers")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Load rasters and extract values ────────────────────────────────────────────

cli::cli_h1("Loading ALS rasters...")

load_site_values <- function(site) {
  mdir    <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")
  sm6_dir <- file.path(paths$output, "intermediate", "sm6", site)
  paths_list <- list(
    dsm_sd = file.path(sm6_dir, "dsm_sd_res_10_m.tif"),
    lai    = file.path(mdir,    "lidarlai_res_10_m.tif"),
    maxh   = file.path(mdir,    "max_res_10_m.tif")
  )
  for (p in paths_list)
    if (!file.exists(p)) stop("Missing file:\n  ", p)

  r <- terra::rast(unlist(paths_list))
  df <- as.data.table(terra::values(r, na.rm = FALSE))
  setnames(df, c("dsm_sd", "lai", "maxh"))
  df <- df[complete.cases(df)]
  df[, site := site]
  cli::cli_alert_success("{site}: {nrow(df)} valid pixels")
  df
}

dt <- data.table::rbindlist(lapply(sites, load_site_values))
dt[, site := factor(site, levels = sites)]

# ── Shared theme ───────────────────────────────────────────────────────────────

theme_pub <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    strip.text       = ggplot2::element_text(face = "bold", size = 11),
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "right"
  )

hex_scale <- ggplot2::scale_fill_viridis_c(
  option    = "viridis",
  trans     = "log10",
  name      = "Count",
  breaks    = c(1, 10, 100, 1000),
  labels    = c("1", "10", "100", "1000")
)

# ── Panel (a): Max Height vs DSM_SD ────────────────────────────────────────────

p_a <- ggplot2::ggplot(dt, ggplot2::aes(x = dsm_sd, y = maxh)) +
  ggplot2::geom_hex(bins = 200) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "blue",
                        linewidth = 0.8) +
  hex_scale +
  ggplot2::facet_wrap(~ site, nrow = 1) +
  ggplot2::labs(
    x = expression(DSM[SD]),
    y = "Max Height (m)"
  ) +
  theme_pub

# ── Panel (b): Max Height vs LAI_ALS ───────────────────────────────────────────

p_b <- ggplot2::ggplot(dt, ggplot2::aes(x = lai, y = maxh)) +
  ggplot2::geom_hex(bins = 200) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "blue",
                        linewidth = 0.8) +
  hex_scale +
  ggplot2::facet_wrap(~ site, nrow = 1) +
  ggplot2::labs(
    x = expression(LAI[ALS]),
    y = "Max Height (m)"
  ) +
  theme_pub

# ── Panel (c): LAI_ALS vs DSM_SD ───────────────────────────────────────────────

p_c <- ggplot2::ggplot(dt, ggplot2::aes(x = dsm_sd, y = lai)) +
  ggplot2::geom_hex(bins = 200) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "blue",
                        linewidth = 0.8) +
  hex_scale +
  ggplot2::facet_wrap(~ site, nrow = 1) +
  ggplot2::labs(
    x = expression(DSM[SD]),
    y = expression(LAI[ALS])
  ) +
  theme_pub

# ── Assemble ────────────────────────────────────────────────────────────────────

fig <- (p_a / p_b / p_c) +
  patchwork::plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)"))
  ) &
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 13, face = "bold"),
    plot.tag.position = c(0.01, 1)
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "als_metrics_scatter.pdf")
out_png <- file.path(out_dir, "als_metrics_scatter.png")

ggplot2::ggsave(out_pdf, fig, width = 12, height = 10, units = "in")
ggplot2::ggsave(out_png, fig, width = 12, height = 10, units = "in",
                device = "png", dpi = 300, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — ALS metrics density scatterplots")
