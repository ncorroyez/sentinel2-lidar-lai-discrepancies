# ---
# title:  28_plot_pad_profiles.R
# desc:   Boxplots of ALS-derived PAD (Plant Area Density) values by canopy
#         depth for Aigoual, Blois, and Mormal, comparing two height-
#         normalisation procedures:
#           CHM_ALS: heights normalised via DSM (DSM_keepTrees2)
#           DTM_ALS: heights normalised via DTM (DTM_keepTrees2)
#         Faceted by site; depth axis reversed (top of canopy at top).
#
#         Source: 02_CODES/LiDAR/pad_boxplots.R
#
#         Reads (read-only):
#           {ext_results}/{site}/Metrics/Deciduous_Only/
#             ladstack_classic.tif   (38 layers: LAD_Layer_2.5 … LAD_Layer_39.5)
#
#         Outputs:
#           output/figures/pad_profiles_boxplot.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("scripts/28_plot_pad_profiles.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(terra)
library(cli)

source(here::here("R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites      <- c("Aigoual", "Blois", "Mormal")
max_height <- 40L

# k rescaling — the ladstack rasters were computed with k_ref = 0.5
# (Bouvier default at PAD precomputation); this study uses k_select = 0.65.
# Beer–Lambert is multiplicative, so LAI rescales by k_ref / k_select.
k_ref      <- 0.5
k_select   <- 0.65
k_scale    <- k_ref / k_select   # ≈ 0.7692

# Two normalisation procedures, each with its own ladstack:
#   CHM (canopy-height-model based) → ladstack_dsm.tif
#   DTM (digital-terrain-model based) → ladstack_dtm.tif
norm_files <- list(
  CHM = "ladstack_dsm.tif",
  DTM = "ladstack_dtm.tif"
)

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Load rasters ───────────────────────────────────────────────────────────────
# For each (site, norm), we read the per-bin LAD raster (38 layers, 1 m bins
# named LAD_Layer_2.5 … LAD_Layer_39.5), reorder layers from canopy top to
# bottom, and compute the cumulative sum so that the column "PAD_value" at
# depth d is the LAI integrated from the canopy top down to depth d.

cli::cli_h1("Loading PAD rasters (CHM + DTM, cumulative from top)...")

load_site_ladstack <- function(site, norm_label, fn) {
  lad_path <- file.path(paths$ext_results, site, "Metrics",
                        "Deciduous_Only", fn)
  if (!file.exists(lad_path)) {
    cli::cli_warn("Missing raster for {site} / {norm_label}: {lad_path}")
    return(NULL)
  }

  lad     <- terra::rast(lad_path)
  heights <- as.numeric(sub("LAD_Layer_", "", names(lad)))
  depths  <- max_height - heights + 0.5   # depth from canopy top (m)

  # Reorder layers so depth grows top → bottom (depth = 1 at top, 38 at bottom)
  ord <- order(depths)
  mat <- terra::values(lad, mat = TRUE)[, ord, drop = FALSE]
  depths_ord <- depths[ord]

  # Cumulative LAI from canopy top down to each depth bin, rescaled to
  # k_select (Beer–Lambert is multiplicative).
  mat[is.na(mat)] <- 0
  mat <- mat * k_scale
  cum_mat <- t(apply(mat, 1L, cumsum))
  # cum_mat[i, j] = LAI integrated from top down to depths_ord[j]

  dt_list <- lapply(seq_along(depths_ord), function(j) {
    vals <- cum_mat[, j]
    vals <- vals[!is.na(vals) & vals > 0]
    data.table(PAD_value = vals, Depth = depths_ord[[j]],
               Norm = norm_label, Site = site)
  })
  data.table::rbindlist(dt_list)
}

dt_list <- list()
for (site in sites) {
  for (norm_label in names(norm_files)) {
    cli::cli_progress_step("Loading {site} / {norm_label}")
    dt_one <- load_site_ladstack(site, norm_label, norm_files[[norm_label]])
    if (!is.null(dt_one)) dt_list[[length(dt_list) + 1L]] <- dt_one
  }
}

dt <- data.table::rbindlist(dt_list)
cli::cli_progress_done()

dt[, Site := factor(Site, levels = sites)]
dt[, Norm := factor(Norm, levels = c("CHM", "DTM"))]

depth_levels <- sort(unique(dt$Depth), decreasing = TRUE)
dt[, Depth_f := factor(Depth, levels = depth_levels)]

cli::cli_alert_success("Loaded {nrow(dt)} PAD pixel-depth records")

# ── Pre-compute boxplot statistics (avoids ggplot2 OOM on 78M rows) ────────────

cli::cli_progress_step("Computing box statistics")

stats_dt <- dt[, {
  q   <- quantile(PAD_value, probs = c(0.25, 0.50, 0.75))
  iqr <- q[3] - q[1]
  .(
    xmin    = max(0,  q[1] - 1.5 * iqr),
    xlower  = q[1],
    xmiddle = q[2],
    xupper  = q[3],
    xmax    = min(15, q[3] + 1.5 * iqr)
  )
}, by = .(Site, Norm, Depth_f)]

cli::cli_progress_done()
rm(dt); gc()

# ── Plot ───────────────────────────────────────────────────────────────────────

norm_colors <- c(CHM = "#F8766D", DTM = "#00BFC4")
norm_labels <- c(CHM = expression(CHM[ALS]), DTM = expression(DTM[ALS]))

pd <- ggplot2::position_dodge(width = 0.75)

p <- ggplot2::ggplot(stats_dt, ggplot2::aes(y = Depth_f, fill = Norm, colour = Norm)) +
  # whiskers (colored by norm)
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = xmin, xmax = xmax),
    width       = 0,
    linewidth   = 0.5,
    orientation = "y",
    position    = pd
  ) +
  # box body + median line (fill = norm color, borders + median = black)
  ggplot2::geom_crossbar(
    ggplot2::aes(x = xmiddle, xmin = xlower, xmax = xupper),
    colour           = "black",
    middle.linewidth = 1.5,
    width            = 0.65,
    position         = pd
  ) +
  ggplot2::facet_wrap(~ Site, scales = "fixed", nrow = 1) +
  ggplot2::scale_y_discrete(drop = FALSE) +
  ggplot2::scale_x_continuous(limits = c(0, 12)) +
  ggplot2::scale_fill_manual(values = norm_colors, labels = norm_labels,
                              name = "Normalization") +
  ggplot2::scale_colour_manual(values = norm_colors, labels = norm_labels,
                                name = "Normalization") +
  ggplot2::labs(x = "PAD value", y = "Depth (m)") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    strip.text       = ggplot2::element_text(face = "bold", size = 13),
    axis.text.y      = ggplot2::element_text(size = 9),
    axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "right"
  )

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "pad_profiles_boxplot.pdf")
out_png <- file.path(out_dir, "pad_profiles_boxplot.png")

ggplot2::ggsave(out_pdf, p, width = 10, height = 10, units = "in", device = cairo_pdf)
ggplot2::ggsave(out_png, p, width = 10, height = 10, units = "in",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — PAD profiles boxplot")
