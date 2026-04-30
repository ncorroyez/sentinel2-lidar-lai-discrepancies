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
#           revision/output/figures/reviewers/pad_profiles_boxplot.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/28_plot_pad_profiles.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(terra)
library(cli)

source(here::here("revision", "R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites      <- c("Aigoual", "Blois", "Mormal")
max_height <- 40L

ladstack_subpath <- file.path("Metrics", "Deciduous_Only", "ladstack_classic.tif")

out_dir <- file.path(paths$output, "figures", "reviewers")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Load rasters ───────────────────────────────────────────────────────────────

cli::cli_h1("Loading PAD rasters...")

load_site_ladstack <- function(site) {
  lad_path <- file.path(paths$ext_results, site, ladstack_subpath)
  if (!file.exists(lad_path)) stop("Missing file:\n  ", lad_path)

  lad     <- terra::rast(lad_path)
  heights <- as.numeric(sub("LAD_Layer_", "", names(lad)))
  depths  <- max_height - heights + 0.5   # depth from canopy top (m)

  mat <- terra::values(lad, mat = TRUE)
  dt_list <- lapply(seq_along(heights), function(i) {
    vals <- mat[, i]
    vals <- vals[!is.na(vals) & vals >= 0]
    data.table(PAD_value = vals, Depth = depths[[i]], Norm = "DSM", Site = site)
  })
  data.table::rbindlist(dt_list)
}

dt_list <- lapply(sites, function(site) {
  cli::cli_progress_step("Loading {site}")
  load_site_ladstack(site)
})

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

norm_colors <- c(DSM = "#F8766D")
norm_labels <- c(DSM = expression(DSM[ALS]))

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
  ggplot2::scale_x_continuous(limits = c(0, 15)) +
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

ggplot2::ggsave(out_pdf, p, width = 10, height = 10, units = "in")
ggplot2::ggsave(out_png, p, width = 10, height = 10, units = "in",
                device = "png", dpi = 300, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — PAD profiles boxplot")
