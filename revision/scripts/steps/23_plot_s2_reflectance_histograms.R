# ---
# title:  23_plot_s2_reflectance_histograms.R
# desc:   Histograms of Sentinel-2 reflectance bands B03, B04, B08 for
#         Aigoual, Blois, and Mormal.
#         Pixels filtered: max_height > 10 m AND fCover > 0.9 (same criteria
#         as the stratified sampling used in the PROSAIL optimisation).
#         Layout: 3 rows (B03 / B04 / B08) × 3 columns (sites),
#         site strip labels shown on first row only.
#
#         Reads (from 03_RESULTS — read-only):
#           {ext_results}/{site}/Metrics/Deciduous_Only/max_res_10_m.tif
#           {ext_results}/{site}/Metrics/Deciduous_Only/fCover_res_10_m.tif
#           03_RESULTS/{site}/{acq_id}/Reflectance/res_10_m/{acq_id}_Refl
#
#         Outputs:
#           revision/output/figures/reviewers/s2_reflectance_histograms.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/23_plot_s2_reflectance_histograms.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(patchwork)
library(terra)
library(cli)

source(here::here("revision", "R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites      <- c("Aigoual", "Blois", "Mormal")
bands_sel  <- c("B03", "B04", "B08")   # bands 2, 3, 7 in the raster
h_min      <- 10        # minimum canopy height (m)
fcover_thr <- 0.9       # fCover threshold (0–1 scale; adjust if stored as 0–100)
n_max      <- 100000L   # max pixels per site (subsampled for speed)

# Acquisition IDs for the summer 2021 images used in the study
acq_ids <- c(
  Aigoual = "L2A_T31TEJ_A031608_20210711T104217",
  Blois   = "L2A_T31TCN_A031222_20210614T105443",
  Mormal  = "L2A_T31UER_A031222_20210614T105443"
)

out_dir <- file.path(paths$output, "figures", "reviewers")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Site colours (consistent with rest of the pipeline)
site_colours <- c(Aigoual = "#1b7837", Blois = "#762a83", Mormal = "#d6604d")

# ── Load reflectance values per site ──────────────────────────────────────────

cli::cli_h1("Loading S2 reflectance bands...")

all_vals <- vector("list", length(sites))
names(all_vals) <- sites

for (site in sites) {
  base_ext  <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only")
  max_path  <- file.path(base_ext, "max_res_10_m.tif")
  fcov_path <- file.path(base_ext, "fCover_res_10_m.tif")
  acq_id    <- acq_ids[[site]]
  refl_path <- file.path(
    paths$ext_results, site, "PROSAIL_Optimization", acq_id,
    "Reflectance",
    paste0(acq_id, "_Refl")
  )

  for (p in c(max_path, fcov_path, refl_path)) {
    if (!file.exists(p)) stop("Missing file:\n  ", p)
  }

  cli::cli_alert_info("{site}: reading rasters...")

  refl  <- terra::rast(refl_path)[[c(2L, 3L, 7L)]]  # B03, B04, B08
  names(refl) <- bands_sel

  max_r  <- terra::rast(max_path)
  fcov_r <- terra::rast(fcov_path)

  # Detect fCover scale (0–1 vs 0–100) from range
  fcov_max <- terra::global(fcov_r, "max", na.rm = TRUE)[[1L]]
  thr <- if (fcov_max > 1.5) fcover_thr * 100 else fcover_thr

  # Build valid-pixel mask
  mask_valid <- (max_r > h_min) & (fcov_r > thr)

  # Apply mask to reflectance stack
  refl_masked <- terra::mask(refl, mask_valid, maskvalues = 0)

  # Extract values
  vals <- as.data.table(terra::values(refl_masked, na.rm = FALSE))
  vals <- vals[complete.cases(vals)]

  # Subsample if needed
  if (nrow(vals) > n_max) vals <- vals[sample(.N, n_max)]

  vals[, site := site]
  all_vals[[site]] <- vals

  cli::cli_alert_success("{site}: {nrow(vals)} pixels retained")
}

dt <- data.table::rbindlist(all_vals, use.names = TRUE)
dt[, site := factor(site, levels = sites)]

# ── Pivot to long format ────────────────────────────────────────────────────────

dt_long <- data.table::melt(
  dt,
  id.vars       = "site",
  measure.vars  = bands_sel,
  variable.name = "band",
  value.name    = "reflectance"
)
dt_long[, band := factor(band, levels = bands_sel)]

# ── Build plots — one per band (style matching legacy exactly) ─────────────────

theme_legacy <- ggplot2::theme_bw(base_size = 20) +
  ggplot2::theme(
    axis.title       = ggplot2::element_text(size = 24),
    axis.text        = ggplot2::element_text(size = 18, color = "black"),
    strip.text       = ggplot2::element_text(size = 24, face = "bold"),
    strip.background = ggplot2::element_rect(fill = "white"),
    legend.position  = "none",
    plot.tag         = ggplot2::element_text(size = 30)
  )

plot_list <- lapply(seq_along(bands_sel), function(i) {
  b   <- bands_sel[[i]]
  sub <- dt_long[band == b]

  p <- ggplot2::ggplot(sub, ggplot2::aes(x = reflectance, fill = site)) +
    ggplot2::geom_histogram(bins = 100, colour = "black", linewidth = 0.1,
                             alpha = 0.8) +
    ggplot2::facet_grid(~ site) +
    ggplot2::scale_fill_viridis_d(option = "viridis", name = "Site") +
    ggplot2::labs(x = "Reflectance", y = paste(b, "Count")) +
    theme_legacy

  # Hide site strip labels on rows 2 and 3
  if (i > 1L) p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())

  p
})

# ── Assemble and tag ────────────────────────────────────────────────────────────

fig <- (plot_list[[1L]] / plot_list[[2L]] / plot_list[[3L]]) +
  patchwork::plot_annotation(
    tag_levels = list(c("(a)", "(b)", "(c)"))
  ) &
  ggplot2::theme(plot.tag.position = c(0.01, 1))

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "s2_reflectance_histograms.pdf")
out_png <- file.path(out_dir, "s2_reflectance_histograms.png")

ggplot2::ggsave(out_pdf, fig, width = 16, height = 18, units = "in",
                bg = "white", device = cairo_pdf)
ggplot2::ggsave(out_png, fig, width = 16, height = 18, units = "in",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — S2 reflectance histograms")
