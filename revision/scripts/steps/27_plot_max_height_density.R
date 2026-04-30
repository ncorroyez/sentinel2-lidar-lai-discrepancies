# ---
# title:  27_plot_max_height_density.R
# desc:   Density histograms of ALS maximum canopy height for Aigoual, Blois,
#         and Mormal, with 99.9th percentile thresholds (Q99.9) shown as
#         dashed vertical lines and annotated labels.
#
#         Reads (read-only):
#           {ext_results}/{site}/Metrics/Deciduous_Only/max_res_10_m.tif
#
#         Outputs:
#           revision/output/figures/reviewers/max_height_density.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/27_plot_max_height_density.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(terra)
library(cli)

source(here::here("revision", "R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites     <- c("Aigoual", "Blois", "Mormal")
q_prob    <- 0.999   # 99.9th percentile

out_dir <- file.path(paths$output, "figures", "reviewers")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Load rasters and extract values ────────────────────────────────────────────

cli::cli_h1("Loading max height rasters...")

dt_list <- lapply(sites, function(site) {
  path <- file.path(paths$ext_results, site, "Metrics", "Deciduous_Only",
                    "max_res_10_m.tif")
  if (!file.exists(path)) stop("Missing: ", path)
  vals <- terra::values(terra::rast(path), mat = FALSE)
  vals <- vals[!is.na(vals)]
  q999 <- quantile(vals, probs = q_prob, na.rm = TRUE)
  cli::cli_alert_success("{site}: {length(vals)} pixels, Q99.9 = {round(q999, 2)} m")
  data.table(max_height = vals, site = site, q999 = q999)
})

dt <- data.table::rbindlist(dt_list)
dt[, site := factor(site, levels = sites)]

# ── Quantile lines and label positions ─────────────────────────────────────────

q_dt <- dt[, .(q999 = first(q999)), by = site]
q_dt[, label := paste0(site, ": Q99.9 = ", round(q999, 2))]

# Y positions: stack labels from top, spaced by 8% of max density
density_max <- max(sapply(sites, function(s) {
  max(density(dt[site == s, max_height])$y)
}))
q_dt[, y := seq(density_max * 0.95,
                by = -density_max * 0.08,
                length.out = .N)]

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggplot2::ggplot(dt, ggplot2::aes(x = max_height,
                                       colour = site, fill = site)) +
  ggplot2::geom_density(alpha = 0.7) +
  ggplot2::geom_vline(
    data      = q_dt,
    ggplot2::aes(xintercept = q999, colour = site),
    linetype  = "dashed", linewidth = 1
  ) +
  ggplot2::geom_text(
    data    = q_dt,
    ggplot2::aes(x = Inf, y = y, label = label, colour = site),
    hjust   = 1.1, size = 5.5, show.legend = FALSE
  ) +
  ggplot2::labs(x = "Max Height (m)", y = "Density") +
  ggplot2::theme_minimal(base_size = 16)

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "max_height_density.pdf")
out_png <- file.path(out_dir, "max_height_density.png")

ggplot2::ggsave(out_pdf, p, width = 10, height = 10, units = "in")
ggplot2::ggsave(out_png, p, width = 10, height = 10, units = "in",
                device = "png", dpi = 300)

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — max height density plot")
