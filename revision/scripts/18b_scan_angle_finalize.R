# ---
# title:  18b_scan_angle_finalize.R
# desc:   Mosaics per-tile cos_theta tifs, writes cos_theta_mean.tif and
#         lai_als_corrected.tif per site, generates PNG figures.
#         Run AFTER 18a has completed for all three sites (all tiles present).
# ---

library("here")
library("terra")
library("data.table")
library("ggplot2")

source(here::here("revision", "R", "paths.R"))

sites <- c("Aigoual", "Blois", "Mormal")

ladstack_paths <- setNames(
  file.path(paths$prosail_lidar, sites, "LiDAR",
            "PAD_Profiles_Classic", "ladstack.tif"),
  sites
)

out_root <- here::here("revision", "output", "intermediate", "scan_angle")
out_figs <- here::here("revision", "output", "figures", "scan_angle")
if (!dir.exists(out_figs)) dir.create(out_figs, recursive = TRUE)

site_stats  <- list()
lai_diag_dt <- list()

for (site in sites) {
  cat("\n══", site, "══\n")

  tmp_dir   <- file.path(out_root, site, "tmp_tiles")
  out_final <- file.path(out_root, site, "cos_theta_mean.tif")

  tile_files <- sort(list.files(tmp_dir, pattern = "_cos\\.tif$",
                                 full.names = TRUE))
  cat("  Tile tifs found:", length(tile_files), "\n")
  if (length(tile_files) == 0L) {
    cat("  SKIP — run 18a first\n"); next
  }

  # Expected tile count from LAS dir (informational)
  uid     <- system("id -u", intern = TRUE)
  las_dir <- file.path(
    paste0("/run/user/", uid, "/gvfs"),
    "smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse",
    "root/_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR_Points_Clouds",
    site, "2-las_utm"
  )
  n_expected <- length(grep("LAS_[0-9]",
                             list.files(las_dir, full.names = FALSE)))
  cat("  Expected:", n_expected, "| Have:", length(tile_files), "\n")
  if (length(tile_files) < n_expected)
    warning(site, ": only ", length(tile_files), " / ", n_expected,
            " tiles present — run 18a again")

  # Mosaic from disk
  cat("  Mosaicking...\n")
  tile_rasts <- lapply(tile_files, terra::rast)
  cos_rast   <- do.call(terra::mosaic, c(tile_rasts, list(fun = "mean")))
  terra::crs(cos_rast) <- "EPSG:32631"
  terra::writeRaster(cos_rast, out_final, overwrite = TRUE,
                     gdal = "COMPRESS=LZW")
  cat("  cos_theta_mean.tif written\n")

  # Apply to ladstack
  ladstack <- terra::rast(ladstack_paths[[site]])
  c_lyr    <- terra::resample(cos_rast, ladstack[[1L]], method = "bilinear")

  lai_orig <- terra::sum(ladstack, na.rm = TRUE)
  lai_corr <- terra::sum(ladstack * c_lyr, na.rm = TRUE)

  terra::writeRaster(lai_corr,
                     file.path(out_root, site, "lai_als_corrected.tif"),
                     overwrite = TRUE, gdal = "COMPRESS=LZW")

  c_v  <- terra::values(c_lyr,   na.rm = TRUE)
  o_v  <- terra::values(lai_orig, na.rm = TRUE)
  cr_v <- terra::values(lai_corr, na.rm = TRUE)

  cat(sprintf("  C mean: %.4f  (reduction: %.1f%%)\n",
              mean(c_v), (1 - mean(c_v)) * 100))
  cat(sprintf("  LAI: %.3f → %.3f\n",
              mean(o_v, na.rm = TRUE), mean(cr_v, na.rm = TRUE)))

  site_stats[[site]] <- data.table(
    site          = site,
    n_tiles       = length(tile_files),
    C_mean        = round(mean(c_v), 4),
    C_sd          = round(sd(c_v),   4),
    C_min         = round(min(c_v),  4),
    reduction_pct = round((1 - mean(c_v)) * 100, 2),
    lai_orig      = round(mean(o_v,  na.rm = TRUE), 3),
    lai_corr      = round(mean(cr_v, na.rm = TRUE), 3)
  )

  lai_diag_dt[[site]] <- data.table(
    site = site,
    type = rep(c("Original", "Corrected"), c(length(o_v), length(cr_v))),
    lai  = c(o_v, cr_v)
  )
}

# ── Summary ───────────────────────────────────────────────────────────────────

if (length(site_stats) == 0L) stop("No sites processed — run 18a first")
sum_dt <- data.table::rbindlist(site_stats)
cat("\n── Summary ──\n")
print(sum_dt)
data.table::fwrite(sum_dt, file.path(out_root, "scan_angle_summary.csv"))

# ── Figure ────────────────────────────────────────────────────────────────────

lai_dt <- data.table::rbindlist(lai_diag_dt)
lai_dt <- lai_dt[is.finite(lai) & lai > 0 & lai < 20]
lai_dt[, site := factor(site, levels = sites)]
lai_dt[, type := factor(type, levels = c("Original", "Corrected"))]

p <- ggplot(lai_dt, aes(x = type, y = lai, fill = type)) +
  geom_violin(alpha = 0.6, linewidth = 0.3, colour = "grey40") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.4,
               fill = "white", colour = "grey20") +
  facet_wrap(~site, ncol = 3) +
  scale_fill_manual(values = c(Original = "#b2b2b2", Corrected = "#4dac26")) +
  labs(x = NULL, y = expression(LAI[ALS]~(m^2~m^{-2})),
       title = "LAI_ALS — scan angle correction (k = 0.6, nadir baseline)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

out_png <- file.path(out_figs, "lai_correction_summary.png")
ggsave(out_png, p, width = 9, height = 4, dpi = 150)
cat("\nFigure:", out_png, "\nDone.\n")
