# ---
# title:  18_scan_angle_correction.R
# desc:   Computes per-pixel scan angle correction rasters C(x,y) = mean(cos|θ|)
#         for Aigoual, Blois, Mormal. Processes LAS tiles one by one via manual
#         data.table aggregation (avoids pixel_metrics issues with SMB catalogs).
#         Applies correction multiplicatively to existing ladstack.tif.
#
# Outputs per site in output/intermediate/scan_angle/{site}/:
#   cos_theta_mean.tif       — correction raster C(x,y)
#   lai_als_corrected.tif    — ladstack sum × C
#
# Figures in output/figures/:
#   lai_correction_summary.png
#
# Run from NC_Full/:  source("scripts/18_scan_angle_correction.R")
# ---

library("here")
library("terra")
library("lidR")
library("data.table")
library("ggplot2")

source(here::here("R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

RES   <- 10L   # raster resolution (m), must match ladstack
sites <- c("Aigoual", "Blois", "Mormal")

uid <- system("id -u", intern = TRUE)
smb_las <- file.path(
  paste0("/run/user/", uid, "/gvfs"),
  "smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse",
  "root/_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR_Points_Clouds"
)

las_dirs <- setNames(file.path(smb_las, sites, "2-las_utm"), sites)

ladstack_paths <- setNames(
  file.path(paths$prosail_lidar, sites, "LiDAR",
            "PAD_Profiles_Classic", "ladstack.tif"),
  sites
)

out_root <- file.path(paths$output, "intermediate", "scan_angle")
out_figs <- file.path(paths$output, "figures")
for (d in c(file.path(out_root, sites), out_figs))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ── Helper: one LAS tile → cos_theta raster ───────────────────────────────────

process_tile <- function(las_path, res = 10L) {
  las <- tryCatch(readLAS(las_path, select = "a"), error = function(e) NULL)
  if (is.null(las) || nrow(las@data) == 0L) return(NULL)

  dt <- data.table::as.data.table(las@data)[, .(X, Y, ScanAngle)]
  dt[, cos_theta := cos(abs(ScanAngle) * pi / 180)]
  dt[, cx := (floor(X / res) + 0.5) * res]
  dt[, cy := (floor(Y / res) + 0.5) * res]

  agg <- dt[, .(
    cos_theta_mean = mean(cos_theta),
    angle_mean_deg = mean(abs(ScanAngle))
  ), by = .(cx, cy)]

  ext_r <- terra::ext(
    min(agg$cx) - res / 2, max(agg$cx) + res / 2,
    min(agg$cy) - res / 2, max(agg$cy) + res / 2
  )
  r <- terra::rast(ext_r, resolution = res, crs = "EPSG:32631", nlyrs = 2L)
  names(r) <- c("cos_theta_mean", "angle_mean_deg")

  cells <- terra::cellFromXY(r, as.matrix(agg[, .(cx, cy)]))
  ok    <- !is.na(cells)
  r[[1L]][cells[ok]] <- agg$cos_theta_mean[ok]
  r[[2L]][cells[ok]] <- agg$angle_mean_deg[ok]

  r
}

# ── Main loop ─────────────────────────────────────────────────────────────────

site_stats  <- list()
lai_diag_dt <- list()

for (site in sites) {

  cat("\n══ Site:", site, "══\n")
  t_site <- proc.time()

  out_dir   <- file.path(out_root, site)
  tmp_dir   <- file.path(out_dir, "tmp_tiles")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir)
  out_final <- file.path(out_dir, "cos_theta_mean.tif")

  if (file.exists(out_final)) {
    cat("  cos_theta_mean.tif exists — loading\n")
    cos_rast <- terra::rast(out_final)
  } else {
    las_files <- sort(list.files(las_dirs[[site]], full.names = TRUE))
    las_files <- las_files[grepl("LAS_[0-9]", basename(las_files))]
    n <- length(las_files)
    cat("  Tiles:", n, "\n")

    tile_rasts <- vector("list", n)
    for (i in seq_len(n)) {
      tmp_tif <- file.path(tmp_dir,
        paste0(tools::file_path_sans_ext(basename(las_files[[i]])), "_cos.tif"))

      if (file.exists(tmp_tif)) {
        tile_rasts[[i]] <- terra::rast(tmp_tif)
      } else {
        r <- process_tile(las_files[[i]], res = RES)
        if (!is.null(r)) {
          terra::writeRaster(r, tmp_tif, overwrite = TRUE,
                             gdal = "COMPRESS=LZW")
          tile_rasts[[i]] <- r
        }
      }
      if (i %% 20 == 0 || i == n)
        cat(sprintf("  %d / %d tiles (%.0f s)\n", i, n,
                    (proc.time() - t_site)[["elapsed"]]))
    }

    tile_rasts <- Filter(Negate(is.null), tile_rasts)
    cat("  Mosaicking", length(tile_rasts), "tiles...\n")
    cos_rast <- do.call(terra::mosaic, c(tile_rasts, list(fun = "mean")))
    terra::crs(cos_rast) <- "EPSG:32631"
    terra::writeRaster(cos_rast[["cos_theta_mean"]], out_final,
                       overwrite = TRUE, gdal = "COMPRESS=LZW")
    cat("  cos_theta_mean.tif written\n")
  }

  # ── Apply to ladstack ────────────────────────────────────────────────────────

  ladstack <- terra::rast(ladstack_paths[[site]])
  c_lyr    <- if (terra::nlyr(cos_rast) > 1L) cos_rast[["cos_theta_mean"]] else cos_rast
  c_lyr    <- terra::resample(c_lyr, ladstack[[1L]], method = "bilinear")

  lai_orig <- terra::sum(ladstack, na.rm = TRUE)
  lai_corr <- terra::sum(ladstack * c_lyr, na.rm = TRUE)

  terra::writeRaster(lai_corr,
                     file.path(out_dir, "lai_als_corrected.tif"),
                     overwrite = TRUE, gdal = "COMPRESS=LZW")

  # ── Stats ────────────────────────────────────────────────────────────────────

  c_v  <- terra::values(c_lyr,  na.rm = TRUE)
  o_v  <- terra::values(lai_orig, na.rm = TRUE)
  cr_v <- terra::values(lai_corr, na.rm = TRUE)

  elapsed <- round((proc.time() - t_site)[["elapsed"]] / 60, 1)
  cat(sprintf("  C mean: %.4f  (reduction: %.1f%%)\n",
              mean(c_v), (1 - mean(c_v)) * 100))
  cat(sprintf("  LAI_ALS: %.3f → %.3f  (Δ = %.1f%%)\n",
              mean(o_v, na.rm = TRUE), mean(cr_v, na.rm = TRUE),
              (mean(cr_v, na.rm = TRUE) / mean(o_v, na.rm = TRUE) - 1) * 100))
  cat(sprintf("  Elapsed: %.1f min\n", elapsed))

  site_stats[[site]] <- data.table(
    site          = site,
    n_pixels      = length(c_v),
    C_mean        = round(mean(c_v), 4),
    C_sd          = round(sd(c_v), 4),
    reduction_pct = round((1 - mean(c_v)) * 100, 2),
    lai_orig      = round(mean(o_v, na.rm = TRUE), 3),
    lai_corr      = round(mean(cr_v, na.rm = TRUE), 3)
  )

  lai_diag_dt[[site]] <- data.table(
    site = site,
    type = rep(c("Original", "Corrected"), c(length(o_v), length(cr_v))),
    lai  = c(o_v, cr_v)
  )
}

# ── Summary ───────────────────────────────────────────────────────────────────

sum_dt <- data.table::rbindlist(site_stats)
cat("\n── Summary ──\n")
print(sum_dt)
data.table::fwrite(sum_dt,
  file.path(out_root, "scan_angle_summary_all_sites.csv"))

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
       title = "LAI_ALS before/after scan angle correction (k = 0.6, θ corrected)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave(file.path(out_figs, "lai_correction_summary.png"),
       p, width = 9, height = 4, dpi = 150)
cat("\nFigure: lai_correction_summary.png\nDone.\n")
