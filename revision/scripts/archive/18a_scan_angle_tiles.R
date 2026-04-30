# ---
# title:  18a_scan_angle_tiles.R
# desc:   Processes LAS tiles for ONE site and writes per-tile cos_theta tifs.
#         Exits cleanly after TIME_LIMIT seconds so foreground runs stay within
#         the 10-min bash timeout. Re-run until all tiles are done (file.exists
#         checkpoint skips already-processed tiles).
#
#         Usage:  Rscript 18a_scan_angle_tiles.R <site>
#         Exit codes:
#           0 — all tiles done
#           2 — time limit reached, re-run needed
# ---

library("here")
library("terra")
library("lidR")
library("data.table")

TIME_LIMIT <- 520L   # seconds; leave buffer before 10-min bash timeout

args <- commandArgs(trailingOnly = TRUE)
site <- if (length(args) >= 1L) args[[1L]] else "Aigoual"
stopifnot(site %in% c("Aigoual", "Blois", "Mormal"))

RES <- 10L

uid     <- system("id -u", intern = TRUE)
las_dir <- file.path(
  paste0("/run/user/", uid, "/gvfs"),
  "smb-share:server=pnas3.stockage.inrae.fr,share=mo-mtd-pulse",
  "root/_PROJETS/2023_2026_These_Nathan_Corroyez/LiDAR_Points_Clouds",
  site, "2-las_utm"
)
if (!dir.exists(las_dir)) stop("LAS dir not found: ", las_dir)

tmp_dir <- file.path(paths$output, "intermediate",
                       "scan_angle", site, "tmp_tiles")
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

las_files <- sort(list.files(las_dir, full.names = TRUE))
las_files <- las_files[grepl("LAS_[0-9]", basename(las_files))]
n <- length(las_files)
cat("Site:", site, "| Tiles:", n, "\n")

process_tile <- function(f, res = 10L) {
  las <- tryCatch(readLAS(f, select = "a"), error = function(e) NULL)
  if (is.null(las) || nrow(las@data) == 0L) return(NULL)
  dt <- data.table::as.data.table(las@data)[, .(X, Y, ScanAngle)]
  dt[, cos_theta := cos(abs(ScanAngle) * pi / 180)]
  dt[, cx := (floor(X / res) + 0.5) * res]
  dt[, cy := (floor(Y / res) + 0.5) * res]
  agg <- dt[, .(cos_theta_mean = mean(cos_theta)), by = .(cx, cy)]
  ext_r <- terra::ext(
    min(agg$cx) - res / 2, max(agg$cx) + res / 2,
    min(agg$cy) - res / 2, max(agg$cy) + res / 2
  )
  r <- terra::rast(ext_r, resolution = res, crs = "EPSG:32631", nlyrs = 1L)
  names(r) <- "cos_theta_mean"
  cells <- terra::cellFromXY(r, as.matrix(agg[, .(cx, cy)]))
  ok    <- !is.na(cells)
  r[[1L]][cells[ok]] <- agg$cos_theta_mean[ok]
  r
}

t0         <- proc.time()
n_done     <- 0L
time_limit <- FALSE

for (i in seq_len(n)) {
  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (elapsed > TIME_LIMIT) {
    time_limit <- TRUE
    cat("Time limit reached at tile", i, "—",
        n_done, "/", n, "done. Re-run needed.\n")
    break
  }

  tmp_tif <- file.path(tmp_dir,
    paste0(tools::file_path_sans_ext(basename(las_files[[i]])), "_cos.tif"))

  if (file.exists(tmp_tif)) {
    n_done <- n_done + 1L
    next
  }

  r <- process_tile(las_files[[i]], res = RES)
  if (!is.null(r)) {
    terra::writeRaster(r, tmp_tif, overwrite = TRUE, gdal = "COMPRESS=LZW")
    system("sync")
    n_done <- n_done + 1L
  }

  if (i %% 20L == 0L)
    cat(sprintf("  %d / %d  (%.0f s)\n", i, n, elapsed))
}

cat(sprintf("Done: %d / %d tiles  (%.1f min)\n",
            n_done, n, (proc.time() - t0)[["elapsed"]] / 60))

quit(status = if (time_limit) 2L else 0L)
