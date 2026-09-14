# ---
# download_s2_fullyear.R — download Sentinel-2 L2A for the 3 sites at the
# user-specified cloud-free monthly dates (2021), to extend the LAI time series
# beyond the summer plateau (enables real leaf-out / senescence phenology).
#
# Output: data/results/<site>/PROSAIL_Optimization (revision data folder,
#         gitignored — NOT 01_DATA / 03_RESULTS).
# Credentials: CDSE_ID / CDSE_SECRET from ~/.Renviron (never the 01_DATA yml).
#
# Env vars:
#   S2_DRY_ONE=1   -> download only the first (site,date) pair (timing/sanity)
#   S2_SITES=Blois -> restrict to a comma-separated subset of sites
# ---
library("here"); library("preprocS2")
source(here::here("R", "paths.R"))
source(here::here("R", "sentinel2_preprocessing.R"))

# cloud-free monthly dates per site (NA months omitted)
dates <- list(
  Aigoual = c("2021-01-07","2021-02-26","2021-03-28","2021-04-07","2021-05-27",
              "2021-06-26","2021-07-11","2021-08-15","2021-09-24","2021-10-14",
              "2021-11-18","2021-12-18"),
  Blois   = c("2021-02-24","2021-03-24","2021-04-23","2021-05-28","2021-06-14",
              "2021-07-19","2021-08-26","2021-09-22","2021-11-09","2021-12-21"),
  Mormal  = c("2021-02-24","2021-03-31","2021-04-25","2021-05-30","2021-06-14",
              "2021-07-19","2021-08-08","2021-09-07","2021-10-22","2021-12-21")
)

sites_env <- Sys.getenv("S2_SITES", "")
if (nzchar(sites_env)) dates <- dates[strsplit(sites_env, ",")[[1]]]
dry_one <- nzchar(Sys.getenv("S2_DRY_ONE", ""))
mode    <- Sys.getenv("S2_MODE", "reflectance")   # "reflectance" | "geometry"

name_vect <- "utm_init.shp"
tiling_grid_kml <- here::here("tileS2_kml", "Sentinel-2_tiling_grid.kml")
stopifnot(file.exists(tiling_grid_kml))
jobs <- do.call(rbind, lapply(names(dates), function(s)
  data.frame(site = s, date = dates[[s]], stringsAsFactors = FALSE)))
if (dry_one) jobs <- jobs[1, , drop = FALSE]
cat(sprintf("[download] mode=%s | %d (site,date) jobs%s\n", mode, nrow(jobs),
            if (dry_one) " [DRY: first only]" else ""))

max_try <- as.integer(Sys.getenv("S2_RETRY", "4"))   # MPC STAC is flaky (text/plain)
# Output: per-date folders under <root>/<site>/Sentinel-2/<date>/, matching the
# existing 01_DATA archive layout. Root = external HDD 01_DATA if mounted, else
# internal 01_DATA. Override with S2_OUTROOT.
hdd_root <- "/media/corroyez/MyPassport/01_DATA"
default_root <- if (dir.exists(hdd_root)) hdd_root else here::here("01_DATA")
out_root <- Sys.getenv("S2_OUTROOT", default_root)
cat(sprintf("[download] output root (per-date): %s/<site>/Sentinel-2/<date>/\n", out_root))
already_done <- function(output_dir, s, d) {
  # reflectance: raster_samples/<site>_*_<date>.tiff ; geometry: geom_acq_S2/sza_<date>.tiff
  pat <- if (mode == "geometry") sprintf("sza_%s\\.tif", d) else sprintf("_%s\\.tiff$", d)
  sub <- if (mode == "geometry") "geom_acq_S2" else "raster_samples"
  length(list.files(file.path(output_dir, sub), pattern = pat)) > 0
}
for (i in seq_len(nrow(jobs))) {
  s <- jobs$site[i]; d <- jobs$date[i]
  aoi_path   <- file.path(paths$raw_data, s, "Geo_Files", name_vect)
  output_dir <- file.path(out_root, s, "Sentinel-2", d)   # per-date folder
  if (already_done(output_dir, s, d)) {
    cat(sprintf("[%d/%d] %s %s -> SKIP (exists)\n", i, nrow(jobs), s, d)); next
  }
  cat(sprintf("[%d/%d] %s %s\n", i, nrow(jobs), s, d))
  ok <- FALSE
  for (attempt in seq_len(max_try)) {
    t0 <- Sys.time()
    ok <- tryCatch({ download_s2_cdse(site = s, date = d, aoi_path = aoi_path,
                                      output_dir = output_dir, mode = mode,
                                      tiling_grid = tiling_grid_kml); TRUE },
                   error = function(e) { cat("  try", attempt, "ERROR:", conditionMessage(e), "\n"); FALSE })
    cat(sprintf("  try %d %s (%.1f min)\n", attempt, if (ok) "done" else "fail",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    if (ok) break
    Sys.sleep(10 * attempt)   # backoff before retry
  }
  cat(sprintf("  => %s\n", if (ok) "DONE" else "GAVE_UP"))
  Sys.sleep(3)                # politeness between dates
}
cat("DOWNLOAD_DONE\n")
