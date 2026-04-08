# ================================================================
# GEDI extractor - Multi-site version (no manual file selection)
# ================================================================

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(rhdf5)
library(data.table)
library(sf)
library(stringr)
library(lubridate)
library(bit64)

# ---------------------
# Parameters
# ---------------------
data_dir <- "../../01_DATA"
sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Aigoual")
dz <- 5  # bin size for pavd_z profiles

# ---------------------
# Helper functions
# ---------------------

# h5file <- "/home/corroyez/Documents/NC_Full/01_DATA/Aigoual/GEDI/GEDI02_B_2021184062047_O14477_03_T10849_02_003_01_V002.h5"
# h5_contents <- h5ls(h5file, recursive = TRUE)
# beam_subset <- h5_contents[grepl("BEAM0000", h5_contents$group), ]
# write.csv(beam_subset, file = "beam_subset.csv", row.names = FALSE)

Hmax_padz <- function(padz, dz = 5) {
  if (length(padz) == 0 || all(is.na(padz))) return(0)
  padz <- as.numeric(padz)
  nz <- which(padz > 0)
  if (length(nz) == 0) return(0)
  last_idx <- max(nz)
  (length(padz) - last_idx) * dz
}

infos_beam <- function(h5file_path, beam_name, dz = 5) {
  fid <- H5Fopen(h5file_path)
  base <- paste0("/", beam_name)
  on.exit(H5Fclose(fid), add = TRUE)
  
  safe_read <- function(dataset, bit64conv = FALSE) {
    path <- paste0(base, dataset)
    if (H5Lexists(fid, path)) {
      if (bit64conv) {
        return(H5Dread(eval(parse(text=paste0("fid&'", path, "'"))), bit64conversion = "bit64"))
      } else {
        return(h5read(fid, path))
      }
    } else {
      return(NULL)
    }
  }
  
  algorithmrun_flag <- safe_read("/algorithmrun_flag")
  degrade_flag      <- safe_read("/geolocation/degrade_flag")
  elev_bin0         <- safe_read("/geolocation/elevation_bin0")
  elev_dem          <- safe_read("/geolocation/digital_elevation_model")
  shot_number       <- safe_read("/geolocation/shot_number", bit64conv = TRUE)
  delta_time        <- safe_read("/geolocation/delta_time")
  lat_lowestmode    <- safe_read("/geolocation/lat_lowestmode")
  lon_lowestmode    <- safe_read("/geolocation/lon_lowestmode")
  rh100             <- safe_read("/rh100")
  pai               <- safe_read("/pai")
  cover             <- safe_read("/cover")
  pavd_z            <- safe_read("/pavd_z")

  if (is.null(algorithmrun_flag) || is.null(elev_bin0) ||
      is.null(elev_dem) || is.null(shot_number)) return(NULL)
  
  diff_alt0 <- elev_bin0 - elev_dem
  subset_idx <- which(
    algorithmrun_flag == 1 &
      diff_alt0 < 100 & diff_alt0 > 0 &
      (if (!is.null(degrade_flag)) degrade_flag == 0 else TRUE)
  )
  if (length(subset_idx) == 0) return(NULL)
  
  dt <- data.table(
    shot_number = as.character(shot_number[subset_idx]),
    delta_time  = as.numeric(delta_time[subset_idx]),
    lat         = as.numeric(lat_lowestmode[subset_idx]),
    lon         = as.numeric(lon_lowestmode[subset_idx]),
    RH100       = as.numeric(rh100[subset_idx]) / 100,
    pai         = as.numeric(pai[subset_idx]),
    cover       = as.numeric(cover[subset_idx]),
    beam        = beam_name
  )
  
  if (!is.null(pavd_z) && is.matrix(pavd_z)) {
    if (ncol(pavd_z) >= max(subset_idx)) {
      Hmax_vals <- apply(pavd_z[, subset_idx, drop = FALSE], 2, Hmax_padz, dz = dz)
      dt[, Hmax := as.numeric(Hmax_vals)]
    } else {
      dt[, Hmax := NA_real_]
    }
  } else {
    dt[, Hmax := NA_real_]
  }
  
  dt
}

infos_allbeams <- function(h5file_path, dz = 5) {
  groups <- h5ls(h5file_path, recursive = FALSE)$name
  groups <- groups[groups != "metadata"]
  result_list <- list()
  
  for (beam in groups) {
    dt_beam <- tryCatch(
      infos_beam(h5file_path, beam, dz = dz),
      error = function(e) {
        message("Error reading beam ", beam, " in ", basename(h5file_path), ": ", e$message)
        NULL
      }
    )
    if (!is.null(dt_beam) && nrow(dt_beam) > 0) result_list[[beam]] <- dt_beam
  }
  
  if (length(result_list) == 0) return(NULL)
  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# ================================================================
# MAIN LOOP OVER SITES
# ================================================================

for (site in sites) {
  message("\n==============================")
  message("Processing site: ", site)
  message("==============================")
  
  gedi_dir <- file.path(data_dir, site, "GEDI")
  shp_path <- file.path(data_dir, site, "Geo_Files", "utm_init.shp")
  
  if (!dir.exists(gedi_dir)) {
    message("GEDI directory not found for ", site, " — skipping.")
    next
  }
  if (!file.exists(shp_path)) {
    message("Shapefile not found for ", site, " — skipping.")
    next
  }
  
  h5_files <- list.files(gedi_dir, pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) == 0) {
    message("No GEDI .h5 files found in ", gedi_dir)
    next
  }
  
  all_infos <- list()
  for (i in seq_along(h5_files)) {
    f <- h5_files[i]
    message(sprintf("Processing file %d / %d : %s", i, length(h5_files), basename(f)))
    
    ok <- tryCatch({
      h5ls(f, recursive = FALSE)
      TRUE
    }, error = function(e) {
      message(sprintf("Cannot open %s : %s", f, e$message))
      FALSE
    })
    if (!ok) next
    
    dt_file <- infos_allbeams(f, dz = dz)
    if (is.null(dt_file)) next
    stop()
    fname  <- basename(f)
    parts  <- str_split(fname, "_")[[1]]
    orbit  <- if (length(parts) >= 4) substr(parts[4], 2, 6) else NA_character_
    dt_file[, orbit := orbit]
    
    all_infos[[basename(f)]] <- dt_file
  }
  
  if (length(all_infos) == 0) {
    message("No valid GEDI data extracted for ", site)
    next
  }
  
  df <- rbindlist(all_infos, use.names = TRUE, fill = TRUE)
  
  # Add time fields
  gedi_epoch <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
  df[, time_utc := gedi_epoch + as.numeric(delta_time)]
  df[, year  := format(time_utc, "%Y")]
  df[, month := format(time_utc, "%m")]
  df[, day   := format(time_utc, "%d")]
  df[, heure := format(time_utc, "%H:%M:%S")]
  
  # Convert to sf and project
  df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  df_sf_utm <- st_transform(df_sf, crs = 32631)
  
  # Filter by shapefile
  shp <- st_read(shp_path, quiet = TRUE)
  shp <- st_transform(shp, crs = 32631)
  df_sf_clipped <- st_intersection(df_sf_utm, shp)
  
  if (nrow(df_sf_clipped) == 0) {
    message("No GEDI shots within shapefile for ", site)
    next
  }
  
  num_cols <- c("lat", "lon", "RH100", "pai", "cover", "Hmax")
  num_cols <- intersect(num_cols, colnames(df_sf_clipped))
  
  df_sf_clipped <- df_sf_clipped %>%
    dplyr::mutate(across(all_of(num_cols), ~ round(.x, 2)))
  
  # Also round coordinates in geometry (optional)
  coords_utm <- st_coordinates(df_sf_clipped)
  df_sf_clipped$X_utm <- round(coords_utm[,1], 2)
  df_sf_clipped$Y_utm <- round(coords_utm[,2], 2)
  
  # Export
  out_base <- file.path(data_dir, site, "Geo_Files")
  if (!dir.exists(out_base)) dir.create(out_base)
  
  gpkg_path <- file.path(out_base, paste0("gedi_shots_", site, ".gpkg"))
  st_write(df_sf_clipped, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
  message("Saved GeoPackage: ", gpkg_path)
  
  coords_utm <- st_coordinates(df_sf_clipped)
  df_out <- cbind(as.data.table(df_sf_clipped)[, !"geometry"],
                  X_utm = coords_utm[,1], Y_utm = coords_utm[,2])
  csv_path <- file.path(out_base, paste0("gedi_shots_", site, ".csv"))
  fwrite(df_out, csv_path)
  message("Saved CSV: ", csv_path)
}

message("\nAll sites processed successfully.")
