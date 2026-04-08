# ---
# title: "0_process_parquet.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-12-08"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(arrow)
library(sf)
library(dplyr)
library(purrr)
library(stringr)
library(readr)

input_folder  <- "/media/corroyez/MyPassport/01_DATA/outputs_meso/france"
output_folder <- "../../01_DATA"
shapefile_paths <- c(
  "Aigoual" = "/media/corroyez/MyPassport/01_DATA/Aigoual/Geo_Files/l93_init.shp",
  "Blois"   = "/media/corroyez/MyPassport/01_DATA/Blois/Geo_Files/l93_init.shp",
  "Mormal"  = "/media/corroyez/MyPassport/01_DATA/Mormal/Geo_Files/l93_init.shp"
)

# GEDI Column names
lat_col <- "lat"
lon_col <- "lon"

if (!dir.exists(output_folder)) dir.create(output_folder)
message("Loading shapefiles...")
shps <- lapply(shapefile_paths, function(x) {
  s <- st_read(x, quiet = TRUE)
  if (!all(st_is_valid(s))) {
    message(paste("Found invalid geometry in:", 
                  basename(x), "- Repairing automatically..."))
    s <- st_make_valid(s)
  }
  st_transform(s, 4326)
})

# Create a list to store matching dataframes for each zone
results_list <- list(
  Aigoual = list(),
  Blois = list(),
  Mormal = list()
)

# Get list of parquet files
files <- list.files(input_folder, pattern = "\\.parquet$", full.names = TRUE)
message(paste("Found", length(files), "files to process."))

# Loop through GEDI files *once* to minimize I/O overhead
for (f in files) {
  
  # Read the file
  tryCatch({
    df <- read_parquet(f)
    
    # Check bounds for EACH shapefile
    for (zone_name in names(shps)) {
      
      zone_sf <- shps[[zone_name]]
      bbox    <- st_bbox(zone_sf)
      
      # Check if the GEDI file is completely disjoint from the zone
      file_min_lat <- min(df[[lat_col]], na.rm = TRUE)
      file_max_lat <- max(df[[lat_col]], na.rm = TRUE)
      file_min_lon <- min(df[[lon_col]], na.rm = TRUE)
      file_max_lon <- max(df[[lon_col]], na.rm = TRUE)
      
      if (file_max_lat < bbox["ymin"] || file_min_lat > bbox["ymax"] ||
          file_max_lon < bbox["xmin"] || file_min_lon > bbox["xmax"]) {
        next
      }
      
      # Keep rows roughly inside the bounding box
      df_rough <- df %>%
        filter(
          .data[[lat_col]] >= bbox["ymin"], .data[[lat_col]] <= bbox["ymax"],
          .data[[lon_col]] >= bbox["xmin"], .data[[lon_col]] <= bbox["xmax"]
        )
      
      if (nrow(df_rough) == 0) next
      
      # Convert to sf object
      df_sf <- st_as_sf(df_rough, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
      
      # Spatial intersection (keeps only points strictly inside the polygon)
      df_precise <- st_filter(df_sf, zone_sf)
      
      # If we have matches, store them in our results list
      if (nrow(df_precise) > 0) {
        # Convert back to regular dataframe (drop geometry column for speed/storage)
        df_clean <- df_precise %>% st_drop_geometry()
        results_list[[zone_name]][[length(results_list[[zone_name]]) + 1]] <- df_clean
      }
    }
    
    if (which(files == f) %% 100 == 0) {
      message(paste("Processed", which(files == f), "files..."))
    }
    
  }, error = function(e) {
    message(paste("Error reading file:", f, "\n", e))
  })
}

message("Processing complete. Merging and saving final files...")
for (zone_name in names(results_list)) {
  data_chunks <- results_list[[zone_name]]
  
  if (length(data_chunks) > 0) {
    # Merge all chunks for this zone
    final_df <- bind_rows(data_chunks)
    
    # Convert shot_number to character (Safe for IDs)
    if ("shot_number" %in% names(final_df)) {
      final_df$shot_number <- as.character(final_df$shot_number)
    }
    
    # Coordinate Transformation (Lambert 93 -> UTM 31N)
    if (all(c("x", "y") %in% names(final_df))) {
      message(paste("Transforming coordinates for", zone_name, "..."))
      
      # Create a temporary spatial object from the existing Lambert 93 columns
      # EPSG:2154 = RGF93 / Lambert-93
      temp_sf <- st_as_sf(final_df, coords = c("x", "y"), crs = 2154, remove = FALSE)
      
      # Transform to UTM Zone 31N
      # EPSG:32631 = WGS 84 / UTM zone 31N
      temp_sf_utm <- st_transform(temp_sf, 32631)
      
      # Extract the new coordinates
      new_coords <- st_coordinates(temp_sf_utm)
      
      # Assign them to new columns
      final_df$x_corrected <- new_coords[, "X"]
      final_df$y_corrected <- new_coords[, "Y"]
      
    } else {
      warning(paste("Skipping coordinate transform for", zone_name, "- 'x' or 'y' columns missing."))
    }
    
    # Save to Parquet
    out_file <- file.path(output_folder, paste0(zone_name, "_GEDI_Combined.parquet"))
    write_parquet(final_df, out_file)
    
    message(paste("Saved:", out_file, "with", nrow(final_df), "rows."))
    
  } else {
    message(paste("No data found for", zone_name))
  }
}
# stop()
# 
# a <- read_parquet("../../01_DATA/Blois_GEDI_Combined.parquet")
# gpkg <- st_read("../../01_DATA/Blois/Geo_Files/gedi_shots_Blois.gpkg")
# joined <- inner_join(a, gpkg, by = "shot_number") %>%
#   select(-ends_with(".y")) %>%
#   rename_with(~ str_remove(., "\\.x$"), ends_with(".x"))


sites <- c("Aigoual", "Blois", "Mormal")
for (site in sites) {
  
  message(paste("Processing site:", site, "..."))
  path_parquet <- file.path("../../01_DATA", paste0(site, "_GEDI_Combined.parquet"))
  path_gpkg_in <- file.path("../../01_DATA", site, "Geo_Files", paste0("gedi_shots_", site, ".gpkg"))
  path_gpkg_out <- file.path("../../01_DATA", site, "Geo_Files", paste0("gedi_joined_", site, ".gpkg"))
  path_csv_out  <- file.path("../../01_DATA", site, "Geo_Files", paste0("gedi_joined_", site, ".csv"))
  
  if (file.exists(path_parquet) && file.exists(path_gpkg_in)) {
    parquet <- read_parquet(path_parquet)
    gpkg <- st_read(path_gpkg_in, quiet = TRUE)
    
    # Variable names handling
    new_names <- tolower(names(parquet))
    new_names <- gsub("\\.", "_", new_names)
    names(parquet) <- make.unique(new_names, sep = "_")
    new_names <- tolower(names(gpkg))
    new_names <- gsub("\\.", "_", new_names)
    names(gpkg) <- make.unique(new_names, sep = "_")
    
    # Join and Clean
    joined <- inner_join(parquet, gpkg, by = "shot_number") %>%
      select(-ends_with(".y")) %>%
      rename_with(~ str_remove(., "\\.x$"), ends_with(".x"))
    
    # Convert back to sf object (Critical for saving as GPKG)
    joined <- st_as_sf(joined)
    st_write(joined, path_gpkg_out, delete_dsn = TRUE, quiet = TRUE)
    write_csv(joined, path_csv_out)
    message(paste("Saved GPKG & CSV for:", site))
    
  } else {
    warning(paste("Skipping", site, "- Input files not found."))
  }
}
message("Processing complete.")
