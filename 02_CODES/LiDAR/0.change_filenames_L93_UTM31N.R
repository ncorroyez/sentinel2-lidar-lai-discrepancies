# ---
# title: "0.change_filenames_L93_UTM31N.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-12-02"
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

library(sf)
library(dplyr)
library(stringr)

las_l93_dir <- "/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/1-las_l93"
las_utm_dir <- "/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/3-las_normalized_utm"
execute_action <- TRUE 

# Ensure output directory exists
if(!dir.exists(las_utm_dir)) dir.create(las_utm_dir, recursive = TRUE)

# 2. CREATE MAPPING DATAFRAME
# ---------------------------------------------------------
files_l93 <- list.files(las_l93_dir, pattern = "^LAS_.*\\.las$", full.names = TRUE)

file_map <- data.frame(old_path = files_l93, stringsAsFactors = FALSE) %>%
  mutate(
    old_name = basename(old_path),
    # Extract coordinates from "LAS_737500_6338000.las"
    coords_raw = str_extract_all(old_name, "\\d{6,7}")
  ) %>%
  rowwise() %>%
  mutate(
    # Assign L93 coordinates
    X_L93 = as.numeric(coords_raw[1]),
    Y_L93 = as.numeric(coords_raw[2])
  ) %>%
  ungroup()

# 3. SPATIAL TRANSFORMATION
# ---------------------------------------------------------
# Convert to sf object (Lambert-93 / EPSG:2154)
points_l93 <- st_as_sf(file_map, coords = c("X_L93", "Y_L93"), crs = 2154)

# Transform to UTM 31N (EPSG:32631)
points_utm <- st_transform(points_l93, crs = 32631)

# Extract new coordinates
coords_utm <- st_coordinates(points_utm)

file_map$X_UTM <- floor(coords_utm[,1])
file_map$Y_UTM <- floor(coords_utm[,2])

# --- MODIFIED SECTION ---
file_map <- file_map %>%
  mutate(
    new_name = paste0("LAS_", X_UTM, "_", Y_UTM, ".las"),
    # HERE IS THE CHANGE: Use las_utm_dir instead of dirname(old_path)
    new_path = file.path(las_utm_dir, new_name)
  )

# 4. REVIEW
# ---------------------------------------------------------
print(head(file_map[, c("old_name", "new_name", "new_path")]))

# 5. EXECUTE COPY
# ---------------------------------------------------------
if (execute_action) {
  message("Renaming ", nrow(file_map), " files to UTM directory...")
  success <- file.rename(from = file_map$old_path, to = file_map$new_path)
  
  if(all(success)) {
    message("Success! All files copied to UTM folder with new coordinates.")
  } else {
    warning(paste(sum(!success), "files failed to copy."))
  }
} else {
  message("Dry run complete. Set 'execute_action <- TRUE' to actually copy files.")
}