# ---
# title: "0_Very_Simple_LiDAR_LAI_Workflow"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "23-10-2023"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("data.table")
library("raster")
library("rgdal")
library("plotly")
library("terra")
library("viridis")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

datapath = "../../01_DATA"

las_dir <- file.path(datapath, "Mormal/00-LAS/3-las_normalize_utm")
cat("lasDpath: ", las_dir, "\n")
ctg <- readLAScatalog(las_dir)
plot(ctg)

opt_chunk_size(ctg) <- 0 # Processing by files
opt_chunk_buffer(ctg) <- 10

set_lidr_threads(0)
get_lidr_threads()

n = 5
ctg_2 = head(ctg, n = n)

# plan(multisession)

# DTM
opt_output_files(ctg) <- paste0("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/dtm/res_1m/{XLEFT}_{YBOTTOM}", overwrite = TRUE)
dtm_kriging <- rasterize_terrain(ctg, res = 1, algorithm = kriging(k = 40)) # overwrite does not seem to work
plot_dtm3d(dtm_kriging, bg = "white") 

# Canopy
opt_output_files(ctg) <- paste0("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/canopy/res_1m/{XLEFT}_{YBOTTOM}", overwrite = TRUE)
chm <- rasterize_canopy(ctg, res = 1, p2r(0.2, na.fill = tin())) # overwrite does not seem to work
plot(chm)

# MNE

mne <- chm - dtm_kriging

chm_charged <- rast("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/canopy/res_1m/rasterize_canopy.vrt")
dtm_charged <- rast("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/dtm/res_1m/rasterize_terrain.vrt")

mne_charged <- chm_charged - dtm_charged

min(mne)
min(mne_charged)

mask_sol <- mne > 2
plot(mask_sol)
plot(mne_charged)

mne_file <- file.path(mne_file, basename('res_1m'))
mne_file <- file.path(las_dir, basename('mne'))
dir.create(path = mne_file,showWarnings = FALSE,recursive = TRUE)
writeRaster(mne_charged, filename = file.path(mne_file, basename("mne.tif")), overwrite=TRUE)
writeRaster(mask_sol, filename = file.path(mne_file, basename("mne_mask_2m.tif")), overwrite=TRUE)

myPAI <- function(z, zmin, k=0.5){
  Nout = sum(z<zmin)
  Nin = length(z)
  PAI = -log(Nout/Nin)/k
}

n = 5
ctg_2 = head(ctg, n = n)

lidR::pixel_metrics(las, ~myPAI(Z, zmin=2))
pai_las <- lidR::pixel_metrics(ctg, ~myPAI(Z, zmin=2))
plot(pai_las, col = height.colors(50))

# Delete inf values
pai_las_cleaned <- pai_las
pai_las_cleaned[pai_las_cleaned == Inf] <- NA

# Get min and max values 
palette_min <- min(values(pai_las_cleaned), na.rm = TRUE)
palette_max <- max(values(pai_las_cleaned), na.rm = TRUE)

# Create viridis colorbar with a fix number of colors (e.g. 50)
num_colors <- 50
custom_palette <- viridis(num_colors)

# Plot the graph with the viridis colorbar
plot(pai_las_cleaned, col = custom_palette)

# Add a legend
legend("topright", legend = c(palette_min, palette_max), fill = custom_palette)

output_file <- "/home/corroyez/Documents/NC_Full/03_RESULTS/pai_las_cleaned.tif"

# Export raster to GeoTIFF
writeRaster(pai_las_cleaned, filename = output_file, overwrite = TRUE)

# Mask 2m + LAI result
lai <- rast("../../03_RESULTS/pai_las_cleaned.tif")
mask2m <- rast("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/mne/res_10m/mne_mask_2m.tif")

final <- mask(lai, mask2m)
plot(final)



