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

chm_charged <- rast("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/canopy/res_1m/rasterize_canopy.vrt")
dtm_charged <- rast("../../01_DATA/Mormal/00-LAS/3-las_normalize_utm/dtm/res_1m/rasterize_terrain.vrt")



