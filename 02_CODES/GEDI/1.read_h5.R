# ---
# title: "1.read_h5.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-11-03"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Libraries
library(sf)
library(dplyr)
library(rGEDI)
library(ggplot2)
library(httr)
library(jsonlite)

# Setup
data_dir <- "../../01_DATA"
results_dir <- "../../03_RESULTS"
sites <- c("Aigoual", "Blois", "Mormal")
site <- c("Aigoual")

gedi_file <- file.path(data_dir, site, "GEDI",
                       "GEDI02_B_2021184062047_O14477_03_T10849_02_003_01_V002.h5")

gedilevel2b <- readLevel2B(gedi_file)

level2BVPM<-getLevel2BVPM(gedilevel2b)
level2BVPM_valid <- level2BVPM[
  !is.na(longitude_lastbin) & !is.na(latitude_lastbin) &
    longitude_lastbin != -9999 & latitude_lastbin != -9999
]

# Converting shot_number as "integer64" to "character"
level2BVPM_valid$shot_number<-as.character(level2BVPM_valid$shot_number)

# Converting GEDI Vegetation Profile Biophysical Variables as data.table to sf
level2BVPM_spdf<-sf::st_as_sf(
  level2BVPM_valid,
  coords = c("longitude_lastbin","latitude_lastbin"),
  crs = "epsg:4326"
)
# Exporting GEDI Vegetation Profile Biophysical Variables as ESRI Shapefile
sf::st_write(level2BVPM_spdf,
             file.path(outdir,"GEDI02_B_2021184062047_O14477_03_T10849_02_003_01_V002.h5_sub_VPM.shp"))

str(level2BVPM)












level2BPAIProfile<-getLevel2BPAIProfile(gedilevel2b)
head(level2BPAIProfile[,c("beam","shot_number","pai_z0_5m","pai_z5_10m")])
level2BPAIProfile$shot_number<-as.character(level2BPAIProfile$shot_number)

level2BPAIProfile_spdf <- sf::st_as_sf(
  level2BPAIProfile,
  coords = c("lon_lowestmode", "lat_lowestmode"),
  crs = "epsg:4326"
)

#specify GEDI beam
beam="BEAM0101"

# Plot Level2B PAI Profile
gPAIprofile<-plotPAIProfile(level2BPAIProfile, beam=beam, elev=TRUE)

l2b_df <- getLevel2BDF(gedilevel2b)
