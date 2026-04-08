# ---
# title: "0.download_gedi_data.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-09-23"
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
sites <- c("Aigoual")

# rGEDI::gediSetLogin(
#   usr = "likho1404",
#   pwd = "P4?Et##@$CYppQTc"
# )

# Define the base URL for the GEDI API
# base_url <- "https://data.science.lu/gedi"

for (site in sites){
  aoi <- st_read(file.path(data_dir, site, "Geo_Files", "utm_init.shp"))
  aoi_wgs <- st_transform(aoi, crs = 4326)
  bbox_wgs <- st_bbox(aoi_wgs)
  
  footprints <- gedifinder(
    ul_lat = bbox_wgs["ymax"],
    ul_lon = bbox_wgs["xmin"],
    lr_lat = bbox_wgs["ymin"],
    lr_lon = bbox_wgs["xmax"],
    version = "002",
    product = "GEDI02_B"
  )
  
  # Study area boundary box coordinates
  ul_lat <- bbox_wgs["ymax"]
  lr_lat <- bbox_wgs["ymin"]
  ul_lon <- bbox_wgs["xmin"]
  lr_lon <- bbox_wgs["xmax"]
  # ul_lat <- 50
  # lr_lat <- 0
  # ul_lon <- 0
  # lr_lon <- 50
  
  # Specifying the date range
  daterange=c("2021-05-01T00:00:00Z","2022-12-31T23:59:59Z")
  daterange <- c("2019-04-18", "2023-03-14")
  
  # Get path to GEDI data
  gLevel1B <- rGEDI::gedifinder(product = "GEDI01_B",
                         ul_lat, ul_lon, lr_lat, lr_lon,
                         version = "002", daterange = daterange)
  gLevel2A <- rGEDI::gedifinder(product = "GEDI02_A",
                         ul_lat, ul_lon, lr_lat, lr_lon,
                         version = "002", daterange = daterange)
  gLevel2B <- rGEDI::gedifinder(product = "GEDI02_B",
                         ul_lat, ul_lon, lr_lat, lr_lon,
                         version = "002", daterange = daterange)
  # gLevel3 <- rGEDI::gedifinder(product = "GEDI03",
  #                               ul_lat, ul_lon, lr_lat, lr_lon,
  #                               version = "002", daterange = daterange)
  stop()
  outdir <- file.path(data_dir, site, "GEDI")
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Create .netrc file
  # netrc = file.path(outdir, ".netrc")
  netrc <- "~/.netrc"
  # netrc_conn <- file(netrc)
  
  writeLines(c("machine urs.earthdata.nasa.gov",
               sprintf("login %s", Sys.getenv("EARTHDATA_USER")),
               sprintf("password %s", Sys.getenv("EARTHDATA_PASS"))
  ), netrc)
  Sys.chmod(netrc, mode = "600")
  
  # close(netrc_conn)
  
  gediDownload(filepath = gLevel2B, outdir = outdir)
}

gedilevel1b <- readLevel1B(level1Bpath = file.path(outdir,
                                                   "GEDI01_B_2019128004537_O02270_02_T04863_02_005_01_V002.h5"))
