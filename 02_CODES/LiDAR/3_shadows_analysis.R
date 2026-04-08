# ---
# title: "1.shadows_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2023-05-28"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# ------------------------------ Libraries -------------------------------------

library("lidR")
library("data.table")
library("raster")
library("plotly")
library("terra")
library("viridis")
library("future")
library("rayshader")

# --------------------------- Import useful functions --------------------------

source("../libraries/functions_shadows_analysis.R")

# --------------------------------- Setup --------------------------------------
datapath <- "../../01_DATA"
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- "Mormal" # Aigoual Blois Mormal
# -------------------------------- Process -------------------------------------
for (site in sites) {
  cat("Processing site:", site, "\n")
  perform_shadows_analysis(datapath, site)
}
