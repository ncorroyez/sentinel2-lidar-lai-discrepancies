# ---
# title: "0_create_chm.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, 
# CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-05-24"
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
library("raster")
library("terra")

# --------------------------- Import useful functions --------------------------

source("../libraries/functions_chm.R")

# --------------------------------- Setup --------------------------------------
data_path <- "../../01_DATA"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Aigoual Blois Mormal
resolution <- 10 # 1
  
# -------------------------------- Process -------------------------------------
for (site in sites) {
  results_path <- file.path("../../03_RESULTS", site, "LiDAR")
  cat("Processing site:", site, "\n")
  create_chm(data_path, site, results_path, resolution)
}
