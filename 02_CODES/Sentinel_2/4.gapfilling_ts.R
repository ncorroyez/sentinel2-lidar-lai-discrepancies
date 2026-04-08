# ---
# title: "4.gapfilling_ts.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-02-13"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
library(terra)
library(zoo)
library(readODS)
library(stringr)
library(signal)
library(ptw)
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_general_tools.R")

# Pre-processing Parameters
data_dir <- '/media/corroyez/MyPassport/01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
results_dir <- '../../03_RESULTS'

# Read S2 Acquisition Dates (1/month)
data <- read_ods(file.path(data_dir, "S2_Dates_Summer2021.ods"))
common_dates <- seq(from = as.Date("2021-06-01"), 
                    to   = as.Date("2021-09-30"), 
                    by   = "10 days")

cat("Grille commune définie :", length(common_dates), "dates.\n")
print(common_dates)

distribs <- c("atbd", "atbd_optim_common")


for (site in sites) {
  print(site)
  res_path <- file.path(results_dir, site, "Metrics/Deciduous_Only")
  dateAcqs <- get_site_full_dates(site, data)
  # stop()
  for (distrib in distribs) {
    print(distrib)
    files <- list.files(res_path, 
                        pattern = paste0(
                          "^s2lai_(",
                          paste(dateAcqs, collapse = "|"),
                          ")_",
                          distrib,
                          "_res_10_m\\.tif$"
                        ), 
                        full.names = TRUE)
    r_stack <- rast(files)
    terra::time(r_stack) <- dateAcqs
    cat("Interpolating missing values (Gap-Filling)...\n")
    r_filled <- terra::approximate(r_stack, method = "linear", rule = 2)
    
    smooth_fun <- function(y) {
      if (all(is.na(y))) return(y)
      return(rollmean(y, k = 3, fill = "extend", align = "center"))
    }
    
    sg_fun <- function(y) {
      if (all(is.na(y))) return(y)
      return(signal::sgolayfilt(y, p = 2, n = 5)) 
    }
    
    whittaker_fun <- function(y) {
      if (all(is.na(y))) return(y)
      return(whit2(y, lambda = 10)) 
    }
    
    cat("Applying smoothing filter (Rolling Mean)...\n")
    # r_smoothed <- app(r_filled, fun = smooth_fun, cores = 1)
    # r_smoothed <- app(r_filled, fun = sg_fun, cores = 1)
    r_smoothed <- app(r_filled, fun = whittaker_fun, cores = 1)
    names(r_smoothed) <- paste0("LAI_Smooth_", dateAcqs)
    
    out_dir <- dirname(files[1])
    out_name <- file.path(out_dir, paste0(site, 
                                          "_LAI_TimeSeries_SmoothedW_",
                                          distrib,
                                          ".tif"))
    writeRaster(r_smoothed, out_name, overwrite = TRUE)
    
    # Harmonization (Resampling to Common Grid)
    # --------------------------------------------
    cat("     Harmonizing to common grid...\n")
    
    # Define function to project local dates to common_dates
    resample_fun <- function(y) {
      if (all(is.na(y))) return(rep(NA, length(common_dates)))
      # Interpolate from current_dates (x) to common_dates (xout)
      return(approx(x = as.numeric(dateAcqs), 
                    y = y, 
                    xout = as.numeric(common_dates), 
                    rule = 2
      )$y
      )
    }
    
    r_harmonized <- app(r_smoothed, fun = resample_fun, cores = 1)
    names(r_harmonized) <- paste0("LAI_", common_dates)
    terra::time(r_harmonized) <- common_dates
    
    # Save Harmonized Time Series (Ready for ML)
    out_name_harmonized <- file.path(res_path, 
                                     paste0(site, 
                                            "_LAI_HarmonizedW_10d_", 
                                            distrib, 
                                            ".tif"))
    writeRaster(r_harmonized, out_name_harmonized, overwrite = TRUE)
    
    cat("     [OK] Saved:", basename(out_name_harmonized), "\n")
  }
}