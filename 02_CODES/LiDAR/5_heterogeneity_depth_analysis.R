# ---
# title: "5.heterogeneity_depth_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-06-12"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library("lidR")
library("plotly")
library("terra")
library("viridis")
library("future")

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import useful functions
source("../libraries/functions_plots.R")
source("../libraries/functions_heterogeneity_depth_analysis.R")

# -------------------------------- Setup ---------------------------------------

# Pre-processing Parameters & Directories
results_dir <- '../../03_RESULTS'
figures_dir <- '../../04_FIGURES'

# Sites
sites <- c("Aigoual", "Blois", "Mormal")
sites <- "Blois" # Mormal Blois Aigoual
# sites <- c("Aigoual", "Blois")

# Heterogeneity metrics
heterogeneity_metrics <- list("mean", 
                              # "max", 
                              # "std", 
                              # "cv", 
                              # "cv_lad",
                              # "gap_fraction",
                              # "slope",
                              # "rumple",
                              "lcv"
                              # "lskew",
                              # "vci"
                              # "vdr",
                              # "shade"
)

# heterogeneity_metrics <- list("cvlad", "cvladnorm"
# )
# 
# heterogeneity_metrics <- list("cvladnorm"
# )


# N_intervals
n_intervals <- c(3, 5, 10)
n_intervals <- c(3)

# Choice_equal_deciles
choice_equal_deciles <- "Deciles"

# -------------------- Vegetation Heterogeneity Analysis -----------------------
for (site in sites){
  cat("Processing site:", site, "\n")
  # Composition masks
  if (site == "Aigoual"){
    composition_masks <- list(
      "Full_Composition",
      "Coniferous_Only",
      "Deciduous_Only",
      "Coniferous_Flex",
      "Deciduous_Flex",
      "Not_Masked"
    )
  } else if (site == "Blois"){
    composition_masks <- list(
      "Full_Composition",
      "Coniferous_Only",
      "Deciduous_Only",
      "Coniferous_Flex",
      "Deciduous_Flex",
      "Not_Masked"
    )
  } else if (site == "Mormal"){
    composition_masks <- list(
      "Full_Composition",
      "Coniferous_Only",
      "Deciduous_Only",
      "Coniferous_Flex",
      "Deciduous_Flex",
      "Not_Masked"
    )
  } else{
    stop("Error: Site is not Aigoual or Blois or Mormal")
  }
  for (composition_mask in composition_masks){
    cat("Processing mask:", composition_mask, "\n")
    for (heterogeneity_metric in heterogeneity_metrics){
      cat("Processing metric:", heterogeneity_metric, "\n")
      for (choice_equal_decile in choice_equal_deciles){
        cat("Processing interval type:", choice_equal_decile, "\n")
        for (n_interval in n_intervals){
          if (n_interval == 5){
            inc <- 0.2
          } 
          else if (n_interval == 10){
            inc <- 0.1
          }
          else if (n_interval == 3){
            inc <- 1/3
          }
          else {
            stop("Error on n_intervals")
          }
          cat("Processing number of intervals:", n_interval, "\n")
          quantiles_dir <- file.path(results_dir,
                                     site,
                                     "LiDAR/Heterogeneity_Masks",
                                     "Quantiles",
                                     composition_mask,
                                     heterogeneity_metric,
                                     choice_equal_decile,
                                     paste0(n_interval, "_Quantiles"))
          figures_path <- file.path(figures_dir, 
                                    site,
                                    "Vegetation_Heterogeneity_Analysis",
                                    composition_mask,
                                    heterogeneity_metric,
                                    choice_equal_decile,
                                    paste0(n_interval, "_Quantiles"))
          if (!dir.exists(figures_path)) {
            dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)
          }
          perform_heterogeneity_analysis_2(site,
                                           composition_mask,
                                           heterogeneity_metric,
                                           quantiles_dir,
                                           n_interval,
                                           inc,
                                           figures_path)
        }
      }
    }
  }
}


# ------------------------ Vegetation Depth Analysis ---------------------------
for (site in sites){
  cat("Processing site:", site, "\n")
  # Composition masks
  if (site == "Aigoual"){
    composition_masks <- list(
      "Full_Composition",
      "Coniferous_Only",
      "Deciduous_Only",
      "Coniferous_Flex",
      "Deciduous_Flex"
      # "Not_Masked"
    )
  } else if (site == "Blois"){
    composition_masks <- list(
      "Full_Composition",
      # "Coniferous_Only",
      "Deciduous_Only",
      # "Coniferous_Flex",
      "Deciduous_Flex"
      # "Not_Masked"
    )
  } else if (site == "Mormal"){
    composition_masks <- list(
      "Full_Composition",
      # "Coniferous_Only",
      "Deciduous_Only",
      # "Coniferous_Flex",
      "Deciduous_Flex"
      # "Not_Masked"
    )
  } else{
    stop("Error: Site is not Aigoual or Blois or Mormal")
  }
  for (composition_mask in composition_masks){
    cat("Processing mask:", composition_mask, "\n")
    pad_profiles_dir <- file.path(results_dir,
                                  site,
                                  "Metrics",
                                  composition_mask,
                                  "PAD_Profiles_updated")
    figures_path <- file.path(figures_dir, 
                              site,
                              "Vegetation_Depth_Analysis",
                              "update")
    chm_10m <- terra::rast(file.path(results_dir,
                                     site,
                                     "Metrics",
                                     composition_mask,
                                     "mean_res_10_m.tif"))
    if (!dir.exists(figures_path)) {
      dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)
    }
    perform_depth_analysis_2(site,
                             composition_mask,
                             chm_10m,
                             pad_profiles_dir,
                             figures_path)
  }
}

# --------------------- Heterogeneity & Depth Analysis -------------------------
for (site in sites){
  cat("Processing site:", site, "\n")
  # Composition masks
  if (site == "Aigoual"){
    composition_masks <- list(
      # "Full_Composition",
      # "Coniferous_Only",
      "Deciduous_Only"
      # "Coniferous_Flex",
      # "Deciduous_Flex",
      # "Not_Masked"
    )
  } else if (site == "Blois"){
    composition_masks <- list(
      # "Full_Composition",
      # "Coniferous_Only",
      "Deciduous_Only"
      # "Coniferous_Flex",
      # "Deciduous_Flex",
      # "Not_Masked"
    )
  } else if (site == "Mormal"){
    composition_masks <- list(
      # "Full_Composition",
      # "Coniferous_Only",
      "Deciduous_Only"
      # "Coniferous_Flex",
      # "Deciduous_Flex",
      # "Not_Masked"
    )
  } else{
    stop("Error: Site is not Aigoual or Blois or Mormal")
  }
  for (composition_mask in composition_masks){
    cat("Processing mask:", composition_mask, "\n")
    pad_params <- c("PAD_Profiles_NA", "PAD_Profiles_updated_modifminz_NA")
    for (pad_param in pad_params){
      pad_profiles_dir <- file.path(results_dir,
                                    site,
                                    "Metrics",
                                    composition_mask,
                                    # "PAD_Profiles_updated_modifminz" 
                                    # "PAD_Profiles_updated"
                                    # "PAD_Profiles"
                                    pad_param)
      for (heterogeneity_metric in heterogeneity_metrics){
        cat("Processing metric:", heterogeneity_metric, "\n")
        for (choice_equal_decile in choice_equal_deciles){
          cat("Processing interval type:", choice_equal_decile, "\n")
          for (n_interval in n_intervals){
            if (n_interval == 5){
              inc <- 0.2
            } 
            else if (n_interval == 10){
              inc <- 0.1
            }
            else if (n_interval == 3){
              inc <- 1/3
            }
            else {
              stop("Error on n_intervals")
            }
            cat("Processing number of intervals:", n_interval, "\n")
            quantiles_dir <- file.path(results_dir,
                                       site,
                                       "LiDAR/Heterogeneity_Masks",
                                       "Quantiles",
                                       composition_mask,
                                       heterogeneity_metric,
                                       choice_equal_decile,
                                       paste0(n_interval, "_Quantiles"))
            figures_path <- file.path(figures_dir, 
                                      site,
                                      "Vegetation_Heterogeneity_Depth_Analysis",
                                      pad_param, # dtm dsm
                                      composition_mask,
                                      heterogeneity_metric,
                                      choice_equal_decile)
            chm_10m <- terra::rast(file.path(results_dir,
                                             site,
                                             "Metrics",
                                             composition_mask,
                                             "mean_res_10_m.tif"))
            if (!dir.exists(figures_path)) {
              dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)
            }
            for (mask_chm_bool in c(TRUE)){
              cat("Processing mask chm bool:", mask_chm_bool, "\n")
              perform_heterogeneity_depth_analysis_2(site,
                                                     composition_mask,
                                                     heterogeneity_metric,
                                                     chm_10m,
                                                     quantiles_dir,
                                                     n_interval,
                                                     inc,
                                                     pad_profiles_dir,
                                                     figures_path,
                                                     mask_chm_bool)
            }
          }
        }
      }
    }
  }
}

# ----------------------------- CV_LAD Analysis --------------------------------
for (site in sites){
  cat("Processing site:", site, "\n")
  # Composition masks
  if (site == "Aigoual"){
    composition_masks <- list(
      "Full_Composition"
      # "Coniferous_Only",
      # "Deciduous_Only",
      # "Beech_Only"
    )
  } else if (site == "Blois"){
    composition_masks <- list(
      "Full_Composition"
      # "Deciduous_Only",
      # "Oak_Only"
    )
  } else if (site == "Mormal"){
    composition_masks <- list(
      "Full_Composition"
      # "Deciduous_Only",
      # "Oak_Only",
      # "Beech_Only"
    )
  } else{
    stop("Error: Site is not Aigoual or Blois or Mormal")
  }
  for (composition_mask in composition_masks){
    cat("Processing mask:", composition_mask, "\n")
    cv_lad_dir <- file.path(results_dir,
                            site,
                            "Metrics",
                            composition_mask,
                            "PAD_Profiles_updated/CV_LAD")
    for (heterogeneity_metric in heterogeneity_metrics){
      cat("Processing metric:", heterogeneity_metric, "\n")
      for (choice_equal_decile in choice_equal_deciles){
        cat("Processing interval type:", choice_equal_decile, "\n")
        n_intervals <- 3
        for (n_interval in n_intervals){
          if (n_interval == 5){
            inc <- 0.2
          } 
          else if (n_interval == 10){
            inc <- 0.1
          }
          else if (n_interval == 3){
            inc <- 1/3
          }
          else {
            stop("Error on n_intervals")
          }
          quantiles_dir <- file.path(results_dir,
                                     site,
                                     "LiDAR/Heterogeneity_Masks",
                                     "Quantiles",
                                     composition_mask,
                                     "CV_LAD",
                                     choice_equal_decile,
                                     paste0(n_interval, "_Quantiles"))
          cat("Processing number of intervals:", n_interval, "\n")
          figures_path <- file.path(figures_dir, 
                                    site,
                                    "Vegetation_Heterogeneity_Depth_Analysis",
                                    "update",
                                    composition_mask,
                                    "CV_LAD",
                                    choice_equal_decile)
          chm_10m <- terra::rast(file.path(results_dir,
                                           site,
                                           "Metrics",
                                           composition_mask,
                                           "mean_res_10_m.tif"))
          if (!dir.exists(figures_path)) {
            dir.create(figures_path, showWarnings = FALSE, recursive = TRUE)
          }
          perform_heterogeneity_depth_analysis_cvlad(site,
                                                     composition_mask,
                                                     chm_10m,
                                                     quantiles_dir,
                                                     n_interval,
                                                     inc,
                                                     figures_path)
        }
      }
    }
  }
}
