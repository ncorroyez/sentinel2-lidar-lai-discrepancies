# ---
# title: "3.train_predict_prosail_ts.R"
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
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")
source("../libraries/functions_general_tools.R")

# Pre-processing Parameters
data_dir <- '/media/corroyez/MyPassport/01_DATA'
# sites <- c("Aigoual", "Blois", "Mormal")
sites <- c('Blois')

# Read S2 Acquisition Dates (1/month)
data <- read_ods(file.path(data_dir, "S2_Dates_Summer2021.ods"))
nb_bands_s2 <- 10
resolutions <- c(10)

for (site in sites){
  cat("\n========================================\n")
  cat("PROCESSING SITE:", site, "\n")
  cat("========================================\n")
  
  # Setup
  results_dir <- '../../03_RESULTS'
  masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics")
  if (!dir.exists(metrics_dir)) {
    dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
  }
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  # Extract all dates for the current site
  dateAcqs <- get_site_full_dates(site, data)
  
  pad_opt_file <- switch(site,
                         "Aigoual" = "PAD_35.5_40.tif",
                         "Blois"   = "PAD_35.5_40.tif",
                         "Mormal"  = "PAD_31.5_40.tif")
  lidar_lai_depth <- values(rast(file.path(metrics_dir,
                                           "Deciduous_Only",
                                           "PAD_Profiles_dsm_keepTrees",
                                           pad_opt_file)), 
                            na.rm = T)
  lidar_lai_depth <- na.omit(lidar_lai_depth)
  
  # Shapefile
  path_vector <- file.path(data_dir, site, "Geo_Files", 
                           "utm_from_lidar.gpkg")
  
  for (resolution in resolutions){
    for (dateAcq in dateAcqs){
      dateAcq <- as.Date(dateAcq)
      print(dateAcq)
      if (resolution == 10){
        cat("Chosen resolution", resolution, "\n")
      } else if (resolution == 20){
        cat("Chosen resolution", resolution, "\n")
      } else{
        stop("Error: Resolution must be 10 or 20.\n")
      }
      # s2_creation_directory <- paste(data_dir,
      #                                site,
      #                                'S2_Images',
      #                                dateAcq,
      #                                sep = "/")
      # dir.create(path = s2_creation_directory,
      #            showWarnings = FALSE,
      #            recursive = TRUE)
      # results_path <- paste(results_dir, site, sep = '/')
      # S2source <- 'SAFE'
      # saveRaw <- TRUE
      # 
      # # S-2 Pre-Processing: Cloud Mask, Reflectance
      # results <- preprocess_S2_no_dl(
      #                                # dateAcq,
      #                                path_vector,
      #                                s2_creation_directory,
      #                                results_path,
      #                                resolution = resolution,
      #                                S2source = S2source,
      #                                saveRaw = TRUE)
      # results_site_path <- results$results_site_path
      # refl_path <- results$refl_path # Reflectance
      # HDR_Refl <- results$HDR_Refl
      # Cloud_File <- results$Cloud_File
      
      refl_path <- file.path(data_dir, site, "Sentinel-2", dateAcq,
                             paste0(site, "_", dateAcq, "_masked.tif"))
      results_site_path <- file.path(results_dir,
                                     site,
                                     "PROSAIL_2026")
      
      # get S2 geometry from preprocS2
      angles_path <- get_s2_angles(
        path_angles = file.path(data_dir, site, "Sentinel-2", dateAcq, "geom_acq_S2"),
        path_bbox = path_vector,
        dateAcq = dateAcq
      )
      # define ranges for zenith angles of acquisition
      GeomAcq <- list('min' = data.frame('tto' = angles_path$MinAngle['vza'],
                                         'tts' = angles_path$MinAngle['sza'], 
                                         'psi' = angles_path$MinAngle['psi']), 
                      'max' = data.frame('tto' = angles_path$MaxAngle['vza'], 
                                         'tts' = angles_path$MaxAngle['sza'], 
                                         'psi' = angles_path$MaxAngle['psi']))
      # stop()
      # Get Sensor Response
      SRF <- GetRadiometry('Sentinel_2')
      
      # Choose Parameters
      nbModels = 10
      nbSamples = 1000
      FigPlot = FALSE
      MultiplyingFactor = 10000
      
      # Define Spectral Bands and Variables to Estimate
      # Parameters
      Parms2Estimate <- c('lai') #lai, CHL, EWT, LMA, fCover, fAPAR, albedo
      
      # Bands
      bands_10m <- c('B3','B4','B8')
      bands_20m <- c('B3','B4','B5','B6','B7','B8A','B11','B12')
      # bands_select <- list(bands_10m, bands_20m)
      bands_select <- list(bands_10m)
      
      # Define Noise Level
      # NULL is AD/MD noise, NoiseLevel is Gaussian
      # NoiseLevel <- define_noise_levels()
      # noises <- list(NULL_value = NULL, NoiseLevel = NoiseLevel)
      noises <- list(NULL_value = NULL)
      
      # Modify one or more variable distributions
      distribs <- c("atbd")
      distribs <- c("atbd", "atbd_optim_common")
      distribs <- c("atbd_optim_common_brownmodif")
      distribs <- c("optim_Blois_brownmodif")
      
      for (bands in bands_select){
        S2BandSelect <- Bands2Select <- list()
        for (parm in Parms2Estimate){
          S2BandSelect[[parm]] <- bands
          Bands2Select[[parm]] <- match(S2BandSelect[[parm]],SRF$Spectral_Bands)
        }
        
        if (all(bands_20m %in% bands)) {
          res <- "bands_3_4_5_6_7_8A_11_12"
          if (resolution == 10){
            ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
          } else if (resolution == 20){
            ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8A','B11','B12')
          } else{
            stop("Error: Resolution must be 10 or 20.\n")
          }
        }
        else if (all(bands_10m %in% bands)) {
          res <- "bands_3_4_8"
          ImgBandNames <- c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12')
        }
        else{
          stop(paste("Error in bands selection. \n"))
        }
        for (noise in noises) {
          if (is.null(noise)) {
            noise_type <- "addmult"
          }
          else if (is.list(noise)) {
            noise_type <- "mult"
          }
          else {
            stop("Error: Unknown noise type.\n")
          }
          
          for (distrib in distribs) {
            
            # Set input distribution
            InputPROSAIL <- set_prosail_distribution(distrib, refl_path,
                                                     lidar_lai_depth = lidar_lai_depth)
            
            # Define where Results are Saved: one directory per experiment
            
            # Chosen S-2 Image Resolution
            cat("S-2 Image Resolution", resolution, "\n")
            cat("Reflectance Bands", ImgBandNames, "\n")
            PROSAIL_ResPath <- file.path(results_site_path,
                                         paste0('res_', resolution, 'm'))
            # dir.create(path = PROSAIL_ResPath,
            #            showWarnings = FALSE,
            #            recursive = TRUE)
            
            # Chosen Distribution
            cat("Input Distribution", distrib, "\n")
            PROSAIL_ResPath <- file.path(PROSAIL_ResPath,
                                         paste0('PRO4SAIL_INVERSION_', distrib))
            # dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
            
            # Chosen Noise
            cat("Noise Type", noise_type, "\n")
            PROSAIL_ResPath <- file.path(PROSAIL_ResPath, noise_type)
            # dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
            
            # Chosen Inversion Bands Resolution
            cat("Inversion Bands Resolution", res, "\n")
            PROSAIL_ResPath <- file.path(PROSAIL_ResPath, res)
            dir.create(path = PROSAIL_ResPath,
                       showWarnings = FALSE,
                       recursive = TRUE)
            # stop()
            # Extract and Write Statistics
            # summary_df <- extract_summary_stats(InputPROSAIL, PROSAIL_ResPath)
            # plot_distributions(InputPROSAIL,
            #                    save_plots = TRUE,
            #                    dirname = PROSAIL_ResPath,
            #                    filename = paste0("hist_", distrib))
            
            # Choose SAILversion
            if (grepl("clumping", distrib)) {
              SAILversion <- "4SAIL2"
            } else {
              SAILversion <- "4SAIL"
            }
            cat("SAILversion", SAILversion, "\n")
            
            # Train
            modelSVR <- train_prosail_inversion(InputPROSAIL = InputPROSAIL,
                                                Parms2Estimate = Parms2Estimate,
                                                Bands2Select = Bands2Select,
                                                NoiseLevel = noise,
                                                SAILversion = SAILversion,
                                                SRF = SRF,
                                                SpecPROSPECT = NULL,
                                                SpecSOIL = NULL,
                                                SpecATM = NULL,
                                                Path_Results = PROSAIL_ResPath,
                                                nbModels = nbModels,
                                                nbSamples = nbSamples,
                                                FigPlot = FigPlot,
                                                GeomAcq = GeomAcq)
            
            # Predict
            Apply_prosail_inversion(raster_path = refl_path,
                                    HybridModel = modelSVR,
                                    PathOut = PROSAIL_ResPath,
                                    SelectedBands = S2BandSelect,
                                    bandname = ImgBandNames,
                                    MultiplyingFactor = MultiplyingFactor)
            
            # Apply forest masks
            lai_s2 <- rast(paste0(PROSAIL_ResPath, "/",
                                  site, "_", dateAcq,
                                  "_masked_lai.tif"))
            lai_s2_file <- paste0("s2lai_", dateAcq,
                                  "_", distrib,
                                  "_res_10_m.tif")
            writeRaster(lai_s2,
                        file.path(raw_dir, lai_s2_file),
                        filetype = "GTiff",
                        overwrite = T)
            apply_and_save_masks(lai_s2,
                                 lai_s2_file,
                                 masks_dir,
                                 metrics_dir)
            print(paste(distrib, "was saved in all compositions"))
          }
        }
      }
    }
  }
}
