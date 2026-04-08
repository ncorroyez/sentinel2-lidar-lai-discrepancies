# ---
# title: "3.train_predict_prosail.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-01-11"
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
source("../libraries/functions_create_masks.R")
source("../libraries/functions_sentinel_2.R")
source("../libraries/functions_plots.R")

# Pre-processing Parameters
# data_dir <- '../../01_DATA'
data_dir <- '/media/corroyez/MyPassport/01_DATA'
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- "Aigoual" # Mormal Blois Aigoual
sites <- c("Blois")
# sites <- c("Hayes", "Reine")

for (site in sites){
  
  pad_opt_file <- switch(site,
                         "Aigoual" = "PAD_35.5_40.tif",
                         "Blois"   = "PAD_35.5_40.tif",
                         "Mormal"  = "PAD_31.5_40.tif",)
  
  # Setup
  results_dir <- '../../03_RESULTS'
  masks_dir <- file.path(results_dir, site, "LiDAR/Heterogeneity_Masks")
  metrics_dir <- file.path(results_dir, site, "Metrics")
  # lidar_lai_optD <- values(rast(file.path(metrics_dir, 
  #                                         "Deciduous_Only/PAD_Profiles_dsm_keepTrees", 
  #                                         pad_opt_file)), na.rm = T)
  
  # Creation of directories
  if (!dir.exists(metrics_dir)) {
    dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
  }
  raw_dir <- file.path(metrics_dir, "Raw")
  if (!dir.exists(raw_dir)) {
    dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # S2 Acquisition Dates: yyyy-mm-dd format mandatory for PROSAIL
  if (site == "Mormal") {
    # dateAcqs <- c('2021-06-14', '2021-12-21')
    # dateAcqs <- c('2021-06-14')
  } else if (site == "Blois") {
    # dateAcqs <- c('2021-06-14', '2021-12-21')
    dateAcqs <- c('2021-06-14')
  } else if (site == "Aigoual") {
    dateAcqs <- c('2021-06-11',
                  '2021-06-16',
                  '2021-06-21',
                  '2021-06-26',
                  '2021-07-01',
                  '2021-07-11',
                  '2021-07-21',
                  '2021-07-26',
                  '2021-08-15',
                  '2021-08-20',
                  '2021-08-25',
                  '2021-08-30',
                  '2021-09-24')
    # dateAcqs <- c('2021-07-11', '2021-12-18')
    # dateAcqs <- c('2021-07-11')
    # } else if (site == "Hayes") {
    #   # dateAcqs <- c('2021-07-11', '2021-12-18')
    #   dateAcqs <- c('2024-07-30')
    # } else if (site == "Reine") {
    #   # dateAcqs <- c('2021-07-11', '2021-12-18')
    #   dateAcqs <- c('2024-07-30')
  }
  else {
    stop("Unknown site. Please provide a valid site name.")
  }
  path_vector <- paste(data_dir, site, "Geo_Files/utm_init.shp", sep = "/")
  # path_vector <- file.path(metrics_dir, "Not_Masked/mean_res_10_m.tif")
  # vector_polygons <- as.polygons(rast(path_vector), dissolve = FALSE, mask = TRUE)
  # raster_data <- rast(path_vector)
  # vector_perimeter <- as.lines(vector_polygons)
  # writeVector(vector_perimeter, file.path(metrics_dir, "utm_final.shp"), overwrite = T, insert = F)
  # path_vector <- file.path(metrics_dir, "utm_final.shp")
  # aaaaaa
  
  # resolutions <- c(10,20)
  resolutions <- c(10) # 10 20
  
  for (resolution in resolutions){
    for (dateAcq in dateAcqs){
      if (resolution == 10){
        cat("Chosen resolution", resolution, "\n")
      } else if (resolution == 20){
        cat("Chosen resolution", resolution, "\n")
      } else{
        stop("Error: Resolution must be 10 or 20.\n")
      }
      # Check month for seasonality: 12 is winter, else is summer
      if (gregexpr("-12-", dateAcq)[[1]] == 5) {
        season <- "winter"
      } else if (gregexpr("-06-|-07-", dateAcq)[[1]] == 5) {
        season <- "summer"
      } else {
        warning("Month is not 06, 07, or 12: An error may be there.\n")
      }
      
      s2_creation_directory <- paste(data_dir,
                                     site,
                                     'Sentinel-2',
                                     dateAcq,
                                     sep = "/")
      dir.create(path = s2_creation_directory,
                 showWarnings = FALSE,
                 recursive = TRUE)
      results_path <- paste(results_dir, site, sep = '/')
      S2source <- 'SAFE'
      saveRaw <- TRUE
      
      # S-2 Pre-Processing: Cloud Mask, Reflectance
      results <- preprocess_S2(
                               dateAcq,
                               path_vector,
                               s2_creation_directory,
                               results_path,
                               resolution = resolution,
                               S2source = 'SAFE',
                               saveRaw = TRUE)
      results_site_path <- results$results_site_path
      refl_path <- results$refl_path # Reflectance
      HDR_Refl <- results$HDR_Refl
      Cloud_File <- results$Cloud_File
      stop()
      # Refl: 10m, 10m no B08, 10m -> 20m, 10m -> 20m no B08
      # 10m
      # refl <- terra::rast(refl_path)
      #
      # # Modifications
      # res_20m <- "_20m"
      # no_b08 <- "_no_B08"
      #
      # # 10m no B08
      # layer_index <- which(names(refl) == "B08 (835.1 nanometers)")
      # refl_10m_no_b08 <- refl[[setdiff(1:nlyr(refl), layer_index)]]
      # refl_10m_no_b08_path <- paste0(refl_path, no_b08, ".envi")
      # writeRaster(refl_10m_no_b08, refl_10m_no_b08_path, overwrite=T)
      #
      # # 20m
      # refl_20m_path <- paste0(refl_path, res_20m, ".envi")
      # writeRaster(terra::aggregate(refl,
      #                              fact=2,
      #                              fun='mean'),
      #             refl_20m_path, overwrite=T)
      #
      # # 20m no B08
      # refl_20m_no_b08_path <- paste0(refl_path, res_20m, no_b08, ".envi")
      # writeRaster(terra::aggregate(refl[[setdiff(1:nlyr(refl), layer_index)]],
      #                              fact=2,
      #                              fun='mean'),
      #             refl_20m_no_b08_path, overwrite=T)
      #
      # # List with all reflectances
      # reflectances_list <- list(refl_path,
      #                           refl_10m_no_b08,
      #                           refl_20m_path,
      #                           refl_20m_no_b08_path)
      
      # Get Sensor Response
      SRF <- get_sensor_response(HDR_Refl)
      
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
      bands_select <- list(bands_10m, bands_20m)
      bands_select <- list(bands_10m)
      
      # Define Noise Level
      # NULL is AD/MD noise, NoiseLevel is Gaussian
      NoiseLevel <- define_noise_levels()
      noises <- list(NULL_value = NULL, NoiseLevel = NoiseLevel)
      noises <- list(NULL_value = NULL)
      
      # Modify one or more variable distributions
      
      # distribs <- c("atbd", "q_zhang_et_al_2005", "brede_et_al_2020",
      #               "hauser_et_al_2021", "sinha_et_al_2020",
      #               "verhoef_and_bach_2007", "shiklomanov_et_al_2016")
      # distribs <- c("sinha_et_al_2020",
      # #               "verhoef_and_bach_2007", "shiklomanov_et_al_2016")
      # distribs <- c(
      #   "atbd", "atbd_S2Geom", "atbd_JBGeom",
      #   "atbd_fixed_psi", "atbd_psi_low", "atbd_psi_high",
      #   "atbd_fixed_tto", "atbd_tto_low", "atbd_tto_high",
      #   "atbd_fixed_tts", "atbd_tts_low", "atbd_tts_high",
      #   "q_zhang_et_al_2005", "brede_et_al_2020",
      #   "hauser_et_al_2021", "sinha_et_al_2020",
      #   "verhoef_and_bach_2007", "shiklomanov_et_al_2016",
      #   "atbd_n_high", "atbd_n_low", "atbd_n_full",
      #   "atbd_chl_high", "atbd_chl_low", "atbd_chl_full",
      #   "atbd_brown_high", "atbd_brown_low", "atbd_brown_full",
      #   "atbd_ewt_high", "atbd_ewt_low", "atbd_ewt_full",
      #   "atbd_lma_high", "atbd_lma_low", "atbd_lma_full",
      #   "atbd_lidfa_high", "atbd_lidfa_low", "atbd_lidfa_full",
      #   "atbd_lai_high", "atbd_lai_low", "atbd_lai_full",
      #   "atbd_q_high", "atbd_q_low", "atbd_q_full",
      #   "atbd_psoil_high", "atbd_psoil_low", "atbd_psoil_full"
      # )
      
      # distribs <- c("atbd_clumping_uniform", "atbd_clumping_gaussian")
      # distribs <- c("atbd_clumping_uniform")
      # distribs <- c("atbd_S2Geom")
      # distribs <- c("atbd_optim_6m_try")
      # distribs <- c("atbd_optim_common")
      # distribs <- paste0("optim_", site)
      # distribs <- paste0("optim_lidarD_", site)
      # distribs <- c("atbd_codist_false")
      distribs <- c("atbd")
      
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
            InputPROSAIL <- set_prosail_distribution(distrib, refl_path, lidar_lai_depth = lidar_lai_optD)
            
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
            
            # Extract and Write Statistics
            summary_df <- extract_summary_stats(InputPROSAIL, PROSAIL_ResPath)
            plot_distributions(InputPROSAIL,
                               save_plots = TRUE,
                               dirname = PROSAIL_ResPath,
                               filename = paste0("hist_", distrib))
            
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
                                                FigPlot = FigPlot)
            
            # Predict
            Apply_prosail_inversion(raster_path = refl_path, #reflectance_path
                                    HybridModel = modelSVR,
                                    PathOut = PROSAIL_ResPath,
                                    SelectedBands = S2BandSelect,
                                    bandname = ImgBandNames,
                                    MultiplyingFactor = MultiplyingFactor)
            
            # Convert S2 Metrics from .envi to .tif (currently LAI)
            if (distrib == "atbd" 
                || distrib == "atbd_S2Geom" 
                || distrib == "atbd_codist_false" 
                || distrib == "atbd_optim_common"
                || distrib == paste0("optim_", site)
                || distrib == paste0("optim_lidarD_", site)){ # Do it for only the best distribution
              lai_files <- list.files(PROSAIL_ResPath,
                                      pattern = "lai\\.tif$",
                                      full.names = TRUE)
              
              # Only one file ends by "lai.envi"
              lai_s2 <- terra::rast(lai_files[1])
              lai_s2_file <- paste0("s2lai", "_", season, "_",
                                    # "depth_study_res_10_m.tif")
                                    # "depth_study_common_res_10_m.tif")
                                    # "best_indiv_lidarD_res_10_m.tif")
                                    # "atbd_codist_false_res_10_m.tif")
                                    "atbd_res_10_m.tif")
              # "best_indiv_res_10_m.tif")
              writeRaster(lai_s2,
                          file.path(raw_dir, lai_s2_file),
                          filetype = "GTiff",
                          overwrite = T)
              apply_and_save_masks(lai_s2,
                                   lai_s2_file,
                                   masks_dir,
                                   metrics_dir)
              
              # Get VIs
              # vi_path <- file.path(results_site_path,
              #                      "SpectralIndices/res_10_m")
              # vi_files <- list.files(vi_path, full.names = TRUE)
              # 
              # ndvi_file <- vi_files[grepl("_NDVI$", vi_files)]
              # evi_file <- vi_files[grepl("_EVI$", vi_files)]
              # savi_file <- vi_files[grepl("_SAVI$", vi_files)]
              
              # NDVI
              # ndvi_s2 <- terra::rast(ndvi_file)
              # ndvi_filename <- paste0("ndvi", "_", season, "_", "res_10_m.tif")
              # writeRaster(ndvi_s2,
              #             file.path(raw_dir, ndvi_filename),
              #             filetype = "GTiff",
              #             overwrite = T)
              # apply_and_save_masks(ndvi_s2,
              #                      ndvi_filename,
              #                      masks_dir,
              #                      metrics_dir)
              
              # EVI
              # evi_s2 <- terra::rast(evi_file)
              # evi_filename <- paste0("evi", "_", season, "_", "res_10_m.tif")
              # writeRaster(evi_s2,
              #             file.path(raw_dir, evi_filename),
              #             filetype = "GTiff",
              #             overwrite = T)
              # apply_and_save_masks(evi_s2,
              #                      evi_filename,
              #                      masks_dir,
              #                      metrics_dir)
              
              # SAVI
              # savi_s2 <- terra::rast(savi_file)
              # savi_filename <- paste0("savi", "_", season, "_", "res_10_m.tif")
              # writeRaster(savi_s2,
              #             file.path(raw_dir, savi_filename),
              #             filetype = "GTiff",
              #             overwrite = T)
              # apply_and_save_masks(savi_s2,
              #                      savi_filename,
              #                      masks_dir,
              #                      metrics_dir)
              
              # Print message showing the save is finished
              print(paste(distrib, "was saved in all compositions"))
            }
            # if (distrib == "atbd_clumping_uniform"){
            #   lai_files <- list.files(PROSAIL_ResPath,
            #                           pattern = "lai\\.envi$",
            #                           full.names = TRUE)
            #
            #   # Only one file ends by "lai.envi"
            #   lai_s2 <- terra::rast(lai_files[1])
            #   lai_s2_file <- "s2laiclumpinguniform_res_10_m.tif"
            #   writeRaster(lai_s2,
            #               file.path(raw_dir, lai_s2_file),
            #               filetype = "GTiff",
            #               overwrite = T)
            #   apply_and_save_masks(lai_s2,
            #                        lai_s2_file,
            #                        masks_dir,
            #                        metrics_dir)
            #   print(paste(distrib, "was saved in all compositions"))
            # }
          }
        }
      }
    }
  }
}
