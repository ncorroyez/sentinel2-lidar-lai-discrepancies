# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(prosail)

# 1- define input & output directories, and 
sites <- c('Aigoual', 'Blois', 'Mormal')
input_dir <- file.path('../01_DATA/sites',sites)
input_bbox <- file.path(input_dir, 'aoi_bbox.GPKG')
bbox <- lapply(X = input_bbox, FUN = sf::st_read)
output_dir <- file.path('../03_RESULTS',sites)
names(bbox) <- names(output_dir) <- sites
dates <- list('Aigoual' = '2021-07-11', 'Blois' = '2021-06-14', 'Mormal' = '2021-06-14')
name_vect <- 'utm_init.shp'
angles <- c('saa', 'sza', 'vaa', 'vza')

# 2- define 
# get sensor response for Sentinel-2
SensorName <- 'Sentinel_2'
SRF <- prosail::GetRadiometry(SensorName)
# define parameters to estimate
Parms2Estimate <- 'lai'
# define spectral bands
S2BandSelect <- Bands2Select <- list()
S2BandSelect <- list('lai' = c('B3','B4','B8'))
Bands2Select <- list('lai' = match(S2BandSelect$lai,SRF$Spectral_Bands))
  
for (site in sites){
  # 1- define image and mask
  image_path <- file.path(output_dir[site], 'raster_samples', 
                          paste(site, '_', dates[[site]], '.tiff', sep = ''))
  mask_path <- file.path(output_dir[site], 'raster_samples', 
                          paste(site, '_', dates[[site]], '_BIN_v2.tiff', sep = ''))
  
  # 2- get geometry of acquisition
  angles_path <- file.path(output_dir[site], 'geomAcq_S2', 
                           paste(angles, '_', dates[[site]], '.tiff', sep = ''))
  names(angles_path) <- angles
  geom <- MeanAngle <- list()
  for (angle in angles){
    geom[[angle]] <- terra::crop(x = terra::rast(angles_path[angle]), y = as(bbox[[site]], "Spatial"))
    MeanAngle[[angle]] <- mean(terra::values(geom[[angle]]))
  }
  MeanAngle$psi <- abs(MeanAngle$vaa - MeanAngle$saa)
  if (MeanAngle$psi>180)
    MeanAngle$psi <- 360-MeanAngle$psi
  
  # 3- apply prosail
  # define output directory where LUTs will be saved
  PROSAIL_ResPath <- file.path(output_dir[[site]],'HybridInversion_1')
  dir.create(path = PROSAIL_ResPath, showWarnings = FALSE,recursive = TRUE)
  
  # get S2 geometry: read metadata file from S2 image

  # define ranges for zenith angles of acquisition
  GeomAcq <- list('min' = data.frame('tto' = max(c(0, MeanAngle$vza-2)),
                                     'tts' = MeanAngle$sza-3, 
                                     'psi' = MeanAngle$psi), 
                  'max' = data.frame('tto' = MeanAngle$vza+2, 
                                     'tts' = MeanAngle$sza+3,
                                     'psi' = MeanAngle$psi))

  # train regression model
  # set method = 'svmRadial' or method = 'svmLinear' to use SVR implemented in the 
  # R package caret instead of liquidSVM
  modelSVR <- train_prosail_inversion(Parms2Estimate = Parms2Estimate,
                                      atbd = TRUE, GeomAcq = GeomAcq, 
                                      SRF = SRF, 
                                      Bands2Select = Bands2Select, 
                                      Path_Results = PROSAIL_ResPath, 
                                      method = 'liquidSVM')
  
  # apply regression model on Sentinel-2 raster data
  Apply_prosail_inversion(raster_path = image_path, 
                          HybridModel = modelSVR, 
                          PathOut = PROSAIL_ResPath,
                          SelectedBands = S2BandSelect, 
                          bandname = SRF$Spectral_Bands, 
                          MaskRaster = mask_path, bigRaster = T,
                          MultiplyingFactor = 10000)
}
