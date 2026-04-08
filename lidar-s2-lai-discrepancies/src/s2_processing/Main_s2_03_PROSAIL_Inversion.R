# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(prosail)
library(cli)
library(progressr)
library(progress)
library(parallel)
library(liquidSVM)
library(future)
library(future.apply)
library(data.table)
library(fst)

# 1- Directories & input data
set.seed(42)
input_dir <- '../../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Blois', 'Mormal')
# sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')

# 2- additional parameters
S2BandSelect <- c('B3','B4','B8')
S2BandSelect_2 <- c('B03','B04','B08')
Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
nbSamplesS2Refl <- 5000
nbSamplesTrain <- 1000
Samplingmethods <- 'stratified_uniform'

# 3- define parameters to be tested during inversion and save them
# name_strategy <- 'LIDFa_only'
# name_strategy <- 'LAI_only'
# name_strategy <- 'LMA_only'
# name_strategy <- 'BROWN_only'
# name_strategy <- 'LIDFa_lai'
# name_strategy <- 'LIDFa_BROWN'
# name_strategy <- 'LIDFa_LAI_BROWN'
# name_strategy <- 'LAI_LMA'
name_strategy <- "LIDFa_lai_LMA_BROWN_N_CHL_psoil_q_Final"
# name_strategy <- 'LIDFa_lai_LMA_BROWN'

# parms2test <- c('LIDFa')
# parms2test <- c('lai')
# parms2test <- c('LMA')
# parms2test <- c('BROWN')
# parms2test <- c('LIDFa', 'lai', 'BROWN')
# parms2test <- c('LIDFa', 'BROWN')
# parms2test <- c('lai', 'LMA')
parms2test <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
# parms2test <- c("LIDFa", "lai", "LMA", "BROWN")

####
# name_strategy <- c("LIDFa_lai_LMA")
# parms2test <- c("LIDFa", "lai", "LMA")
# name_strategy <- c("BROWN_CHL_N_psoil_q")
# parms2test <- c("BROWN", "CHL", "N", "psoil", "q")

simulations <- define_parm_combinations(parms2test = parms2test, 
                                        output_dir = output_dir, 
                                        name_strategy = name_strategy, 
                                        overwrite = T)

prosail_models_dir <- file.path(output_dir, sites, 
                                'PROSAIL_Optimization', 'PROSAIL_Models',
                                name_strategy)
names(prosail_models_dir) <- sites
lapply(X = prosail_models_dir, dir.create, showWarnings = F, recursive = T)

# 4- get geometry of acquisition
geom_s2 <- list()
for (site in sites){
  path_angles <- file.path(output_dir, site, 'PROSAIL_Optimization', 'geomAcq_S2')
  path_bbox <- file.path(input_dir, site, 'Geo_Files', 'aoi_bbox.GPKG')
  geom_s2[[site]] <- get_s2_angles(path_angles = path_angles, 
                                   dateAcq = datesAcq[[site]], 
                                   path_bbox = path_bbox)
}

# define combinations of labels
combination_labels <- apply(simulations$grid_simu, 1, function(row) {
  paste(names(row), row, sep = "=", collapse = "_")
})

# 5- run PROSAIL simulations and apply on reflectance data
for (site in sites){
  tictoc::tic()
  print(paste('training PROSAIL for site' , site))
  filename_SVR <- file.path(prosail_models_dir[site], paste0(combination_labels, '.rds'))
  lidar_lai_path <- file.path(input_dir, site, 'LiDAR/PAD_Profiles_Classic',
                              'ladstack.tif')
  lidar_lai_rast <- sum(rast(lidar_lai_path), na.rm = T)
  
  lidar_lai_site_rast <- file.path(input_dir, site, 'LiDAR/lidar_lai_best_site_depth.tif')
  lidar_lai_common_rast <- file.path(input_dir, site, 'LiDAR/lidar_lai_best_common_depth.tif')
  lidar_lai_3m_rast <- file.path(input_dir, site, 'LiDAR/lidar_lai_3m.tif')
  
  # train prosail 
  PROSAIL_train_sensitivity(simulations = simulations,
                            geom_s2 = geom_s2[[site]],
                            lidar_lai_rast = lidar_lai_rast,
                            lidar_lai_site_rast = lidar_lai_site_rast,
                            lidar_lai_common_rast = lidar_lai_common_rast,
                            lidar_lai_3m_rast = lidar_lai_3m_rast,
                            nbSamples = nbSamplesTrain,
                            filename = filename_SVR,
                            S2BandSelect = S2BandSelect)
  
  for (Samplingmethod in Samplingmethods){
    
    # load reflectance data
    filename <- file.path(output_dir, site, 'PROSAIL_Optimization', 'sampling',
                          paste0('S2_reflectance_', Samplingmethod,
                                 '_nbSamples_', nbSamplesS2Refl, '.csv'))
    Refl <- readr::read_delim(file = filename, delim = '\t', show_col_types = F)
    selbands <- match(x = S2BandSelect_2, 
                      table = sub(" .*", "", names(Refl))
                      # table = names(Refl)
    ) 
    Refl <- Refl / 10000
    
    # apply models
    pb <- progress_bar$new(
      format = "PROSAIL inversion [:bar] :percent in :elapsed",
      total = length(combination_labels), clear = FALSE, width = 100)
    LAI_Mean <- LAI_SD <- as.data.frame(matrix(data = NA, 
                                               nrow = nrow(Refl), 
                                               ncol = length(combination_labels)))
    for (i in seq_len(length(combination_labels))){
      pb$tick()
      filename <- filename_SVR[i]
      modelSVR <- lapply(readRDS(file = filename), 
                         liquidSVM::unserialize.liquidSVM)
      LAIest <- PROSAIL_Hybrid_Apply(RegressionModels = modelSVR,
                                     Refl = Refl[,selbands]) 
      LAI_Mean[,i] <- LAIest$MeanEstimate
      LAI_SD[,i] <- LAIest$StdEstimate
    }
    # rename combinations from simulations$grid_simu
    names(LAI_Mean) <- combination_labels
    names(LAI_SD) <- combination_labels
    
    # save estimates
    filename <- file.path(prosail_models_dir[site], 
                          paste0('LAI_estimated_', Samplingmethod,
                                 '_nbSamples_', nbSamplesS2Refl, '.csv'))
    readr::write_delim(x = round(x = LAI_Mean, digits = 4), 
                       file = filename, 
                       delim = '\t')
    filename <- file.path(prosail_models_dir[site], 
                          paste0('LAI_estimated_', Samplingmethod,
                                 '_nbSamples_', nbSamplesS2Refl, '_SD.csv'))
    readr::write_delim(x = round(x = LAI_SD, digits = 4), 
                       file = filename, 
                       delim = '\t')
    rm(LAI_Mean)
    rm(LAI_SD)
    gc()
  }
  tictoc::toc()
}