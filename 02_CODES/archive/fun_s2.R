# Function to generate the summary table
generate_summary_table <- function(data) {
  var_names <- colnames(data)
  n_vars <- length(var_names)
  
  summary_table <- data.frame(matrix(ncol = 7, nrow = n_vars))
  colnames(summary_table) <- c("Variable", "Min", "Max", "Mean", "Std", "Val_u", "Distribution")
  
  for (i in seq_along(var_names)) {
    variable <- var_names[i]
    distribution_info <- get_distribution_info(data[[variable]])
    
    summary_table[i, ] <- c(variable, distribution_info)
  }
  
  return(summary_table)
}

train_prosail_inversion_updated <- function(nbSamples = 2000, 
                                            nbSamplesPerRun = 100, 
                                            nbModels = 20,
                                            SAILversion='4SAIL',
                                            Parms2Estimate = 'lai', 
                                            Bands2Select = NULL, 
                                            NoiseLevel = NULL,
                                            SpecPROSPECT = NULL, 
                                            SpecSOIL = NULL, 
                                            SpecATM = NULL,
                                            SRF = NULL,
                                            Path_Results = './', 
                                            FigPlot = FALSE, 
                                            Force4LowLAI = TRUE,
                                            method = 'liquidSVM', 
                                            verbose = FALSE, 
                                            atdb = FALSE, 
                                            InputPROSAIL = NULL, 
                                            GeomAcq = GeomAcq,
                                            BRF_LUT = BRF_LUT,
                                            BRF_LUT_Noise = BRF_LUT_Noise){
  
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###            PRODUCE A LUT TO TRAIN THE HYBRID INVERSION            ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  
  # generate LUT of BRF corresponding to InputPROSAIL, for a sensor
  # str(SpecPROSPECT)
  # str(SpecSOIL)
  # str(SpecATM)
  # 
  # InputPROSAIL <- get_InputPROSAIL(atbd = TRUE, 
  #                                  GeomAcq = GeomAcq)
  # 
  # str(InputPROSAIL)
  # print(1)
  # res <- Generate_LUT_PROSAIL(SAILversion = '4SAIL',
  #                             InputPROSAIL = InputPROSAIL,
  #                             SpecPROSPECT = SpecPROSPECT,
  #                             SpecSOIL = SpecSOIL,
  #                             SpecATM = SpecATM)
  # 
  # print(2)
  
  # BRF_LUT_1nm <- BRF_LUT$BRF
  # InputPROSAIL$fCover <- BRF_LUT$fCover
  # InputPROSAIL$fAPAR <- BRF_LUT$fAPAR
  # InputPROSAIL$albedo <- BRF_LUT$albedo
  # 
  # str(SpecPROSPECT$lambda)
  # str(SRF)
  # str(BRF_LUT_1nm)
  # 
  # BRF_LUT <- applySensorCharacteristics(wvl = SpecPROSPECT$lambda, 
  #                                     SRF = SRF, 
  #                                     InRefl = BRF_LUT_1nm)
  # # identify spectral bands in LUT
  # rownames(BRF_LUT) <- SRF$Spectral_Bands
  # 
  # BRF_LUT_Noise <- list()
  # for (parm in Parms2Estimate){
  #   BRF_LUT_Noise[[parm]] <- apply_noise_atbd(BRF_LUT)
  # }
  # 
  
  
  
  # write parameters LUT
  output <- matrix(unlist(InputPROSAIL), ncol = length(InputPROSAIL), 
                   byrow = FALSE)
  filename <- file.path(Path_Results,'PROSAIL_LUT_InputParms.txt')
  write.table(x = format(output, digits=3),file = filename,append = F, 
              quote = F, col.names = names(InputPROSAIL),
              row.names = F,sep = '\t')
  a <- as.character(SpecPROSPECT$lambda)
  str(a)
  # Write BRF LUT corresponding to parameters LUT
  filename <- file.path(Path_Results,'PROSAIL_LUT_Reflectance.txt')
  write.table(x = format(t(BRF_LUT), digits=5),file = filename,append = F, 
              quote = F, col.names = a, 
              row.names = F,sep = '\t')
  
  # Which bands will be used for inversion?
  if (is.null(Bands2Select)){
    Bands2Select <- list()
    for (parm in Parms2Estimate){
      Bands2Select[[parm]] <- seq(1,length(a))
    }
  }
  # # Add gaussian noise to reflectance LUT: one specific LUT per parameter
  # if (is.null(NoiseLevel)){
  #   NoiseLevel <- list()
  #   for (parm in Parms2Estimate){
  #     NoiseLevel[[parm]] <- 0.05
  #   }
  # }
  # 
  # # produce LIT with noise
  # BRF_LUT_Noise <- list()
  # for (parm in Parms2Estimate){
  #   BRF_LUT_Noise[[parm]] <- BRF_LUT[Bands2Select[[parm]],]+BRF_LUT[Bands2Select[[parm]],]*matrix(rnorm(nrow(BRF_LUT[Bands2Select[[parm]],])*ncol(BRF_LUT[Bands2Select[[parm]],]),
  #                                                                                                       0,NoiseLevel[[parm]]),nrow = nrow(BRF_LUT[Bands2Select[[parm]],]))
  # }
  # 
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  ###                     PERFORM HYBRID INVERSION                      ###
  ### == == == == == == == == == == == == == == == == == == == == == == ###
  # train SVR for each variable and each run
  modelSVR = list()
  for (parm in Parms2Estimate){
    ColParm <- which(parm==names(InputPROSAIL))
    InputVar <- InputPROSAIL[[ColParm]]
    modelSVR[[parm]] <- PROSAIL_Hybrid_Train(BRF_LUT = BRF_LUT_Noise[[parm]],
                                             InputVar = InputVar,
                                             FigPlot = FigPlot,
                                             nbEnsemble = nbModels,
                                             WithReplacement = Replacement,
                                             method = method, 
                                             verbose = verbose)
  }
  return(modelSVR)
}






# Call the function to extract summary statistics
# summary_statistics <- extract_summary_stats(InputPROSAIL) 

# Print the summary statistics
# print(summary_statistics)

# set_initial_atbd_distribution <- function(atdb){
#   if(atdb){
#     SAILversion <- '4SAIL'
#     
#     # Define distribution (uniform / gaussian)
#     TypeDistrib <- data.frame('lai'      = 'Gaussian',
#                               'LIDFa'    = 'Gaussian',
#                               'q'        = 'Gaussian',
#                               'N'        = 'Gaussian',
#                               'CHL'      = 'Gaussian',
#                               'LMA'      = 'Gaussian',
#                               'Cw_rel'   = 'Gaussian',
#                               'BROWN'    = 'Gaussian',
#                               'psoil'    = 'Gaussian',
#                               'tto'      = 'Uniform',
#                               'tts'      = 'Uniform',
#                               'psi'      = 'Uniform',
#                               'PROT'     = 'Unique',
#                               'CBC'      = 'Unique',
#                               'ANT'      = 'Unique',
#                               'alpha'    = 'Unique',
#                               'TypeLidf' = 'Unique',
#                               'LIDFb'    = 'Unique')
#     
#     # Define min value to set gaussian distribution
#     minval <- data.frame('lai'      = 0.0,
#                          'LIDFa'    = 30.0,
#                          'q'        = 0.1,
#                          'N'        = 1.2,
#                          'CHL'      = 20.0,
#                          'LMA'      = 0.003,
#                          'Cw_rel'   = 0.6,
#                          'BROWN'    = 0.0,
#                          'psoil'    = 0.0,
#                          'tto'      = 0.0,
#                          'tts'      = 20.0,
#                          'psi'      = 0.0)
#     
#     # Define max value to set gaussian distribution
#     maxval <- data.frame('lai'      = 15.0,
#                          'LIDFa'    = 80.0,
#                          'q'        = 0.5,
#                          'N'        = 1.8,
#                          'CHL'      = 90.0,
#                          'LMA'      = 0.011,
#                          'Cw_rel'   = 0.85,
#                          'BROWN'    = 2.0,
#                          'psoil'    = 1.0,
#                          'tto'      = 10.0,
#                          'tts'      = 30.0,
#                          'psi'      = 360.0)
#     
#     # Define mean value to set gaussian distribution
#     Mean <- data.frame('lai'    = 2.0,
#                        'LIDFa'  = 60.0,
#                        'q'      = 0.2,
#                        'N'      = 1.5,
#                        'CHL'    = 45.0,
#                        'LMA'    = 0.005,
#                        'Cw_rel' = 0.75,
#                        'BROWN'  = 0.0,
#                        'psoil'  = 0.25)
#     Std <- data.frame('lai'    = 3.0,
#                       'LIDFa'  = 30.0,
#                       'q'      = 0.5,
#                       'N'      = 0.3,
#                       'CHL'    = 30.0,
#                       'LMA'    = 0.005,
#                       'Cw_rel' = 0.8,
#                       'BROWN'  = 0.3,
#                       'psoil'  = 0.6)
#     GaussianDistrib <- list('Mean' = Mean, 'Std' = Std)
#     
#     # Define unique values
#     uniqueval <- data.frame('PROT'     = 0.0,
#                             'CBC'      = 0.0,
#                             'ANT'      = 0.0,
#                             'alpha'    = 40.0,
#                             'TypeLidf' = 2.0,
#                             'LIDFb'    = 0.0)
#     
#     # Generate truncated Gaussian samples
#     set.seed(42)
#     generate_samples <- function(mean, std, min_val, max_val, size) {
#       samples <- rtruncnorm(n = size, a = min_val, b = max_val, 
#                             mean = mean, sd = std)
#       return(samples)
#     }
#     
#     sample_size <- 2000
#     
#     # Get the unique values of `typedistrib`
#     unique_types <- unique(TypeDistrib)
#     
#     # Create an empty list to store the input samples
#     InputPROSAIL <- list()
#     
#     # Create an empty dataframe to store the summary information
#     summary_df <- data.frame(
#       var_name = character(),
#       min = numeric(),
#       max = numeric(),
#       mean = numeric(),
#       std = numeric(),
#       unique_val = numeric(),
#       distrib = character(),
#       stringsAsFactors = FALSE
#     )
#     
#     # Loop over the unique types
#     for (type in unique_types) {
#       # Get the parameter names for the current type
#       param_names <- names(TypeDistrib)[TypeDistrib == type]
#       
#       # Check if it's a Gaussian distribution
#       if (type == 'Gaussian') {
#         # Get the corresponding mean, std, minval, and maxval dataframes
#         mean_df <- Mean[param_names]
#         std_df <- Std[param_names]
#         minval_df <- minval[param_names]
#         maxval_df <- maxval[param_names]
#         
#         # Generate samples for each parameter
#         for (param in param_names) {
#           samples <- generate_samples(mean_df[[param]], std_df[[param]], 
#                                       minval_df[[param]], maxval_df[[param]], 
#                                       sample_size)
#           InputPROSAIL[[param]] <- samples
#         }
#       } else if (type == 'Uniform') {
#         # Get the corresponding minval_uni and maxval_uni dataframes
#         minval_df <- minval[param_names]
#         maxval_df <- maxval[param_names]
#         uniqueval_df <- NA
#         
#         # Generate samples using runif for each parameter
#         for (param in param_names) {
#           samples <- runif(sample_size, minval_df[[param]], maxval_df[[param]])
#           InputPROSAIL[[param]] <- samples
#         }
#       } else if (type == 'Unique') {
#         # Get the corresponding unique_values dataframe
#         uniqueval_df <- uniqueval[param_names]
#         
#         # Generate samples using runif for each parameter
#         for (param in param_names) {
#           samples <- runif(sample_size, uniqueval_df[[param]])
#           InputPROSAIL[[param]] <- samples
#         }
#       }
#     }
#     
#     # Add remaining parameters that are not in TypeDistrib
#     InputPROSAIL$EWT <- ((InputPROSAIL$LMA)/(1-InputPROSAIL$Cw_rel))-InputPROSAIL$LMA
#     InputPROSAIL$CAR <- 0.25*InputPROSAIL$CHL
#     
#     complete_TypeDistrib <- data.frame(
#       'EWT' = 'Gaussian',
#       'CAR' = 'Gaussian'
#     )
#     
#     TypeDistrib <- cbind(TypeDistrib, complete_TypeDistrib)
#     
#     # remaining_params <- setdiff(names(minval), unlist(names(InputPROSAIL)))
#     type_summary <- data.frame(
#       var_name = "EWT",
#       min = min(InputPROSAIL$EWT),
#       max = max(InputPROSAIL$EWT),
#       mean = mean(InputPROSAIL$EWT),
#       std = sd(InputPROSAIL$EWT),
#       unique_val = NA,
#       distrib = "Gaussian"
#     )
#     summary_df <- rbind(summary_df, type_summary)
#     type_summary <- data.frame(
#       var_name = "CAR",
#       min = min(InputPROSAIL$CAR),
#       max = max(InputPROSAIL$CAR),
#       mean = mean(InputPROSAIL$CAR),
#       std = sd(InputPROSAIL$CAR),
#       unique_val = NA,
#       distrib = "Gaussian"
#     )
#     summary_df <- rbind(summary_df, type_summary)
#   }
#   return(list(InputPROSAIL = data.frame(InputPROSAIL), Summary = summary_df))
# }

# "user_defined" = {
#   cat("User-defined Training\n")
#   PROSAIL_ResPath <- file.path(PROSAIL_ResPath, format(Sys.time(), 
#                                                        "%Y_%m_%d_%H:%M:%S"))
#   dir.create(path = PROSAIL_ResPath, showWarnings = FALSE, recursive = TRUE)
#   
#   # Train
#   modelSVR_user <- train_prosail_inversion(minval = minval, maxval = maxval,
#                                            TypeDistrib = TypeDistrib,
#                                            GaussianDistrib = GaussianDistrib,
#                                            Parms2Estimate = Parms2Estimate,
#                                            Bands2Select = Bands2Select,
#                                            NoiseLevel = NoiseLevel,
#                                            SAILversion = SAILversion,
#                                            SpecPROSPECT = SpecPROSPECT_Sensor,
#                                            SpecSOIL = SpecSOIL_Sensor,
#                                            SpecATM = SpecATM_Sensor,
#                                            Path_Results = PROSAIL_ResPath,
#                                            nbModels = nbModels,
#                                            nbSamples = nbSamples,
#                                            FigPlot = FigPlot, method = method)
#   
#   # Predict
#   pred_user <- Apply_prosail_inversion(raster_path = refl_path,
#                                        HybridModel = modelSVR_user,
#                                        PathOut = PROSAIL_ResPath,
#                                        SelectedBands = S2BandSel,
#                                        bandname = ImgBandNames,
#                                        MaskRaster = Cloud_File,
#                                        MultiplyingFactor = MultiplyingFactor,
#                                        method = method)
# },
# "n_a" = {
#   cat("No Training\n")
# }

# Step 1: Define Distributions (for User-Defined Training)
# S2Geom = get_s2_geometry(refl_path)
# min_max_values <- define_min_max_values(S2Geom)
# minval <- min_max_values$minval
# maxval <- min_max_values$maxval
# 
# distribs <- define_distributions()
# TypeDistrib <- distribs$TypeDistrib
# GaussianDistrib <- distribs$GaussianDistrib