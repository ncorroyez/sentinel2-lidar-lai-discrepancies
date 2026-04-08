# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(ggplot2)
library(dplyr)
library(stringr)

# 1- Directories & input data
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Aigoual', 'Blois')
# sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')

# name_strategy <- 'LIDFa_only'
# name_strategy <- 'LAI_only'
# name_strategy <- 'LMA_only'
# name_strategy <- 'BROWN_only'
# name_strategy <- 'LIDFa_BROWN'
# name_strategy <- 'LAI_LMA'
# name_strategy <- c("LIDFa_lai_LMA_q_BROWN_N_CHL_psoil")
# name_strategy <- c("LIDFa_lai_LMA")
# name_strategy <- 'LAI_LMA'
# name_strategy <- 'LIDFa_LAI_BROWN'
# parms2test <- c("LIDFa", "lai", "BROWN")
# parms2test <- c("LIDFa")

####
# name_strategy <- "LIDFa_lai_LMA"
# parms2test <- c("LIDFa", "lai", "LMA")
# name_strategy <- c("BROWN_CHL_N_psoil_q")
# parms2test <- c("BROWN", "CHL", "N", "psoil", "q")
name_strategy <- c("LIDFa_lai_LMA_BROWN_N_CHL_psoil_q")
parms2test <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
nbSamplesS2Refl <- 5000

strategies <- readRDS(file.path(output_dir, "Simulations_Strategy",
                                name_strategy, "simulation_strategy.rds"))
prosail_models_dir <- file.path(output_dir, sites, 
                                'PROSAIL_Optimization', 'PROSAIL_Models',
                                name_strategy)
names(prosail_models_dir) <- sites
lapply(X = prosail_models_dir, dir.create, showWarnings = F, recursive = T)
Samplingmethods <- c('random', 'stratified', 'stratified_uniform')
# Samplingmethods <- 'random'

# Initialize storage for metrics and optimal parameters
all_metrics <- list()
optimal_params_r2 <- data.frame()  # For R²
optimal_params_rmse <- data.frame()  # For RMSE
optimal_params_bias <- data.frame()  # For Bias

for (site in sites) {
  for (Samplingmethod in Samplingmethods) {
    site_dir <- file.path(output_dir, site, "Plots", 
                          name_strategy, Samplingmethod)
    ensure_dir(site_dir)
    
    # Read data files
    s2_csv_path <- file.path(prosail_models_dir[site], 
                             paste0('LAI_estimated_', Samplingmethod,
                                    '_nbSamples_', nbSamplesS2Refl, '.csv'))
    s2_lai <- read.csv(s2_csv_path, header = TRUE, sep = "\t")
    
    lidar_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                'sampling',
                                paste0('LiDAR_LAI_Samples_', Samplingmethod,
                                       '_nbSamples_', nbSamplesS2Refl, '.csv'))
    lidar_lai <- read.csv(lidar_csv_path, header = TRUE, sep = "\t")
    
    # Initialize metrics storage
    metrics <- data.frame()
    
    # Process each S2 LAI column
    for (col in colnames(s2_lai)) {
      s2_values <- s2_lai[[col]]
      
      # Check if this is the ATBD case
      is_atbd <- str_count(col, "\\.1(\\D|$)") == length(parms2test)
      
      # Calculate metrics
      valid <- complete.cases(lidar_lai$lidar_values, s2_values)
      lidar_valid <- lidar_lai$lidar_values[valid]
      s2_valid <- s2_values[valid]
      
      r <- cor(lidar_valid, s2_valid)
      r2 <- r^2
      nrmse <- sqrt(mean((lidar_valid - s2_valid)^2)) / IQR(lidar_valid, na.rm = TRUE)
      bias <- mean(lidar_valid - s2_valid, na.rm = TRUE)
      
      metrics <- rbind(metrics, data.frame(
        Site = site,
        Method = Samplingmethod,
        Parameter = col,
        R = r,
        R2 = r2,
        NRMSE = nrmse,
        Bias = bias,
        Bias_abs = abs(bias),
        Is_ATBD = is_atbd
      ))
      
      # Create ATBD plot if applicable
      if (is_atbd) {
        create_scatterplot(s2_valid, lidar_valid,
                           paste("ATBD vs LiDAR -", site, "(", col, ")"),
                           file.path(site_dir,
                                     paste0(col, "_ATBD_vs_LiDAR_", 
                                            site, ".png")))
      }
    }
    
    # Store metrics
    all_metrics[[paste(site, Samplingmethod)]] <- metrics
    
    # Identify best parameters for this site/method (excluding ATBD)
    best_local_r2 <- metrics %>% 
      filter(!Is_ATBD) %>% 
      slice_max(R, n = 1)
    
    best_local_rmse <- metrics %>% 
      filter(!Is_ATBD) %>% 
      slice_min(NRMSE, n = 1)
    
    best_local_bias <- metrics %>% 
      filter(!Is_ATBD) %>% 
      slice_min(Bias_abs, n = 1)
    
    optimal_params_r2 <- rbind(optimal_params_r2, best_local_r2)
    optimal_params_rmse <- rbind(optimal_params_rmse, best_local_rmse)
    optimal_params_bias <- rbind(optimal_params_bias, best_local_bias)
  }
}

# Find global optimal parameters (excluding ATBD)
global_optimal_r2 <- optimal_params_r2 %>%
  group_by(Parameter) %>%
  summarise(Mean_R2 = mean(R2)) %>%
  slice_max(Mean_R2, n = 1)

global_optimal_rmse <- optimal_params_rmse %>%
  group_by(Parameter) %>%
  summarise(Mean_RMSE = mean(NRMSE)) %>%
  slice_min(Mean_RMSE, n = 1)

global_optimal_bias <- optimal_params_bias %>%
  group_by(Parameter) %>%
  summarise(Mean_Bias = mean(Bias)) %>%
  slice_min(Mean_Bias, n = 1)

# Find the best site-wide parameter (across all methods) for each site
site_optimal_params_r2 <- optimal_params_r2 %>%
  group_by(Site) %>%
  slice_max(R2, n = 1) %>%
  ungroup()

site_optimal_params_rmse <- optimal_params_rmse %>%
  group_by(Site) %>%
  slice_min(NRMSE, n = 1) %>%
  ungroup()

site_optimal_params_bias <- optimal_params_bias %>%
  group_by(Site) %>%
  slice_min(abs(Bias), n = 1) %>%
  ungroup()

# Find the best parameter per sampling method for each site
best_per_method_r2 <- optimal_params_r2 %>%
  group_by(Site, Method) %>%
  slice_max(R2, n = 1) %>%
  ungroup()

best_per_method_rmse <- optimal_params_rmse %>%
  group_by(Site, Method) %>%
  slice_min(NRMSE, n = 1) %>%
  ungroup()

best_per_method_bias <- optimal_params_bias %>%
  group_by(Site, Method) %>%
  slice_min(abs(Bias), n = 1) %>%
  ungroup()

# Generate plots for R² and NRMSE
for (site in sites) {
  for (Samplingmethod in Samplingmethods) {
    # Define directory paths
    site_dir <- file.path(output_dir, site, "Plots", 
                          name_strategy, Samplingmethod)
    
    # Get site-specific optimal parameters (excluding ATBD)
    site_optimal_r2 <- site_optimal_params_r2 %>%
      filter(Site == site)

    site_optimal_rmse <- site_optimal_params_rmse %>%
      filter(Site == site)
    
    site_optimal_bias <- site_optimal_params_bias %>%
      filter(Site == site)
    
    # Get best config for a given sampling method at a given site
    best_method_config_r2 <- best_per_method_r2 %>%
      filter(Site == site, Method == Samplingmethod)

    best_method_config_rmse <- best_per_method_rmse %>%
      filter(Site == site, Method == Samplingmethod)
    
    best_method_config_bias <- best_per_method_bias %>%
      filter(Site == site, Method == Samplingmethod)
    
    # Get global optimal parameter & method
    global_param_r2 <- global_optimal_r2$Parameter[1]
    global_param_rmse <- global_optimal_rmse$Parameter[1]
    global_param_bias <- global_optimal_bias$Parameter[1]
    
    # Site-specific optimal plot (R²)
    site_optimal_plot_path_r2 <- file.path(site_dir, paste0("Site_Optimal_R2_", site_optimal_r2$Parameter, ".png"))
    generate_optimal_plot(site_optimal_r2$Parameter, Samplingmethod,
                          nbSamplesS2Refl,
                          paste("Site Optimal (R²) -", site_optimal_r2$Parameter, "(", Samplingmethod, ")"),
                          site_optimal_plot_path_r2)

    # Site-specific optimal plot (RMSE)
    site_optimal_plot_path_rmse <- file.path(site_dir, paste0("Site_Optimal_NRMSE_", site_optimal_rmse$Parameter, ".png"))
    generate_optimal_plot(site_optimal_rmse$Parameter, Samplingmethod,
                          nbSamplesS2Refl,
                          paste("Site Optimal (NRMSE) -", site_optimal_rmse$Parameter, "(", Samplingmethod, ")"),
                          site_optimal_plot_path_rmse)
    
    # Site-specific optimal plot (Bias)
    site_optimal_plot_path_bias <- file.path(site_dir, paste0("Site_Optimal_Bias_", site_optimal_bias$Parameter, ".png"))
    generate_optimal_plot(site_optimal_bias$Parameter, Samplingmethod, 
                          nbSamplesS2Refl,
                          paste("Site Optimal (Bias) -", site_optimal_bias$Parameter, "(", Samplingmethod, ")"), 
                          site_optimal_plot_path_bias)
    
    # Global optimal plot for this sampling method (R²)
    global_plot_path_r2 <- file.path(site_dir, paste0("Global_Optimal_R2_", global_param_r2, ".png"))
    generate_optimal_plot(global_param_r2, Samplingmethod,
                          nbSamplesS2Refl,
                          paste("Global Optimal (R²) -", global_param_r2, "(", Samplingmethod, ")"),
                          global_plot_path_r2)

    # Global optimal plot for this sampling method (RMSE)
    global_plot_path_rmse <- file.path(site_dir, paste0("Global_Optimal_NRMSE_", global_param_rmse, ".png"))
    generate_optimal_plot(global_param_rmse, Samplingmethod,
                          nbSamplesS2Refl,
                          paste("Global Optimal (NRMSE) -", global_param_rmse, "(", Samplingmethod, ")"),
                          global_plot_path_rmse)
    
    # Global optimal plot for this sampling method (Bias)
    global_plot_path_bias <- file.path(site_dir, paste0("Global_Optimal_Bias_", global_param_bias, ".png"))
    generate_optimal_plot(global_param_bias, Samplingmethod, 
                          nbSamplesS2Refl,
                          paste("Global Optimal (Bias) -", global_param_bias, "(", Samplingmethod, ")"), 
                          global_plot_path_bias)
    
    # Best configuration for this sampling method at this site (R²)
    if (nrow(best_method_config_r2) > 0) {
      best_method_plot_path_r2 <- file.path(site_dir,
                                            paste0("Best_R2_", Samplingmethod, "_",
                                                   best_method_config_r2$Parameter,
                                                   ".png"))
      generate_optimal_plot(best_method_config_r2$Parameter,
                            best_method_config_r2$Method,
                            nbSamplesS2Refl,
                            paste("Best Config (R²) -", best_method_config_r2$Parameter,
                                  "(", Samplingmethod, " at ", site, ")"),
                            best_method_plot_path_r2)
    }

    # Best configuration for this sampling method at this site (RMSE)
    if (nrow(best_method_config_rmse) > 0) {
      best_method_plot_path_rmse <- file.path(site_dir,
                                              paste0("Best_NRMSE_", Samplingmethod, "_",
                                                     best_method_config_rmse$Parameter,
                                                     ".png"))
      generate_optimal_plot(best_method_config_rmse$Parameter,
                            best_method_config_rmse$Method,
                            nbSamplesS2Refl,
                            paste("Best Config (NRMSE) -", best_method_config_rmse$Parameter,
                                  "(", Samplingmethod, " at ", site, ")"),
                            best_method_plot_path_rmse)
    }
    
    # Best configuration for this sampling method at this site (Bias)
    if (nrow(best_method_config_bias) > 0) {
      best_method_plot_path_bias <- file.path(site_dir, 
                                              paste0("Best_Bias_", Samplingmethod, "_", 
                                                     best_method_config_bias$Parameter,
                                                     ".png"))
      generate_optimal_plot(best_method_config_bias$Parameter, 
                            best_method_config_bias$Method, 
                            nbSamplesS2Refl,
                            paste("Best Config (Bias) -", best_method_config_bias$Parameter,
                                  "(", Samplingmethod, " at ", site, ")"), 
                            best_method_plot_path_bias)
    }
  }
}