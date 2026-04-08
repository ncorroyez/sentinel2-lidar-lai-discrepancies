# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = '../libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(dplyr)
library(ggplot2)
library(fmsb)
library(tidyr)
library(rlang)
library(progressr)
library(progress)
library(data.table)
library(emoa)
library(stringr)
library(patchwork)

# 1- Directories & input data
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Aigoual', 'Blois')
sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')

# name_strategy <- 'LIDFa_only'
# name_strategy <- 'LIDFa_BROWN'
# name_strategy <- 'LAI_LMA'
# name_strategy <- 'LIDFa_lai'
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

sampling_methods <- c('random', 'stratified', 'stratified_uniform')
norm_methods <- c('DTM', 'DSM')
depths <- 1:30
# depths <- c(2, 5, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30)
# depths <- c(2, 5)

# Initialize combined results container
combined_results <- data.frame()

# Loop through sites, normalization methods, depths, and sampling methods
for (site in sites) {
  cat(paste('Processing site:', site, "\n"))
  for (norm in norm_methods) {
    for (depth in depths) {
      
      depth_results <- data.frame()
      
      for (sampling_method in sampling_methods) {
        
        # Define file paths
        lidar_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                    'sampling', paste0("PAD_", norm, "_Depth_", 
                                                       depth, "_Samples_", sampling_method,
                                                       '_nbSamples_', nbSamplesS2Refl, ".csv"))
        s2_csv_path <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                 'PROSAIL_Models', name_strategy, 
                                 paste0('LAI_estimated_', sampling_method,
                                        '_nbSamples_', nbSamplesS2Refl, '.csv'))
        
        if (!file.exists(lidar_csv_path) || !file.exists(s2_csv_path)) next
        
        lidar_lai <- fread(lidar_csv_path, sep = "\t")
        s2_lai <- as.data.frame(fread(s2_csv_path, sep = "\t", 
                                      header = TRUE, colClasses = "character"))
        
        if (nrow(lidar_lai) != nrow(s2_lai)) next
        
        result_dir <- file.path(output_dir, site, 'PROSAIL_Optimization', 
                                'sampling_results', name_strategy, norm, 
                                paste0("Depth_", depth))
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        
        baseline_config <- paste(paste0(parms2test, "=1"), collapse = "_")
        s2_vals <- as.numeric(s2_lai[[baseline_config]])
        valid_idx <- complete.cases(lidar_lai$lidar_values, s2_vals)
        lidar_vals <- lidar_lai$lidar_values[valid_idx]
        s2_vals <- s2_vals[valid_idx]
        
        lm_model <- lm(lidar_vals ~ s2_vals)
        metrics_baseline <- data.frame(Column = baseline_config,
                                       Parameter = "baseline",
                                       Level = "1",
                                       R = round(cor(lidar_vals, s2_vals), 2),
                                       R2 = round(summary(lm_model)$r.squared, 2),
                                       NRMSE = round(sqrt(mean((lidar_vals - s2_vals)^2, 
                                                               na.rm = TRUE)) / IQR(lidar_vals, na.rm = TRUE), 2),
                                       Slope = round(coef(lm_model)[2], 2),
                                       Intercept = round(coef(lm_model)[1], 2),
                                       Bias = round(mean(lidar_vals - s2_vals, na.rm = TRUE), 2),
                                       ATBD = TRUE,
                                       stringsAsFactors = FALSE)
        
        depth_results <- rbind(depth_results, 
                               cbind(Site = site, Norm = norm, 
                                     Method = sampling_method, Depth = depth, metrics_baseline))
        
        # Test other parameter levels
        for (param in parms2test) {
          available_levels <- extract_levels(param, colnames(s2_lai))
          for (lev in available_levels) {
            config_values <- sapply(parms2test, function(x) if (x == param) lev else 1)
            config_str <- paste(paste0(names(config_values), "=", config_values), collapse = "_")
            
            # if (!(config_str %in% colnames(s2_lai))) next
            
            s2_vals <- as.numeric(s2_lai[[config_str]])
            valid_idx <- complete.cases(lidar_lai$lidar_values, s2_vals)
            lidar_vals <- lidar_lai$lidar_values[valid_idx]
            s2_vals <- s2_vals[valid_idx]
            
            # if (length(lidar_vals) == 0) next
            
            lm_model <- lm(lidar_vals ~ s2_vals)
            metrics_var <- data.frame(Column = config_str,
                                      Parameter = param,
                                      Level = lev,
                                      R = round(cor(lidar_vals, s2_vals), 2),
                                      R2 = round(summary(lm_model)$r.squared, 2),
                                      NRMSE = round(sqrt(mean((lidar_vals - s2_vals)^2, 
                                                              na.rm = TRUE)) / IQR(lidar_lai$lidar_values, na.rm = TRUE), 2),
                                      Slope = round(coef(lm_model)[2], 2),
                                      Intercept = round(coef(lm_model)[1], 2),
                                      Bias = round(mean(lidar_lai$lidar_values - s2_vals, na.rm = TRUE), 2),
                                      ATBD = FALSE,
                                      stringsAsFactors = FALSE)
            depth_results <- rbind(depth_results, 
                                   cbind(Site = site, Norm = norm, 
                                         Method = sampling_method, Depth = depth, metrics_var))
          }
        }
      }
      
      # Save combined CSV for this depth
      write.csv(depth_results, file.path(result_dir, 
                                         paste0("metrics_depth_", depth, ".csv")), row.names = FALSE)
      combined_results <- rbind(combined_results, depth_results)
    }
  }
}

# Save full results
write.csv(combined_results, file.path(output_dir,
                                      paste0("all_results_combined_", 
                                             name_strategy, "_param_variations.csv")), row.names = FALSE)

# --------------------- Parameters Variation Analysis --------------------------
# _param variations only file (fixing all parameters to ATBD except one)
csv <- read.csv(file.path(output_dir, paste0("all_results_combined_",
                                             name_strategy, "_param_variations.csv")))
site <- "Aigoual" # Aigoual Blois Mormal
site <- c("Aigoual", "Blois", "Mormal")
norm <- "DSM" # DTM DSM c("DTM", "DSM")
params <- c("baseline", "LIDFa") # baseline LIDFa lai LMA BROWN N CHL psoil q

# Filter data for the given Site, Norm, and Parameter
filtered_data <- csv %>%
  filter(Site == site, Norm == norm, Parameter == params) %>%
  # filter(Depth == 8) %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column)

# Normalize objectives (scale from 0 to 1)
# normalize <- function(x) (x - min(x)) / (max(x) - min(x))
normalize <- function(x) (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
filtered_data <- filtered_data %>%
  mutate(
    R_norm = normalize(R),  # Higher is better
    NRMSE_norm = 1 - normalize(NRMSE),  # Lower is better
    Bias_norm = 1 - normalize(abs(Bias)),  # Lower |Bias| is better
    Slope_norm = 1 - normalize(abs(Slope - 1)),  # Slope closer to 1 is better
    Intercept_norm = 1 - normalize(abs(Intercept))  # Intercept closer to 0 is better
  )

# Compute ranking score
ranked_data <- filtered_data %>%
  # mutate(score = NRMSE_norm) %>%
  mutate(score = R_norm + NRMSE_norm + Slope_norm) %>%
  # mutate(score = R_norm + NRMSE_norm + Bias_norm + Slope_norm + Intercept_norm) %>%
  arrange(desc(score)) %>%  # Higher score is better
  slice(1:30)

# Print the best configurations
ranked_data %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) %>%
  print()

# --------------------- All/Top Combinations Analysis --------------------------
# All combinations
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "3.csv")))

# top performers parameters
# csv_top <- read.csv(file.path(output_dir, paste0("top_performers_log_",
#                                                  name_strategy, ".csv")))

site <- "Mormal" # Aigoual Blois Mormal
# norm <- "DSM" # DTM DSM c("DTM", "DSM")
method <- "stratified_uniform"

csv_to_study <- csv_all # csv_all csv_top

# Filter data for the given Site, Norm, and Parameter
filtered_data <- csv_to_study %>%
  filter(Site == site, Norm == norm, Method == method) %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) 

# Normalize objectives (scale from 0 to 1)
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
# normalize <- function(x) (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
filtered_data <- filtered_data %>%
  mutate(
    R_norm = normalize(R),
    NRMSE_norm = -normalize(NRMSE),
    Bias_norm = -normalize(abs(Bias)),
    Slope_norm = -normalize(abs(Slope - 1)),
    Intercept_norm = -normalize(abs(Intercept))
  )

# Compute ranking score
ranked_data <- filtered_data %>%
  # mutate(score = NRMSE_norm) %>%
  # mutate(score = R_norm + NRMSE_norm) %>%
  # mutate(score = R_norm + NRMSE_norm + Slope_norm + Intercept_norm) %>%
  mutate(score = R_norm + NRMSE_norm + Bias_norm + Slope_norm + Intercept_norm) %>%
  arrange(desc(score)) %>%  # Higher score is better
  slice(1:100)

# Print the best configurations
ranked_data %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) %>%
  print()

# ------------------- One Config Between Sites Analysis ------------------------

compute_bounds <- function(x, remove_high = TRUE) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_value <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_value
  upper_bound <- Q3 + 1.5 * IQR_value
  
  if (remove_high) {
    return(x <= upper_bound)  # Remove only high values
  } else {
    return(x >= lower_bound)  # Remove only low values
  }
}

# All combinations
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "3.csv")))

# top performers parameters
# csv_top <- read.csv(file.path(output_dir, paste0("top_performers_log_",
#                                                  name_strategy, ".csv")))

site <- "Blois" # Aigoual Blois Mormal
# norm <- "DSM" # DTM DSM c("DTM", "DSM")
norm <- "DSM_keepTrees" # DSM DSM_Above20 DSM_keepTrees
method <- "stratified_uniform"
# column <- "LIDFa.2_lai.5_LMA.1_BROWN.3_N.1_CHL.2_psoil.2_q.1"
column <- "LIDFa.2_lai.5_LMA.1_BROWN.3_N.1_CHL.2_psoil.2_q.1"

# Filter data for the given Site, Norm, and Parameter
init_filtered_data <- csv_to_study %>%
  # filter(
  #   # Column == column,
  #   # Norm == norm, 
  #   # Method == method,
  #   Site == site
  #   # Depth == 4
  #   # Depth >= 5
  #   # Depth <= 10
  # ) %>%
  filter(
    compute_bounds(R, remove_high = FALSE),
    compute_bounds(NRMSE, remove_high = TRUE),
    compute_bounds(abs(Bias), remove_high = TRUE),
    compute_bounds(abs(Slope - 1), remove_high = TRUE),
    compute_bounds(abs(Intercept), remove_high = TRUE)
  ) %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) 

# Normalize objectives (scale from 0 to 1)
# normalize <- function(x) (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
normalize <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
# normalize <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)

filtered_data <- init_filtered_data %>%
  mutate(
    R_norm = normalize(R),
    NRMSE_norm = 1 - normalize(NRMSE),
    Bias_norm = 1 - normalize(abs(Bias)),
    Slope_norm = 1 - normalize(abs(Slope - 1)),
    Intercept_norm = 1 - normalize(abs(Intercept)),
    score = (R_norm + NRMSE_norm) / 2
    # score = (R_norm + NRMSE_norm + Slope_norm) / 3
    # score = (R_norm + NRMSE_norm + Slope_norm + Intercept_norm) / 4
  ) %>%
  group_by(Site, Norm) %>%
  # filter(score >= quantile(score, 0.995, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(score)) %>%
  # Top 1000
  # group_by(Site, Norm) %>%
  group_by(Depth) %>%
  # slice_head(n = 1000) %>%
  ungroup() %>%
  filter (Site == "Blois")
# filter (Norm == "DSM")
# top 1000 / profondeur

# ------------------------------------------------------------------------------
# Boxplot
# Extract parameters from Column
parameter_counts <- filtered_data %>%
  separate_rows(Column, sep = "_") %>%
  count(Site, Norm, Depth, Column, name = "count")

# Keep only relevant parameters
parameter_counts <- parameter_counts %>%
  filter(str_detect(Column, "LIDFa|lai|CHL|LMA|BROWN|N|psoil|q")) %>%
  mutate(category = case_when(
    str_detect(Column, "LIDFa") ~ "LIDFa",
    str_detect(Column, "lai") ~ "lai",
    str_detect(Column, "CHL") ~ "CHL",
    str_detect(Column, "LMA") ~ "LMA",
    str_detect(Column, "BROWN") ~ "BROWN",
    str_detect(Column, "N") ~ "N",
    str_detect(Column, "psoil") ~ "psoil",
    str_detect(Column, "q") ~ "q",
    TRUE ~ "Other"
  ))  %>%
  arrange(category, desc(count))

# One plot per Site
unique_sites <- unique(parameter_counts$Site)
for (site in unique_sites) {
  site_data <- parameter_counts %>% filter(Site == site) %>%
    group_by(Norm, Column) %>%
    summarise(total_count = sum(count), .groups = "drop")
  depth_data <- filtered_data %>% filter(Site == site)
  depth_data$Depth <- as.numeric(as.character(depth_data$Depth))
  depth_data <- depth_data %>% filter(!is.na(Depth))
  if (all(!is.na(depth_data$Depth))) {
    
    # Depth distribution plot (filled by Norm)
    depth_plot <- ggplot(depth_data, aes(x = Depth, fill = Norm)) +
      geom_histogram(bins = 30, alpha = 0.5) +
      labs(title = "Depth Distribution", x = "Depth", y = "Count") +
      theme_minimal()
    
    # Score distribution plot (filled by Norm)
    score_plot <- ggplot(filtered_data, aes(x = score, fill = Norm)) +
      geom_histogram(bins = 30, alpha = 0.5) +
      labs(title = "Score Distribution", x = "Score", y = "Count") +
      theme_minimal()
    
    # Main plot (parameter counts filled by Norm)
    main_plot <- ggplot() +
      geom_bar(data = site_data, aes(x = Column, y = total_count, fill = Norm), 
               stat = "identity", position = "dodge") +
      labs(title = paste("Parameter Counts for Site:", site),
           x = "Parameter",
           y = "Count") +
      theme_minimal() +
      theme(legend.position = "top", 
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Arrange all plots in a grid (2 top, 1 bottom)
    grid.arrange(
      depth_plot, score_plot, main_plot,
      ncol = 2, nrow = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3)),
      heights = c(1, 2)
    )
  } else {
    message("Depth is not numeric for site: ", site)
  }
}

# One plot per Norm
unique_norms <- unique(parameter_counts$Norm)
for (norm in unique_norms) {
  norm_data <- parameter_counts %>% filter(Norm == norm) %>%
    group_by(Site, Column) %>%
    summarise(total_count = sum(count), .groups = "drop")
  
  # Depth data filtered by Norm
  depth_data <- filtered_data %>% filter(Norm == norm)
  depth_data$Depth <- as.numeric(as.character(depth_data$Depth))
  depth_data <- depth_data %>% filter(!is.na(Depth))
  
  # Ensure valid Depth for plotting
  if (all(!is.na(depth_data$Depth))) {
    depth_plot <- ggplot(depth_data, aes(x = Depth, fill = Site)) +
      geom_histogram(bins = 30, alpha = 0.5) +
      labs(title = "Depth Distribution", x = "Depth", y = "Count") +
      theme_minimal()
    
    score_plot <- ggplot(filtered_data, aes(x = score, fill = Site)) +
      geom_histogram(bins = 30, alpha = 0.5) +
      labs(title = "Score Distribution", x = "Score", y = "Count") +
      theme_minimal()
    
    main_plot <- ggplot() +
      geom_bar(data = norm_data, aes(x = Column, y = total_count, fill = Site),
               stat = "identity", position = "dodge") +
      labs(title = paste("Parameter Counts for Norm:", norm),
           x = "Parameter",
           y = "Count") +
      theme_minimal() +
      theme(legend.position = "top", 
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Arrange all plots in a grid (2 top, 1 bottom)
    grid.arrange(
      depth_plot, score_plot, main_plot,
      ncol = 2, nrow = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3)),
      heights = c(1, 2)
    )
    
  } else {
    message("Depth is not numeric for norm: ", norm)
  }
}


# # Facet by Norm (bars grouped by Site)
# ggplot(parameter_counts, aes(x = Column, y = count, fill = Site)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   facet_wrap(~ Norm, scales = "free_x") +  # One plot per Norm
#   labs(title = "Parameter Counts Grouped by Site (Faceted by Norm)",
#        x = "Parameter",
#        y = "Count") +
#   theme_minimal() +
#   theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Facet by Site (bars grouped by Norm)
# ggplot(parameter_counts, aes(x = Column, y = count, fill = Norm)) +
#   geom_bar(stat = "identity", position = "dodge") + 
#   facet_wrap(~ Site, scales = "free_x") +  # One plot per Site
#   labs(title = "Parameter Counts Grouped by Norm (Faceted by Site)",
#        x = "Parameter",
#        y = "Count") +
#   theme_minimal() +
#   theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

# Plot: Bar Chart (Recommended for clarity)
ggplot(parameter_counts, aes(x = Column, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  # coord_flip() +
  labs(title = paste(filtered_data$Site[1], "Parameter Counts in Column (Grouped by Type)"),
       x = "Parameter",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "top")
# Combine per site


##
best_columns_per_depth <- filtered_data %>%
  filter(Site == "Aigoual",
         Norm == "DSM") %>%
  group_by(Depth) %>%
  arrange(desc(score)) %>%
  slice(1) %>%  # Keep only the best row per depth
  ungroup()

# Filter data to retain only the best-performing Columns across all depths
best_configs <- filtered_data %>%
  filter(Column %in% best_columns_per_depth$Column)

# Plot score evolution across depths
ggplot(best_configs, aes(x = Depth, y = score, color = Column, group = Column)) +
  geom_line(size = 1) +   # Connect scores with lines
  geom_point(size = 2) +  # Highlight individual points
  labs(title = paste(best_configs$Site[1], 
                     "Performance of Best Configurations Across Depths"),
       x = "Depth (m)",
       y = "Score",
       color = "Best Column") +
  theme_minimal() +
  theme(legend.position = "right")  # Legend for column names

average_scores <- filtered_data %>%
  group_by(Column) %>%
  summarise(avg_score = mean(score, na.rm = TRUE)) %>%
  arrange(desc(avg_score))
best_column <- average_scores %>%
  slice(1)
print(best_column)
print(head(average_scores, 5))

best_columns_per_depth <- filtered_data %>%
  group_by(Depth) %>%
  arrange(desc(score)) %>%
  slice(1) %>%  # Keep only the best row per depth
  ungroup()

# Count occurrences of each best column and list associated depths
best_column_depths <- best_columns_per_depth %>%
  group_by(Column) %>%
  summarise(
    depths = paste(sort(unique(Depth)), collapse = ", "),  # Concatenate depths
    count = n(),  # Count occurrences across depths
    avg_score = mean(score, na.rm = TRUE)  # Compute average score for reference
  ) %>%
  arrange(desc(avg_score))  # Sort by frequency and performance

# Print the best-performing columns and their associated depths
print(best_column_depths)

filtered_data <- init_filtered_data %>%
  mutate(
    R_norm = normalize(R),
    NRMSE_norm = 1 - normalize(NRMSE),
    Bias_norm = 1 - normalize(abs(Bias)),
    Slope_norm = 1 - normalize(abs(Slope - 1)),
    Intercept_norm = 1 - normalize(abs(Intercept)),
    # score = (R_norm + NRMSE_norm) / 2
    score = (R_norm + NRMSE_norm + Slope_norm) / 3
    # score = (R_norm + NRMSE_norm + Slope_norm + Intercept_norm) / 4
    # score = (R_norm + NRMSE_norm + Slope_norm + Intercept_norm + Bias_norm) / 5
  ) %>%
  filter(Site == "Blois", # Aigoual Blois Mormal
         Norm == "DSM_Above20" # DSM DSM_keepTrees DSM_Above20
  )

ggplot(filtered_data, aes(x = factor(Depth), y = score, fill = factor(Depth))) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  scale_fill_viridis_d() +
  labs(title = paste(filtered_data$Site[1], "Score Distribution per Depth"),
       x = "Depth (m)", y = "Score") +
  ylim(c(0,1)) +
  theme_bw()

# ------------------------------------------------------------------------------

# Compute ranking score
ranked_data <- filtered_data %>%
  # mutate(score = NRMSE_norm) %>%
  # mutate(score = R_norm) %>%
  # mutate(score = R_norm + NRMSE_norm) %>%
  # mutate(score = R_norm + Slope_norm + NRMSE_norm) %>%
  # mutate(score = R_norm + NRMSE_norm + Bias_norm) %>%
  mutate(score = R_norm + NRMSE_norm + Slope_norm) %>%
  # mutate(score = R_norm + NRMSE_norm + Bias_norm + Slope_norm + Intercept_norm) %>%
  arrange(desc(score)) %>%  # Higher score is better
  slice(1:10)

# Print the best configurations
ranked_data %>%
  dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) %>%
  print()

# Count

# Split the "Column" variable into individual parameter components
setDT(ranked_data)
params_list <- strsplit(ranked_data$Column, "_")
param_counts_list <- lapply(params_list, function(x) as.data.table(table(x)))
param_counts <- rbindlist(param_counts_list, fill = TRUE)
final_counts <- param_counts[, .(Count = sum(N)), by = x]
setnames(final_counts, "x", "Parameter")
final_counts <- final_counts[order(Parameter)]
print(final_counts)

# Best per depth
for (d in 3:15) {
  best_config <- csv_to_study %>%
    filter(
      Norm == norm, 
      Method == method,
      Site == site,
      Depth == d
    ) %>%
    dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) %>%
    mutate(
      R_norm = normalize(R),
      NRMSE_norm = 1 - normalize(NRMSE),  
      Bias_norm = 1 - normalize(abs(Bias)),
      Slope_norm = 1 - normalize(abs(Slope - 1)),
      Intercept_norm = 1 - normalize(abs(Intercept))
    ) %>%
    mutate(score = R_norm + NRMSE_norm + Slope_norm + Intercept_norm) %>%
    arrange(desc(score)) %>%
    slice(1)  # Take the best one
  
  cat("Depth:", d, 
      "\n  Best Column:", best_config$Column, 
      "\n  R:", best_config$R, 
      "\n  NRMSE:", best_config$NRMSE, 
      "\n  Slope:", best_config$Slope, 
      "\n  Intercept:", best_config$Intercept, "\n\n")
}



# Initialize dataframe
best_per_depth <- data.frame()
sites <- c("Aigoual", "Blois", "Mormal")

# Loop through sites and depths to extract best column
for (site in sites) {
  for (d in 4:10) {
    best_config <- csv_to_study %>%
      filter(
        Norm == norm,
        Method == method,
        Site == site,
        Depth == d
      ) %>%
      dplyr::select(Site, Norm, Depth, Method, R, NRMSE, Bias, Slope, Intercept, Column) %>%
      mutate(
        R_norm = normalize(R),
        NRMSE_norm = 1 - normalize(NRMSE),
        Bias_norm = 1 - normalize(abs(Bias)),
        Slope_norm = 1 - normalize(abs(Slope - 1)),
        Intercept_norm = 1 - normalize(abs(Intercept))
      ) %>%
      # mutate(score = R_norm + NRMSE_norm + Slope_norm + Intercept_norm) %>%
      mutate(score = R_norm + NRMSE_norm + Bias_norm) %>%
      arrange(desc(score)) %>%
      slice(1)  # Take the best one
    
    best_per_depth <- rbind(best_per_depth, best_config)
  }
}

# Extract parameter numbers from Column names
param_names <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
param_data <- best_per_depth %>%
  mutate(Column_split = str_extract_all(Column, paste0("(", paste(param_names, collapse = "|"), ")\\.\\d+"))) %>%  # Extract patterns
  unnest(cols = c(Column_split)) %>%  # Unnest to get one value per row
  mutate(
    param = str_extract(Column_split, "^[^.]+"),  # Extract the parameter (e.g., LIDFa, lai)
    number = as.numeric(str_extract(Column_split, "\\d+$"))  # Extract the number after the dot
  ) %>%
  dplyr::select(Site, Depth, param, number)

# Plot the evolution of parameter numbers
ggplot(param_data, aes(x = Depth, y = number, color = Site, group = Site)) +
  geom_line(aes(linetype = Site), size = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Evolution of Parameter Numbers with Depth for All Sites",
       x = "Depth",
       y = "Parameter Number",
       color = "Site") +
  theme_bw() +
  facet_wrap(~ param, scales = "free_y")  # Facet by parameter

# ------------------------------------------------------------------------------
# Optim prof w/ ATBD & DSM_keepTreesAbove20
plot_norm <- function(data, title) {
  ggplot(data, aes(x = Norm, y = score, fill = Norm)) +
    geom_boxplot() +
    labs(title = title, x = "Normalization Technique", y = "Score") +
    theme_bw() +
    # ylim(c(0,1)) +
    theme(legend.position = "none")
}

csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "final_study.csv")))

csv_to_study <- csv_all %>%
  # filter(!(Norm == "DSM_Above20" & Depth > 20)) %>%
  filter(
    compute_bounds(R, remove_high = FALSE),
    compute_bounds(NRMSE, remove_high = TRUE),
    compute_bounds(abs(Bias), remove_high = TRUE),
    compute_bounds(abs(Slope - 1), remove_high = TRUE),
    compute_bounds(abs(Intercept), remove_high = TRUE)
  )

data_all <- csv_to_study %>%
  mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3)
data_Aigoual <- data_all %>%
  filter(Site == "Aigoual") %>%
  group_by(Norm) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_Blois <- data_all %>%
  filter(Site == "Blois") %>%
  group_by(Norm) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_Mormal <- data_all %>%
  filter(Site == "Mormal") %>%
  group_by(Norm) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_all_sites <- rbind(data_Aigoual, data_Blois, data_Mormal)

# Build individual plots
p_Aigoual <- plot_norm(data_Aigoual, "Aigoual")
p_Blois   <- plot_norm(data_Blois, "Blois")
p_Mormal  <- plot_norm(data_Mormal, "Mormal")
p_All     <- plot_norm(data_all_sites, "All Sites")

# Arrange them in a 2x2 grid (first row: two panels; second row: two panels)
# For example, we can arrange as:
layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_All)
print(layout)

# Combine data for all sites with an additional column for site name
data_combined <- bind_rows(
  data_Aigoual %>% mutate(Site = "Aigoual"),
  data_Blois %>% mutate(Site = "Blois"),
  data_Mormal %>% mutate(Site = "Mormal"),
  data_all_sites %>% mutate(Site = "3 sites")
)

# Create a single plot where boxplots for each site are side by side for each Norm
ggplot(data_combined, aes(x = Norm, y = score, fill = Site)) +
  geom_boxplot() +
  labs(title = "Comparison of Normalization Techniques Across Sites",
       x = "Normalization Technique", 
       y = "Score") +
  theme_bw() +
  theme(legend.position = "top")

# ------------------------------------------------------------------------------
plot_depth <- function(data, title) {
  ggplot(data, aes(x = factor(Depth), y = score, fill = factor(Depth))) +
    geom_violin(width = 1.5, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    # geom_violin() +
    # geom_boxplot(outlier.shape = NA) +
    scale_fill_viridis_d() +
    labs(title = title, x = "Depth (m)", y = "Score") +
    ylim(0.79, 1) +
    theme_bw() +
    theme(legend.position = "none")
}

# csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              # name_strategy, "5.csv")))
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "AigoualOnly.csv")))
csv_to_study <- csv_all
# csv_to_study <- csv_all %>%
# filter(!(Norm == "DSM_Above20" & Depth > 20)) # csv_all csv_top
# norm <- "DSM_keepTreesAbove20"
norm <- "DSM_keepTrees"

csv_test <- csv_all %>%
  mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3) %>%
  filter(Column == "LIDFa=1_lai=1_LMA=1_BROWN=1_N=1_CHL=1_psoil=1_q=1")

data_all <- csv_to_study %>%
  # mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3) %>%
  # mutate(score = (normalize(R) + (1 - normalize(NRMSE))) / 2) %>%
  # mutate(score = (1 - normalize(NRMSE))) %>%
  mutate(
    norm_R = (R - min(R)) / (max(R) - min(R)),
    norm_NRMSE = 1 - (NRMSE - min(NRMSE)) / (max(NRMSE) - min(NRMSE)),
    norm_slope = 1 - (abs(Slope - 1) - min(abs(Slope - 1))) / 
      (max(abs(Slope - 1)) - min(abs(Slope - 1))),
    # score = (norm_R + norm_NRMSE) / 2,
    score = (norm_R + norm_NRMSE + norm_slope) / 3) %>%
  filter(Norm == norm)
data_Aigoual <- data_all %>%
  filter(Site == "Aigoual") %>%
  # filter(Depth == 6) %>%
  group_by(Depth) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_Blois <- data_all %>%
  filter(Site == "Blois") %>%
  group_by(Depth) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_Mormal <- data_all %>%
  filter(Site == "Mormal") %>%
  group_by(Depth) %>%
  arrange(desc(score)) %>%
  slice_head(n = 1000) %>%
  ungroup()
data_all_sites <- rbind(data_Aigoual, data_Blois, data_Mormal)

p_Aigoual <- plot_depth(data_Aigoual, "Aigoual")
p_Blois   <- plot_depth(data_Blois, "Blois")
p_Mormal  <- plot_depth(data_Mormal, "Mormal")
p_All     <- plot_depth(data_all_sites, paste("All Sites", norm))

layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_All)
print(layout)

# ------------------------------------------------------------------------------
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "AigoualOnly.csv")))
csv_to_study <- csv_all
TOP_N <- 1000  

# Filter for Norm == "DSM" and Depth == 6, and compute score
data_all <- csv_to_study %>%
  # mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3) %>%
  mutate(score = (normalize(R) + (1 - normalize(NRMSE))) / 2) %>%
  # filter(Norm == "DSM_keepTrees")
  filter(Depth == 6)

# Identify top N best configurations based on median score
top_configs <- data_all %>%
  group_by(Column) %>%
  summarise(median_score = median(score, na.rm = TRUE)) %>%
  arrange(desc(median_score)) %>%
  slice_head(n = TOP_N) %>%
  pull(Column)

# Keep only the best configurations in the dataset
data_filtered <- data_all %>%
  filter(Column %in% top_configs) %>%
  mutate(Column = factor(Column, levels = top_configs))  # Ensure ordered factor

# Subset data per site
data_Aigoual <- data_filtered %>% filter(Site == "Aigoual")
data_Blois   <- data_filtered %>% filter(Site == "Blois")
data_Mormal  <- data_filtered %>% filter(Site == "Mormal")

# Define a plotting function for the top configurations
plot_config <- function(data, title) {
  ggplot(data, aes(x = Column, y = score, fill = Column)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.3, color = "black") +
    coord_flip() +
    labs(title = title, x = "Configuration", y = "Score") +
    scale_fill_viridis_d() +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8))
}

# Create plots for each site and the combined dataset
p_Aigoual <- plot_config(data_Aigoual, "Aigoual")
p_Blois   <- plot_config(data_Blois, "Blois")
p_Mormal  <- plot_config(data_Mormal, "Mormal")
p_All     <- plot_config(data_filtered, "All Sites")

# Arrange plots in a 2x2 grid
layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_All)
print(layout)

# ------------------------------------------------------------------------------
median_scores_per_site <- data_all %>%
  group_by(Site, Column) %>%
  summarise(median_score = median(score, na.rm = TRUE), .groups = "drop")

# Extract best configuration for all sites combined
best_config_all <- median_scores_per_site %>%
  group_by(Column) %>%
  summarise(avg_median_score = mean(median_score, na.rm = TRUE)) %>%
  arrange(desc(avg_median_score)) %>%
  slice_head(n = 1) %>%
  pull(Column)

# Extract best configuration for each site separately
best_config_per_site <- median_scores_per_site %>%
  group_by(Site) %>%
  slice_max(order_by = median_score, n = 2, with_ties = FALSE) %>%  # Slice top 2 to handle ties
  arrange(Site, desc(median_score)) %>%  # Sort by median score
  group_by(Site) %>%
  slice_head(n = 1) %>%  # Select the top configuration
  select(Site, Column) %>%
  pivot_wider(names_from = Site, values_from = Column)

# Store site-specific best configurations
best_config_Aigoual <- best_config_per_site$Aigoual
best_config_Blois <- best_config_per_site$Blois
best_config_Mormal <- best_config_per_site$Mormal

# Initial reference configuration (ATBD)
initial_config_column <- "LIDFa.1_lai.1_LMA.1_BROWN.1_N.1_CHL.1_psoil.1_q.1"

# Filter data for relevant configurations
data_config <- csv_to_study %>%
  filter(Norm == "DSM_keepTrees") %>%
  mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3)

# Extract data for each configuration
data_initial_config <- data_config %>% filter(Column == initial_config_column)
data_best_all <- data_config %>% filter(Column == best_config_all)
data_best_Aigoual <- data_config %>% filter(Column == best_config_Aigoual)
# data_best_Blois <- data_config %>% filter(Column == best_config_Blois)
# data_best_Mormal <- data_config %>% filter(Column == best_config_Mormal)

# Function to prepare data for plotting
prepare_plot_data <- function(...) {
  bind_rows(...) %>%
    mutate(Config = case_when(
      Column == initial_config_column ~ "ATBD",
      Column == best_config_all ~ "Averaged All Sites",
      Column == best_config_Aigoual ~ "Best Aigoual",
      # Column == best_config_Blois ~ "Best Blois",
      # Column == best_config_Mormal ~ "Best Mormal",
      TRUE ~ Column  # Just in case
    ))
}

# Prepare data for each site
data_Aigoual <- prepare_plot_data(data_initial_config, data_best_all, data_best_Aigoual) %>% filter(Site == "Aigoual")
# data_Blois <- prepare_plot_data(data_initial_config, data_best_all, data_best_Blois) %>% filter(Site == "Blois")
# data_Mormal <- prepare_plot_data(data_initial_config, data_best_all, data_best_Mormal) %>% filter(Site == "Mormal") 

# Prepare data for all sites combined
# data_all_sites <- prepare_plot_data(data_initial_config, data_best_all, data_best_Aigoual, data_best_Blois, data_best_Mormal) %>%
  # group_by(Depth, Config) %>%
  # summarise(score = mean(score, na.rm = TRUE), .groups = "drop")  # Average across sites

# Function to plot score vs depth
plot_score_vs_depth <- function(data, title) {
  ggplot(data, aes(x = Depth, y = score, color = Config, linetype = Config)) +
    geom_line(size = 1) +  
    labs(title = title, x = "Depth", y = "Score") +
    ylim(c(0.7, 1)) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank())
}

# Generate plots
p_Aigoual <- plot_score_vs_depth(data_Aigoual, "Aigoual")
print(p_Aigoual)
# p_Blois <- plot_score_vs_depth(data_Blois, "Blois")
# p_Mormal <- plot_score_vs_depth(data_Mormal, "Mormal")
# p_all_sites <- plot_score_vs_depth(data_all_sites, "All Sites (Averaged)")

# Arrange the plots in a single 2x2 layout
layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_all_sites)

# Display the layout
print(layout)

# ------------------------------------------------------------------------------
# Optim depth w/ ATBD
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "5.csv")))

# csv_to_study <- csv_all %>%
# filter(Norm == "DSM_keepTreesAbove20" & Depth < 21) %>%
# filter(Column == "LIDFa.1_lai.1_LMA.1_BROWN.1_N.1_CHL.1_psoil.1_q.1") %>%
# filter(
#   compute_bounds(R, remove_high = FALSE),
#   compute_bounds(NRMSE, remove_high = TRUE),
#   compute_bounds(abs(Bias), remove_high = TRUE),
#   compute_bounds(abs(Slope - 1), remove_high = TRUE),
#   compute_bounds(abs(Intercept), remove_high = TRUE)
# )

csv_to_study <- csv_all %>%
  # mutate(score = (normalize(R) + (1 - normalize(NRMSE)) + (1 - normalize(abs(Slope - 1)))) / 3) %>%
  mutate(score = (normalize(R) + (1 - normalize(NRMSE))) / 2) %>%
  filter(Norm == "DSM_keepTreesAbove20" & Depth < 21) %>%
  filter(Column == "LIDFa.1_lai.1_LMA.1_BROWN.1_N.1_CHL.1_psoil.1_q.1")

avg_score <- csv_to_study %>%
  group_by(Depth) %>%
  summarise(score = mean(score, na.rm = TRUE)) %>%
  mutate(Site = "Average")

# Combine the data
plot_data <- bind_rows(csv_to_study, avg_score)

# Plot
ggplot(plot_data, aes(x = Depth, y = score, color = Site)) +
  geom_line() +
  geom_point() +
  labs(title = "Score Evolution with Depth",
       x = "Depth",
       y = "Score",
       color = "Site") +
  theme_bw()
# ------------------------------------------------------------------------------
# S-2 = f(PAD LiDAR / depth)
depths <- 1:10
combined_df <- data.frame()
# norm <- "DSM_keepTreesAbove20"
norm <- "DTM"

# Loop over sites and depths
for (site in sites){
  for (d in depths){  
    print(paste("Processing site:", site, "depth:", d))
    
    # Load the LiDAR data
    csv_dsm <- fread(file.path(output_dir, site, "PROSAIL_Optimization/sampling",
                               paste0("PAD_", norm, "_Depth_", d,
                                      "_Samples_stratified_uniform_nbSamples_5000.csv")))
    lidar_values <- csv_dsm$lidar_values
    
    # Load the Sentinel-2 data
    csv_s2 <- fread(file.path(output_dir, site, "PROSAIL_Optimization",
                              "PROSAIL_Models",
                              "LIDFa_lai_LMA_BROWN_N_CHL_psoil_q",
                              "LAI_estimated_stratified_uniform_nbSamples_5000.csv"))
    s2_values <- csv_s2[[1]]  # Extract first column (assuming it contains LAI values)
    
    # Create a data frame with site and depth info
    new_data <- data.frame(
      site = site, 
      depth = d,
      lidar_values = lidar_values, 
      s2_values = s2_values
    )
    
    # Append to the combined data frame
    combined_df <- rbind(combined_df, new_data)
  }
}

# Determine plot limits
min_val <- min(combined_df$lidar_values, combined_df$s2_values, na.rm = TRUE)
max_val <- max(combined_df$lidar_values, combined_df$s2_values, na.rm = TRUE)

# Plot with faceting by site, including the averaged panel
ggplot(combined_df, aes(x = lidar_values, y = s2_values, color = factor(depth))) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 1:1 line
  scale_x_continuous(limits = c(min_val, max_val)) +
  scale_y_continuous(limits = c(min_val, max_val)) +
  labs(title = "S-2 vs LiDAR Values for Different Depths & Sites", 
       x = "LiDAR Values", 
       y = "S-2 Values", 
       color = "Depth") +
  facet_wrap(~site, ncol = 1) +  # Facet by site, including the averaged site
  theme_bw()




for (site in unique(combined_df$site)) {
  site_data <- combined_df %>% filter(site == !!site)
  
  # Determine plot limits
  min_val <- min(site_data$lidar_values, site_data$s2_values, na.rm = TRUE)
  max_val <- max(site_data$lidar_values, site_data$s2_values, na.rm = TRUE)
  
  p <- ggplot(site_data, aes(x = lidar_values, y = s2_values, color = factor(depth))) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(min_val, max_val)) +
    scale_y_continuous(limits = c(min_val, max_val)) +
    labs(title = paste("S-2 vs LiDAR Values -", site), 
         x = "LiDAR Values", 
         y = "S-2 Values", 
         color = "Depth") +
    facet_wrap(~depth, ncol = 2) +  # Facet by depth
    theme_bw()
  
  print(p)
}