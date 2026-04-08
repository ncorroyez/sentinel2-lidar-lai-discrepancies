# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = 'libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(dplyr)
library(ggplot2)
library(stringr)
library(pdp)
library(tidyr)
library(dplyr)
library(fmsb)
library(rlang)
library(data.table)
library(rPref)
library(zoo)
library(randomForest)
library(pdp)
library(patchwork)
library(stringr)
library(grid)
library(RColorBrewer)
library(broom)
library(Metrics)

# 1- Directories & input data
output_dir <- '../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Aigoual', 'Blois')
# sites <- c("Aigoual")
datesAcq <- list('Aigoual' = '2021-07-11', 
                 'Blois' = '2021-06-14', 
                 'Mormal' = '2021-06-14')
name_strategy <- c("LIDFa_lai_LMA_BROWN_N_CHL_psoil_q")
set.seed(42)

# ------------------------------------------------------------------------------
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "_Final13_05.csv")))
csv_to_study <- csv_all
norm <- "DSM_keepTrees"

numeric_cols <- names(csv_to_study)[sapply(csv_to_study, is.numeric)]
csv_to_study[, (numeric_cols) := lapply(.SD, round, 2), .SDcols = numeric_cols]

data_all1 <- csv_to_study %>%
  filter(Norm == norm)


data_Aigoual <- data_all1 %>%
  filter(Site == "Aigoual") %>%
  filter(Depth == 5)
data_Aigoual <- unique(data_Aigoual)
# data_Aigoual <- data_Aigoual %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

data_Blois <- data_all1 %>%
  filter(Site == "Blois") %>%
  filter(Depth == 5)
data_Blois <- unique(data_Blois)
# data_Blois <- data_Blois %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

data_Mormal <- data_all1 %>%
  filter(Site == "Mormal") %>%
  filter(Depth == 9)
data_Mormal <- unique(data_Mormal)
# data_Mormal <- data_Mormal %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

data_all_sites <- data_all1 %>%
  filter(Site == "All_sites") %>%
  filter(Depth == 3) %>%
  head(-1)
# data_all_sites <- data_all_sites %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

my_data <- list(
  Aigoual = data_Aigoual,
  Blois   = data_Blois,
  Mormal  = data_Mormal
)
p_R <- plot_metric_density2(my_data, "R")
p_nrmse <- plot_metric_density2(my_data, "RMSE")
p_slope <- plot_metric_density2(my_data, "Slope")
p_bias <- plot_metric_density2(my_data, "Bias")
# print(p_R)
# print(p_nrmse)
# print(p_slope)

# met_layout <- p_R + p_nrmse + p_slope
met_layout <- (p_R + p_nrmse) / (p_bias + p_slope)  +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
print(met_layout)
# ggsave("metrics_layout.png", width = 10, height = 10, dpi = 300)
# ggsave("metrics_layout2.png", width = 6.5, height = 6.5, units = "in")
# stop()

p_Aigoual <- plot_depth(data_Aigoual, score = R, title = "Aigoual")
p_Blois   <- plot_depth(data_Blois, score = R, title = "Blois")
p_Mormal  <- plot_depth(data_Mormal, score = R, title = "Mormal")
p_All     <- plot_depth(data_all_sites, score = R, title = paste("All Sites", norm))

layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_All)
print(layout)

data_Aigoual <- extract_parameters(data_Aigoual)
data_Blois <- extract_parameters(data_Blois)
data_Mormal <- extract_parameters(data_Mormal)
data_all_sites <- extract_parameters(data_all_sites)

# Frequency
params <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
param_value <- manual_param_config()
combined <- bind_rows(data_Aigoual, data_Blois, data_Mormal, data_all_sites) %>%
  dplyr::select(Site, all_of(params)) %>%
  pivot_longer(
    cols = -Site,
    names_to = "Parameter",
    values_to = "Value"
  )
# combined$Site <- factor(combined$Site, levels = c("Aigoual", "Blois", "Mormal", "All_sites"))
freq_tbl <- combined %>%
  count(Site, Parameter, Value)
all_sites <- unique(combined$Site)
full_grid <- expand_grid(Site = all_sites, param_value)

# Join with actual frequencies and replace missing with 0
freq_tbl_full <- full_grid %>%
  left_join(freq_tbl, by = c("Site", "Parameter", "Value")) %>%
  mutate(n = replace_na(n, 0))

freq_tbl_norm <- freq_tbl_full %>%
  group_by(Site, Parameter) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()
freq_tbl_norm$Site <- factor(freq_tbl_norm$Site, 
                             levels = c("Aigoual", "Blois", "Mormal", "All_sites"),
                             labels = c("Aigoual", "Blois", "Mormal", "Synthesis of 3 Sites"))

reference_lines <- param_value %>%
  group_by(Parameter) %>%
  summarise(
    n_configs = n_distinct(Value),
    uniform_freq = 1 / n_configs,
    .groups = "drop"
  )

p_norm <- ggplot(freq_tbl_norm, 
                 aes(x = factor(Value), y = freq, fill = Site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(data = reference_lines, 
             aes(yintercept = uniform_freq), 
             linetype = "dashed", 
             color = "black", 
             size = 1) +
  facet_wrap(~ Parameter, 
             ncol = 4, 
             scales = "free_x",
             strip.position = "top") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Parameter Value",
    y = "Proportion of Parameter Value in Models",
    fill = "Site"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    strip.background  = element_rect(fill = "grey95", color = NA),
    # strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    # axis.title        = element_text(face = "bold"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none", # none bottom
    legend.key.size   = unit(0.6, "lines"),
    legend.text       = element_text(size = 14),
    legend.title      = element_text(face = "bold")
  )

print(p_norm)

# ------------------------------------------------------------------------------
library(rPref)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(scales)
# stop()
# Function to perform Pareto Front analysis
perform_pareto_analysis <- function(data, site_name) {
  cat("Performing Pareto analysis for", site_name, "\n")
  cat("Original data size:", nrow(data), "\n")
  
  # Create objectives for Pareto optimization
  # Note: rPref uses high() for maximization and low() for minimization
  pareto_pref <- high(R) * low(RMSE) * low(abs(Bias)) * low(abs(Slope - 1))
  
  # Get Pareto front
  pareto_result <- psel(data, pareto_pref)
  
  # Add front levels
  pareto_result$front_level <- 1
  
  # For additional fronts, we can iteratively remove front 1 and find next front
  remaining_data <- data[!rownames(data) %in% rownames(pareto_result), ]
  all_fronts <- pareto_result
  
  front_num <- 2
  # while(nrow(remaining_data) > 0 && front_num <= 10) {  # Limit to 3 fronts for visualization
  while(nrow(remaining_data) > 0) {  # Limit to 3 fronts for visualization
    next_front <- psel(remaining_data, pareto_pref)
    if(nrow(next_front) == 0) break
    
    next_front$front_level <- front_num
    all_fronts <- rbind(all_fronts, next_front)
    
    remaining_data <- remaining_data[!rownames(remaining_data) %in% rownames(next_front), ]
    front_num <- front_num + 1
  }
  
  cat("Front 1 size:", sum(all_fronts$front_level == 1), "\n")
  cat("Front 2 size:", sum(all_fronts$front_level == 2), "\n")
  cat("Front 3 size:", sum(all_fronts$front_level == 3), "\n\n")
  
  return(list(
    front_1 = pareto_result,
    all_fronts = all_fronts,
    site = site_name
  ))
}

# Function to extract parameters from configuration strings
extract_parameters <- function(data) {
  if(!"Config" %in% names(data)) {
    stop("Config column not found in data")
  }
  
  # Extract parameters from Config column
  # Assuming format like "LIDFa_X_lai_Y_LMA_Z_BROWN_A_N_B_CHL_C_psoil_D_q_E"
  data <- data %>%
    mutate(
      LIDFa = as.numeric(str_extract(Config, "LIDFa_([0-9.-]+)", group = 1)),
      lai = as.numeric(str_extract(Config, "lai_([0-9.-]+)", group = 1)),
      LMA = as.numeric(str_extract(Config, "LMA_([0-9.-]+)", group = 1)),
      BROWN = as.numeric(str_extract(Config, "BROWN_([0-9.-]+)", group = 1)),
      N = as.numeric(str_extract(Config, "N_([0-9.-]+)", group = 1)),
      CHL = as.numeric(str_extract(Config, "CHL_([0-9.-]+)", group = 1)),
      psoil = as.numeric(str_extract(Config, "psoil_([0-9.-]+)", group = 1)),
      q = as.numeric(str_extract(Config, "q_([0-9.-]+)", group = 1))
    )
  
  return(data)
}

# Function to create parameter frequency plots
plot_parameter_frequencies <- function(pareto_data, site_name) {
  params <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
  
  # Extract parameters if not already done
  if(!all(params %in% names(pareto_data))) {
    pareto_data <- extract_parameters(pareto_data)
  }
  
  # Create long format for plotting
  long_data <- pareto_data %>%
    dplyr::select(all_of(params)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "Value"
    ) %>%
    filter(!is.na(Value))
  
  # Calculate frequencies
  freq_data <- long_data %>%
    count(Parameter, Value) %>%
    group_by(Parameter) %>%
    mutate(
      freq = n / sum(n),
      total_configs = sum(n)
    ) %>%
    ungroup()
  
  # Create the plot
  p <- ggplot(freq_data, aes(x = factor(Value), y = freq)) +
    geom_col(fill = "steelblue", alpha = 0.7, width = 0.7) +
    geom_text(aes(label = n), 
              vjust = -0.5, 
              size = 3,
              color = "black") +
    facet_wrap(~ Parameter, 
               scales = "free_x",
               ncol = 4,
               strip.position = "top") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = paste("Parameter Frequencies - Pareto Front 1:", site_name),
      subtitle = paste("Total configurations in Front 1:", nrow(pareto_data)),
      x = "Parameter Value",
      y = "Frequency (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  return(p)
}

# Function to create Pareto front visualization
plot_pareto_fronts <- function(pareto_results, site_name) {
  data <- pareto_results$all_fronts
  
  # Create a 2x2 plot showing different objective combinations
  p1 <- ggplot(data, aes(x = RMSE, y = R, color = factor(front_level))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Front") +
    labs(title = "R vs RMSE", x = "RMSE", y = "R") +
    theme_bw()
  
  p2 <- ggplot(data, aes(x = abs(Bias), y = R, color = factor(front_level))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Front") +
    labs(title = "R vs |Bias|", x = "|Bias|", y = "R") +
    theme_bw()
  
  p3 <- ggplot(data, aes(x = abs(Slope - 1), y = R, color = factor(front_level))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Front") +
    labs(title = "R vs |Slope - 1|", x = "|Slope - 1|", y = "R") +
    theme_bw()
  
  p4 <- ggplot(data, aes(x = RMSE, y = abs(Bias), color = factor(front_level))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Front") +
    labs(title = "RMSE vs |Bias|", x = "RMSE", y = "|Bias|") +
    theme_bw()
  
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste("Pareto Fronts:", site_name),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  
  return(combined_plot)
}

# Main analysis function
run_pareto_analysis <- function() {
  # Perform Pareto analysis for each dataset
  cat("=== PARETO FRONT ANALYSIS ===\n\n")
  
  pareto_aigoual <- perform_pareto_analysis(data_Aigoual, "Aigoual")
  pareto_blois <- perform_pareto_analysis(data_Blois, "Blois")
  pareto_mormal <- perform_pareto_analysis(data_Mormal, "Mormal")
  pareto_all_sites <- perform_pareto_analysis(data_all_sites, "All Sites")
  
  # Store results
  pareto_results <- list(
    Aigoual = pareto_aigoual,
    Blois = pareto_blois,
    Mormal = pareto_mormal,
    All_Sites = pareto_all_sites
  )
  
  # Create parameter frequency plots for Front 1
  cat("Creating parameter frequency plots for Front 1...\n")
  
  freq_plots <- list()
  pareto_plots <- list()
  
  for(site in names(pareto_results)) {
    site_data <- pareto_results[[site]]
    
    # Parameter frequencies for Front 1
    freq_plots[[site]] <- plot_parameter_frequencies(site_data$front_1, site)
    
    # Pareto front visualization
    pareto_plots[[site]] <- plot_pareto_fronts(site_data, site)
  }
  
  # Display all plots
  cat("Displaying plots...\n")
  
  # Show parameter frequency plots
  for(site in names(freq_plots)) {
    print(freq_plots[[site]])
  }
  
  # Show Pareto front plots
  # for(site in names(pareto_plots)) {
  #   print(pareto_plots[[site]])
  # }
  
  # Print summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  for(site in names(pareto_results)) {
    cat("\n", site, ":\n")
    front_1 <- pareto_results[[site]]$front_1
    
    cat("Front 1 - Performance Metrics:\n")
    cat("  R: ", sprintf("%.3f ± %.3f (%.3f - %.3f)", 
                         mean(front_1$R), sd(front_1$R), 
                         min(front_1$R), max(front_1$R)), "\n")
    cat("  RMSE: ", sprintf("%.3f ± %.3f (%.3f - %.3f)", 
                            mean(front_1$RMSE), sd(front_1$RMSE), 
                            min(front_1$RMSE), max(front_1$RMSE)), "\n")
    cat("  |Bias|: ", sprintf("%.3f ± %.3f (%.3f - %.3f)", 
                              mean(abs(front_1$Bias)), sd(abs(front_1$Bias)), 
                              min(abs(front_1$Bias)), max(abs(front_1$Bias))), "\n")
    cat("  |Slope-1|: ", sprintf("%.3f ± %.3f (%.3f - %.3f)", 
                                 mean(abs(front_1$Slope - 1)), sd(abs(front_1$Slope - 1)), 
                                 min(abs(front_1$Slope - 1)), max(abs(front_1$Slope - 1))), "\n")
  }
  
  return(pareto_results)
}

# Run the analysis
pareto_results <- run_pareto_analysis()

pareto_Aigoual <- pareto_results$Aigoual$front_1
pareto_Blois <- pareto_results$Blois$front_1
pareto_Mormal <- pareto_results$Mormal$front_1
pareto_AllSites <- pareto_results$All_Sites$front_1

# stop()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 1. Bind the three Front 1 data.frames together
front1_all <- bind_rows(
  pareto_Aigoual  %>% mutate(site = "Aigoual"),
  pareto_Blois    %>% mutate(site = "Blois"),
  pareto_Mormal   %>% mutate(site = "Mormal"),
  pareto_AllSites %>% mutate(site = "Sites combined")
)

# 2. Pivot to long form for the four metrics
metrics_long <- front1_all %>%
  dplyr::select(site, R, RMSE, Bias, Slope, ATBD) %>%
  pivot_longer(
    cols      = c(R, RMSE, Bias, Slope),
    names_to  = "Metric",
    values_to = "Value"
  )

# 3. Find the ATBD config in all fronts for each site
atbd_all <- bind_rows(
  pareto_results$Aigoual$all_fronts %>%
    filter(if_all(c(LIDFa, lai, LMA, BROWN, N, CHL, psoil, q), ~ . == 1)) %>%
    mutate(site = "Aigoual"),
  pareto_results$Blois$all_fronts %>%
    filter(if_all(c(LIDFa, lai, LMA, BROWN, N, CHL, psoil, q), ~ . == 1)) %>%
    mutate(site = "Blois"),
  pareto_results$Mormal$all_fronts %>%
    filter(if_all(c(LIDFa, lai, LMA, BROWN, N, CHL, psoil, q), ~ . == 1)) %>%
    mutate(site = "Mormal"),
  pareto_results$All_Sites$all_fronts %>%
    filter(if_all(c(LIDFa, lai, LMA, BROWN, N, CHL, psoil, q), ~ . == 1)) %>%
    mutate(site = "Sites combined")
) %>% as_tibble()

# 4. Pivot ATBD to long form
atbd_long <- atbd_all %>%
  dplyr::select(site, R, RMSE, Bias, Slope) %>%
  pivot_longer(
    cols      = c(R, RMSE, Bias, Slope),
    names_to  = "Metric",
    values_to = "Value"
  )

# 5. Build one ggplot per metric, overlaying the ATBD star
plots <- metrics_long %>%
  group_split(Metric) %>%
  purrr::map(~{
    df <- .
    metric_name <- unique(df$Metric)
    
    # Set y-axis limits based on the metric
    y_limits <- switch(metric_name,
                       "Bias"  = c(-1, NA),
                       "R"     = c(0.4, NA),
                       "RMSE"  = c(0, NA),
                       "Slope" = c(0, NA),
                       c(NA, NA))  # default: no limit
    
    # Set formatted label
    metric_label <- switch(metric_name,
                           "R"     = "r",
                           "RMSE"  = "RMSE",
                           "Bias"  = "Bias",
                           "Slope" = "Slope",
                           metric_name)
    
    ggplot(df, aes(x = site, y = Value)) +
      geom_boxplot(outlier.shape = NA) +
      # geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
      # overlay ATBD (★) from all fronts
      geom_point(
        data = atbd_long %>% filter(Metric == metric_name),
        aes(x = site, y = Value),
        shape = 8, size = 3, stroke = 1.2, color = "red"
      ) +
      labs(title = metric_label, x = NULL, y = metric_label) +
      coord_cartesian(ylim = y_limits) +
      theme_bw(base_size = 12) +
      theme(
        plot.title   = element_text(face = "bold", hjust = 0.5),
        axis.text.x  = element_text(size = 10)
      )
  })

# 6. Combine into a 2×2 panel
combined <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) +
  plot_annotation(
    title = "Front 1 Pareto Performance by Site  (★ = ATBD config from all fronts)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(combined)


# 1) Define your parameter names and the full grid of possible values:
params      <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
param_value <- manual_param_config()   # same as before: a data.frame with columns Parameter + Value

# 2) Pull only Front 1 from each site and bind them:
front1_combined <- bind_rows(
  pareto_results$Aigoual$front_1 %>%   dplyr::select(Site, all_of(params)),
  pareto_results$Blois$front_1   %>%   dplyr::select(Site, all_of(params)),
  pareto_results$Mormal$front_1  %>%   dplyr::select(Site, all_of(params)),
  pareto_results$All_Sites$front_1 %>% dplyr::select(Site, all_of(params))
)

# 3) Pivot to long form and count frequencies:
freq_tbl_front1 <- front1_combined %>%
  pivot_longer(cols = -Site,
               names_to  = "Parameter",
               values_to = "Value") %>%
  count(Site, Parameter, Value)

# 4) Build the full (Site × Parameter × Value) grid and fill zeros:
all_sites <- unique(freq_tbl_front1$Site)
full_grid <- expand_grid(Site = all_sites, param_value)

freq_tbl_full <- full_grid %>%
  left_join(freq_tbl_front1, by = c("Site", "Parameter", "Value")) %>%
  mutate(n = replace_na(n, 0))

# 5) Normalize within each Site+Parameter:
freq_tbl_norm <- freq_tbl_full %>%
  group_by(Site, Parameter) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(Parameter = recode(Parameter, "LIDFa" = "ALA")) %>%
  ungroup()

# 6) Compute the uniform (baseline) frequency for each Parameter:
reference_lines <- param_value %>%
  group_by(Parameter) %>%
  summarize(
    n_configs    = n_distinct(Value),
    uniform_freq = 1 / n_configs,
    .groups      = "drop"
  ) %>% mutate(Parameter = recode(Parameter, "LIDFa" = "ALA"))

# 7) Make Site a factor in your desired order/labels:
freq_tbl_norm$Site <- factor(
  freq_tbl_norm$Site,
  levels = c("Aigoual", "Blois", "Mormal", "All_sites"),
  labels = c("Aigoual", "Blois", "Mormal", "Synthesis of 3 Sites")
)

# 8) Plot:
p_front1 <- ggplot(freq_tbl_norm, 
                   aes(x = factor(Value), y = freq, fill = Site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(data = reference_lines, 
             aes(yintercept = uniform_freq),
             linetype = "dashed", size = 0.8) +
  facet_wrap(~ Parameter, 
             ncol = 4, 
             scales = "free_x",
             strip.position = "top") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    # title = "Parameter Frequencies – Pareto Front 1",
    x     = "Parameter Value",
    y     = "Proportion of Configurations",
    fill  = "Site"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background  = element_rect(fill = "grey95", color = NA),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom"
  )

print(p_front1)

stop()










# ------------------------------------------------------------------------------
df_list <- list()
sites <- c("Aigoual", "Blois", "Mormal")
for (site in sites){
  lidar_lai_full <- rast(file.path(paste0("../01_DATA/", site, "/LiDAR/lidarlai_res_10_m.tif")))
  lidar_lai_opt <- rast(file.path(paste0("../01_DATA/", site, "/LiDAR/lidar_lai_best_site_depth.tif")))
  s2_lai_raw <- rast(file.path(paste0("../01_DATA/", site, "/Sentinel-2/s2lai_summer_atbd_res_10_m.tif")))
  s2_lai_opt <- rast(file.path(paste0("../01_DATA/", site, "/Sentinel-2/s2lai_summer_best_indiv_res_10_m.tif")))
  vci <- rast(file.path(paste0("../01_DATA/", site, "/LiDAR/vci_res_10_m.tif")))
  mean <- rast(file.path(paste0("../01_DATA/", site, "/LiDAR/mean_res_10_m.tif")))
  lcv <- rast(file.path(paste0("../01_DATA/", site, "/LiDAR/lcv_no_ground.tif")))
  # lcv_no_ground.tif vci_res_10_m.tif
  
  
  mns <- rast(paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/", site, "/LiDAR/dsm/res_1_m/rasterize_canopy.vrt"))
  mns_sd <- aggregate(mns, fact = 10, fun = "sd", na.rm = TRUE)
  mns_mean <- aggregate(mns, fact = 10, fun = "mean", na.rm = TRUE)
  mns_cv <- mns_sd / mns_mean
  mns_cv <- project(mns_cv, lidar_lai_opt)
  mns_cv <- mask(mns_cv, lidar_lai_opt)
  mns_cv_values <- values(mns_cv, na.rm = T)
  
  
  s <- c(lidar_lai_full, lidar_lai_opt, s2_lai_raw, s2_lai_opt, vci, mean, lcv, mns_cv)
  names(s) <- c("lidar_lai_full", "lidar_lai_opt", "s2_lai_raw", "s2_lai_opt", "vci", "mean", "lcv", "mns_cv")
  
  # Remove rows with NA in any layer
  s_clean <- mask(s, sum(is.na(s)) == 0)
  
  # Coords
  xy <- as.data.frame(crds(s_clean))
  values_df <- as.data.frame(s_clean, na.rm = TRUE)
  df <- cbind(xy[complete.cases(values_df), ], values_df)
  df$site <- site
  
  # Create bins from 2 to max lidar_lai_full (rounded up)
  p95 <- quantile(df$lidar_lai_full, 0.95, na.rm = TRUE)
  max_val <- ceiling(p95)
  breaks <- seq(2, max_val, by = 1)
  df$bin <- cut(df$lidar_lai_full, breaks = c(breaks, Inf), include.lowest = TRUE)
  
  # Target total: 5000 samples
  nsamp <- 5000
  nbins <- length(unique(na.omit(df$bin)))
  samples_per_bin <- floor(nsamp / nbins)
  
  # Sample uniformly across bins
  df_sampled <- df %>%
    group_by(bin) %>%
    sample_n(size = min(n(), samples_per_bin), replace = FALSE) %>%
    ungroup()
  
  df_list[[site]] <- df_sampled
  
  # Save to list
  # df_list[[site]] <- df
  # df_list[[site]] <- df[sample(nrow(df), 5000), ]
}

df_Aigoual <- df_list$Aigoual
df_Blois <- df_list$Blois
df_Mormal <- df_list$Mormal




model1 <- lm(lidar_lai_opt ~ s2_lai_raw, data = df_Aigoual)
model2 <- lm(lidar_lai_opt ~ s2_lai_opt, data = df_Aigoual)
model2b <- lm(lidar_lai_opt ~ mns_cv, data = df_Aigoual) # vci mns_cv
model3_add <- lm(lidar_lai_opt ~ s2_lai_opt + mns_cv, data = df_Aigoual)
model3_int <- lm(lidar_lai_opt ~ s2_lai_opt * mns_cv, data = df_Aigoual)

anova(model1, model2, model2b, model3_add, model3_int)
anova(model3_add, model3_int)



























set.seed(42)
library(randomForest)
library(Metrics)
library(dplyr)

data <- df_Mormal # df_Aigoual df_Blois df_Mormal
# 1. Split des données
n <- nrow(data)
train_index <- sample(seq_len(n), size = 0.7 * n)
train <- data[train_index, ]
test  <- data[-train_index, ]

# 2. Modèles RF
rf_full <- randomForest(lidar_lai_full ~ s2_lai_raw + mns_cv + vci, data = train, mtry = 3) # vci mns_cv
rf_opt <- randomForest(lidar_lai_opt ~ s2_lai_opt + mns_cv + vci, data = train, mtry = 3)
rf_raw <- randomForest(lidar_lai_opt ~ s2_lai_raw + mns_cv + vci, data = train, mtry = 3)

# 3. Prédictions sur test
test$pred_full <- predict(rf_full, newdata = test)
test$pred_opt <- predict(rf_opt, newdata = test)
test$pred_raw <- predict(rf_raw, newdata = test)

# 4. Fonction d’évaluation
eval_model <- function(truth, pred) {
  fit <- lm(pred ~ truth)
  data.frame(
    bias = mean(pred - truth),
    rmse = rmse(truth, pred),
    mae  = mae(truth, pred),
    slope = coef(fit)[2],
    r2 = summary(fit)$r.squared
  )
}

# 5. Comparaison des performances
res_full <- eval_model(test$lidar_lai_full, test$pred_full) %>% mutate(model = "full")
res_opt <- eval_model(test$lidar_lai_opt, test$pred_opt) %>% mutate(model = "s2_lai_opt")
res_raw <- eval_model(test$lidar_lai_opt, test$pred_raw) %>% mutate(model = "s2_lai_raw")
results <- bind_rows(res_full, res_opt, res_raw)
print(results)




library(spaMM)

# 1. Ajouter coordonnées
coords <- terra::xyFromCell(rast(file.path("../01_DATA/Mormal/LiDAR/lidarlai_res_10_m.tif")), 
                            cell = as.integer(rownames(df_Mormal)))
df_Mormal$x <- coords[, 1]
df_Mormal$y <- coords[, 2]

# 2. Modèle sans spatial (référence)
mod_nospat <- fitme(lidar_lai_opt ~ s2_lai_opt + vci + mns_cv, 
                    data = df_Mormal, method = "REML")

# 3. Modèle spatial avec effet aléatoire corrélé selon distance
mod_spat <- fitme(lidar_lai_opt ~ s2_lai_opt + vci + mns_cv + Matern(1 | x + y), data = df_Mormal, method = "REML")

# 4. Comparaison des modèles
AIC(mod_nospat, mod_spat)


















mns <- rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/LiDAR/dsm/res_1_m/rasterize_canopy.vrt")
plot(mns)
mns_sd <- aggregate(mns, fact = 10, fun = "sd", na.rm = TRUE)
mns_mean <- aggregate(mns, fact = 10, fun = "mean", na.rm = TRUE)
mns_cv <- mns_sd / mns_mean
mns_cv <- project(mns_cv, lidar_lai_opt)
mns_cv <- mask(mns_cv, lidar_lai_opt)
mns_cv_values <- values(mns_cv, na.rm = T)
