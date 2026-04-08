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
# csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
# name_strategy, "5.csv")))
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              # name_strategy, "_Final29_04.csv")))
                                              name_strategy, "_Final13_05.csv")))
csv_to_study <- csv_all
# csv_to_study[, (names(csv_to_study)) := lapply(.SD, function(col) {
#   if (is.numeric(col)) round(col, 2) else col
# })]
# csv_to_study <- csv_all %>%
# filter(!(Norm == "DSM_Above20" & Depth > 20)) # csv_all csv_top
# norm <- "DSM"
norm <- "DSM_keepTrees"

data_all1 <- csv_to_study %>%
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

# plot_score <- data_all$NRMSE

data_Aigoual <- data_all1 %>%
  filter(Site == "Aigoual") %>%
  # filter(Depth == 6) %>%
  filter(Depth == 5)
# group_by(Depth) %>%
# arrange(desc(R)) %>%
# arrange(desc(R))
# slice_head(n = 1000) %>%
# ungroup()
data_Aigoual <- unique(data_Aigoual)
# data_Aigoual <- data_Aigoual %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

data_Blois <- data_all1 %>%
  filter(Site == "Blois") %>%
  # filter(Depth == 6) %>%
  filter(Depth == 5)
# group_by(Depth) %>%
# arrange(desc(R)) %>%
# arrange(desc(R))
# slice_head(n = 1000) %>%
# ungroup()
data_Blois <- unique(data_Blois)
# data_Blois <- data_Blois %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

data_Mormal <- data_all1 %>%
  filter(Site == "Mormal") %>%
  # filter(Depth == 9) %>%
  filter(Depth == 9)
# group_by(Depth) %>%
# arrange(desc(R)) %>%
# slice_head(n = 1000) %>%
# arrange(desc(R))
# ungroup()
data_Mormal <- unique(data_Mormal)
# data_Mormal <- data_Mormal %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

# data_all_sites <- rbind(data_Aigoual, data_Blois, data_Mormal)

data_all_sites <- data_all1 %>%
  filter(Site == "All_sites") %>%
  filter(Depth == 3) %>%
  # filter(Depth == 4) %>%
  head(-1)
# data_all_sites <- data_all_sites %>%
#   filter(R >= quantile(R, 0.8),
#          RMSE <= quantile(RMSE, 0.2),
#          abs(Slope - 1) <= quantile(abs(Slope - 1), 0.2),
#          abs(Bias) <= quantile(abs(Bias), 0.2)
#   )

# stop()

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
# stop()
# justifier qt
# f rel
# lai_t = f(lai6m + met)

# ------------------------------------------------------------------------------
data_Aigoual <- extract_parameters(data_Aigoual)
data_Blois <- extract_parameters(data_Blois)
data_Mormal <- extract_parameters(data_Mormal)
data_all_sites <- extract_parameters(data_all_sites)

# Frequency
params <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")
param_value <- manual_param_config()
combined <- bind_rows(data_Aigoual, data_Blois, data_Mormal, data_all_sites) %>%
  select(Site, all_of(params)) %>%
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
stop()
# ggsave("count.png", width = 10, height = 6, unit = "in")

# Step 2 — Get common columns
# common_columns <- Reduce(intersect, list(data_Aigoual$Column, 
#                                          data_Blois$Column, 
#                                          data_Mormal$Column))
common_columns <- data_all_sites$Column

# Step 3 — Filter datasets to only common columns
data_Aigoual_common <- data_Aigoual[Column %in% common_columns]
data_Blois_common   <- data_Blois[Column %in% common_columns]
data_Mormal_common  <- data_Mormal[Column %in% common_columns]

data_common_all <- rbindlist(list(
  cbind(Site = "Aigoual", data_Aigoual_common),
  cbind(Site = "Blois", data_Blois_common),
  cbind(Site = "Mormal", data_Mormal_common)
))

# One ranking row per site/column
data_all_sites[, `:=`(
  R_rank     = frank(-R, ties.method = "min"),
  RMSE_rank  = frank(RMSE, ties.method = "min"),
  Slope_rank = frank(abs(Slope - 1), ties.method = "min"),
  Bias_rank  = frank(abs(Bias), ties.method = "min")
), by = Site]

# Sum ranks across sites per column
common_column_ranking <- data_all_sites[, .(
  total_rank = sum(R_rank + RMSE_rank + Slope_rank + Bias_rank)
  # total_rank = sum(R_rank + RMSE_rank)
), by = Column]

# Pick the best one
# best_common_column <- common_column_ranking[which.min(total_rank), Column]
best_common_column <- common_column_ranking[order(total_rank)][1:5, Column]
common_stats <- data_all_sites[Column %in% best_common_column, .(Site, R, RMSE, Slope, Bias)]

# Step 5 — Apply ranking to all data for individual best
best_indiv_Aigoual <- rank_parameters(data_Aigoual) %>%
  arrange(rank) %>%
  slice_head(n = 5)
best_indiv_Blois   <- rank_parameters(data_Blois) %>%
  arrange(rank) %>%
  slice_head(n = 5)
best_indiv_Mormal  <- rank_parameters(data_Mormal) %>%
  arrange(rank) %>%
  slice_head(n = 5)

# Step 6 — Final table
result <- data.table(
  Site = c("Aigoual", "Blois", "Mormal"),
  
  # Best common column info
  Best_Common_Column = best_common_column,
  R_common     = common_stats$R,
  RMSE_common = common_stats$RMSE,
  Slope_common = common_stats$Slope,
  Bias_common = common_stats$Bias,
  
  # Best individual column info
  Best_Indiv_Column = c(best_indiv_Aigoual$Column,
                        best_indiv_Blois$Column,
                        best_indiv_Mormal$Column),
  R_indiv     = c(best_indiv_Aigoual$R,
                  best_indiv_Blois$R,
                  best_indiv_Mormal$R),
  RMSE_indiv = c(best_indiv_Aigoual$RMSE,
                 best_indiv_Blois$RMSE,
                 best_indiv_Mormal$RMSE),
  Slope_indiv = c(best_indiv_Aigoual$Slope,
                  best_indiv_Blois$Slope,
                  best_indiv_Mormal$Slope),
  Bias_indiv = c(best_indiv_Aigoual$Bias,
                 best_indiv_Blois$Bias,
                 best_indiv_Mormal$Bias)
)
print(result)
stop()
# all_data <- rbindlist(list(data_Aigoual_common, 
#                            data_Blois_common, 
#                            data_Mormal_common), 
#                       use.names = TRUE)

all_data <- rbindlist(list(data_Aigoual, 
                           data_Blois, 
                           data_Mormal,
                           data_all_sites), 
                      use.names = TRUE, fill = TRUE)

# Select parameters to analyze
params <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")

# Melt to long format for faceting
long_data <- melt(all_data,
                  measure.vars = params,
                  variable.name = "Parameter",
                  value.name = "Value")

# Plot RMSE ~ Value for each Parameter, with boxplots only
ggplot(long_data, aes(x = factor(Value), y = RMSE, fill = Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ Parameter, scales = "free_x", nrow = 4) +
  labs(title = "RMSE by Parameter Values (All Sites)",
       x = "Parameter Value",
       y = "RMSE") +
  theme_bw() +
  theme(legend.position = "bottom")

# Plot Slope ~ Value for each Parameter, with boxplots only
ggplot(long_data, aes(x = factor(Value), y = Slope, fill = Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ylim(c(0.5, 1)) +
  facet_wrap(~ Parameter, scales = "free_x", nrow = 4) +
  labs(title = "Slope by Parameter Values (All Sites)",
       x = "Parameter Value",
       y = "Slope") +
  theme_bw() +
  theme(legend.position = "bottom")

# Plot R ~ Value for each Parameter, with boxplots only
ggplot(long_data, aes(x = factor(Value), y = R, fill = Site)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ylim(c(0.3, 1)) +
  facet_wrap(~ Parameter, scales = "free_x", nrow = 4) +
  labs(title = "R by Parameter Values (All Sites)",
       x = "Parameter Value",
       y = "R") +
  theme_bw() +
  theme(legend.position = "bottom")

# Assuming long_data is a data.frame or data.table with: Site, Parameter, Value, RMSE, Slope
# Find value(s) with lowest RMSE
best_rmse <- long_data %>%
  group_by(Site, Parameter) %>%
  slice_min(order_by = RMSE, with_ties = FALSE) %>%
  mutate(Metric = "Best RMSE") %>%
  ungroup()

# Find value(s) with slope closest to 1
best_slope <- long_data %>%
  group_by(Site, Parameter) %>%
  slice_min(order_by = abs(Slope - 1), with_ties = FALSE) %>%
  mutate(Metric = "Best Slope") %>%
  ungroup()

# Combine both results
best_values <- bind_rows(best_rmse, best_slope) %>%
  select(Site, Parameter, Value, RMSE, Slope, Metric)
print(best_values)

data_Aigoual_filtered <- data_Aigoual[Column %in% c(best_common_column, best_indiv_Aigoual$Column)]
data_Blois_filtered   <- data_Blois[Column %in% c(best_common_column, best_indiv_Blois$Column)]
data_Mormal_filtered  <- data_Mormal[Column %in% c(best_common_column, best_indiv_Mormal$Column)]

stop()

# ------------------------------------------------------------------------------
atbd <- "s2lai_summer_atbd_res_10_m.tif"
common <- "s2lai_summer_depth_study_common_res_10_m.tif"
best_indiv <- "s2lai_summer_best_indiv_res_10_m.tif"
vci <- "vci_res_10_m.tif"

# sites <- "Blois"
list_plots <- list()
# placeholders for plots
plots1  <- list()  # will hold init1, c, e
plots2  <- list()  # will hold init1_opt, d, f
init1_list     <- list()
c_list         <- list()
e_list         <- list()
init1_opt_list <- list()
d_list         <- list()
f_list         <- list()

for (site in sites){
  dir <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/", 
                site, "/Metrics/Deciduous_Only")
  dsm <- "testPADs/PAD_Profiles_DSM_keepTrees"
  data_dir <- paste0("/home/corroyez/Documents/NC_Full/PROSAIL-Optimization/01_DATA/", site, "/LiDAR/")
  pai_dir <- paste0(data_dir, dsm)
  gpkg_dir <- paste0("/home/corroyez/Documents/NC_Full/PROSAIL-Optimization/03_RESULTS/", 
                     site, "/PROSAIL_Optimization/sampling")
  atbd_r <- rast(file.path(dir, atbd))
  common_r <- rast(file.path(dir, common))
  best_indiv_r <- rast(file.path(dir, best_indiv))
  lidar_r <- sum(rast(file.path(data_dir,
                                'PAD_Profiles_Classic',
                                'ladstack.tif')), na.rm = T)
  h <- switch(site, "Aigoual" = "35.5", "Blois" = "35.5", "Mormal" = "31.5")
  pai_r <- rast(file.path(pai_dir, paste0("PAD_", h, "_40.tif")))
  max_r <- rast(file.path(data_dir, "max_res_10_m.tif"))
  
  atbd_v <- values(atbd_r)
  common_v <- values(common_r)
  best_indiv_v <- values(best_indiv_r)
  lidar_v <- values(lidar_r)
  pai_v <- values(pai_r)
  max_v <- values(max_r)
  
  # Full site
  r_stack <- cbind(atbd_v, common_v, best_indiv_v, lidar_v, pai_v, max_v)
  valid_idx <- complete.cases(r_stack)
  atbd_r_clean <- atbd_v[valid_idx]
  common_r_clean <- common_v[valid_idx]
  best_indiv_r_clean <- best_indiv_v[valid_idx]
  lidar_r_clean <- lidar_v[valid_idx]
  pai_clean <- pai_v[valid_idx]
  max_clean <- max_v[valid_idx]
  
  # Output filenames
  # out_dir <- file.path("output_figures_prosail", site, "density_scatter")
  out_dir <- file.path("output_figures_prosail", site, "scatter")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # create_scatterplot(atbd_r_clean, lidar_r_clean,
  #                    title = paste("ATBD vs LiDAR -", site),
  #                    filename = file.path(out_dir, paste0(site, "_scatter_atbd.png")))
  # create_scatterplot(common_r_clean, lidar_r_clean,
  #                    title = paste("Common vs LiDAR -", site),
  #                    filename = file.path(out_dir, paste0(site, "_scatter_common.png")))
  # create_scatterplot(best_indiv_r_clean, lidar_r_clean,
  #                    title = paste("Best Indiv vs LiDAR -", site),
  #                    filename = file.path(out_dir, paste0(site, "_scatter_best_indiv.png")))
  # 
  # create_scatterplot(atbd_r_clean, pai_clean,
  #                    title = paste("ATBD vs PAI 6m -", site),
  #                    filename = file.path(out_dir, paste0("pai_6m", site, "_scatter_atbd.png")))
  # create_scatterplot(common_r_clean, pai_clean,
  #                    title = paste("Common vs PAI 6m -", site),
  #                    filename = file.path(out_dir, paste0("pai_6m", site, "_scatter_common.png")))
  # create_scatterplot(best_indiv_r_clean, pai_clean,
  #                    title = paste("Best Indiv vs PAI 6m -", site),
  #                    filename = file.path(out_dir, paste0("pai_6m", site, "_scatter_best_indiv.png")))
  
  a <- plot_lai_scatter(lidar_r_clean, atbd_r_clean, site,
                        xlab = expression(LAI[ALS]), ylab = expression(LAI[S2_ATBD]),
                        xlim = c(0, 15), ylim = c(0, 15), title = site,
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = T, remove_y = T)
  e <- plot_lai_scatter(lidar_r_clean, best_indiv_r_clean, site,
                        xlab = expression(LAI[ALS]), ylab = expression(LAI[S2_opt]),
                        xlim = c(0, 15), ylim = c(0, 15),
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = F, remove_y = T)
  
  b <- plot_lai_scatter(pai_clean, atbd_r_clean, site,
                        xlab = expression(LAI[ALS_dopt]), ylab = expression(LAI[S2_ATBD]),
                        xlim = c(0, 15), ylim = c(0, 15), title = site,
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = T, remove_y = T)
  f <- plot_lai_scatter(pai_clean, best_indiv_r_clean, site,
                        xlab = expression(LAI[ALS_dopt]), ylab = expression(LAI[S2_opt]),
                        xlim = c(0, 15), ylim = c(0, 15),
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = F, remove_y = T)
  
  c <- plot_lai_scatter(lidar_r_clean, common_r_clean, site,
                        xlab = expression(LAI[ALS]), ylab = expression(LAI[S2_opt3s]),
                        xlim = c(0, 15), ylim = c(0, 15),
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = T, remove_y = T)
  d <- plot_lai_scatter(pai_clean, common_r_clean, site,
                        xlab = expression(LAI[ALS_dopt]), ylab = expression(LAI[S2_opt3s]),
                        xlim = c(0, 15), ylim = c(0, 15),
                        add_annotation = TRUE, annotation_pos = "topleft",
                        remove_x = T, remove_y = T)
  # final <- (a | b) / (c | d)
  # stop()
  final <- (a | c | e) / (b | d |f)
  print(final)
  # ggsave(paste0(site, ".png"), width = 10, height = 10, dpi = 300)
  ggsave(paste0(site, ".png"), width = 12, height = 8, unit = "in")
  
  init1 <- plot_lai_scatter(lidar_r_clean, atbd_r_clean, site,
                            xlab = expression(LAI[ALS]), ylab = expression(LAI[S2_ATBD]),
                            xlim = c(0, 15), ylim = c(0, 15), title = site,
                            add_annotation = TRUE, annotation_pos = "topleft",
                            remove_y = F)
  init1_opt <- plot_lai_scatter(pai_clean, atbd_r_clean, site,
                                xlab = expression(LAI[ALS_dopt]), ylab = expression(LAI[S2_ATBD]),
                                xlim = c(0, 15), ylim = c(0, 15), title = site,
                                add_annotation = TRUE, annotation_pos = "topleft",
                                remove_y = F)
  init2 <- plot_lai_scatter(lidar_r_clean, atbd_r_clean, site,
                            xlab = expression(LAI[ALS]), ylab = expression(LAI[S2_ATBD]),
                            xlim = c(0, 15), ylim = c(0, 15), title = site,
                            add_annotation = TRUE, annotation_pos = "topleft",
                            remove_y = TRUE)
  i <- switch(site, "Aigoual" = init1, "Blois" = init2, "Mormal" = init2)
  list_plots[[site]] <- i
  # stop()
  
  # Sampled
  # gpkg <- vect(file.path(gpkg_dir, "Above20_Sampling_stratified_uniform_nbSamples_5000.GPKG"))
  # r_stack <- c(atbd_r, common_r, best_indiv_r, lidar_r, pai_r, max_r)
  # names(r_stack) <- c("ATBD", "Common", "Best_Indiv", "LiDAR", "PAI", "MaxH")
  # sampled_vals <- terra::extract(r_stack, gpkg)
  # sampled_vals <- sampled_vals[complete.cases(sampled_vals), -1]
  # 
  # # Plot using sampled values
  # create_scatterplot(sampled_vals$ATBD, sampled_vals$LiDAR,
  #                            title = paste("ATBD vs LiDAR -", site),
  #                            filename = file.path(out_dir, paste0(site, "_scatter_atbd_sampled.png")))
  # 
  # create_scatterplot(sampled_vals$Common, sampled_vals$LiDAR,
  #                            title = paste("Common vs LiDAR -", site),
  #                            filename = file.path(out_dir, paste0(site, "_scatter_common_sampled.png")))
  # 
  # create_scatterplot(sampled_vals$Best_Indiv, sampled_vals$LiDAR,
  #                            title = paste("Best Indiv vs LiDAR -", site),
  #                            filename = file.path(out_dir, paste0(site, "_scatter_best_indiv_sampled.png")))
  # 
  # create_scatterplot(sampled_vals$ATBD, sampled_vals$PAI,
  #                            title = paste("ATBD vs PAI 6m -", site),
  #                            filename = file.path(out_dir, paste0("pai_6m_", site, "_scatter_atbd_sampled.png")))
  # 
  # create_scatterplot(sampled_vals$Common, sampled_vals$PAI,
  #                            title = paste("Common vs PAI 6m -", site),
  #                            filename = file.path(out_dir, paste0("pai_6m_", site, "_scatter_common_sampled.png")))
  # 
  # create_scatterplot(sampled_vals$Best_Indiv, sampled_vals$PAI,
  #                            title = paste("Best Indiv vs PAI 6m -", site),
  #                            filename = file.path(out_dir, paste0("pai_6m_", site, "_scatter_best_indiv_sampled.png")))
  
  # Filtered height
  # height_bins <- c(0, 6, 10, 15, 20, 30)
  # for (i in 1:(length(height_bins) - 1)) {
  #   h_min <- height_bins[i]
  #   h <- height_bins[i + 1]
  #   
  #   r_vals <- values(r_stack)
  #   valid_idx <- complete.cases(r_vals) &
  #     r_vals[, "MaxH"] > h_min &
  #     r_vals[, "MaxH"] <= h
  #   sampled_vals <- as.data.frame(r_vals[valid_idx, ])
  #   
  #   out_dir <- file.path("output_figures_prosail", site, "height_study", 
  #                        paste0(h, "m"))
  #   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  #   
  #   create_density_scatterplot(sampled_vals$ATBD, sampled_vals$PAI,
  #                              title = paste("ATBD vs PAI", h, "m -", site),
  #                              filename = file.path(out_dir, paste0("pai_6m_",
  #                                                                   "max_", h, "m_",
  #                                                                   site, "_scatter_atbd_sampled.png")))
  #   
  #   create_density_scatterplot(sampled_vals$Common, sampled_vals$PAI,
  #                              title = paste("Common vs PAI", h, "m -", site),
  #                              filename = file.path(out_dir, paste0("pai_6m_",
  #                                                                   "max_", h, "m_",
  #                                                                   site, "_scatter_common_sampled.png")))
  #   
  #   create_density_scatterplot(sampled_vals$Best_Indiv, sampled_vals$PAI,
  #                              title = paste("Best Indiv vs PAI", h, "m -", site),
  #                              filename = file.path(out_dir, paste0("pai_6m_",
  #                                                                   "max_", h, "m_",
  #                                                                   site, "_scatter_best_indiv_sampled.png")))
  # }
  # store
  plots1 <- c(plots1, list(init1, c, e))
  plots2 <- c(plots2, list(init1_opt, d, f))
  init1_list[[site]]     <- a
  c_list[[site]]         <- c
  e_list[[site]]         <- e
  init1_opt_list[[site]] <- b
  d_list[[site]]         <- d
  f_list[[site]]         <- f
}
init <- (list_plots$Aigoual + theme(axis.title.x = element_blank())) | 
  (list_plots$Blois) | 
  (list_plots$Mormal + theme(axis.title.x = element_blank()))
print(init)
ggsave("init.png", width = 12, height = 4, unit = "in")
# ggsave("init.png", width = 10, height = 4, unit = "in")

fig1 <- (init1_list$Aigoual | init1_list$Blois | init1_list$Mormal) /
  (c_list$Aigoual     | c_list$Blois     | c_list$Mormal) /
  (e_list$Aigoual + theme(axis.title.x = element_blank()) | e_list$Blois | e_list$Mormal + theme(axis.title.x = element_blank()))
fig2 <- (init1_opt_list$Aigoual | init1_opt_list$Blois | init1_opt_list$Mormal) /
  (d_list$Aigoual     | d_list$Blois     | d_list$Mormal) /
  (f_list$Aigoual + theme(axis.title.x = element_blank()) | f_list$Blois | f_list$Mormal + theme(axis.title.x = element_blank()))
print(fig1)
print(fig2)

# ------------------------------------------------------------------------------
for (site in sites){
  dir <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/", 
                site, "/Metrics/Deciduous_Only")
  lcv <- file.path(dir, 'lcv_no_ground.tif')
  
}

# ------------------------------------------------------------------------------
csv_all <- fread(file.path(output_dir, paste0("all_results_combined_",
                                              name_strategy, "16_04.csv")))
csv_to_study <- csv_all
TOP_N <- 100

# Filter for Norm == "DSM" and Depth == 6, and compute score
data_all <- csv_to_study %>%
  mutate(score_r_rmse_slope = (normalize(R) + (1 - normalize(RMSE)) + (1 - normalize(abs(Slope - 1)))) / 3) %>%
  mutate(score_r_rmse = (normalize(R) + (1 - normalize(RMSE))) / 2) %>%
  # mutate(score = R) %>%
  # filter(Norm == "DSM_keepTrees")
  filter(Depth == 6)

# Top N per score for each site
top_configs_r_rmse_slope <- data_all %>%
  group_by(Site) %>%
  arrange(Site, desc(score_r_rmse_slope)) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  mutate(criteria = "R_RMSE_Slope")

top_configs_r_rmse <- data_all %>%
  group_by(Site) %>%
  arrange(Site, desc(score_r_rmse)) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  mutate(criteria = "R_RMSE")

top_configs_R <- data_all %>%
  group_by(Site) %>%
  arrange(Site, desc(R)) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  mutate(criteria = "R")

top_configs_RMSE <- data_all %>%
  group_by(Site) %>%
  arrange(Site, RMSE) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  mutate(criteria = "RMSE")

top_configs_Slope <- data_all %>%
  group_by(Site) %>%
  arrange(Site, abs(Slope - 1)) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  mutate(criteria = "Slope")

# Combine all top N results into one dataset
data_range <- bind_rows(top_configs_r_rmse_slope, top_configs_r_rmse,
                        top_configs_R, top_configs_RMSE, top_configs_Slope)

# Create the plots for each score type
p_score_r_rmse <- plot_score_ranges(top_configs_r_rmse, 
                                    "score_r_rmse", "Weighted Norm (R, RMSE)")
p_score_r_rmse_slope <- plot_score_ranges(top_configs_r_rmse_slope, 
                                          "score_r_rmse_slope", "Weighted Norm (R, RMSE, Slope)")
p_R <- plot_score_ranges(top_configs_R, "R", "R")
p_RMSE <- plot_score_ranges(top_configs_RMSE, "RMSE", "RMSE")
p_Slope <- plot_score_ranges(top_configs_Slope, "Slope", "Slope")

# Arrange the plots in a grid
layout <- (p_score_r_rmse + p_score_r_rmse_slope) / (p_R + p_RMSE + p_Slope)
print(layout)

# Find top configurations per site and all
# top_configs_Aigoual <- data_all %>%
#   filter(Site == "Aigoual") %>%
#   arrange(desc(R)) %>%
#   slice_head(n = TOP_N) %>%
#   pull(Column)
# 
# top_configs_Blois <- data_all %>%
#   filter(Site == "Blois") %>%
#   arrange(desc(R)) %>%
#   slice_head(n = TOP_N) %>%
#   pull(Column)
# 
# top_configs_Mormal <- data_all %>%
#   filter(Site == "Mormal") %>%
#   arrange(desc(R)) %>%
#   slice_head(n = TOP_N) %>%
#   pull(Column)

# top_configs_overall <- c(top_configs_Aigoual, top_configs_Blois, top_configs_Mormal)

# Filter data for each case
# data_Aigoual <- data_all %>% filter(Column %in% top_configs_Aigoual & Site == "Aigoual")
# data_Blois   <- data_all %>% filter(Column %in% top_configs_Blois & Site == "Blois")
# data_Mormal  <- data_all %>% filter(Column %in% top_configs_Mormal & Site == "Mormal")
# data_overall <- data_all %>% filter(Column %in% top_configs_overall)
data_Aigoual <- data_range %>%
  filter(Site == "Aigoual")
data_Blois <- data_range %>%
  filter(Site == "Blois")
data_Mormal <- data_range %>%
  filter(Site == "Mormal")
data_overall <- data_range

# Create plots for each site and the combined dataset
p_Aigoual <- plot_config(data_Aigoual, "Aigoual")
p_Blois   <- plot_config(data_Blois, "Blois")
p_Mormal  <- plot_config(data_Mormal, "Mormal")
p_All     <- plot_config(data_overall, "All Sites")

# Arrange plots in a 2x2 grid
# layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_All)
# print(layout)
# stop()
# ------------------------------------------------------------------------------

# Identify non-unique Column values
# non_unique_columns <- data_overall[duplicated(data_overall$Column) | duplicated(data_overall$Column, fromLast = TRUE), ]

# Get all possible Site values
# all_sites <- unique(data_overall$Site)

# For each non-unique Column, check if all Site values are present
# valid_columns <- non_unique_columns %>%
#   group_by(Column) %>%
#   filter(all(all_sites %in% Site)) %>%
#   ungroup()

# ------------------------------------------------------------------------------
score_columns <- unique(data_range$criteria)

# Initialize an empty list to store the parameter frequencies
param_freq_list <- list()

# Compute parameter frequencies for each score column
for (score_value in score_columns) {
  param_freq_list[[score_value]] <- list(
    Aigoual = compute_param_freq_by_score(data_Aigoual, score_value, all_params),
    Blois = compute_param_freq_by_score(data_Blois, score_value, all_params),
    Mormal = compute_param_freq_by_score(data_Mormal, score_value, all_params),
    Overall = compute_param_freq_by_score(data_overall, score_value, all_params)
  )
}

# Generate plots for each score column
plot_list <- list()
for (score_value in score_columns) {
  param_freq_Aigoual <- param_freq_list[[score_value]]$Aigoual
  param_freq_Blois   <- param_freq_list[[score_value]]$Blois
  param_freq_Mormal  <- param_freq_list[[score_value]]$Mormal
  param_freq_overall <- param_freq_list[[score_value]]$Overall
  
  p_Aigoual <- plot_param_freq(param_freq_Aigoual, param_colors, paste("Aigoual -", score_value))
  p_Blois   <- plot_param_freq(param_freq_Blois, param_colors, paste("Blois -", score_value))
  p_Mormal  <- plot_param_freq(param_freq_Mormal, param_colors, paste("Mormal -", score_value))
  p_Overall <- plot_param_freq(param_freq_overall, param_colors, paste("Overall -", score_value))
  
  layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_Overall)
  plot_list[[score_value]] <- layout
}

# Print all plots
for (plot_name in names(plot_list)) {
  p <- plot_list[[plot_name]]
  print(p)
}

# Create an empty data frame to hold the stacked data
param_freq_combined <- bind_rows(
  lapply(score_columns, function(score_value) {
    bind_rows(
      mutate(param_freq_list[[score_value]]$Aigoual, Site = "Aigoual", Score = score_value),
      mutate(param_freq_list[[score_value]]$Blois, Site = "Blois", Score = score_value),
      mutate(param_freq_list[[score_value]]$Mormal, Site = "Mormal", Score = score_value),
      mutate(param_freq_list[[score_value]]$Overall, Site = "Overall", Score = score_value)
    )
  })
)

# Plot
ggplot(param_freq_combined, aes(x = param_list, y = N, fill = param_list)) +
  geom_bar(stat = "identity", position = "dodge") +  # Position bars side by side
  facet_grid(Site ~ Score, scales = "free_y") +  # Facet by Site and Score combinations
  scale_fill_manual(values = param_colors) +  # Apply custom colors
  theme_bw() +
  labs(title = "Parameter Frequencies by Site and Score", x = "Parameter", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",  # Hide legend
        strip.text = element_text(size = 10))  # Customize facet labels

duplicated_columns <- data_range %>%
  group_by(Column) %>%
  filter(n_distinct(Site) == 3) %>%
  ungroup()
duplicated_columns
unique(duplicated_columns$Column)

# Keep unique rows excluding the 'criteria' variable
unique_rows <- duplicated_columns %>%
  select(-criteria) %>%
  distinct()
unique_rows

avg <- unique_rows %>%
  group_by(Column) %>%
  mutate(avg_R = mean(R, na.rm = TRUE),
         avg_RMSE = mean(RMSE, na.rm = TRUE),
         avg_Slope = mean(Slope, na.rm = TRUE),
         avg_R_RMSE = mean(score_r_rmse, na.rm = TRUE),
         avg_R_RMSE_Slope = mean(score_r_rmse_slope, na.rm = TRUE)) %>%
  select(Column, avg_R, avg_RMSE, avg_Slope, avg_R_RMSE, avg_R_RMSE_Slope) %>%
  distinct()
boxplot(avg$avg_R_RMSE)

# Bests - setup
order_descending <- c(TRUE, FALSE, TRUE, TRUE, TRUE)

# Best per criteria per site
criteria_columns <- c("R", "RMSE", "Slope", "R_RMSE", "R_RMSE_Slope")
best_columns_per_site <- list()
for (site in unique(data_range$Site)) {
  site_data <- data_range %>% filter(Site == site)
  site_best_columns <- list()
  for (i in 1:length(criteria_columns)) {
    column_name <- criteria_columns[i]
    is_descending <- order_descending[i]
    site_criteria_data <- site_data %>% filter(criteria == column_name)
    if (is_descending) {
      site_best_columns[[column_name]] <- site_criteria_data[order(-site_criteria_data[[column_name]]), ][1, ]
    } else {
      site_best_columns[[column_name]] <- site_criteria_data[order(site_criteria_data[[column_name]]), ][1, ]
    }
  }
  best_columns_per_site[[site]] <- site_best_columns
}
best_columns_per_site

# Best per criteria overall
columns <- c("avg_R", "avg_RMSE", "avg_Slope", "avg_R_RMSE", "avg_R_RMSE_Slope")
best_columns <- list()
for (i in 1:length(columns)) {
  column_name <- columns[i]
  is_descending <- order_descending[i]
  if (is_descending) {
    best_columns[[column_name]] <- avg[order(-avg[[column_name]]), ][1, ]
  } else {
    best_columns[[column_name]] <- avg[order(avg[[column_name]]), ][1, ]
  }
}
best_columns

# ------------------------------------------------------------------------------
# all_params <- unique(unlist(str_split(data_overall$Column, "_")))
all_params <- c("LIDFa=1", "LIDFa=2", "LIDFa=3", "LIDFa=4", "LIDFa=5",
                "lai=1", "lai=2", "lai=3", "lai=4", "lai=5", "lai=6",
                "BROWN=1", "BROWN=2", "BROWN=3",
                "LMA=1", "LMA=2", "LMA=3", "LMA=4", "LMA=5",
                "N=1", "N=2",
                "CHL=1", "CHL=2",
                "psoil=1", "psoil=2",
                "q=1", "q=2", "q=3", "q=4")
param_freq_Aigoual <- compute_param_freq(data_Aigoual, all_params)
param_freq_Blois   <- compute_param_freq(data_Blois, all_params)
param_freq_Mormal  <- compute_param_freq(data_Mormal, all_params)
param_freq_overall <- compute_param_freq(data_overall, all_params)

param_colors <- generate_param_colors(all_params)

# Generate plots
p_Aigoual <- plot_param_freq(param_freq_Aigoual, param_colors, "Aigoual")
p_Blois   <- plot_param_freq(param_freq_Blois, param_colors, "Blois")
p_Mormal  <- plot_param_freq(param_freq_Mormal, param_colors, "Mormal")
p_Overall <- plot_param_freq(param_freq_overall, param_colors, "Overall")
layout <- (p_Aigoual + p_Blois) / (p_Mormal + p_Overall)
print(layout)

# ------------------------------------------------------------------------------
# predictors <- as.formula(paste(param_names, collapse = " + "))
rf_R <- randomForest(R ~ LIDFa + lai + LMA + BROWN + N + CHL + q + psoil, 
                     data = data_extended, importance = TRUE)
rf_NRMSE <- randomForest(NRMSE ~ LIDFa + lai + LMA + BROWN + N + CHL + q + psoil, 
                         data = data_extended, importance = TRUE)
rf_Slope <- randomForest(Slope ~ LIDFa + lai + LMA + BROWN + N + CHL + q + psoil, 
                         data = data_extended, importance = TRUE)
rf_Intercept <- randomForest(Intercept ~ LIDFa + lai + LMA + BROWN + N + CHL + q + psoil,
                             data = data_extended, importance = TRUE)
# rf_Bias <- randomForest(Bias ~ LIDFa + lai + LMA + BROWN + N + CHL + q + psoil, 
#                         data = data_extended, importance = TRUE)
importance(rf_R)

pdp_LIDFa <- partial(rf_NRMSE, pred.var = "BROWN")
plotPartial(pdp_LIDFa)

# ------------------------------------------------------------------------------
# Generate partial dependence plots for each metric and parameter
pdp_list <- list()

# For each metric (R, NRMSE, Slope, Intercept), generate PDPs for each parameter
metrics <- c("R", "NRMSE", "Slope", "Intercept")
param_names <- c("LIDFa", "lai", "LMA", "BROWN", "N", "CHL", "psoil", "q")

# Loop to generate PDPs for all parameters and metrics
for (metric in metrics) {
  print(metric)
  for (param in param_names) {
    rf_model <- randomForest(as.formula(paste(metric, "~", paste(param_names, collapse = " + "))),
                             data = data_extended, importance = TRUE)
    
    pdp_result <- partial(rf_model, pred.var = param, train = data_extended)
    pdp_result$Metric <- metric
    pdp_result$Parameter <- param
    pdp_list[[paste(metric, param, sep = "_")]] <- pdp_result
  }
}

# Combine all PDP results into one data frame
pdp_combined <- bind_rows(pdp_list)

# Plot: Facet 2x4 (or 4x2) grid where each plot is a parameter
ggplot(pdp_combined, aes(x = x, y = y, color = Metric)) +
  geom_line() +
  facet_wrap(~ Parameter, scales = "free", ncol = 4) +
  labs(title = "Impact of Each Parameter on Metrics", x = "Parameter Value", y = "Metric Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------------------------------------------------------------
# Make sure your data frame is in long format
long_df <- pdp_combined %>%
  filter(!Metric == "Intercept") %>%
  pivot_longer(cols = c(LIDFa, lai, LMA, BROWN, N, CHL, psoil, q),
               names_to = "Parameter_name",
               values_to = "Parameter_value")

# Plot 1: Facet 2x4 where each plot is a parameter
ggplot(long_df, aes(x = Parameter_value, y = yhat, color = Metric)) +
  geom_line() +
  facet_wrap(~ Parameter_name, scales = "free", ncol = 4) +
  theme_minimal() +
  ggtitle("Metrics as a function of Parameter Number")

# Plot 2: Facet 2x2 where each plot is a metric
ggplot(long_df, aes(x = Parameter_value, y = yhat, color = Parameter_name)) +
  geom_line() +
  facet_wrap(~ Metric, scales = "free") +
  theme_minimal() +
  ggtitle("Parameters as a function of Metric")
































# Pareto Front Analysis for Aigoual Dataset
# Using R, RMSE, Slope, and Bias as criteria

# Load required libraries
library(data.table)
library(ggplot2)
library(plotly)
library(dplyr)

# Function to determine if a point is dominated by another
is_dominated <- function(point, other_points, maximize = c(), minimize = c()) {
  # Check if 'point' is dominated by any point in 'other_points'
  for (i in 1:nrow(other_points)) {
    other_point <- other_points[i, ]
    
    # Check domination conditions
    dominates <- TRUE
    at_least_one_better <- FALSE
    
    # For maximization criteria (higher is better)
    for (col in maximize) {
      if (point[[col]] > other_point[[col]]) {
        dominates <- FALSE
        break
      } else if (point[[col]] < other_point[[col]]) {
        at_least_one_better <- TRUE
      }
    }
    
    # For minimization criteria (lower is better)
    for (col in minimize) {
      if (point[[col]] < other_point[[col]]) {
        dominates <- FALSE
        break
      } else if (point[[col]] > other_point[[col]]) {
        at_least_one_better <- TRUE
      }
    }
    
    # If this point dominates and is better in at least one criterion
    if (dominates && at_least_one_better) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Function to find Pareto Front
find_pareto_front <- function(data, maximize = c(), minimize = c()) {
  pareto_indices <- c()
  
  for (i in 1:nrow(data)) {
    current_point <- data[i, ]
    other_points <- data[-i, ]
    
    if (!is_dominated(current_point, other_points, maximize, minimize)) {
      pareto_indices <- c(pareto_indices, i)
    }
  }
  
  return(data[pareto_indices, ])
}

# Define optimization criteria
# R: maximize (higher correlation is better)
# RMSE: minimize (lower error is better)
# Bias: minimize absolute value (closer to 0 is better)
# Slope: minimize absolute deviation from 1 (closer to 1 is better)

# Prepare data for Pareto analysis
data_for_pareto <- data_Aigoual[, .(
  R = R,
  RMSE = RMSE,
  Slope = Slope,
  Bias = Bias,
  Bias_abs = abs(Bias),
  Slope_dev = abs(Slope - 1),  # Deviation from ideal slope of 1
  # Keep other relevant columns for identification
  Site, Norm, Depth, Method, Column, LIDFa, lai, LMA, BROWN, N, CHL, psoil, q,
  # Keep all original metrics
  R2, NRMSE, Intercept, MAE, ATBD, norm_R, norm_NRMSE, norm_slope, score
)]

# Find Pareto Front
pareto_front <- find_pareto_front(
  data_for_pareto,
  maximize = c("R"),
  minimize = c("RMSE", "Bias_abs", "Slope_dev")
)

# Print results
cat("Pareto Front Analysis Results:\n")
cat("=============================\n")
cat("Total number of solutions:", nrow(data_for_pareto), "\n")
cat("Number of Pareto optimal solutions:", nrow(pareto_front), "\n")
cat("Percentage of Pareto optimal solutions:", round(nrow(pareto_front)/nrow(data_for_pareto)*100, 2), "%\n\n")

# Display Pareto optimal solutions
cat("Pareto Optimal Solutions:\n")
cat("========================\n")
pareto_summary <- pareto_front[, .(
  Configuration = Column,
  R = round(R, 4),
  RMSE = round(RMSE, 4),
  Slope = round(Slope, 4),
  Slope_dev = round(Slope_dev, 4),
  Bias = round(Bias, 4),
  Bias_abs = round(Bias_abs, 4),
  LIDFa, lai, LMA, BROWN, N, CHL, psoil, q
)]

print(pareto_summary)

# Create visualization function
create_pareto_plots <- function(data, pareto_data) {
  # Basic 2D plots for pairwise comparisons
  
  # 1. R vs RMSE
  p1 <- ggplot(data, aes(x = RMSE, y = R)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = RMSE, y = R), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: R vs RMSE",
         x = "RMSE (minimize)", y = "R (maximize)") +
    theme_minimal()
  
  # 2. R vs |Bias|
  p2 <- ggplot(data, aes(x = Bias_abs, y = R)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = Bias_abs, y = R), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: R vs |Bias|",
         x = "|Bias| (minimize)", y = "R (maximize)") +
    theme_minimal()
  
  # 3. R vs |Slope - 1|
  p3 <- ggplot(data, aes(x = Slope_dev, y = R)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = Slope_dev, y = R), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: R vs |Slope - 1|",
         x = "|Slope - 1| (minimize)", y = "R (maximize)") +
    theme_minimal()
  
  # 4. RMSE vs |Bias|
  p4 <- ggplot(data, aes(x = Bias_abs, y = RMSE)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = Bias_abs, y = RMSE), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: RMSE vs |Bias|",
         x = "|Bias| (minimize)", y = "RMSE (minimize)") +
    theme_minimal()
  
  # 5. RMSE vs |Slope - 1|
  p5 <- ggplot(data, aes(x = Slope_dev, y = RMSE)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = Slope_dev, y = RMSE), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: RMSE vs |Slope - 1|",
         x = "|Slope - 1| (minimize)", y = "RMSE (minimize)") +
    theme_minimal()
  
  # 6. |Bias| vs |Slope - 1|
  p6 <- ggplot(data, aes(x = Bias_abs, y = Slope_dev)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = Bias_abs, y = Slope_dev), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: |Bias| vs |Slope - 1|",
         x = "|Bias| (minimize)", y = "|Slope - 1| (minimize)") +
    theme_minimal()
  
  # 7. Composite Score Plot (R vs combined error metrics)
  data$combined_error <- sqrt(data$RMSE^2 + data$Bias_abs^2 + data$Slope_dev^2)
  pareto_data$combined_error <- sqrt(pareto_data$RMSE^2 + pareto_data$Bias_abs^2 + pareto_data$Slope_dev^2)
  
  p7 <- ggplot(data, aes(x = combined_error, y = R)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_point(data = pareto_data, aes(x = combined_error, y = R), 
               color = "red", size = 3, alpha = 0.8) +
    labs(title = "Pareto Front: R vs Combined Error",
         x = "Combined Error (minimize)", y = "R (maximize)") +
    theme_minimal()
  
  return(list(p1, p2, p3, p4, p5, p6, p7))
}

# Create plots
plots <- create_pareto_plots(data_for_pareto, pareto_front)

# Display plots
print(plots[[1]])  # R vs RMSE
print(plots[[2]])  # R vs |Bias|
print(plots[[3]])  # R vs |Slope - 1|
print(plots[[4]])  # RMSE vs |Bias|
print(plots[[5]])  # RMSE vs |Slope - 1|
print(plots[[6]])  # |Bias| vs |Slope - 1|
print(plots[[7]])  # R vs Combined Error

# Additional advanced plots
create_advanced_plots <- function(data, pareto_data) {
  # 8. Parallel Coordinates Plot
  library(GGally)
  
  # Prepare data for parallel coordinates
  plot_data <- rbind(
    data[, .(R, RMSE, Bias_abs, Slope_dev, Type = "All Solutions")],
    pareto_data[, .(R, RMSE, Bias_abs, Slope_dev, Type = "Pareto Front")]
  )
  
  p8 <- ggparcoord(plot_data, columns = 1:4, groupColumn = "Type",
                   scale = "globalminmax", alphaLines = 0.3) +
    scale_color_manual(values = c("All Solutions" = "lightblue", "Pareto Front" = "red")) +
    labs(title = "Parallel Coordinates: Pareto Front Analysis",
         x = "Criteria", y = "Normalized Values") +
    theme_minimal()
  
  # 9. Radar Chart for top Pareto solutions
  library(fmsb)
  
  # Select top 5 Pareto solutions by R value
  top_pareto <- pareto_data[order(-R)][1:min(5, nrow(pareto_data))]
  
  # 10. Heatmap of criteria correlations
  criteria_data <- data[, .(R, RMSE, Bias_abs, Slope_dev)]
  correlation_matrix <- cor(criteria_data)
  
  library(reshape2)
  corr_melted <- melt(correlation_matrix)
  
  p10 <- ggplot(corr_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation") +
    labs(title = "Correlation Matrix of Optimization Criteria") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  return(list(p8, p10))
}

# Create advanced plots
advanced_plots <- create_advanced_plots(data_for_pareto, pareto_front)
print(advanced_plots[[1]])  # Parallel coordinates
print(advanced_plots[[2]])  # Correlation heatmap

# Alternative approach: Rank-based Pareto Front
# This approach ranks each criterion and finds solutions that are not dominated
rank_based_pareto <- function(data) {
  # Create ranks (1 = best)
  data_ranked <- data[, .(
    rank_R = rank(-R),           # Higher R is better, so negative rank
    rank_RMSE = rank(RMSE),      # Lower RMSE is better
    rank_Bias = rank(Bias_abs),  # Lower |Bias| is better
    rank_Slope = rank(Slope_dev), # Lower |Slope - 1| is better
    # Keep original data
    R, RMSE, Slope, Bias, Bias_abs, Slope_dev,
    Site, Norm, Depth, Method, Column, LIDFa, lai, LMA, BROWN, N, CHL, psoil, q
  )]
  
  # Find Pareto front using ranks
  pareto_rank <- find_pareto_front(
    data_ranked,
    maximize = c(),
    minimize = c("rank_R", "rank_RMSE", "rank_Bias", "rank_Slope")
  )
  
  return(pareto_rank)
}

# Calculate rank-based Pareto front
cat("\n\nRank-based Pareto Front Analysis:\n")
cat("==================================\n")
pareto_rank <- rank_based_pareto(data_for_pareto)
cat("Number of rank-based Pareto optimal solutions:", nrow(pareto_rank), "\n")

# Summary statistics of Pareto front
cat("\n\nSummary Statistics of Pareto Front:\n")
cat("===================================\n")
summary_stats <- pareto_front[, .(
  R_mean = mean(R),
  R_range = paste(round(min(R), 4), "-", round(max(R), 4)),
  RMSE_mean = mean(RMSE),
  RMSE_range = paste(round(min(RMSE), 4), "-", round(max(RMSE), 4)),
  Slope_mean = mean(Slope),
  Slope_range = paste(round(min(Slope), 4), "-", round(max(Slope), 4)),
  Slope_dev_mean = mean(Slope_dev),
  Slope_dev_range = paste(round(min(Slope_dev), 4), "-", round(max(Slope_dev), 4)),
  Bias_mean = mean(abs(Bias)),
  Bias_range = paste(round(min(abs(Bias)), 4), "-", round(max(abs(Bias)), 4))
)]

print(summary_stats)

# Parameter analysis of Pareto solutions
cat("\n\nParameter Analysis of Pareto Solutions:\n")
cat("=======================================\n")
param_analysis <- pareto_front[, .(
  LIDFa_freq = .N
), by = .(LIDFa)][order(-LIDFa_freq)]

cat("LIDFa distribution in Pareto front:\n")
print(param_analysis)

lai_analysis <- pareto_front[, .(
  lai_freq = .N
), by = .(lai)][order(-lai_freq)]

cat("\nLAI distribution in Pareto front:\n")
print(lai_analysis)

# Save results
# write.csv(pareto_front, "pareto_front_solutions.csv", row.names = FALSE)
# cat("\nPareto front solutions saved to 'pareto_front_solutions.csv'\n")

brown_analysis <- pareto_front[, .(
  brown_freq = .N
), by = .(BROWN)][order(-brown_freq)]

cat("\nBROWN distribution in Pareto front:\n")
print(brown_analysis)

q_analysis <- pareto_front[, .(
  q_freq = .N
), by = .(q)][order(-q_freq)]

cat("\nq distribution in Pareto front:\n")
print(q_analysis)

# boxplot cumD / site

# front 1 et 2 ?