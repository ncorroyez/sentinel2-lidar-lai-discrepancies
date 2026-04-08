# clean environment
rm(list = ls(all=TRUE)); gc()
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))
lapply(X = as.list(list.files(path = 'libraries', full.names = T)), FUN = source, verbose = F)
library(terra)
library(dplyr)
library(ggplot2)
library(fmsb)
library(tidyr)
library(rlang)
library(data.table)
library(zoo)
library(purrr)

# 1- Directories & input data
set.seed(42)
data_dir <- '../01_DATA'
output_dir <- '../../03_RESULTS'
sites <- c('Aigoual', 'Blois', 'Mormal')
# sites <- c('Aigoual', 'Blois')
# sites <- c("Blois")
# sites <- c("Mormal")

results_list <- list()
reflectance_data_list <- list()
all_sampled_points <- all_final_df <- data.frame()
norms <- c("DTM", "DSM")
# norms <- c("DSM")
max_h_values <- c(10, 15, 20)
# max_h_values <- 0
n_samp <- 5000
normalize <- function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
# normalize_iqr <- function(x, q1, q3) {
#   (x - q1) / (q3 - q1)
# }

for (site in sites) {
  max_rast <- rast(file.path(data_dir, site, "LiDAR",
                             "max_res_10_m.tif"))
  mean_rast <- rast(file.path(data_dir, site, "LiDAR",
                              "mean_res_10_m.tif"))
  lcv_rast <- rast(file.path(data_dir, site, "LiDAR",
                             "lcv_no_ground.tif"))
  lskew_rast <- rast(file.path(data_dir, site, "LiDAR",
                               "lskew_res_10_m.tif"))
  vci_rast <- rast(file.path(data_dir, site, "LiDAR",
                             "vci_res_10_m.tif"))
  rumple_rast <- rast(file.path(data_dir, site, "LiDAR",
                                "rumple_res_10_m.tif"))
  fCover_rast <- rast(file.path(data_dir, site, "LiDAR",
                                "fCover_res_10_m.tif"))
  lidar_lai_rast <- rast(file.path(data_dir, site, 'LiDAR/PAD_Profiles_Classic',
                                   'ladstack.tif'))
  lidar_lai_rast <- sum(lidar_lai_rast, na.rm = T)
  s2_lai_rast <- rast(file.path(data_dir, site, "Sentinel-2",
                                "s2lai_summer_depth_study_res_10_m.tif"))
  
  refle <- switch(site,
                  "Aigoual" = rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Aigoual/L2A_T31TEJ_A031608_20210711T104217/Reflectance/res_10_m/L2A_T31TEJ_A031608_20210711T104217_Refl"),
                  "Blois" = rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31TCN_A031222_20210614T105443_Refl"),
                  "Mormal" = rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/L2A_T31UER_A031222_20210614T105443/Reflectance/res_10_m/L2A_T31UER_A031222_20210614T105443_Refl")
  )
  refl <- refle[[c(2, 3, 7)]]
  names(refl) <- c("B03", "B04", "B08")
  # stop()
  # lidar_values <- values(lidar_lai_rast)
  # s2_values <- values(s2_lai_rast)
  # lidar_values[is.nan(lidar_values)] <- NA
  # s2_values[is.nan(s2_values)] <- NA
  # aligned_values <- data.frame(lidar = lidar_values[, 1], s2 = s2_values[, 1])
  # aligned_values <- na.omit(aligned_values)
  # lidar_values <- aligned_values$lidar
  # s2_values <- aligned_values$s2
  
  # stop()
  
  # Stack aligned rasters for filtering
  combined_stack <- c(lidar_lai_rast, s2_lai_rast, max_rast, lcv_rast, mean_rast,
                      lskew_rast, vci_rast, rumple_rast, fCover_rast)
  names(combined_stack) <- c("lidar_lai", "s2_lai", "max_height", "lcv", "mean",
                             "lskew", "vci", "rumple", "fCover")
  full_data_mask <- !app(combined_stack, fun = function(...) any(is.na(c(...))))
  
  # Apply mask
  combined_stack <- mask(combined_stack, full_data_mask, maskvalues = FALSE)
  cor_LAI_total <- cor(values(combined_stack$lidar_lai, na.rm = T), 
                       values(combined_stack$s2_lai, na.rm = T))
  # stop()
  
  for (max_h in max_h_values) {
    # couper depth à max_h
    # p1 que 1 metric: A/ norm B/ non norm
    # norm sur ensemble h ?
    # scatterplot s2_lai = f(pai_depth) -> density scatt
    
    # Convert to data frame & filter where max > 20
    data_df <- as.data.frame(combined_stack, xy = TRUE, na.rm = TRUE)
    filtered_df <- data_df %>% dplyr::filter(max_height > max_h)
    # hist(filtered_df$lidar_lai, 100, xlim = c(0, 20))
    # stop()
    
    # Define breaks within the 5th-95th percentile range
    lai_5th <- floor(quantile(filtered_df$lidar_lai, probs = 0.02, na.rm = TRUE))
    lai_95th <- floor(quantile(filtered_df$lidar_lai, probs = 0.98, na.rm = TRUE))
    breaks <- seq(lai_5th, lai_95th, length.out = lai_95th - lai_5th + 1)
    
    sampled_points <- filtered_df %>%
      dplyr::filter(lidar_lai >= lai_5th & lidar_lai <= lai_95th) %>%
      mutate(lai_bin = cut(lidar_lai, breaks = breaks, include.lowest = TRUE)) %>%
      group_by(lai_bin) %>%
      mutate(bin_id = cur_group_id()) %>%
      ungroup()
    
    # Stratified sampling parameters
    bin_ids <- unique(sampled_points$bin_id)
    num_bins <- length(bin_ids)
    base_samples <- floor(n_samp / num_bins)
    extra_samples <- n_samp %% num_bins
    extra_bins <- sample(bin_ids, extra_samples)
    
    # Create sampling plan
    sampling_plan <- tibble(
      bin_id = bin_ids,
      n_samples = base_samples + as.integer(bin_ids %in% extra_bins)
    )
    
    # Join plan to points
    sampled_points <- sampled_points %>%
      left_join(sampling_plan, by = "bin_id") # %>%
    # group_by(bin_id) %>%
    # slice_sample(n = 833) %>%
    # ungroup()
    
    # Sample per bin using row_number()
    sampled_points <- sampled_points %>%
      group_by(bin_id) %>%
      arrange(runif(n())) %>%
      mutate(row_id = row_number()) %>%
      dplyr::filter(row_id <= first(n_samples)) %>%
      ungroup() %>%
      dplyr::select(-row_id, -n_samples) %>%
      mutate(max_h = max_h, site = site, cor_LAI_total = cor_LAI_total)
    
    # sampled_points <- sample_n(filtered_df, 5000)
    
    cor_LAI_max_h <- cor(sampled_points$lidar_lai, sampled_points$s2_lai)
    sampled_points <- sampled_points %>%
      mutate(cor_LAI_max_h = cor_LAI_max_h)
    
    # hist(sampled_points$max_height)
    
    # sampled_points <- data_df %>%
    #   sample_n(size = 89502) %>%
    #   mutate(max_h = max_h, site = site)
    # stop()
    
    refl_vals <- terra::extract(refl, sampled_points[, c("x", "y")])
    refl_vals <- refl_vals[, -1]  # remove ID column
    # stop()
    
    sampled_points <- bind_cols(sampled_points, refl_vals)
    all_sampled_points <- bind_rows(all_sampled_points, sampled_points)
    # print(hist(sampled_points$lidar_lai))
    
    # stop()
    # Extract LiDAR PAI values at sampled locations
    for (norm_type in norms) {
      lidar_pai_values <- list()
      
      for (d in seq(39.5, 2.5, by = -1)) {
        lidar_pai_rast <- rast(file.path(data_dir, site, 
                                         paste0("LiDAR/testPADs/PAD_Profiles_", 
                                                norm_type,
                                                "_keepTrees"),
                                         # norm_type),
                                         paste0("PAD_", d, "_40.tif")))
        
        # NA test
        # na_mask <- is.na(lidar_pai_rast) | is.na(vci_rast)
        # lidar_pai_rast <- mask(lidar_pai_rast, na_mask, maskvalues = TRUE)
        
        # Extract PAI values at sampled locations
        sampled_pai <- terra::extract(lidar_pai_rast, sampled_points[, c("x", "y")])[, 2]
        
        depth <- 40 - d + 0.5
        
        lidar_pai_values[[paste0("pai_", depth)]] <- sampled_pai
      }
      
      # Combine PAI values into data frame
      sampled_pai_df <- as.data.frame(lidar_pai_values)
      final_df <- cbind(sampled_points, sampled_pai_df)
      final_df <- final_df %>% mutate(site = site, norm_type = norm_type, max_height = max_h)
      all_final_df <- bind_rows(all_final_df, final_df)
      
      # Identify PAI columns
      pai_cols <- grep("^pai_", names(final_df), value = TRUE)
      
      # Fill missing values row-wise with the last valid value
      final_df <- final_df %>%
        rowwise() %>%
        mutate(across(all_of(pai_cols), ~ ifelse(is.na(.x), tail(na.omit(c_across(all_of(pai_cols))), 1), .x))) %>%
        ungroup()
      final_df <- final_df %>%
        mutate(mean_refl = rowMeans(across(c(B03, B04, B08)), na.rm = TRUE))
      
      if (norm_type == "DSM"){
        # Scatterplot of LiDAR LAI vs S2 LAI
        # scat <- ggplot(sampled_points, aes(x = lidar_lai, y = s2_lai)) +
        scat <- ggplot(final_df, aes(x = pai_5, y = s2_lai)) +
          # geom_point(alpha = 0.3, color = "#0072B2") +
          geom_hex(bins = 100) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
          labs(
            title = paste0("LiDAR LAI vs Sentinel-2 LAI\nSite: ", site, ", Max Height > ", max_h),
            x = "LiDAR LAI",
            y = "Sentinel-2 LAI"
          ) +
          xlim(0, 10) +
          ylim(0, 10) +
          theme_bw()
        plot(scat)
        
        # Histogram for each reflectance band
        # Pivot reflectance columns to long format
        refl_long <- final_df %>%
          dplyr::select(B03, B04, B08) %>%
          pivot_longer(
            cols = everything(),
            names_to = "band",
            values_to = "reflectance"
          )
        
        # Plot all bands in the same histogram
        hist <- ggplot(refl_long, aes(x = reflectance, fill = band)) +
          geom_histogram(position = "identity", bins = 100, alpha = 0.7) +
          labs(
            title = paste0("Histogram of Reflectance Bands (Site: ", site,
                           ", Hmax: ", max_h, " m)"),
            x = "Reflectance",
            y = "Count"
          ) +
          theme_bw() +
          scale_fill_manual(values = c("B03" = "green", "B04" = "red", "B08" = "gray"))
        plot(hist)
        
        mean_hist <- ggplot(final_df, aes(x = mean_refl)) +
          geom_histogram(bins = 100, fill = "gray40", alpha = 0.7) +
          labs(
            title = paste0("Histogram of Mean Reflectance (Bands B03, B04, B08) (Site: ", site,
                           ", Hmax: ", max_h, " m)"),
            x = "Mean Reflectance",
            y = "Count"
          ) +
          theme_bw()
        plot(mean_hist)
      }
      
      if (norm_type == "DSM" & max_h == 20) {
        band_pairs <- list(c("B03", "B04"), c("B03", "B08"), c("B04", "B08"))
        
        for (pair in band_pairs) {
          p <- ggplot(final_df, aes_string(x = pair[1], y = pair[2])) +
            geom_hex(bins = 100) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
            labs(
              title = paste0("Scatterplot ", pair[1], " vs ", pair[2], " (Site: ", site, ")"),
              x = paste0(pair[1], " Reflectance"),
              y = paste0(pair[2], " Reflectance")
            ) +
            theme_bw()
          plot(p)
        }
      }
      
      if (norm_type == "DSM" & max_h == 20) {
        reflectance_subset <- final_df %>%
          dplyr::select(B03, B04, B08) %>%
          mutate(site = site)
        
        reflectance_data_list[[site]] <- reflectance_subset
      }
      
      # plot_long_df <- final_df %>%
      #   pivot_longer(
      #     cols = starts_with("pai_"),
      #     names_to = "depth",
      #     names_prefix = "pai_",
      #     names_transform = list(depth = as.integer),
      #     values_to = "pai"
      #   ) %>%
      #   mutate(depth = factor(depth, levels = 1:20)) %>%
      #   drop_na(depth)
      
      # Create the scatterplot (one file per norm_type and max_h for clarity)
      # p <- ggplot(plot_long_df, aes(x = pai, y = s2_lai)) +
      #   geom_hex(bins = 400) +
      #   scale_fill_viridis_c() +
      #   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      #   facet_wrap(~ depth, scales = "fixed", ncol = 2) +
      #   labs(
      #     x = "PAI", y = "S2 LAI",
      #     title = paste0(
      #       "S2 LAI vs PAI (Depths 1–20)\nSite: ", site,
      #       ", Max H > ", max_h,
      #       ", Norm: ", norm_type
      #     )
      #   ) +
      #   xlim(c(0,10)) +
      #   ylim(c(0,10)) +
      #   theme_bw()
      
      # hist(plot_long_df$max_height)
      # hist(plot_long_df$lidar_lai)
      # stop()
      
      # Save plot to file
      # ggsave(
      #   filename = paste0("Plots/scatter_s2lai_pai_depths1to10_", site, "_", norm_type, "_maxh", max_h, ".png"),
      #   plot = p,
      #   width = 10, height = 10
      # )
      
      # Compute metrics for all pai_x columns
      # metrics_df <- bind_rows(lapply(pai_cols, compute_metrics)) %>%
      #   mutate(site = site, 
      #          norm_type = norm_type, 
      #          max_height = max_h, 
      #          depth = as.numeric(gsub("pai_", "", pai_col)),
      #          cor_LAI_total = cor_LAI_total,
      #          cor_LAI_max_h = cor_LAI_max_h)
      
      pai_cols <- names(final_df)[grepl("^pai_", names(final_df))]
      metrics_df <- map_dfr(pai_cols, compute_metrics)
      metrics_df <- metrics_df %>%
        mutate(depth = as.numeric(gsub("pai_", "", pai_col)))
      
      results_list[[paste0(site, "_", norm_type, "_", max_h)]] <- metrics_df # metrics_score_df
    }
  }
}

# Combine reflectance data from all sites
reflectance_all_sites <- bind_rows(reflectance_data_list)

# Pivot to long format for plotting
refl_long_all <- reflectance_all_sites %>%
  pivot_longer(
    cols = c(B03, B04, B08),
    names_to = "band",
    values_to = "reflectance"
  )

# Generate one histogram per band (with data from all sites)
for (b in unique(refl_long_all$band)) {
  p <- ggplot(refl_long_all %>% filter(band == b),
              aes(x = reflectance, fill = site)) +
    geom_histogram(position = "identity", bins = 100, alpha = 0.5) +
    labs(
      title = paste("Histogram of", b, "Reflectance Across Sites"),
      x = "Reflectance",
      y = "Count"
    ) +
    theme_bw()
  
  print(p)
  
  # Optional: save the plot
  # ggsave(paste0("Plots/histogram_", b, "_reflectance_across_sites.png"), p, width = 8, height = 6)
}

stop()
# Combine results into a single data frame
final_results <- bind_rows(results_list)

# Normalize metrics
final_results_norm <- final_results %>%
  mutate(
    norm_R = normalize(R),
    norm_NRMSE = 1 - normalize(NRMSE),
    norm_RMSE = 1 - normalize(RMSE),
    norm_slope = 1 - normalize(abs(slope - 1)),
    score_r_rmse = (norm_R + norm_RMSE) / 2,
    score_r_nrmse = (norm_R + norm_NRMSE) / 2,
    score = (norm_R + norm_NRMSE + norm_slope) / 3
  )

final_results_ranked <- final_results %>%
  mutate(
    rank_R     = min_rank(desc(R)),
    rank_NRMSE = min_rank(NRMSE),
    rank_RMSE  = min_rank(RMSE),
    rank_slope = min_rank(abs(slope - 1)),
    score_r_rmse  = (rank_R + rank_RMSE) / 2,
    score_r_nrmse = (rank_R + rank_NRMSE) / 2,
    score = (rank_R + rank_NRMSE + rank_slope) / 3
  )

# ==== REAGGREGATE RAW VALUES TO GET TRUE "3 SITES" STATS ====
# stop()
# Combine all sampled points
# all_sampled_points <- bind_rows(all_sampled_points)
# Step 1: Define PAI columns
pai_cols <- grep("^pai_", names(all_final_df), value = TRUE)

# Step 2: Long format the data
long_df <- all_final_df %>%
  pivot_longer(cols = all_of(pai_cols), names_to = "depth", values_to = "pai") %>%
  mutate(
    depth = as.numeric(gsub("pai_", "", depth))
  ) %>%
  filter(!is.na(pai))

# Step 3: Metrics per site
site_metrics <- long_df %>%
  group_by(site, depth, norm_type, max_height) %>%
  summarise(
    cor_LAI_total = unique(cor_LAI_total),
    cor_LAI_max_h = unique(cor_LAI_max_h),
    R = cor(pai, s2_lai, use = "complete.obs"),
    RMSE = sqrt(mean((pai - s2_lai)^2, na.rm = TRUE)),
    NRMSE = RMSE / (max(s2_lai, na.rm = TRUE) - min(s2_lai, na.rm = TRUE)),
    slope = coef(lm(s2_lai ~ pai))[2],
    .groups = "drop"
  )

# Step 4: Synthesis across all sites
synthesis_metrics <- long_df %>%
  group_by(depth, norm_type, max_height) %>%
  summarise(
    cor_LAI_total = mean(cor_LAI_total),
    cor_LAI_max_h = cor(s2_lai, lidar_lai),
    R = cor(pai, s2_lai, use = "complete.obs"),
    RMSE = sqrt(mean((pai - s2_lai)^2, na.rm = TRUE)),
    NRMSE = RMSE / (max(s2_lai, na.rm = TRUE) - min(s2_lai, na.rm = TRUE)),
    slope = coef(lm(s2_lai ~ pai))[2],
    .groups = "drop"
  ) %>%
  mutate(site = "Sites combined")

# Step 5: Combine both
synthesis_df <- bind_rows(site_metrics, synthesis_metrics)

# Normalize synthesis separately
synthesis_norm <- synthesis_df %>%
  mutate(
    norm_R = normalize(R),
    norm_NRMSE = 1 - normalize(NRMSE),
    norm_RMSE = 1 - normalize(RMSE),
    norm_slope = 1 - normalize(abs(slope - 1)),
    score_r_rmse = (norm_R + norm_RMSE) / 2,
    score_r_nrmse = (norm_R + norm_NRMSE) / 2,
    score = (norm_R + norm_RMSE + norm_slope) / 3
  )

synthesis_ranked <- synthesis_df %>%
  mutate(
    rank_R     = min_rank(desc(R)),
    rank_NRMSE = min_rank(NRMSE),
    rank_RMSE  = min_rank(RMSE),
    rank_slope = min_rank(abs(slope - 1)),
    score_r_rmse  = (rank_R + rank_RMSE) / 2,
    score_r_rmse_norm = 1 - ((score_r_rmse - min(score_r_rmse)) 
                             / (max(score_r_rmse) - min(score_r_rmse))),
    score_r_nrmse = (rank_R + rank_NRMSE) / 2,
    score = (rank_R + rank_NRMSE + rank_slope) / 3,
    score_normalized = 1 - ((score - min(score)) / (max(score) - min(score)))
  )

# Append synthesis to results
# final_results_norm <- bind_rows(final_results_norm, synthesis_norm)

# Plot 1: Faceted by site (including "3 sites")
# R log(NRMSE) slope norm_R norm_NRMSE norm_slope score_r_nrmse score

# ajouter R total
# couper
# hist ech LAI, max
# synthesis -> aggregate & recal =/= avg
hline_data <- synthesis_norm %>%
  distinct(site, max_height, cor_LAI_total, cor_LAI_max_h)
max_R_data <- synthesis_norm %>%
  # Convert max_height to numeric for comparison
  mutate(max_height_num = as.numeric(as.character(max_height))) %>%
  # Keep only rows where depth <= max_height
  filter(depth <= max_height_num) %>%
  group_by(site, max_height) %>%
  filter(R == max(R)) %>%
  slice(1) %>%  # in case of ties
  ungroup()
# p1 <- ggplot(synthesis_norm, aes(x = depth, y = R, color = norm_type)) +
#   geom_point(alpha = 0.5) +
#   geom_line() + 
#   geom_hline(data = hline_data, aes(yintercept = cor_LAI_total, 
#                                     linetype = "Total LAI Corr")) +
#   geom_hline(data = hline_data, aes(yintercept = cor_LAI_max_h, 
#                                     linetype = "Sampled LAI Corr")) +
#   geom_vline(data = max_R_data, aes(xintercept = depth, linetype = "Optimal Depth"), color = "darkred") +
#   labs(x = "Depth (m)", y = "Normalized Score",
#        title = paste("Score as a Function of Depth (by Site)")) +
#   scale_linetype_manual(values = c("Total LAI Corr" = "solid", 
#                                    "Sampled LAI Corr" = "twodash",
#                                    "Optimal Depth" = "solid")) +
#   theme_bw() +
#   facet_wrap(~ site + max_height, ncol = 3)
# print(p1)

# stop()

# p1 <- ggplot(synthesis_norm, aes(x = depth, y = R, color = norm_type)) +
#   geom_point(alpha = 0.4, size = 1.5) +
#   geom_line(size = 1) + 
#   geom_hline(data = hline_data, 
#              aes(yintercept = cor_LAI_total, linetype = "Total LAI Corr"), 
#              color = "black", linewidth = 0.8) +
#   geom_hline(data = hline_data, 
#              aes(yintercept = cor_LAI_max_h, linetype = "Sampled LAI Corr"), 
#              color = "black", linewidth = 0.8) +
#   geom_vline(data = max_R_data, 
#              aes(xintercept = depth, linetype = "Optimal Depth"), 
#              color = "darkred", linewidth = 0.8) +
#   labs(
#     x = "Depth (m)",
#     y = "Pearson Correlation",
#     # title = "Normalized Score as a Function of Depth",
#     color = "Normalization",
#     linetype = NULL
#   ) +
#   scale_color_brewer(palette = "Set1") +
#   scale_linetype_manual(values = c(
#     "Total LAI Corr" = "dashed", 
#     "Sampled LAI Corr" = "dotted",
#     "Optimal Depth" = "solid")
#   ) +
#   theme_bw(base_size = 14) +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     # plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#     strip.text = element_text(size = 14, face = "bold"),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 14),
#     legend.text = element_text(size = 14)
#   ) +
#   facet_wrap(~ site + max_height, ncol = 3)


# synthesis_norm1 <- synthesis_norm %>%
#   mutate(
#     site = factor(site, levels = c("Aigoual", "Blois", "Mormal", "Sites combined")),
#     max_height = factor(as.character(max_height),
#                         levels = c("10", "15", "20"))
#   )
synthesis_norm <- synthesis_norm %>%
  mutate(
    site = factor(site, levels = c("Aigoual", "Blois", "Mormal", "Sites combined")),
    max_height = factor(max_height)
  )

# Create a dataframe for shaded areas
shade_df <- synthesis_norm %>%
  mutate(max_height_num = as.numeric(as.character(max_height))) %>%
  distinct(site, max_height, max_height_num) %>%
  mutate(xmin = max_height_num,
         xmax = Inf,
         ymin = -Inf,
         ymax = Inf)

p1 <- ggplot(synthesis_norm, aes(x = depth, y = R, color = norm_type)) +
  # Shaded background for depth > max_height
  geom_rect(data = shade_df, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "lightgray", alpha = 0.3) +
  
  geom_point(alpha = 0.4, size = 1.5) +
  geom_line(size = 1) + 
  geom_hline(data = hline_data, 
             aes(yintercept = cor_LAI_total, linetype = "Total LAI Corr"), 
             color = "black", linewidth = 0.8) +
  geom_hline(data = hline_data, 
             aes(yintercept = cor_LAI_max_h, linetype = "Sampled LAI Corr"), 
             color = "black", linewidth = 0.8) +
  geom_vline(data = max_R_data, 
             aes(xintercept = depth, linetype = "Optimal Depth"), 
             color = "darkred", linewidth = 0.8) +
  labs(
    x = "Depth (m)",
    y = "Pearson Correlation",
    color = "Normalization",
    linetype = NULL
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c(
    "Total LAI Corr" = "dashed", 
    "Sampled LAI Corr" = "dotted",
    "Optimal Depth" = "solid")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  facet_grid(max_height ~ site,
             labeller = labeller(
               max_height = c("10" = "10m", "15" = "15m", "20" = "20m")))
print(p1)
ggsave("depth_rbis.png", width = 10, height = 10, unit = "in")
stop()

p11 <- ggplot(synthesis_norm, aes(x = depth, y = score, color = norm_type)) +
  geom_point(alpha = 0.5) +
  geom_line() + 
  labs(x = "Depth (m)", y = "Normalized Score",
       title = paste("Score as a Function of Depth (by Site)")) +
  theme_bw() +
  facet_wrap(~ site + max_height, ncol = 3)
print(p11)
stop()

max_score_depths <- synthesis_ranked %>%
  group_by(site, max_height) %>%
  filter(score_r_rmse_norm == max(score_r_rmse_norm)) %>% # score_r_rmse_norm score_normalized
  summarize(max_score_depth = unique(depth))
synthesis_rankedd <- synthesis_ranked %>%
  left_join(max_score_depths, by = c("site", "max_height"))

p12 <- ggplot(synthesis_rankedd, aes(x = depth, y = score_r_rmse_norm, color = norm_type)) +
  geom_point(alpha = 0.5) +
  geom_line() + 
  geom_vline(aes(xintercept = max_score_depth), color = "red", linetype = "dashed") +
  labs(x = "Depth (m)", y = "Normalized Score",
       title = paste("Score as a Function of Depth (by Site)")) +
  theme_bw() +
  facet_wrap(~ site + max_height, ncol = 3)
print(p12)
stop()

# Plot 2: Faceted by normalization type (including "3 sites")
# p2 <- ggplot(final_results_norm, aes(x = depth, y = R, color = site)) +
#   geom_point(alpha = 0.5) +
#   geom_line() + 
#   labs(x = "Depth (m)", y = "Normalized Score", 
#        title = paste("Score as a Function of Depth (by Normalization) for Trees >", max_h)) +
#   theme_bw() +
#   facet_wrap(~ norm_type + max_height, ncol = 3)
# print(p2)

# Plot 3: Compare score as a function of depth for each site, with curves for different max_h values
# p3 <- ggplot(final_results_norm, aes(x = depth, y = R, color = factor(max_height))) +
#   geom_point(alpha = 0.5) +
#   geom_line() + 
#   labs(x = "Depth (m)", y = "Normalized Score", 
#        title = paste("Score as a Function of Depth by Site and Max Height")) +
#   theme_bw() +
#   facet_wrap(~ site + norm_type, ncol = 2)
# print(p3)
# 
p4 <- ggplot(all_sampled_points, aes(x = lidar_lai, fill = as.factor(max_h))) +
  geom_histogram(binwidth = 0.5, alpha = 0.5, position = "identity") +
  labs(title = "Histogram of lidar_lai by Site and Max Height",
       x = "LiDAR LAI", y = "Frequency", fill = "Max Height") +
  theme_bw() +
  facet_wrap(~ site, scales = "free_y") +  # Facet by site
  theme(legend.position = "bottom") +  # Adjust legend position
  guides(fill = guide_legend(title = "Max Height"))  # Legend for max_height
print(p4)



# Print plots
# print(p1)
# print(p2)
# print(p3)
# print(p4)
# ggsave(filename = file.path("Plots", "plot_1_score_vs_depth_site_max_height_normR.png"),
#        plot = p1, width = 10, height = 10, dpi = 300)
# ggsave(filename = file.path("Plots", "plot_2_score_vs_depth_norm_type_max_height.png"), 
#        plot = p2, width = 10, height = 10, dpi = 300)
# ggsave(filename = file.path("Plots", "plot_3_score_vs_depth_site_max_height.png"), 
#        plot = p3, width = 10, height = 10, dpi = 300)
# ggsave(filename = file.path("Plots", "plot_4_histogram_lidar_lai_site_max_height.png"), 
#        plot = p4, width = 10, height = 10, dpi = 300)
# stop()

# Mean LAD = f(depth)
mean_pai_results_list <- pai_df_list <- list()
for (site in sites) {
  for (norm_type in c("DTM", "DSM")) {  # Loop over DTM and DSM normalization
    lidar_pai_values <- lidar_pai_rasters <- list()
    s2_lai_rast <- rast(file.path(data_dir, site, "Sentinel-2",
                                  "s2lai_summer_depth_study_res_10_m.tif"))
    
    for (d in seq(39.5, 2.5, by = -1)) {
      lidar_pai_rast <- rast(file.path(data_dir, site, 
                                       # paste0("LiDAR/PAD_Profiles_", norm_type, "_keepTrees"), 
                                       paste0("LiDAR/testPADs/PAD_Profiles_", norm_type),
                                       paste0("PAD_", d, "_40.tif")))
      
      # Extract PAI values over the entire raster (not just sampled points)
      lidar_pai_values[[paste0("pai_", 40 - d + 0.5)]] <- values(lidar_pai_rast)
      lidar_pai_rasters[[paste0("pai_", 40 - d + 0.5)]] <- lidar_pai_rast
    }
    # Convert to data frame
    lidar_pai_stack <- rast(lidar_pai_rasters)
    pai_values <- as.data.frame(values(lidar_pai_stack))
    
    s2_lai_values <- as.data.frame(values(s2_lai_rast))
    
    # pai_df <- as.data.frame(lidar_pai_values)
    pai_df <- cbind(s2_lai_values, pai_values)
    names(pai_df) <- c("s2_lai", pai_cols)
    
    # Remove rows where all values are NA
    pai_df <- pai_df %>% filter(rowSums(!is.na(.)) > 0)
    
    # Remove rows where s2_lai is NA
    pai_df <- pai_df %>% filter(!is.na(s2_lai))
    
    # Apply forward fill row-wise for PAI columns
    pai_df[, -1] <- t(apply(pai_df[, -1], 1, fill_last_seen))
    
    # pai_df <- pai_df %>%
    #   mutate(across(everything(), ~ ifelse(is.nan(.), NA, .))) %>%
    #   filter(rowSums(!is.na(.)) > 0) %>%
    #   filter(rowSums(. != 0, na.rm = TRUE) > 0) %>%
    #   filter(rowSums(is.na(select(., starts_with("pai_")))) < ncol(select(., starts_with("pai_"))))
    
    pai_df$site <- site
    pai_df$norm_type <- norm_type
    
    # Compute mean PAI per depth
    mean_pai_df <- pivot_longer(pai_df, cols = starts_with("pai_"), 
                                names_to = "depth", values_to = "PAI") %>%
      mutate(depth = as.numeric(gsub("pai_", "", depth))) %>%
      group_by(site, norm_type, depth) %>%
      summarise(mean_PAI = mean(PAI, na.rm = TRUE), .groups = "drop")
    
    mean_pai_results_list[[paste0(site, "_", norm_type)]] <- mean_pai_df
    pai_df_list[[paste0(site, "_", norm_type)]] <- pai_df
  }
}

# Combine all sites into a single data frame
final_pai_results <- bind_rows(mean_pai_results_list)
pai_df <- bind_rows(pai_df_list)

# Plot mean PAI per depth for each site separately
p5 <- ggplot(final_pai_results, aes(x = depth, y = mean_PAI, color = norm_type)) +
  geom_point(alpha = 0.5) + 
  geom_line() + 
  labs(x = "Depth (m)", y = "Mean PAI", title = "Mean PAI as a Function of Depth") +
  theme_bw() +
  facet_wrap(~ site, ncol = 3)

# Print plot
# print(p5)
ggsave(filename = file.path("Plots", "plot_5_mean_pai.png"),
       plot = p5, width = 10, height = 8, dpi = 300)

pai_long <- pai_df %>%
  select(s2_lai, starts_with("pai_"), site, norm_type) %>%
  pivot_longer(cols = starts_with("pai_"), names_to = "depth", values_to = "pai") %>%
  filter(depth %in% paste0("pai_", 1:10)) %>%
  mutate(depth = as.numeric(gsub("pai_", "", depth)), max_h = max_h)

stop()
# Create a plot for each site / norm_type combination
# ensemble / > 10 / > 15 / > 20
for (s in unique(pai_long$site)) {
  for (n in unique(pai_long$norm_type)) {
    plot_data <- pai_long %>% filter(site == s, norm_type == n)
    
    p <- ggplot(plot_data, aes(x = pai, y = s2_lai)) +
      geom_hex(bins = 400) +
      scale_fill_viridis_c() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~ depth, scales = "fixed", ncol = 2) +
      labs(
        x = "PAI", y = "S2 LAI",
        title = paste("Scatterplots: S2 LAI vs PAI (Depth 1-10)", "\nSite:", s, "| Norm:", n)
      ) +
      xlim(c(0,10)) +
      ylim(c(0,10)) +
      theme_bw()
    
    ggsave(filename = paste0("scatter_s2lai_pai_", s, "_", n, ".png"), plot = p, width = 10, height = 10)
  }
}
stop()
set.seed(42)
# lai_dsm <- values(rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only/ladstack_dsm.tif"), na.rm=T)
# lai_dtm <- values(rast("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only/ladstack_dtm.tif"), na.rm=T)

lai_dsm <- values(rast("/home/corroyez/Documents/NC_Full/PROSAIL-Optimization/01_DATA/Blois/LiDAR/testPADs/PAD_Profiles_DSM_keepTrees/PAD_35.5_40.tif"))
lai_dtm <- values(rast("/home/corroyez/Documents/NC_Full/PROSAIL-Optimization/01_DATA/Blois/LiDAR/testPADs/PAD_Profiles_DTM_keepTrees/PAD_35.5_40.tif"))

samp <- sample(length(lai_dsm), 10000)
lai_dsm_samp <- lai_dsm[samp]
lai_dtm_samp <- lai_dtm[samp]

lai_df <- data.frame(dtm = lai_dtm_samp, dsm = lai_dsm_samp)
ggplot(lai_df, aes(x = lai_dtm_samp, y = lai_dsm_samp)) +
  geom_point(alpha = 0.5)
