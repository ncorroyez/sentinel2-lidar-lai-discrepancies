# ---
# title: "three_factors_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-07-11"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------
rm(list=ls(all=TRUE)) 
gc() 

# --- Define working directory as the directory where the script is located ----
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# --------------------------------- Libraries ----------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(relaimpo)
library(ranger)
library(rpart)
library(rpart.plot)

# --------------------------------- Setup --------------------------------------
results_dir    <- "../03_RESULTS"
metrics_dir    <- "Metrics/Deciduous_Only"
lidar_lai_file <- "lidarlai_res_10_m.tif"
max_file       <- "max_res_10_m.tif"
mns_sd_file    <- "dsm_sd_res_10_m.tif"
s2_files       <- list(
  ATBD      = "s2lai_summer_atbd_res_10_m.tif",
  Optimized = "s2lai_summer_best_indiv_res_10_m.tif"
)
sites <- c("Aigoual", "Blois", "Mormal")
set.seed(42)

# ---------------------- 1. Load & Assemble Data ------------------------------
analyze_three_factor_relationships <- function() {
  all_data <- list()
  
  for (site in sites) {
    optimal_fn <- switch(site,
                         "Aigoual" = "PAD_35.5_40.tif",
                         "Blois"   = "PAD_35.5_40.tif",
                         "Mormal"  = "PAD_31.5_40.tif")
                         # "Aigoual" = "PAD_34.5_40.tif",
                         # "Blois"   = "PAD_20.5_40.tif",
                         # "Mormal"  = "PAD_26.5_40.tif")
    base      <- file.path(results_dir, site, metrics_dir)
    
    # load rasters
    rast_opt  <- rast(file.path(base, "PAD_Profiles_dsm_keepTrees", optimal_fn))
    rast_full <- rast(file.path(base, lidar_lai_file))
    rast_het  <- rast(file.path(base, mns_sd_file))
    rast_maxh <- rast(file.path(base, max_file))
    
    print(sum(!is.na(values(rast_opt))))
    print(sum(!is.na(values(rast_full))))
    
    # get global maximum LAI if you really need an explicit upper bound:
    max_lai_value <- max(values(rast_full), na.rm = TRUE)
    
    # ----- Build the reference mask (using FULL lidar) -----
    ref_mask <- rast_full
    ref_mask[] <- is.na(rast_full[]) |
      is.na(rast_het[])  |
      is.na(rast_maxh[]) |
      (rast_full[] < 2)  |   # <-- threshold based on FULL
      (rast_full[] > max_lai_value) |
      (rast_maxh[] <= 10)
    
    for (cfg in names(s2_files)) {
      rast_s2 <- rast(file.path(base, s2_files[[cfg]]))
      
      for (depth_tag in c("Optimal", "Full")) {
        rast_lid <- if (depth_tag == "Optimal") rast_opt else rast_full
        
        # build NA‐mask (also excludes height NA)
        # na_mask <- rast_lid
        # na_mask[] <- is.na(rast_lid[]) |
        #   is.na(rast_s2[])  |
        #   is.na(rast_het[]) |
        #   is.na(rast_maxh[])
        
        # add S2 into mask
        # print(ref_mask)
        # print(rast_s2)
        na_mask <- ref_mask | is.na(rast_s2)
        
        # apply mask
        lid   <- mask(rast_lid,  na_mask, maskvalue=TRUE)
        s2    <- mask(rast_s2,   na_mask, maskvalue=TRUE)
        het   <- mask(rast_het,  na_mask, maskvalue=TRUE)
        maxh  <- mask(rast_maxh, na_mask, maskvalue=TRUE)
        
        # stack into a data.frame
        df <- data.frame(
          lidar   = as.numeric(values(lid)),
          s2      = as.numeric(values(s2)),
          mns_sd  = as.numeric(values(het)),
          height  = as.numeric(values(maxh))
        ) %>%
          drop_na() # %>%
        # filter to lidar ≥2, height >10, and (optionally) ≤ max_lai_value
        # filter(
        #   lidar >= 2,
        #   height > 10,
        #   lidar <= max_lai_value
        # )
        
        # skip if too few pixels
        # if (nrow(df) < 10) next
        
        # uniform sample up to 5000 points
        # if (nrow(df) > 50000) {
        #   df <- df %>% slice_sample(n = 50000)
        # }
        
        # annotate
        df <- df %>%
          mutate(
            site         = site,
            s2_config    = cfg,
            lidar_config = depth_tag
          )
        
        all_data[[length(all_data) + 1]] <- df
        print(str(df))
      }
    }
  }
  
  bind_rows(all_data)
}

# ---------------------- 2. Compute Binned Metrics -----------------------------
# compute_metrics_by_factors <- function(df, n_bins = 15) {
#   df %>%
#     dplyr::group_by(site, s2_config, lidar_config) %>%
#     dplyr::filter(dplyr::n() > 50) %>%
#     dplyr::mutate(
#       mns_sd_bin = cut(
#         mns_sd,
#         breaks = unique(quantile(mns_sd, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)), # manuel ?
#         include.lowest = TRUE
#       )
#     ) %>%
#     dplyr::group_by(site, s2_config, lidar_config, mns_sd_bin) %>%
#     dplyr::filter(dplyr::n() >= 10) %>%
#     dplyr::summarise(
#       mns_sd_center = mean(mns_sd, na.rm = TRUE),
#       R             = cor(lidar, s2, use = "complete.obs"),
#       RMSE          = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
#       Bias          = mean(lidar - s2, na.rm = TRUE),
#       Slope         = coef(lm(lidar ~ s2))[2],
#       .groups       = "drop"
#     ) %>%
#     dplyr::ungroup()
# }

# compute_metrics_by_factors <- function(df) {
# 
#   # Step 1: Create heterogeneity classes and group by site and mns_sd_class
#   # This is the key change to ensure consistent sampling across all configs.
#   grouped_df <- df %>%
#     dplyr::mutate(
#       mns_sd_class = cut(
#         mns_sd,
#         breaks = c(-Inf, 2.5, 5, Inf),
#         labels = c("1: < 2.5", "2: 2.5-5", "3: > 5"),
#         right = FALSE,
#         include.lowest = TRUE
#       )
#     ) %>%
#     dplyr::group_by(site, mns_sd_class)
# 
#   # Count points in each site × mns_sd_class group BEFORE sampling
#   cat("\n-- Nombre de points par classe d'hétérogénéité (avant échantillonnage) --\n")
#   print(
#     grouped_df %>%
#       dplyr::ungroup() %>%
#       dplyr::count(site, mns_sd_class),
#     n = 36
#   )
#   cat("--------------------------------------------------\n\n")
# 
#   # Step 2: Sample down to a maximum of 5000 points per class.
#   # This happens once per site/mns_sd_class group.
#   sampled_df <- grouped_df %>%
#     dplyr::group_modify(~ {
#       if(nrow(.x) > 5000) {
#         .x %>% slice_sample(n = 5000)
#       } else {
#         .x
#       }
#     })
# 
#   # Step 3: Print the number of points in each final class
#   cat("\n-- Nombre de points par classe d'hétérogénéité --\n")
#   print(sampled_df %>%
#           dplyr::ungroup() %>%
#           dplyr::count(site, mns_sd_class), n = 36
#   )
#   cat("--------------------------------------------------\n\n")
# 
#   # Step 4: Recalculate the metrics for each final group
#   sampled_df %>%
#     dplyr::ungroup() %>%
#     dplyr::group_by(site, s2_config, lidar_config, mns_sd_class) %>%
#     dplyr::filter(dplyr::n() >= 10) %>%
#     dplyr::summarise(
#       mns_sd_center = mean(mns_sd, na.rm = TRUE),
#       R             = cor(lidar, s2, use = "complete.obs"),
#       RMSE          = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
#       Bias          = mean(lidar - s2, na.rm = TRUE),
#       Slope         = coef(lm(lidar ~ s2))[2],
#       .groups       = "drop"
#     ) %>%
#     dplyr::ungroup()
# }

compute_metrics_by_factors <- function(df) {
  # Step 1: Create heterogeneity classes and a unique ID for each
  # group of four (lidar_config, s2_config) combinations.
  df_with_id <- df %>%
    dplyr::mutate(
      mns_sd_class = cut(
        mns_sd,
        breaks = c(-Inf, 2.5, 5, Inf),
        labels = c("1: < 2.5", "2: 2.5-5", "3: > 5"),
        right = FALSE,
        include.lowest = TRUE
      )
    ) %>%
    # Create a unique ID for each logical "quadruplet" of configurations
    dplyr::group_by(site, mns_sd_class, lidar_config, s2_config) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  print(df_with_id %>% dplyr::count(site, mns_sd_class, lidar_config, s2_config), n = 36)
  
  # Step 2: Sample down to a maximum of 5000 IDs per class.
  # This ensures all four configurations are consistently represented.
  set.seed(42)
  sampled_ids <- df_with_id %>%
    dplyr::group_by(site, mns_sd_class) %>%
    dplyr::distinct(row_id) %>% # Get a unique list of IDs to sample from
    dplyr::slice_sample(n = 5000) %>% # Sample the IDs
    dplyr::ungroup()
  
  # Step 3: Filter the original data frame using the sampled IDs
  sampled_df <- df_with_id %>%
    dplyr::inner_join(sampled_ids, by = c("site", "mns_sd_class", "row_id"))
  
  # Step 4: Recalculate the metrics for each final group
  sampled_df %>%
    dplyr::group_by(site, s2_config, lidar_config, mns_sd_class) %>%
    dplyr::filter(dplyr::n() >= 10) %>%
    dplyr::summarise(
      mns_sd_center = mean(mns_sd, na.rm = TRUE),
      R = cor(lidar, s2, use = "complete.obs"),
      RMSE = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
      Bias = mean(lidar - s2, na.rm = TRUE),
      Slope = coef(lm(lidar ~ s2))[2],
      .groups = "drop"
    ) %>%
    dplyr::ungroup()
}

# ---------------- 3. Compute Relative Importance per Metric (Random Forest) --
compute_weights_per_metric_rf <- function(metrics_df) {
  metrics_long <- metrics_df %>%
    tidyr::pivot_longer(
      cols = c(R, RMSE, Bias, Slope),
      names_to = "metric",
      values_to = "value"
    )
  
  out <- lapply(unique(metrics_long$metric), function(met) {
    lapply(unique(metrics_long$site), function(st) {
      df_sub <- metrics_long %>%
        dplyr::filter(metric == met, site == st) %>%
        # dplyr::select(value, mns_sd_center) %>%
        dplyr::rename(
          Y = value,
          PROSAIL = s2_config,
          Depth = lidar_config,
          Het = mns_sd_class
        ) %>%
        tidyr::drop_na(Y, Het) %>%
        dplyr::filter(is.finite(Y), is.finite(Het))
      
      if (nrow(df_sub) < 5) return(NULL)
      
      mod <- ranger::ranger(
        Y ~ PROSAIL + Depth + Het,
        data = df_sub,
        importance = "permutation",
        num.trees = 500,
        mtry = 3,
        seed = 42
      )
      
      print(paste("Metric:", met, ", Site:", st))
      # Récupération de l'importance brute
      imp_vals <- ranger::importance(mod)
      print(imp_vals)
      # Forcer les valeurs négatives à zéro
      imp_vals[imp_vals < 0] <- 0
      
      # Normalisation
      total_importance <- sum(imp_vals, na.rm = TRUE)
      if (total_importance == 0) {
        rel_weights <- rep(0, length(imp_vals))
      } else {
        rel_weights <- imp_vals / total_importance
      }
      
      # # Get varImp
      # imp_vals <- ranger::importance(mod)
      # # print(imp_vals)
      # # normImp
      # total_importance <- sum(imp_vals, na.rm = TRUE)
      # rel_weights <- imp_vals / total_importance
      
      # Calculate explained and residual variance
      explained_r2 <- mod$r.squared
      total_variance <- var(df_sub$Y, na.rm = TRUE)
      explained_variance <- total_variance * explained_r2
      residual_variance <- total_variance * (1 - explained_r2)
      
      print(explained_r2)
      print(total_variance)
      print(explained_variance)
      print(residual_variance)
      cat("\n")
      
      data.frame(
        metric = met,
        site = st,
        factor = names(rel_weights),
        weight = as.numeric(rel_weights),
        stringsAsFactors = FALSE
      )
    })
  })
  
  weights <- do.call(rbind, unlist(out, recursive = FALSE))
  
  list(
    weights = weights,
    metrics_long = metrics_long
  )
}

# ----------------------- 4. Plot Pie Charts -----------------------------------
plot_weights_pies <- function(weights_df) {
  df <- weights_df %>%
    dplyr::mutate(
      factor = factor(
        factor,
        levels = c("Depth", "PROSAIL", "Het"),
        labels = c(
          "Canopy Depth",
          "PROSAIL Configuration",
          "Het")
      )
    ) %>%
    dplyr::group_by(site, metric) %>%
    dplyr::arrange(factor) %>%
    dplyr::mutate(
      cum  = cumsum(weight),
      ypos = cum - weight / 2
    ) %>%
    dplyr::ungroup()
  
  ggplot(df, aes(x = 1, y = weight, fill = factor)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    facet_grid(
      rows = vars(site),
      cols = vars(metric),
      labeller = labeller(
        site = label_value,
        metric = c(
          "R"     = "r",
          "RMSE"  = "RMSE",
          "Bias"  = "Bias",
          "Slope" = "Slope"
        )
      ),
      switch = "both"
    ) +
    scale_fill_brewer(
      palette = "Set2",
      name = NULL,
      labels = c(
        "Canopy Depth",
        "PROSAIL Configuration",
        expression("Heterogeneity (CHM"[SD]*")")
      )
    ) +
    theme_void(base_size = 14) +
    theme(
      legend.text     = element_text(size = 14),
      legend.position = "bottom",
      strip.text      = element_text(face = "bold", size = 14),
      plot.title      = element_blank()
    )
}

# ---------------------- 5. Run the Full Pipeline -----------------------------
raw      <- analyze_three_factor_relationships()
metrics  <- compute_metrics_by_factors(raw)

weights  <- compute_weights_per_metric_rf(metrics)
metrics_long <- weights$metrics_long
pie_chart <- plot_weights_pies(weights$weights)

# Display
print(pie_chart)

# Optionally save to file:
# ggsave("Weights_Pie_Charts.png", pie_chart, width = 8, height = 6)

# ---------------------- 6. Compute Summary Tables -----------------------------

# 6a) Mean weights per metric (averaged across all sites)
mean_by_metric <- weights$weights %>%
  dplyr::group_by(metric, factor) %>%
  dplyr::summarise(
    mean_weight = mean(weight, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    pct = round(mean_weight * 100, 1)
  ) %>%
  # dplyr::arrange(metric, desc(mean_weight))
  dplyr::arrange(metric, factor)

# 6b) Mean weights per site (averaged across all metrics)
mean_by_site <- weights$weights %>%
  dplyr::group_by(site, factor) %>%
  dplyr::summarise(
    mean_weight = mean(weight, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    pct = round(mean_weight * 100, 1)
  ) %>%
  # dplyr::arrange(site, desc(mean_weight))
  dplyr::arrange(site, factor)

# Print tables
print("=== Mean Relative Importance per Metric ===")
print(mean_by_metric)

print("=== Mean Relative Importance per Site ===")
print(mean_by_site)

# distribution uniforme
stop()





# ------------------------------------------------------------------------------
check_linearity <- function(df, met) {
  mod <- lm(value ~ mns_sd_center + s2_config + lidar_config, data = df)
  plot(mod, which = 1)  # Residuals vs fitted: should be flat cloud
}

# Example for metric "RMSE" at site "Aigoual"
df_sub <- metrics_long %>% filter(metric == "Slope", site == "Mormal")
check_linearity(df_sub, "Slope")

lm_performance <- metrics_long %>%
  group_by(site, metric) %>%
  group_split() %>%
  purrr::map_df(~{
    df <- .
    mod <- lm(value ~ PROSAIL + Depth + mns_sd_center, data = df)
    tibble(
      site   = unique(df$site),
      metric = unique(df$metric),
      R2     = summary(mod)$r.squared,
      RMSE   = sqrt(mean(residuals(mod)^2)),
      Bias   = mean(residuals(mod))
    )
  })

print(lm_performance)