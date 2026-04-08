# ---
# title: "three_factors_analysis.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-07-11"
# ---

# ----------- (Optional) Clear the environment and free memory -----------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

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
    base     <- file.path(results_dir, site, metrics_dir)
    rast_opt <- terra::rast(file.path(base, "PAD_Profiles_dsm_keepTrees", optimal_fn))
    rast_full<- terra::rast(file.path(base, lidar_lai_file))
    rast_het <- terra::rast(file.path(base, mns_sd_file))
    
    for (cfg in names(s2_files)) {
      rast_s2 <- terra::rast(file.path(base, s2_files[[cfg]]))
      for (depth_tag in c("Optimal", "Full")) {
        rast_lid <- if (depth_tag == "Optimal") rast_opt else rast_full
        
        na_mask <- terra::countNA(rast_lid) + terra::countNA(rast_s2) + terra::countNA(rast_het)
        lid  <- rast_lid;  lid[na_mask > 0] <- NA
        s2   <- rast_s2;   s2[na_mask > 0] <- NA
        het  <- rast_het;  het[na_mask > 0] <- NA
        
        df <- data.frame(
          lidar  = as.numeric(values(lid)),
          s2     = as.numeric(values(s2)),
          mns_sd = as.numeric(values(het))
        ) %>% tidyr::drop_na()
        
        if (nrow(df) < 10) next
        if (nrow(df) > 50000) df <- df[sample(nrow(df), 50000), ]
        
        df <- df %>%
          mutate(
            site         = site,
            s2_config    = cfg,
            lidar_config = depth_tag
          )
        
        all_data[[length(all_data) + 1]] <- df
      }
    }
  }
  dplyr::bind_rows(all_data)
}

analyze_three_factor_relationships <- function() {
  all_data <- list()
  
  for (site in sites) {
    optimal_fn <- switch(site,
                         "Aigoual" = "PAD_35.5_40.tif",
                         "Blois"   = "PAD_35.5_40.tif",
                         "Mormal"  = "PAD_31.5_40.tif")
    base       <- file.path(results_dir, site, metrics_dir)
    
    # load rasters
    rast_opt  <- rast(file.path(base, "PAD_Profiles_dsm_keepTrees", optimal_fn))
    rast_full <- rast(file.path(base, lidar_lai_file))
    rast_het  <- rast(file.path(base, mns_sd_file))
    rast_maxh <- rast(file.path(base, max_file))
    
    # get global maximum LAI if you really need an explicit upper bound:
    max_lai_value <- max(values(rast_full), na.rm = TRUE)
    
    for (cfg in names(s2_files)) {
      rast_s2 <- rast(file.path(base, s2_files[[cfg]]))
      
      for (depth_tag in c("Optimal", "Full")) {
        rast_lid <- if (depth_tag == "Optimal") rast_opt else rast_full
        
        # build NA‐mask (also excludes height NA)
        na_mask <- rast_lid
        na_mask[] <- is.na(rast_lid[]) |
          is.na(rast_s2[])  |
          is.na(rast_het[]) |
          is.na(rast_maxh[])
        
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
          drop_na() %>%
          # filter to lidar ≥2, height >10, and (optionally) ≤ max_lai_value
          filter(
            lidar >= 2,
            height > 10,
            lidar <= max_lai_value
          )
        
        # skip if too few pixels
        if (nrow(df) < 10) next
        
        # uniform sample up to 5000 points
        if (nrow(df) > 5000) {
          df <- df %>% slice_sample(n = 5000)
        }
        
        # annotate
        df <- df %>%
          mutate(
            site         = site,
            s2_config    = cfg,
            lidar_config = depth_tag
          )
        
        all_data[[length(all_data) + 1]] <- df
      }
    }
  }
  
  bind_rows(all_data)
}

# ---------------------- 2. Compute Binned Metrics -----------------------------
compute_metrics_by_factors <- function(df, n_bins = 15) {
  df %>%
    dplyr::group_by(site, s2_config, lidar_config) %>%
    dplyr::filter(dplyr::n() > 50) %>%
    dplyr::mutate(
      mns_sd_bin = cut(
        mns_sd,
        breaks = unique(quantile(mns_sd, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)), # manuel ?
        include.lowest = TRUE
      )
    ) %>%
    dplyr::group_by(site, s2_config, lidar_config, mns_sd_bin) %>%
    dplyr::filter(dplyr::n() >= 10) %>%
    dplyr::summarise(
      mns_sd_center = mean(mns_sd, na.rm = TRUE),
      R             = cor(lidar, s2, use = "complete.obs"),
      RMSE          = sqrt(mean((lidar - s2)^2, na.rm = TRUE)),
      Bias          = mean(lidar - s2, na.rm = TRUE),
      Slope         = coef(lm(lidar ~ s2))[2],
      .groups       = "drop"
    ) %>%
    dplyr::ungroup()
}

# ---------------- 3. Compute Relative Importance per Metric ------------------
compute_weights_per_metric <- function(metrics_df) {
  metrics_long <- metrics_df %>%
    tidyr::pivot_longer(
      cols      = c(R, RMSE, Bias, Slope),
      names_to  = "metric",
      values_to = "value"
    )
  
  out <- lapply(unique(metrics_long$metric), function(met) {
    lapply(unique(metrics_long$site), function(st) {
      df_sub <- metrics_long %>%
        dplyr::filter(metric == met, site == st) %>%
        dplyr::select(value, s2_config, lidar_config, mns_sd_center) %>%
        dplyr::rename(
          Y       = value,
          PROSAIL = s2_config,
          Depth   = lidar_config,
          Het     = mns_sd_center
        ) %>%
        dplyr::mutate(
          PROSAIL = factor(PROSAIL),
          Depth   = factor(Depth)
        ) %>%
        tidyr::drop_na(Y, Het) %>%
        dplyr::filter(is.finite(Y), is.finite(Het))
      
      if (nrow(df_sub) < 5) return(NULL)
      
      mod <- lm(Y ~ PROSAIL + Depth + Het, data = df_sub)
      rel <- calc.relimp(mod, type = "lmg", rela = TRUE)
      
      # Calculate explained and residual variance
      explained_r2 <- summary(mod)$r.squared
      total_variance <- var(df_sub$Y, na.rm = TRUE)
      explained_variance <- total_variance * explained_r2
      residual_variance <- total_variance * (1 - explained_r2)
      
      # Create two data frames and combine them to avoid dimension mismatch
      df_lmg <- data.frame(
        metric = met,
        site = st,
        factor = names(rel$lmg),
        weight = as.numeric(rel$lmg),
        stringsAsFactors = FALSE
      )
      
      df_var <- data.frame(
        metric = met,
        site = st,
        factor = c("Explained_R2", "Explained_Total", "Residual"),
        weight = c(explained_r2, explained_variance, residual_variance),
        stringsAsFactors = FALSE
      )
      
      rbind(df_lmg, df_var)
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
    dplyr::filter(factor %in% c("PROSAIL", "Depth", "Het")) %>%
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
      # panel.spacing = unit(0.5, "lines")
    )
}

# ---------------------- 5. Run the Full Pipeline -----------------------------
# raw      <- analyze_three_factor_relationships()
# metrics  <- compute_metrics_by_factors(raw, n_bins = 10)
# weights  <- compute_weights_per_metric(metrics)
# 
# metrics_long <- weights$metrics_long
# pie_chart <- plot_weights_pies(weights$weights)

# Display
# print(pie_chart)

# Optionally save to file:
# ggsave("Weights_Pie_Charts.png", pie_chart, width = 8, height = 6)

# ---------------------- 5. Run the Full Pipeline -----------------------------
raw      <- analyze_three_factor_relationships()
metrics  <- compute_metrics_by_factors(raw, n_bins = 10)
weights  <- compute_weights_per_metric(metrics)

# Separate the weight data from the explained/residual variance data for plotting
pie_weights <- weights$weights %>%
  dplyr::filter(!factor %in% c("Explained_R2", "Explained_Total", "Residual"))

# Create and display the plot
pie_chart <- plot_weights_pies(pie_weights)
print(pie_chart)

# Display a table of the explained and residual variance, and the R2
explained_residual_df <- weights$weights %>%
  dplyr::filter(factor %in% c("Explained_R2", "Explained_Total", "Residual")) %>%
  tidyr::pivot_wider(names_from = factor, values_from = weight)

cat("\n--- Explained R2, Explained Variance, and Residual Variance ---\n")
print(explained_residual_df)
cat("---------------------------------------\n")
stop()

# ---------------------- 6. Compute Summary Tables -----------------------------

# 6a) Mean weights per metric (averaged across all sites)
mean_by_metric <- weights %>%
  dplyr::group_by(metric, factor) %>%
  dplyr::summarise(
    mean_weight = mean(weight, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    pct = round(mean_weight * 100, 1)
  ) %>%
  dplyr::arrange(metric, desc(mean_weight))

# 6b) Mean weights per site (averaged across all metrics)
mean_by_site <- weights %>%
  dplyr::group_by(site, factor) %>%
  dplyr::summarise(
    mean_weight = mean(weight, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    pct = round(mean_weight * 100, 1)
  ) %>%
  dplyr::arrange(site, desc(mean_weight))

# Print tables
print("=== Mean Relative Importance per Metric ===")
print(mean_by_metric)

print("=== Mean Relative Importance per Site ===")
print(mean_by_site)

# distribution uniforme






# ------------------------------------------------------------------------------
check_linearity <- function(df, met) {
  mod <- lm(value ~ mns_sd_center + s2_config + lidar_config, data = df)
  plot(mod, which = 1)  # Residuals vs fitted: should be flat cloud
}

# Example for metric "RMSE" at site "Aigoual"
df_sub <- metrics_long %>% filter(metric == "R", site == "Mormal")
check_linearity(df_sub, "R")

lm_performance <- metrics_long %>%
  group_by(site, metric) %>%
  group_split() %>%
  purrr::map_df(~{
    df <- .
    mod <- lm(value ~ s2_config + lidar_config + mns_sd_center, data = df)
    tibble(
      site   = unique(df$site),
      metric = unique(df$metric),
      R2     = summary(mod)$r.squared,
      RMSE   = sqrt(mean(residuals(mod)^2)),
      Bias   = mean(residuals(mod))
    )
  })

print(lm_performance)
