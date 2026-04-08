# ---
# title: "5.plot_ts.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2026-01-14"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}

# Import libraries
library(terra)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)
library(viridis)

# PARAMETERS ---------------------------------------------------------------
data_dir <- '/media/corroyez/MyPassport/01_DATA'
results_dir <- '../../03_RESULTS'

sites <- c("Aigoual", "Blois", "Mormal")
sites <- c("Blois")
distribs <- c("atbd", "atbd_optim_common")

# Target Dates 
common_dates <- seq(from = as.Date("2021-06-05"), 
                    to   = as.Date("2021-09-30"), 
                    by   = "10 days")

# Sampling 
sample_size <- 1000000

# HELPER FUNCTION ----------------------------------------------------------
plot_comparison <- function(data, x_var, y_var, site_name, title_suffix, out_dir, filename) {
  
  # Check if data exists
  if(nrow(data) == 0) {
    cat("    [WARN] No data for plot:", filename, "(df empty)\n")
    return(NULL)
  }
  
  cat("    Creating plot:", filename, "...\n")
  
  # ----------------------------------------------------------------
  # 1. Calculate Statistics per Date
  # ----------------------------------------------------------------
  # We use data.table syntax for speed
  dt_stats <- as.data.table(data)
  
  stats_labels <- dt_stats[, .(
    R     = cor(get(x_var), get(y_var), use = "complete.obs"),
    RMSE  = sqrt(mean((get(y_var) - get(x_var))^2, na.rm = TRUE)),
    Bias  = mean(get(y_var) - get(x_var), na.rm = TRUE),
    # Calculate linear regression slope
    Slope = tryCatch(coef(lm(get(y_var) ~ get(x_var)))[2], error = function(e) NA)
  ), by = Date]
  
  # Create a formatted label string
  stats_labels[, label_text := paste0(
    "R = ", sprintf("%.2f", R), "\n",
    "RMSE = ", sprintf("%.2f", RMSE), "\n",
    "Bias = ", sprintf("%.2f", Bias), "\n",
    "Slope = ", sprintf("%.2f", Slope)
  )]
  
  # ----------------------------------------------------------------
  # 2. Plotting
  # ----------------------------------------------------------------
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_bin2d(bins = 100) + 
    scale_fill_viridis_c(option = "magma", name = "Count") +
    
    # 1:1 Line
    geom_abline(intercept = 0, slope = 1, color = "white", linetype = "dashed", linewidth = 0.8) +
    
    # Add Statistics Text (Top-Left)
    # x = 0.5, y = 14.5 fits well within the (0, 15) limits
    geom_text(data = stats_labels, aes(label = label_text), 
              x = 0.5, y = 14.5, 
              hjust = 0, vjust = 1, 
              size = 3, color = "black", fontface = "bold", inherit.aes = FALSE) +
    
    facet_wrap(~Date, ncol = 4) +
    
    labs(title = paste0(site_name, ": ", title_suffix),
         # subtitle = paste0("N = ", nrow(data), " pixels"),
         x = x_var,
         y = y_var) +
    
    # Limits and Theme
    xlim(c(0, 15)) +
    ylim(c(0, 15))
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray90"),
          axis.text = element_text(size = 8))
  
  # Save
  ggsave(filename = file.path(out_dir, filename), plot = p, width = 8, height = 8, dpi = 300)
}

# MAIN LOOP ----------------------------------------------------------------
for (site in sites) {
  cat("\n========================================\n")
  cat("PLOTTING FOR SITE:", site, "\n")
  
  res_path <- file.path(results_dir, site, "Metrics/Deciduous_Only")
  plot_dir <- file.path(res_path, "Plots_Comparison")
  if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # 1. LOAD DATA
  # ----------------------------------------------------------------
  path_lidar <- file.path(res_path, "lidarlai_res_10_m.tif")
  if(!file.exists(path_lidar)) { cat("[WARN] LiDAR not found for", site, "\n"); next }
  r_lidar <- rast(path_lidar)
  names(r_lidar) <- "LiDAR"
  
  load_ts_raster <- function(distrib, type) {
    fname <- file.path(res_path, "Smooth_TS", paste0("LAI_", type, "_", distrib, ".tif"))
    if(!file.exists(fname)) return(NULL)
    r <- rast(fname)
    names(r) <- paste0(type, "_", distrib, "_", common_dates)
    return(r)
  }
  
  r_gf_atbd    <- load_ts_raster("atbd", "Gapfilled_10d")
  r_corr_atbd  <- load_ts_raster("atbd", "Corrected")
  r_gf_optim   <- load_ts_raster("atbd_optim_common", "Gapfilled_10d")
  r_corr_optim <- load_ts_raster("atbd_optim_common", "Corrected")
  
  if(any(sapply(list(r_gf_atbd, r_corr_atbd, r_gf_optim, r_corr_optim), is.null))) {
    cat("[WARN] Missing S2 rasters for", site, "\n"); next
  }
  
  # 2. PREPARE DATAFRAME
  # ----------------------------------------------------------------
  cat("  -> Stacking Data...\n")
  s_stack <- c(r_lidar, r_gf_atbd, r_corr_atbd, r_gf_optim, r_corr_optim)
  
  cat("  -> Sampling & Removing NAs...\n")
  if (!is.null(sample_size)) {
    # 'regular' method with na.rm=TRUE in terra usually keeps rows if at least one layer is not NA
    # We will clean strictly in the dataframe next
    val_df <- spatSample(s_stack, size = sample_size, method = "regular", na.rm = TRUE, xy = FALSE)
  } else {
    val_df <- as.data.frame(s_stack, na.rm = TRUE)
  }
  
  # Strict NA removal immediately
  val_df <- na.omit(val_df)
  
  if(nrow(val_df) == 0) {
    cat("[ERROR] No common valid pixels found after NA removal. Check masks.\n")
    next
  }
  
  dt <- as.data.table(val_df)
  dt[, id := 1:.N] 
  
  # Melt
  cols_ts <- setdiff(names(dt), c("LiDAR", "id"))
  dt_long <- melt(dt, id.vars = c("id", "LiDAR"), measure.vars = cols_ts, variable.name = "layer", value.name = "LAI")
  
  # Extract Date
  dt_long[, layer := as.character(layer)]
  dt_long[, Date_Str := str_extract(layer, "\\d{4}\\.\\d{2}\\.\\d{2}")]
  dt_long[, Date := as.Date(Date_Str, format = "%Y.%m.%d")]
  
  # Categorize
  dt_long[grepl("Gapfilled_10d_atbd_", layer),      Var := "Gapfilled_atbd"]
  dt_long[grepl("Corrected_atbd_", layer),          Var := "Corr_atbd"]
  dt_long[grepl("Gapfilled_10d_atbd_optim", layer), Var := "Gapfilled_optim"]
  dt_long[grepl("Corrected_atbd_optim", layer),     Var := "Corr_optim"]
  
  # Reshape to Wide
  dt_final <- dcast(dt_long, id + LiDAR + Date ~ Var, value.var = "LAI")
  
  # ----------------------------------------------------------------
  # CRITICAL: FINAL NA CHECK 
  # Ensure no NAs were introduced during reshaping or missed
  # ----------------------------------------------------------------
  dt_final <- na.omit(dt_final)
  
  # Clean memory
  rm(dt, dt_long, val_df); gc()
  
  # 3. PLOTS
  # ----------------------------------------------------------------
  cat("  -> Generating Plots...\n")
  
  # 1. Corrected (atbd) vs Gapfilled (atbd)
  plot_comparison(dt_final, "Corr_atbd", "Gapfilled_atbd", site, "Gapfilled vs Corrected (atbd)", plot_dir, "1_GF_vs_Corr_atbd.png")
  
  # 2. Corrected (optim) vs Gapfilled (optim)
  plot_comparison(dt_final, "Corr_optim", "Gapfilled_optim", site, "Gapfilled vs Corrected (optim)", plot_dir, "2_GF_vs_Corr_optim.png")
  
  # 3. Gapfilled (optim) vs Gapfilled (atbd)
  plot_comparison(dt_final, "Gapfilled_atbd", "Gapfilled_optim", site, "Gapfilled: Optim vs Atbd", plot_dir, "3_GF_optim_vs_atbd.png")
  
  # 4. Corrected (optim) vs Corrected (atbd)
  plot_comparison(dt_final, "Corr_atbd", "Corr_optim", site, "Corrected: Optim vs Atbd", plot_dir, "4_Corr_optim_vs_atbd.png")
  
  # 5. Corrected (atbd) vs LiDAR
  plot_comparison(dt_final, "LiDAR", "Corr_atbd", site, "S2 Corrected (atbd) vs LiDAR", plot_dir, "5_Corr_atbd_vs_LiDAR.png")
  
  # 6. Corrected (optim) vs LiDAR
  plot_comparison(dt_final, "LiDAR", "Corr_optim", site, "S2 Corrected (optim) vs LiDAR", plot_dir, "6_Corr_optim_vs_LiDAR.png")
  
  # 7. Gapfilled (atbd) vs LiDAR
  plot_comparison(dt_final, "LiDAR", "Gapfilled_atbd", site, "S2 Gapfilled (atbd) vs LiDAR", plot_dir, "7_GF_atbd_vs_LiDAR.png")
  
  # 8. Gapfilled (optim) vs LiDAR
  plot_comparison(dt_final, "LiDAR", "Gapfilled_optim", site, "S2 Gapfilled (optim) vs LiDAR", plot_dir, "8_GF_optim_vs_LiDAR.png")
  
  cat("  -> Saved all plots.\n")
}