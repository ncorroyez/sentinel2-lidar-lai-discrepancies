rm(list=ls(all=TRUE)) # Clear environment
gc() 

# --- 1. SETUP & LIBRARIES ---
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(terra)
library(sf)
library(exactextractr)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. INPUT PARAMETERS ---
radius <- 25
results_dir <- './03_RESULTS'
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# File Paths
path_json <- "data_Blois_utm31n.geojson"
path_lai_stack <- "ladstack_classic_no_na.tif" # ALS Data

# Simulation parameters
corr_types <- c("gap", "corr", "RF/LAI_ALS")
distribs   <- c('atbd', 
                'atbd_optim_common', 
                'atbd_optim_common_brownmodif', 
                'optim_Blois')

# Target simulation dates (Interpolation grid)
sim_days <- seq(from = as.Date("2021-06-05"), 
                to   = as.Date("2021-09-30"), 
                by   = "10 days")
als_date <- as.Date("2021-06-17")

# Check critical inputs
if (!file.exists(path_json)) stop("GeoJSON not found.")
if (!file.exists(path_lai_stack)) stop("ALS Raster not found.")


# --- 3. PREPARE GEOMETRY & ALS REFERENCE ---
print("--- Step 3: Loading Geometry & ALS Reference ---")

# A. Load & Filter Points
json <- st_read(path_json, quiet = TRUE)
ids_remove <- c("41_13", "41_14", "41_20", "41_41", "41_50", "41_51", "41_53")
json_clean <- json %>% filter(!id_plot %in% ids_remove)

# B. Create Buffer Geometry (Vectorized for speed)
# We create one SF object containing all buffered plots
plots_sf <- json_clean %>%
  st_geometry() %>%
  st_buffer(dist = radius) %>%
  st_as_sf()
plots_sf$ID <- json_clean$id_plot

# C. Calculate ALS Reference (LAD Sum)
# Extract all layers at plot centers (no buffer needed for pure profile sum, 
# but using centers is standard for vertical profile checks)
coords <- st_coordinates(json_clean)
r_lad <- terra::rast(path_lai_stack)

# Extract raw values
lad_vals <- terra::extract(r_lad, coords)
df_lad <- cbind(ID = json_clean$id_plot, lad_vals)

# Compute Sum (LAI) ignoring NAs
# Identify LAD columns (assuming they start with "LAD_")
lad_cols <- grep("LAD_", names(df_lad))

df_lai_ref <- df_lad
df_lai_ref$LAI_ALS <- rowSums(df_lad[, lad_cols], na.rm = TRUE)

# Calculate the Global Mean ALS LAI (Reference for Plot)
mean_als_ref <- mean(df_lai_ref$LAI_ALS, na.rm = TRUE)
print(paste("Global Mean ALS LAI calculated:", round(mean_als_ref, 2)))


# --- 4. S2 DATA PROCESSING LOOP ---
print("--- Step 4: Processing S2 Time Series ---")

all_results_list <- list()

for (ct in corr_types) {
  for (dist in distribs) {
    
    # Construct dynamic path based on loop variables
    # Structure: .../smooth/{corr_type}/LAI_Gapfilled_10d_{distrib}.tif
    f_name <- paste0("LAI_", dist, ".tif")
    
    path_s2 <- file.path("../03_RESULTS",
                         "Blois/Metrics/Not_Masked/Smooth_TS/smooth",
                         ct,
                         f_name)
    
    if (!file.exists(path_s2)) {
      warning(paste("Skipping missing file:", path_s2))
      next
    }
    
    message(paste("Processing:", ct, "-", dist))
    
    # Load Raster
    r_stack <- terra::rast(path_s2)
    dates_s2 <- terra::time(r_stack)
    
    if (is.null(dates_s2)) {
      warning(paste("No dates found in raster:", f_name))
      next
    }
    
    # Extract mean values for ALL plots at once (Vectorized)
    # Returns a list of dataframes (one per polygon) or a matrix depending on settings.
    # 'exact_extract' returns a dataframe with one column per layer per feature.
    vals_matrix <- exact_extract(r_stack, plots_sf, fun = 'mean', progress = FALSE)
    
    # Process each plot's time series
    # We apply interpolation row by row
    plot_interpolated <- apply(vals_matrix, 1, function(ts_values) {
      # Remove NAs just in case
      valid <- !is.na(ts_values)
      if (sum(valid) < 2) return(rep(NA, length(sim_days)))
      
      # Interpolate
      approx(x = dates_s2[valid], y = ts_values[valid], xout = sim_days, rule = 2)$y
    })
    
    # Transpose result: Rows = Dates, Cols = Plots -> Convert to Long format
    df_temp <- as.data.frame(plot_interpolated) # Col names are V1, V2...
    colnames(df_temp) <- plots_sf$ID
    df_temp$Date <- sim_days
    
    # Pivot longer to aggregate
    df_long <- df_temp %>%
      tidyr::pivot_longer(cols = -Date, names_to = "ID", values_to = "LAI")
    
    # Aggregate Mean for this specific Correction/Distribution combination
    df_summary <- df_long %>%
      group_by(Date) %>%
      summarise(Mean_LAI = mean(LAI, na.rm = TRUE), .groups = "drop") %>%
      mutate(Corr_Type = ct, Distribution = dist)
    
    all_results_list[[paste(ct, dist, sep = "_")]] <- df_summary
  }
}

# Bind all results
final_data <- bind_rows(all_results_list)
final_data <- final_data %>%
  mutate(
    # Clean Correction Types
    Corr_Type = recode(Corr_Type,
                       "gap"         = "None",
                       "corr"        = "ALS Prorata",
                       "RF/LAI_ALS"  = "Random Forest"),
    
    # Clean Distributions
    Distribution = recode(Distribution,
                          "atbd"                         = "ATBD",
                          "atbd_optim_common"            = "3-sites Opt",
                          "atbd_optim_common_brownmodif" = "3-sites Opt + Brown Modif",
                          "optim_Blois"                  = "Blois-only Opt")
  ) %>%
  # Convert to Factor to control the order in the plot
  mutate(
    Corr_Type = factor(Corr_Type, 
                       levels = c("None", 
                                  "ALS Prorata", 
                                  "Random Forest")
                       ),
    Distribution = factor(Distribution, 
                          levels = c("ATBD", 
                                     "3-sites Opt",
                                     "3-sites Opt + Brown Modif",
                                     "Blois-only Opt")
                          )
  )

# --- 5. PLOTTING ---
print("--- Step 5: Generating Plot ---")

p_compare <- ggplot(final_data, aes(x = Date, y = Mean_LAI, color = Distribution)) +
  # Lines for S2 Distributions
  geom_line(linewidth = 0.8, alpha = 0.8) +
  
  # Reference Line (ALS)
  geom_hline(aes(yintercept = mean_als_ref, linetype = "ALS Ref (Value/Date)"), 
             color = "black", linewidth = 0.7) +
  geom_vline(aes(xintercept = as.numeric(als_date), linetype = "ALS Ref (Value/Date)"), 
             color = "firebrick", linewidth = 0.6, alpha = 0.8) +
  
  # Facet by Correction Type to clean up the view
  facet_wrap(~Corr_Type, ncol = 1) +
  
  # Aesthetics
  scale_linetype_manual(name = "Field Reference", 
                        values = c("ALS Ref (Value/Date)" = "dashed")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "Comparison of S2 Field LAIs",
    x = "Date (2021)",
    y = "Mean LAI",
    color = "Model Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

# Display Plot
print(p_compare)

p_compare2 <- ggplot(final_data, aes(x = Date, y = Mean_LAI, color = Corr_Type)) +
  # Lines for S2 Distributions
  geom_line(linewidth = 0.8, alpha = 0.8) +
  
  # Reference Line (ALS)
  geom_hline(aes(yintercept = mean_als_ref, linetype = "ALS Ref (Value/Date)"), 
             color = "black", linewidth = 0.7) +
  geom_vline(aes(xintercept = as.numeric(als_date), linetype = "ALS Ref (Value/Date)"), 
             color = "firebrick", linewidth = 0.6, alpha = 0.8) +
  
  # Facet by Correction Type to clean up the view
  facet_wrap(~Distribution, ncol = 1) +
  
  # Aesthetics
  scale_linetype_manual(name = "Field Reference", 
                        values = c("ALS Ref (Value/Date)" = "dashed")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(
    title = "Comparison of S2 Field LAIs",
    x = "Date (2021)",
    y = "Mean LAI",
    color = "Model Correction"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

# Display Plot
print(p_compare2)

# Save Plot
# out_file <- file.path(results_dir, "Compare_LAI_Models_vs_ALS.png")
# ggsave(out_file, plot = p_compare, width = 10, height = 8, bg = "white")
# print(paste("Plot saved to:", out_file))