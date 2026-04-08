# ================================================================
# GEDI LAI Analysis
# ================================================================

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)

data_dir <- "../../01_DATA"
metrics_dir <- "Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois")  # uncomment to test only one site

all_sites_data <- list()

for (site in sites) {
  message("Processing site: ", site)
  
  # Load GEDI CSV
  csv_path <- file.path(data_dir, site, "Geo_Files",
                        paste0("gedi_shots_", site, ".csv"))
  gedi <- read.csv(csv_path)
  
  # Create 'date' column
  gedi <- gedi %>%
    mutate(
      date = make_date(year, month, day),
      lai_gedi = pai
    ) %>%
    dplyr::select(X_utm, Y_utm, lai_gedi, date, RH100)
  
  # Define raster directory and load rasters
  results_dir <- file.path("../../03_RESULTS", site, "Metrics", metrics_dir)
  rasters <- list(
    lai_als = rast(file.path(results_dir, "lidarlai_res_10_m.tif")),
    lai_s2 = rast(file.path(results_dir, "s2lai_summer_v2_res_10_m.tif")),
    h_max = rast(file.path(results_dir, "max_res_10_m.tif")),
    dsm_sd = rast(file.path(results_dir, "dsm_sd_res_10_m.tif"))
  )
  
  # Combine rasters into a single SpatRaster
  rstack <- c(rasters$lai_als, rasters$lai_s2, rasters$h_max, rasters$dsm_sd)
  names(rstack) <- c("lai_als", "lai_s2", "h_max", "dsm_sd")
  
  # Extract raster values at GEDI coordinates
  coords <- gedi[, c("X_utm", "Y_utm")]
  extracted <- terra::extract(rstack, coords)
  df_extracted <- cbind(coords, extracted) %>% dplyr::select(-ID)
  
  # Combine results
  site_data <- gedi %>%
    left_join(df_extracted, by = c("X_utm", "Y_utm")) %>%
    mutate(site = site)
  site_data <- site_data[complete.cases(site_data), ]
  
  all_sites_data[[site]] <- site_data
  
  # Optionally, save per-site dataset
  # out_csv <- file.path(results_dir, paste0("combined_metrics_", site, ".csv"))
  # write.csv(site_data, out_csv, row.names = FALSE)
}

# Combine all sites
combined_data <- bind_rows(all_sites_data)

# 1️⃣ Filter for May–September
df_filtered <- combined_data %>%
  filter(format(date, "%m") %in% c("05", "06", "07", "08", "09"))

# 2️⃣ For each site, compute correlations and simple linear models
site_stats <- df_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    cor_lai_gedi_lai_als = cor(lai_gedi, lai_als, use = "complete.obs"),
    cor_lai_gedi_lai_s2  = cor(lai_gedi, lai_s2, use = "complete.obs"),
    slope_lai_als = coef(lm(lai_gedi ~ lai_als, na.action = na.omit))[2],
    slope_lai_s2  = coef(lm(lai_gedi ~ lai_s2, na.action = na.omit))[2],
    intercept_lai_als = coef(lm(lai_gedi ~ lai_als, na.action = na.omit))[1],
    intercept_lai_s2  = coef(lm(lai_gedi ~ lai_s2, na.action = na.omit))[1],
    cor_rh100_hmaxrh9 = cor(RH100, h_max, use = "complete.obs")
  )

print(site_stats)

# Plot GEDI vs ALS and GEDI vs S2 per site
ggplot(df_filtered, aes(x = lai_als, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~site) +
  labs(
    x = "ALS LAI",
    y = "GEDI LAI",
    title = "GEDI vs ALS LAI (May–September)"
  ) +
  theme_bw()

ggplot(df_filtered, aes(x = lai_s2, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "Sentinel-2 LAI",
    y = "GEDI LAI",
    title = "GEDI vs Sentinel-2 LAI (May–September)"
  ) +
  theme_bw()

ggplot(df_filtered, aes(x = h_max, y = RH100)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "ALS Hmean",
    y = "GEDI RH100",
    title = "GEDI RH100 vs ALS Hmean (May–September)"
  ) +
  theme_bw()
