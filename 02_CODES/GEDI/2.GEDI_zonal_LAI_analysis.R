# ---
# title: "zonal_rasters_to_25m.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-11-13"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

data_dir <- "../../01_DATA"
metrics_dir <- "Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
# sites <- c("Blois")  # uncomment to test only one site

all_sites_data <- list()

# ------------------------------------------------------------------
# Define GEDI footprint radius (in meters)
# Default: 12.5 m (diameter ≈ 25 m)
# ------------------------------------------------------------------
gedi_radius <- 12.5  

for (site in sites) {
  message("Processing site: ", site)
  
  # --------------------------------------------------------------
  # 1) Load GEDI
  # --------------------------------------------------------------
  csv_path <- file.path(data_dir, site, "Geo_Files",
                        paste0("gedi_shots_", site, ".csv"))
  gedi <- read.csv(csv_path)
  
  # convert to sf with UTM coordinates
  gedi_sf <- st_as_sf(gedi, coords = c("X_utm", "Y_utm"), crs = 32631)
  
  gedi_sf <- gedi_sf %>%
    mutate(
      date = make_date(year, month, day),
      lai_gedi = pai
    ) %>%
    dplyr::select(lai_gedi, date, RH100)
  
  # --------------------------------------------------------------
  # 2) Create GEDI circular footprints
  # --------------------------------------------------------------
  gedi_buf <- st_buffer(gedi_sf, dist = gedi_radius)
  
  # --------------------------------------------------------------
  # 3) Load ALS / S2 rasters (10 m rasters)
  # --------------------------------------------------------------
  results_dir <- file.path("../../03_RESULTS", site, "Metrics", metrics_dir)
  
  rasters <- list(
    lai_als = rast(file.path(results_dir, "lidarlai_res_10_m.tif")),
    lai_s2  = rast(file.path(results_dir, "s2lai_summer_v2_res_10_m.tif")),
    h_max  = rast(file.path(results_dir, "max_res_10_m.tif")),
    dsm_sd  = rast(file.path(results_dir, "dsm_sd_res_10_m.tif"))
  )
  
  rstack <- c(rasters$lai_als, rasters$lai_s2, rasters$h_max, rasters$dsm_sd)
  names(rstack) <- c("lai_als", "lai_s2", "h_max", "dsm_sd")
  
  
  # --------------------------------------------------------------
  # 4) Zonal statistics for each GEDI footprint
  # exactextractr is fast + handles partial overlaps correctly
  # --------------------------------------------------------------
  # stats_fun <- c("mean", "median", "stdev")
  stats_fun <- c("mean")
  
  zonal <- exact_extract(rstack, gedi_buf, fun = stats_fun)
  
  zonal_df <- as.data.frame(zonal)
  names(zonal_df) <- gsub("\\.", "_", names(zonal_df))  # clean names
  
  
  # --------------------------------------------------------------
  # 5) Combine GEDI + zonal metrics
  # --------------------------------------------------------------
  site_data <- cbind(
    st_coordinates(gedi_sf),
    gedi_sf |> st_drop_geometry(),
    zonal_df
  )
  
  site_data$site <- site
  site_data <- site_data[complete.cases(site_data), ]
  
  all_sites_data[[site]] <- site_data
}

# --------------------------------------------------------------
# 6) Combine all sites into final dataset
# --------------------------------------------------------------
combined_data <- bind_rows(all_sites_data)
months <- c("05", "06", "07", "08", "09")
# months <- c("06", "07", "08")

# 1️⃣ Filter for May–September
df_filtered <- combined_data %>%
  filter(format(date, "%m") %in% months)

# 2️⃣ For each site, compute correlations and simple linear models
site_stats <- df_filtered %>%
  group_by(site) %>%
  summarise(
    n = n(),
    cor_lai_gedi_lai_als = cor(lai_gedi, mean_lai_als, use = "complete.obs"),
    cor_lai_gedi_lai_s2  = cor(lai_gedi, mean_lai_s2, use = "complete.obs"),
    slope_lai_als = coef(lm(lai_gedi ~ mean_lai_als, na.action = na.omit))[2],
    slope_lai_s2  = coef(lm(lai_gedi ~ mean_lai_s2, na.action = na.omit))[2],
    intercept_lai_als = coef(lm(lai_gedi ~ mean_lai_als, na.action = na.omit))[1],
    intercept_lai_s2  = coef(lm(lai_gedi ~ mean_lai_s2, na.action = na.omit))[1],
    cor_rh100_hmax = cor(RH100, mean_h_max, use = "complete.obs"),
    rmse_lai_als          = Metrics::rmse(lai_gedi, mean_lai_als),
    rmse_lai_s2           = Metrics::rmse(lai_gedi, mean_lai_s2),
    rmse_hmax             = Metrics::rmse(RH100, mean_h_max)
  )

print(site_stats)

# Plot GEDI vs ALS and GEDI vs S2 per site
ggplot(df_filtered, aes(x = mean_lai_als, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~site) +
  labs(
    x = "ALS LAI",
    y = "GEDI LAI",
    title = "GEDI vs ALS LAI (May–September)"
  ) +
  xlim(c(0, 15)) + ylim(c(0, 15)) +
  theme_bw()

ggplot(df_filtered, aes(x = mean_lai_s2, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "Sentinel-2 LAI",
    y = "GEDI LAI",
    title = "GEDI vs Sentinel-2 LAI (May–September)"
  ) +
  xlim(c(0, 15)) + ylim(c(0, 15)) +
  theme_bw()

ggplot(df_filtered, aes(x = mean_h_max, y = RH100)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~site) +
  labs(
    x = "ALS Hmean",
    y = "GEDI RH100",
    title = "GEDI RH100 vs ALS Hmean (May–September)"
  ) +
  theme_bw()

# Create classes of mean_h_max
df_filtered_classes <- df_filtered %>%
  mutate(h_class = cut(
    mean_h_max,
    breaks = c(0, 10, 20, 30, 40, Inf),
    labels = c("0–10 m", "10–20 m", "20–30 m", "30–40 m", ">40 m"),
    right = FALSE
  ))

cor_stats_hclass <- df_filtered_classes %>%
  group_by(site, h_class) %>%
  summarise(
    n = n(),
    cor_lai_s2_gedi = cor(mean_lai_s2, lai_gedi, use = "complete.obs")
  ) %>%
  ungroup()

print(cor_stats_hclass)

ggplot(df_filtered_classes, aes(x = mean_lai_s2, y = lai_gedi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", color = "blue") +
  facet_grid(h_class ~ site) +
  labs(
    x = "Sentinel-2 LAI",
    y = "GEDI LAI",
    title = "GEDI vs Sentinel-2 LAI (May–September)"
  ) +
  theme_bw()
