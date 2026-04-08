# ---
# title: "max_hists.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2025-08-25"
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
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)

# Settings
results_dir <- "../../03_RESULTS"
metrics_dir <- "Metrics/Deciduous_Only"
sites <- c("Aigoual", "Blois", "Mormal")
max_path <- "max_res_10_m.tif"

# Load and process all data
max_df_list <- map(sites, function(site) {
  path <- file.path(results_dir, site, metrics_dir, max_path)
  max_rast <- rast(path)
  max_vals <- values(max_rast)
  max_vals <- max_vals[!is.na(max_vals)]
  
  q999 <- quantile(max_vals, probs = 0.998, na.rm = TRUE)
  
  tibble(
    max_height = max_vals,
    site = site,
    q999 = q999
  )
})

# Combine all data
max_df <- bind_rows(max_df_list)

# Extract unique quantiles per site for vlines and labels
q999_lines <- max_df %>%
  group_by(site) %>%
  summarise(q999 = first(q999)) %>%
  mutate(label = paste0(site, ": Q99.9 = ", round(q999, 2)))

# Compute max density (to scale Y for placing labels)
density_max <- max_df %>%
  group_split(site) %>%
  map_dbl(~ max(density(.$max_height)$y)) %>%
  max(na.rm = TRUE)

# Define Y positions for labels (spacing from the top)
n_labels <- nrow(q999_lines)
y_positions <- seq(from = density_max * 0.95, 
                   by = -density_max * 0.08, 
                   length.out = n_labels)

# Add label positions
q999_lines <- q999_lines %>%
  arrange(site) %>%
  mutate(x = Inf, y = y_positions)

# Plot
p <- ggplot(max_df, aes(x = max_height, color = site, fill = site)) +
  geom_density(alpha = 0.7) +
  geom_vline(data = q999_lines, aes(xintercept = q999, color = site),
             linetype = "dashed", linewidth = 1) +
  geom_text(
    data = q999_lines,
    aes(x = x, y = y, label = label, color = site),
    hjust = 1.1, size = 5.5, show.legend = FALSE
  ) +
  labs(
    x = "Max Height (m)",
    y = "Density"
  ) +
  theme_minimal(base_size = 16)
ggsave(
  filename = "hist_max.png",
  plot = p,
  width = 10,
  height = 10,
  dpi = 300
)
