# ---
# title: "3.calculate_lidar_metrics.R"
# author: Nathan CORROYEZ, UMR TETIS, Univ Montpellier, AgroParisTech, CIRAD, CNRS, INRAE, F-34196, Montpellier, France
# output: html_document
# last_update: "2024-07-22"
# ---

# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(terra)
library(lidR)
library(sf)
library(rgl)
library(viridis)
library(viridisLite)

site <- "Mormal"
path <- paste0("../../03_RESULTS/", site, "/LiDAR/Heterogeneity_Masks/")

# Load the rasters
meanh <- rast(paste0(path, "mnc_mean_heights_res_10_m.envi"))
std <- rast(paste0(path, "mnc_std_res_10_m.envi"))
deciduous <- st_read(paste0(path, "site_deciduous_only.shp"))

# Make sure CRS match
# if (st_crs(deciduous)$epsg != crs(meanh)@epsg) {
#   deciduous <- st_transform(deciduous, crs = crs(meanh))
# }

# Mask
meanh_decid <- mask(meanh, deciduous)
std_decid <- mask(std, deciduous)
# Aigoual: 27 5 Mormal: 27 6 / 27 5.2
mask_cond <- (meanh_decid >= 28) & (std_decid <= 4.9) # mean 20/30 and std 3/5
meanh_masked <- mask(meanh_decid, mask_cond, maskvalue = FALSE)
std_masked <- mask(std_decid, mask_cond, maskvalue = FALSE)
max_std <- max(values(std_masked), na.rm = TRUE)
pixels_max_std <- which(values(std_masked) == max_std)

# Convert cell numbers to xy coordinates
coords <- xyFromCell(std_masked, pixels_max_std)
print(paste("Maximum std among pixels with meanh >= 20 is:", max_std))
print("Coordinates (x,y) of pixels with max std:")
print(coords)

data_dir <- "/media/corroyez/MyPassport/01_DATA"
las_utm_dir <- file.path(data_dir,
                         site,
                         "LiDAR",
                         "2-las_utm"
)
las_utm_files <- list.files(las_utm_dir, pattern = "\\.las$", full.names = TRUE)
ctg <- readLAScatalog(las_utm_files)

# Convert coords to sf points (assuming same CRS as your LAScatalog)
points_sf <- st_as_sf(data.frame(coords), coords = c("x", "y"), crs = crs(ctg))

# Get tile boundaries as an sf object
ctg_tiles <- st_as_sf(ctg)

# Find which tile(s) contain the point(s)
intersection <- st_intersects(points_sf, ctg_tiles)

# Print result
for (i in seq_along(intersection)) {
  if (length(intersection[[i]]) > 0) {
    cat("Point", i, "is in tile(s):", intersection[[i]], "\n")
    print(ctg_tiles[intersection[[i]], ])
  } else {
    cat("Point", i, "is not inside any tile.\n")
  }
}

las_load <- readLAS(las_utm_files[intersection[[i]]])
# summary(las_load$Z)
# las <- clip_circle(las_load, xcenter = coords[1,1], ycenter = coords[1,2],
#                    radius = 20)
las <- clip_rectangle(
  las_load,
  xleft   = coords[1,1] - 5,
  ybottom = coords[1,2] - 5,
  xright  = coords[1,1] + 5,
  ytop    = coords[1,2] + 5
)

thr <- 2
# DTM
dtm <- rasterize_terrain(las, algorithm = tin())
las_dtm <- normalize_height(las, algo = rasterize_terrain(las, res = 1, pkg = "terra", algorithm = tin()))
las_dtm <- filter_poi(las_dtm, Z <= 40)
las_dtm <- classify_poi(las_dtm, class = 2L, poi = ~las_dtm@data$Z < thr)
las_dtm$Z[las_dtm$Classification == 2] <- 0
summary(las_dtm$Z)

# plot(las_dtm, color = "Classification")

# DSM
dsm <- rasterize_canopy(las_dtm,
                        res = 1,
                        pkg = "terra",
                        algorithm = pitfree()
)
las_dsm <- normalize_height(las_dtm, algo = dsm)
las_dsm$Z[las_dsm$Classification == 4 & las_dsm$Z > 0] <- 0
las_dsm$Z[las_dsm$Classification == 4] <- las_dsm$Z[las_dsm$Classification == 4] + 40 - max(las_dsm$Z[las_dsm$Classification == 4])
las_dsm <- classify_poi(las_dsm, class = 2L, poi = ~las_dsm@data$Z < thr)
las_dsm$Z[las_dsm$Classification == 2] <- 0
summary(las_dsm$Z)
# summary(las_dsm$Z[las_dsm$Classification == 4])

# plot(las_dsm, color = "Classification")

chm <- dsm - dtm
# # rgl.snapshot("dsm.png")
# 
# # DSM via DTM
# # dsmm2 <- rasterize_canopy(las_dtm,
# #                           res = 1,
# #                           pkg = "terra",
# #                           algorithm = p2r()
# # )
# # dsm_filled <- terra::focal(dsmm2, w = matrix(1, 3, 3), fun = fill.na)
# # las_dsm_via_dtm <- normalize_height(las_dtm, algo = dsm_filled)
# # las_dsm_via_dtm$Z <- las_dsm_via_dtm$Z + 40
# # las_dsm_via_dtm$Z[las_dsm_via_dtm$Classification == 2] <- 0
# # las_dsm_via_dtm2 <- end_dsm_norm(las_dtm2)
# 
# summary(las_dtm$Z)
# # Calculation
# las_dtm$Z[las_dtm$Z > 40] <- 40
las_dtm$Z[las_dtm$Classification == 4] <- las_dtm$Z[las_dtm$Classification == 4] + 40 - max(las_dtm$Z[las_dtm$Classification == 4])
# las_dtm2$Z[las_dtm2$Classification == 4] <- las_dtm2$Z[las_dtm2$Classification == 4] + 40 - max(las_dtm2$Z[las_dtm2$Classification == 4])
# las_dtm$Z[las_dtm$Classification == 2] <- 0
# lad_dtm <- LAD(las_dtm@data$Z) #;lad_dtm$lad[is.na(lad_dtm$lad)] <- 0
lad_dtm <- LAD(las_dtm@data$Z)
lad_dsm <- LAD(las_dsm@data$Z)
# lad_dsm_via_dtm <- LAD(las_dsm_via_dtm@data$Z)
# lad_dsm_via_dtm2 <- LAD(las_dsm_via_dtm2@data$Z)

# lad_dtm2 <- cloud_metrics(las_dtm2, ~myPAD_dtm_new(Z, Classification))
# lad_dsm <- cloud_metrics(las_dsm, ~myPAD_dsm_new(Z))
# lad_dsm_via_dtm2 <- cloud_metrics(las_dsm_via_dtm2, ~myPAD_dsm_new(Z))

plot(lad_dtm$lad, lad_dtm$z, type = 'l',
     xlim = c(0, 0.85),
     xlab = expression("LAD (m"^2*"/m"^3*")"),
     ylab = "Height (m)",
     lty = 1,
     col = "black",
     cex.lab = 1.4,  # Increase axis label size
     cex.axis = 1.2  # Increase axis tick size
)
# points(lad_dsm$lad, lad_dsm$z, type = 'l', col = 'blue')
# points(lad_dsm_via_dtm$lad, lad_dsm_via_dtm$z, type = 'l', col = 'green')
points(lad_dsm$lad, lad_dsm$z, type = 'l', col = 'black', lty = 2)  # dashed line
legend("bottomright",
       legend = c(expression(DTM[ALS]), expression(CHM[ALS])),
       col = c("black", "black"),
       lty = c(1, 2),
       bty = "n",
       cex = 1.4  # Increase legend text size
)

sum(lad_dsm$lad)
sum(lad_dtm$lad)

lad_dtm <- rbind(
  data.frame(z = 0,  lad = 0),
  lad_dtm,
  data.frame(z = 40, lad = 0)
)

lad_dsm <- rbind(
  data.frame(z = 0,  lad = 0),
  lad_dsm,
  data.frame(z = 40, lad = 0)
)

# plot(las, color = "Classification")
# plot(las_dtm, color = "Classification")
# rgl.snapshot("dtm.png")

las2_shifted <- las_dsm
las2_shifted@data$X <- las2_shifted@data$X + 140

las_dtm_shifted <- las_dtm
las_dtm_shifted@data$X <- las_dtm_shifted@data$X + 120

combined <- rbind(las_dtm_shifted, las2_shifted)
plot(combined, bg = "white", legend = TRUE, nbreaks = 10, size = 1)
stop()

# 2️⃣ Define a color palette for Z
z_range <- range(combined$Z, na.rm = TRUE)
colors <- viridis(100)
z_col <- colors[cut(combined$Z, 100)]

# 3️⃣ Plot the 3D point cloud
plot3d(combined@data$X, combined@data$Y, combined@data$Z,
       col = z_col, size = 1, zlab = "Z", box = F)
bg3d("white")
# axes3d(FALSE)
# box3d(FALSE)
# rgl.bbox(color = NA)  # removes any remaining box background
# rgl.pop("axes")       # remove any default axes that may persist
rgl.pop("bboxdeco")   # remove box decorations (grid lines, labels)

# 4️⃣ Add LAD profiles as lines to the side of the cloud
x0 <- min(combined@data$X) - 10   # shift LAD lines to the left
y0_dtm <- min(combined@data$Y)
y0_dsm <- y0_dtm + 5                   # small offset to separate the two
center_x <- mean(combined@data$X, na.rm = TRUE)
center_y <- mean(combined@data$Y, na.rm = TRUE)

# Scale LAD to make them visible in the same Z range (optional)
scale_factor <- 10  # try different values
lines3d(x = center_x - 20 + lad_dtm$lad * scale_factor,
        y = rep(y0_dtm, length(lad_dtm$z)),
        z = lad_dtm$z,
        col = "black", lwd = 3, lty = 1)

lines3d(x = center_x + 20 + lad_dsm$lad * scale_factor,
        y = rep(y0_dtm, length(lad_dsm$z)),
        z = lad_dsm$z,
        col = "black", lwd = 3, lty = 2)

# 5️⃣ Add labels
text3d(x0 + scale_factor/2, y0_dtm - 2, z = max(lad_dtm$z),
       texts = "LAD Profiles", color = "black", cex = 1.2)
legend3d("topright", legend = c("DTM (ALS)", "CHM (ALS)"),
         col = c("black", "black"), lty = c(1, 2), lwd = 3, bty = "n")
















# fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
# 
# 
# # 1. Get Z range
# z_range <- range(combined$Z)
# 
# # 2. Define a color palette (e.g., from viridisLite)
# n_colors <- 100
# colors <- viridisLite::viridis(n_colors)
# 
# # 3. Create a function to add the color bar to the plot
# add_colorbar_rgl <- function(z_range, colors, x_pos, y_pos, z_pos, length, width) {
#   # Generate the x, y, z coordinates for the color bar rectangles
#   z_coords <- seq(z_range[1], z_range[2], length.out = length(colors) + 1)
#   
#   for (i in 1:length(colors)) {
#     quads3d(
#       x = c(x_pos, x_pos + width, x_pos + width, x_pos),
#       y = c(y_pos, y_pos, y_pos + length, y_pos + length),
#       z = c(z_coords[i], z_coords[i], z_coords[i+1], z_coords[i+1]),
#       color = colors[i]
#     )
#   }
#   
#   # Add labels
#   text3d(x = x_pos + width * 1.5, y = y_pos, z = z_range[1], texts = round(z_range[1], 1))
#   text3d(x = x_pos + width * 1.5, y = y_pos, z = z_range[2], texts = round(z_range[2], 1))
#   text3d(x = x_pos + width * 1.5, y = y_pos, z = mean(z_range), texts = "Z (meters)", adj = 0.5, cex = 1.2)
# }
# 
# # 4. Plot the point cloud and then add the color bar
# plot(combined, engine = "rgl", bg = 'white')
# 
# # Call the function to add the color bar
# add_colorbar_rgl(
#   z_range = z_range,
#   colors = colors,
#   x_pos = 20,    # Adjust these values to position the color bar
#   y_pos = 20,
#   z_pos = 0,
#   length = 50,
#   width = 5
# )

# --- Colorize point cloud by height ---
z_range <- range(combined$Z, na.rm = TRUE)
colors <- viridis(100)
z_col <- colors[cut(combined$Z, 100)]

# --- Plot point cloud ---
plot3d(combined@data$X, combined@data$Y, combined@data$Z,
       col = z_col, size = 2, zlab = "Z", box = FALSE)
bg3d("white")
rgl.pop("bboxdeco")  # remove box/grid

# --- Cloud center ---
center_x <- mean(combined@data$X, na.rm = TRUE)
center_y <- mean(combined@data$Y, na.rm = TRUE)
scale_factor <- 10   # horizontal scaling for LAD lines

# --- Offsets for the two profiles ---
offset_dtm <- -30
offset_dsm <- 20
y_offset <- 0  # both on same Y-plane; small adjustments possible

# --- 1️⃣ Plot LAD lines ---
# DTM
lines3d(
  x = center_x + offset_dtm + lad_dtm$lad * scale_factor,
  y = rep(center_y + y_offset, length(lad_dtm$z)),
  z = lad_dtm$z,
  col = "black", lwd = 3, lty = 1
)
# CHM
lines3d(
  x = center_x + offset_dsm + lad_dsm$lad * scale_factor,
  y = rep(center_y + y_offset, length(lad_dsm$z)),
  z = lad_dsm$z,
  col = "black", lwd = 3, lty = 2
)

max_lad <- max(c(lad_dtm$lad, lad_dsm$lad))  # take overall max
lad_ticks <- seq(0, max_lad, length.out = 5) # same tick positions for both

# --- 2️⃣ Add X/Z axes for DTM profile ---
# X-axis (horizontal)
# lines3d(
#   x = c(min(lad_dtm$lad), max(lad_dtm$lad)) * scale_factor + center_x + offset_dtm,
#   y = rep(center_y + y_offset, 2),
#   z = rep(min(lad_dtm$z), 2),
#   col = "black", lwd = 1
# )
lines3d(
  x = center_x + offset_dtm + lad_ticks * scale_factor,   # just line
  y = rep(center_y + y_offset, length(lad_ticks)),
  z = rep(min(lad_dtm$z), length(lad_ticks)),
  col = "black", lwd = 1
)
# Z-axis (vertical)
lines3d(
  x = rep(min(lad_dtm$lad) * scale_factor + center_x + offset_dtm, 2),
  y = rep(center_y + y_offset, 2),
  z = range(lad_dtm$z),
  col = "black", lwd = 1
)
# DTM ticks and labels
# lad_ticks <- seq(0, max(lad_dtm$lad), length.out = 5)
# for (i in lad_ticks) {
#   text3d(
#     x = center_x + offset_dtm + i * scale_factor,
#     y = center_y + y_offset - 1,
#     z = min(lad_dtm$z),
#     texts = sprintf("%.2f", i),
#     adj = c(0.5, 1), cex = 1
#   )
# }
for (i in lad_ticks) {
  text3d(
    x = center_x + offset_dtm + i * scale_factor,
    y = center_y + y_offset - 1,
    z = min(lad_dtm$z),
    texts = sprintf("%.2f", i),
    adj = c(0.5, 1), cex = 1
  )
}
z_ticks <- seq(min(lad_dtm$z), max(lad_dtm$z), length.out = 5)
for (i in z_ticks) {
  text3d(
    x = center_x + offset_dtm - 1,
    y = center_y + y_offset,
    z = i,
    texts = round(i, 0),
    adj = c(1, 0.5), cex = 1
  )
}
# Axis labels
text3d(
  x = center_x + offset_dtm + max(lad_dtm$lad) * scale_factor / 2,
  y = center_y + y_offset - 1,
  z = min(lad_dtm$z) - 3,
  texts = "LAD (m²/m³)", cex = 1.2
)
text3d(
  x = center_x + offset_dtm - 8,
  y = center_y + y_offset + 5,
  z = mean(lad_dtm$z),
  texts = "H (m)", cex = 1.2, srt = 90
)

# --- 3️⃣ Add X/Z axes for CHM profile ---
# X-axis
# lines3d(
#   x = c(min(lad_dsm$lad), max(lad_dsm$lad)) * scale_factor + center_x + offset_dsm,
#   y = rep(center_y + y_offset, 2),
#   z = rep(min(lad_dsm$z), 2),
#   col = "black", lwd = 1
# )
lines3d(
  x = center_x + offset_dsm + lad_ticks * scale_factor, 
  y = rep(center_y + y_offset, length(lad_ticks)),
  z = rep(min(lad_dsm$z), length(lad_ticks)),
  col = "black", lwd = 1
)
# Z-axis
lines3d(
  x = rep(min(lad_dsm$lad) * scale_factor + center_x + offset_dsm, 2),
  y = rep(center_y + y_offset, 2),
  z = range(lad_dsm$z),
  col = "black", lwd = 1
)
# CHM ticks and labels
# lad_ticks_chm <- seq(0, max(lad_dsm$lad), length.out = 5)
# for (i in lad_ticks_chm) {
#   text3d(
#     x = center_x + offset_dsm + i * scale_factor,
#     y = center_y + y_offset - 1,
#     z = min(lad_dsm$z),
#     texts = sprintf("%.2f", i),
#     adj = c(0.5, 1), cex = 1
#   )
# }
for (i in lad_ticks) {
  text3d(
    x = center_x + offset_dsm + i * scale_factor,
    y = center_y + y_offset - 1,
    z = min(lad_dsm$z),
    texts = sprintf("%.2f", i),
    adj = c(0.5, 1), cex = 1
  )
}
z_ticks_chm <- seq(min(lad_dsm$z), max(lad_dsm$z), length.out = 5)
for (i in z_ticks_chm) {
  text3d(
    x = center_x + offset_dsm - 1,
    y = center_y + y_offset,
    z = i,
    texts = round(i, 0),
    adj = c(1, 0.5), cex = 1
  )
}
# Axis labels
text3d(
  x = center_x + offset_dsm + max(lad_dsm$lad) * scale_factor / 2,
  y = center_y + y_offset - 1,
  z = min(lad_dsm$z) - 3,
  texts = "LAD (m²/m³)", cex = 1.2
)
text3d(
  x = center_x + offset_dsm,
  y = center_y + y_offset + 20,
  z = mean(lad_dsm$z),
  texts = "H (m)", cex = 1.2, srt = 90
)





# --- 1. Setup & Data ---
z_range <- range(combined$Z, na.rm = TRUE)
colors <- viridis(100)
z_col <- colors[cut(combined$Z, 100)]

# --- 2. Configuration ---
scale_factor <- 10        # Exaggerate LAD curve width
prof_y       <- mean(combined@data$Y, na.rm = TRUE) 
offset_dtm   <- -35       
offset_dsm   <- 25        
tick_len     <- 0.5       

# Calculate shared Limits 
max_lad_val <- max(c(lad_dtm$lad, lad_dsm$lad), na.rm = TRUE)
lad_ticks_vals <- pretty(c(0, max_lad_val), n = 4)

# --- 3. Open Device ---
open3d(windowRect = c(50, 50, 1200, 800)) 
bg3d("white")

# Plot Point Cloud
plot3d(combined@data$X, combined@data$Y, combined@data$Z,
       col = z_col, size = 2, aspect = "iso",
       xlab = "", ylab = "", zlab = "", box = FALSE, axes = FALSE)

# --- 4. Profile Drawing Function ---
draw_technical_profile <- function(lad_df, center_x, offset_x) {
  base_x <- center_x + offset_x
  min_z  <- min(lad_df$z)
  max_z  <- max(lad_df$z)
  
  # A. LAD Curve
  lines3d(x = base_x + lad_df$lad * scale_factor,
          y = rep(prof_y, length(lad_df$z)),
          z = lad_df$z, col = "black", lwd = 4)
  
  # B. Axis Lines
  lines3d(x = c(base_x, base_x), y = c(prof_y, prof_y), 
          z = c(min_z, max_z), col = "black", lwd = 2) 
  lines3d(x = c(base_x, base_x + max(lad_ticks_vals) * scale_factor), 
          y = c(prof_y, prof_y), z = c(min_z, min_z), col = "black", lwd = 2) 
  
  # C. Height Ticks (Bigger Text)
  z_ticks <- pretty(range(lad_df$z), n = 5)
  z_ticks <- z_ticks[z_ticks >= min_z & z_ticks <= max_z]
  for (zt in z_ticks) {
    segments3d(x = c(base_x, base_x - tick_len), y = c(prof_y, prof_y), 
               z = c(zt, zt), col = "black", lwd = 1)
    text3d(x = base_x - 2.5, y = prof_y, z = zt, 
           texts = as.character(zt), cex = 1.1, color = "black", adj = 1) # cex 1.2
  }
  
  # D. LAD Ticks (Bigger Text)
  for (lt in lad_ticks_vals) {
    tx = base_x + lt * scale_factor
    segments3d(x = c(tx, tx), y = c(prof_y, prof_y), 
               z = c(min_z, min_z - tick_len), col = "black", lwd = 1)
    text3d(x = tx, y = prof_y, z = min_z - 3.0, 
           texts = sprintf("%.1f", lt), cex = 1.2, color = "black", adj = 0.5) # cex 1.2
  }
  
  # E. Labels (Much Bigger Text)
  text3d(x = base_x, y = prof_y, z = max_z + 2, 
         texts = "H (m)", cex = 1.5, color = "black", adj = 0) # cex 1.5
  
  text3d(x = base_x + (max(lad_ticks_vals)*scale_factor)/2, 
         y = prof_y, z = min_z - 5, 
         texts = "LAD (m\u00B2/m\u00B3)", cex = 1.5, color = "black") # cex 1.5
}

# --- 5. Legend Drawing Function ---
draw_height_legend <- function(center_x, z_min_scene, width = 50) {
  n_segs <- 100
  h_bar  <- 1.5
  z_base <- z_min_scene - 5
  cols   <- viridis(n_segs)
  x_seq <- seq(center_x - width/2, center_x + width/2, length.out = n_segs + 1)
  
  vx <- as.vector(rbind(x_seq[1:n_segs], x_seq[2:(n_segs+1)], 
                        x_seq[2:(n_segs+1)], x_seq[1:n_segs]))
  vy <- rep(prof_y, 4 * n_segs)
  vz <- as.vector(rbind(rep(z_base, n_segs), rep(z_base, n_segs), 
                        rep(z_base + h_bar, n_segs), rep(z_base + h_bar, n_segs)))
  vcol <- rep(cols, each = 4)
  
  quads3d(vx, vy, vz, col = vcol, lit = FALSE)
  lines3d(x = c(min(x_seq), max(x_seq), max(x_seq), min(x_seq), min(x_seq)),
          y = rep(prof_y, 5),
          z = c(z_base, z_base, z_base + h_bar, z_base + h_bar, z_base),
          col = "black", lwd = 1)
  
  # Legend Numbers (Bigger)
  label_vals <- pretty(z_range, n = 5)
  for(val in label_vals) {
    if(val >= z_range[1] && val <= z_range[2]) {
      prop <- (val - z_range[1]) / (z_range[2] - z_range[1])
      x_loc <- (center_x - width/2) + prop * width
      segments3d(x = c(x_loc, x_loc), y = c(prof_y, prof_y),
                 z = c(z_base, z_base - 0.5), col = "black")
      text3d(x = x_loc, y = prof_y, z = z_base - 3.0,
             texts = as.character(val), cex = 1, color = "black")
    }
  }
  
  # Legend Title (Much Bigger)
  text3d(x = center_x, y = prof_y, z = z_base - 5,
         texts = "Height (m)", cex = 1.5, color = "black")
}

# --- 6. Execute ---
center_x <- mean(combined@data$X, na.rm = TRUE)
min_z_cloud <- min(combined@data$Z, na.rm = TRUE)

draw_technical_profile(lad_dtm, center_x, offset_dtm)
draw_technical_profile(lad_dsm, center_x, offset_dsm)
draw_height_legend(center_x, min_z_cloud, width = 40)

# --- 7. Save High Res Helper ---
save_high_res <- function(filename = "lidar_profile_300dpi.png") {
  old_size <- par3d("windowRect")
  par3d(windowRect = c(0, 0, 3000, 2000))
  Sys.sleep(0.5)
  rgl.snapshot(filename, fmt = "png")
  par3d(windowRect = old_size)
  message(paste("Saved high-res image to:", filename))
}

# Call this function to save
save_high_res()







































































library(rgl)
library(viridis)

# --- 1. Setup & Data ---
z_range <- range(combined$Z, na.rm = TRUE)
colors <- viridis(100)
z_col <- colors[cut(combined$Z, 100)]

# --- 2. Configuration ---
scale_factor <- 10        # Width of the LAD curve
prof_y       <- mean(combined@data$Y, na.rm = TRUE) 
tick_len     <- 0.5       

# --- SUPERPOSITION OFFSETS ---
# Adjust these to align the curves perfectly on your specific trees
offset_dtm_unified <- -15   # Left tree
offset_chm_unified <- 5    # Right tree

# Shared Limits
max_lad_val <- max(c(lad_dtm$lad, lad_dsm$lad), na.rm = TRUE)
lad_ticks_vals <- pretty(c(0, max_lad_val), n = 4)

# --- 3. Open Device ---
# open3d(windowRect = c(50, 50, 1200, 800)) 
open3d(windowRect = c(0, 0, 3000, 2000)) 
bg3d("white")

# --- 4. Plot Point Cloud ---
plot3d(combined@data$X, combined@data$Y, combined@data$Z,
       col = z_col, size = 2, aspect = "iso",
       xlab = "", ylab = "", zlab = "", box = FALSE, axes = FALSE)

# --- 5. Profile Function (Tighter Ticks) ---
draw_superposed_profile <- function(lad_df, center_x, offset_x, title_text) {
  
  base_x <- center_x + offset_x
  min_z  <- min(lad_df$z)
  max_z  <- max(lad_df$z)
  
  # --- A. Draw Curve & Axes ---
  # Curve
  lines3d(x = base_x + lad_df$lad * scale_factor,
          y = rep(prof_y, length(lad_df$z)),
          z = lad_df$z, col = "black", lwd = 4)
  # Vertical Spine
  lines3d(x = c(base_x, base_x), y = c(prof_y, prof_y), 
          z = c(min_z, max_z), col = "black", lwd = 2) 
  # Horizontal Base
  lines3d(x = c(base_x, base_x + max(lad_ticks_vals) * scale_factor), 
          y = c(prof_y, prof_y), z = c(min_z, min_z), col = "black", lwd = 2) 
  
  # --- B. Height Ticks (Side) ---
  z_ticks <- pretty(range(lad_df$z), n = 5)
  z_ticks <- z_ticks[z_ticks >= min_z & z_ticks <= max_z]
  for (zt in z_ticks) {
    segments3d(x = c(base_x, base_x - tick_len), y = c(prof_y, prof_y), 
               z = c(zt, zt), col = "black", lwd = 1)
    text3d(x = base_x - 3.0, y = prof_y, z = zt, 
           texts = as.character(zt), cex = 1.2, color = "black", adj = 1)
  }
  
  # --- C. LAD Ticks (Bottom - CLOSER & STRAIGHT) ---
  for (lt in lad_ticks_vals) {
    tx = base_x + lt * scale_factor
    
    # Draw Tick Mark
    segments3d(x = c(tx, tx), y = c(prof_y, prof_y), 
               z = c(min_z, min_z - tick_len), col = "black", lwd = 1)
    
    # Draw Label (Moved from -3.0 to -1.5 for tightness)
    text3d(x = tx, y = prof_y, z = min_z - 1.5, 
           texts = sprintf("%.1f", lt), cex = 1.2, color = "black", adj = c(0.5, 1))
  }
  
  # --- D. Axis Titles ---
  # H(m)
  text3d(x = base_x, y = prof_y, z = max_z + 2, 
         texts = "Height (m)", cex = 1.5, color = "black", adj = 0)
  
  # LAD Title (Moved up from -7 to -5 to match closer ticks)
  text3d(x = base_x + (max(lad_ticks_vals)*scale_factor)/2, 
         y = prof_y, z = min_z - 5.0, 
         texts = "LAD (m\u00B2/m\u00B3)", cex = 1.5, color = "black")
}

# --- 6. Legend Function ---
draw_height_legend <- function(center_x, z_min_scene, width = 50) {
  n_segs <- 100
  h_bar  <- 1.5
  z_base <- z_min_scene - 10
  cols   <- viridis(n_segs)
  x_seq <- seq(center_x - width/2, center_x + width/2, length.out = n_segs + 1)
  
  vx <- as.vector(rbind(x_seq[1:n_segs], x_seq[2:(n_segs+1)], 
                        x_seq[2:(n_segs+1)], x_seq[1:n_segs]))
  vy <- rep(prof_y, 4 * n_segs)
  vz <- as.vector(rbind(rep(z_base, n_segs), rep(z_base, n_segs), 
                        rep(z_base + h_bar, n_segs), rep(z_base + h_bar, n_segs)))
  vcol <- rep(cols, each = 4)
  
  quads3d(vx, vy, vz, col = vcol, lit = FALSE)
  lines3d(x = c(min(x_seq), max(x_seq), max(x_seq), min(x_seq), min(x_seq)),
          y = rep(prof_y, 5),
          z = c(z_base, z_base, z_base + h_bar, z_base + h_bar, z_base),
          col = "black", lwd = 1)
  
  label_vals <- pretty(z_range, n = 5)
  for(val in label_vals) {
    if(val >= z_range[1] && val <= z_range[2]) {
      prop <- (val - z_range[1]) / (z_range[2] - z_range[1])
      x_loc <- (center_x - width/2) + prop * width
      segments3d(x = c(x_loc, x_loc), y = c(prof_y, prof_y),
                 z = c(z_base, z_base - 0.5), col = "black")
      text3d(x = x_loc, y = prof_y, z = z_base - 2.5, # Adjusted for consistency
             texts = as.character(val), cex = 1.2, color = "black")
    }
  }
  text3d(x = center_x, y = prof_y, z = z_base - 5,
         texts = "Height (m)", cex = 1.5, color = "black")
}

# --- 7. Execute ---
center_x <- mean(combined@data$X, na.rm = TRUE)
min_z_cloud <- min(combined@data$Z, na.rm = TRUE)

# Draw Profiles
draw_superposed_profile(lad_dtm, center_x, offset_dtm_unified, "DTM")
draw_superposed_profile(lad_dsm, center_x, offset_chm_unified, "CHM")

# Draw Legend
draw_height_legend(center_x, min_z_cloud, width = 40)

# --- 8. Save Function ---
save_high_res <- function(filename = "lidar_final_tight_ticks.png") {
  old_size <- par3d("windowRect")
  par3d(windowRect = c(0, 0, 3000, 2000))
  Sys.sleep(1) 
  rgl.snapshot(filename, fmt = "png")
  par3d(windowRect = old_size)
  message(paste("Saved high-res image to:", filename))
}