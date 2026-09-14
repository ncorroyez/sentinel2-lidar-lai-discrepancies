# ==============================================================================
# c4_gedi_prototype.R — Chapter 4 prototype: replicate the Chapter 3 Part-A
# analyses (saturation, magnitude ratios, effective depth, magnitude recovery)
# with GEDI added as a third sensor next to ALS and Sentinel-2.
#
# Input: 01_DATA/<site>/GEDI/<site>_GEDI_LiDAR_Metrics.rds (footprint-level
# matchup built by 02_CODES/GEDI/1.ALS_LADs_on_GEDI_footprints.R at k = 0.5),
# plus the Ch3 canonical S2 LAI rasters (summer ATBD, 10 m).
#
# Outputs: manuscripts/ch4/tables/*.csv, manuscripts/ch4/figures/*.png,
# output/intermediate/c4/gedi_matchup_3site.rds
# ==============================================================================
suppressMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(here)
})
source(here::here("scripts", "_article_style.R"))

# --- parameters ---------------------------------------------------------------
k_ref        <- 0.5    # k used when LAI_lidar_* was computed (functions_GEDI.R)
k_select     <- 0.65   # study convention (Ch2/Ch3)
k_scale      <- k_ref / k_select
sens_min     <- 0.90   # GEDI quality filter (0.98 reported as robustness)
leafon_month <- 6:9    # summer identity of Ch3 (May excluded: montane leaf-out)
lai_als_min  <- 2.0    # canonical Ch3 pool filter (applied on raw k=0.5 LAI)
depths_m     <- c(5, 10, 15, 20, 40)  # top-down cumulative depths (40 = full)

sites <- c("Aigoual", "Blois", "Mormal")
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")
dir_int  <- here::here("output", "intermediate", "c4")
for (d in c(dir_tabs, dir_figs, dir_int)) dir.create(d, recursive = TRUE,
                                                     showWarnings = FALSE)

# --- load and harmonise the three matchup datasets ----------------------------
load_site <- function(site) {
  f <- here::here("01_DATA", site, "GEDI", paste0(site, "_GEDI_LiDAR_Metrics.rds"))
  d <- readRDS(f)
  d <- as.data.table(sf::st_drop_geometry(d))
  # geolocation-corrected footprint centres (UTM 31N); fall back to raw UTM
  d[, `:=`(site = site,
           x = fifelse(is.na(x_corrected), x_utm, x_corrected),
           y = fifelse(is.na(y_corrected), y_utm, y_corrected))]
  d
}
dt <- rbindlist(lapply(sites, load_site), fill = TRUE)
cat("Loaded footprints:", nrow(dt), "\n")
print(dt[, table(site, useNA = "ifany")])
cat("\nBDForet types (llib_frt):\n")
print(dt[, .N, by = .(site, llib_frt)][order(site, -N)])

# --- quality + leaf-on + canonical pool filters -------------------------------
dt[, month_num := as.integer(as.character(month))]
dt_all <- copy(dt)
dt <- dt[!is.na(pai) & pai > 0 & !is.na(LAI_lidar_cor) &
           sensitivity >= sens_min &
           month_num %in% leafon_month &
           LAI_lidar_cor >= lai_als_min]

# llib_frt only carries the forest name, not the species: the deciduous mask is
# inherited from the masked S2 raster below (footprints with no valid deciduous
# S2 pixel drop out of the 3-sensor statistics)
cat("\nAfter quality/leaf-on/LAI>=2 filter:", nrow(dt), "\n")

# k harmonisation: express ALS LAI at k = 0.65 (Ch3 convention)
dt[, lai_als := LAI_lidar_cor * k_scale]
dt[, pai_gedi := pai]

# --- extract canonical S2 LAI (summer ATBD) at footprint centres --------------
for (s in sites) {
  r <- terra::rast(here::here("output", "intermediate", "sm6", s,
                              "s2lai_summer_atbd_T_res_10_m.tif"))
  idx <- which(dt$site == s)
  pts <- terra::vect(dt[idx, .(x, y)], geom = c("x", "y"), crs = "EPSG:32631")
  # mean over a 12.5 m buffer approximates the 25 m GEDI footprint on the 10 m grid
  buf <- terra::buffer(pts, 12.5)
  dt[idx, lai_s2 := terra::extract(r, buf, fun = mean, na.rm = TRUE, ID = FALSE)[, 1]]
}
cat("\nS2 extracted, non-NA:", dt[, sum(!is.na(lai_s2))], "/", nrow(dt), "\n")

saveRDS(dt, file.path(dir_int, "gedi_matchup_3site.rds"))

# ==============================================================================
# 1. SATURATION — dynamic range and plateau, three sensors vs ALS full
# ==============================================================================
sat_stats <- function(d) {
  d <- d[!is.na(lai_s2)]
  data.table(
    n              = nrow(d),
    sd_als         = sd(d$lai_als),
    sd_gedi        = sd(d$pai_gedi),
    sd_s2          = sd(d$lai_s2),
    cor_gedi_als   = cor(d$pai_gedi, d$lai_als),
    cor_s2_als     = cor(d$lai_s2, d$lai_als),
    r2_gedi        = cor(d$pai_gedi, d$lai_als)^2,
    r2_s2          = cor(d$lai_s2, d$lai_als)^2,
    gedi_below4    = d[lai_als <= 4, mean(pai_gedi)],
    gedi_above4    = d[lai_als > 4, mean(pai_gedi)],
    s2_below4      = d[lai_als <= 4, mean(lai_s2)],
    s2_above4      = d[lai_als > 4, mean(lai_s2)],
    cor_gedi_hisens = d[sensitivity >= 0.98, cor(pai_gedi, lai_als)],
    n_hisens        = d[, sum(sensitivity >= 0.98)]
  )
}
tab0 <- rbind(dt[, sat_stats(.SD), by = site],
              cbind(site = "POOLED", sat_stats(dt)))
fwrite(tab0, file.path(dir_tabs, "Table0_c4_saturation_3sensor.csv"))
print(tab0, digits = 3)

# Fig1 — scatter GEDI PAI and S2 LAI vs ALS full LAI
dl <- melt(dt[!is.na(lai_s2)],
           id.vars = c("site", "lai_als"),
           measure.vars = c("pai_gedi", "lai_s2"),
           variable.name = "sensor", value.name = "lai")
dl[, sensor := factor(sensor, c("pai_gedi", "lai_s2"),
                      c("GEDI PAI (L2B)", "Sentinel-2 LAI (ATBD)"))]
p1 <- ggplot(dl, aes(lai_als, lai)) +
  geom_point(alpha = 0.15, size = 0.6, colour = "grey30") +
  geom_smooth(aes(colour = sensor), method = "gam",
              formula = y ~ s(x, k = 5), se = FALSE) +
  geom_abline(linetype = 2, colour = "grey55") +
  facet_grid(sensor ~ site) +
  scale_colour_manual(values = c("#7570b3", "#d95f02"), guide = "none") +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 8)) +
  labs(x = expression("ALS LAI, full canopy, k = 0.65 (" * m^2 * "/" * m^2 * ")"),
       y = "Optical / waveform estimate",
       title = "Do GEDI and Sentinel-2 track full-canopy ALS LAI?",
       subtitle = paste("Leaf-on footprints (Jun-Sep), sensitivity ≥ 0.90, PAI > 0,",
                        "LAI_ALS ≥ 2, deciduous mask inherited from the S2 raster")) +
  theme_article()
ggsave(file.path(dir_figs, "Fig1_c4_saturation_scatter.png"), p1,
       width = 9.5, height = 6.5, dpi = 200)

# ==============================================================================
# 2. MAGNITUDE RATIOS — per-site LAI_ALS / sensor, transferability (CV)
# ==============================================================================
tabr <- dt[!is.na(lai_s2),
           .(n = .N,
             mean_als  = mean(lai_als),
             mean_gedi = mean(pai_gedi),
             mean_s2   = mean(lai_s2),
             ratio_gedi = mean(lai_als) / mean(pai_gedi),
             ratio_s2   = mean(lai_als) / mean(lai_s2)),
           by = site]
tabr_cv <- data.table(
  cv_ratio_gedi = tabr[, sd(ratio_gedi) / mean(ratio_gedi)],
  cv_ratio_s2   = tabr[, sd(ratio_s2) / mean(ratio_s2)])
fwrite(tabr,    file.path(dir_tabs, "Table1_c4_ratios_per_site.csv"))
fwrite(tabr_cv, file.path(dir_tabs, "Table1b_c4_ratio_cv.csv"))
print(tabr, digits = 3); print(tabr_cv, digits = 3)

# ==============================================================================
# 3. EFFECTIVE DEPTH — analogue of d_opt: correlation of each sensor with the
#    top-down cumulative ALS LAI, as a function of depth into the canopy
# ==============================================================================
lad_cols <- paste0("LAD_", seq(0, 35, 5), "_", seq(5, 40, 5), "_cor")
stopifnot(all(lad_cols %in% names(dt)))
# LAD_* bins are per-metre densities (m2/m3): x5 m gives LAI per bin (k-harmonised)
lad_mat <- as.matrix(dt[, ..lad_cols]) * k_scale * 5
lad_mat[is.na(lad_mat)] <- 0

# top-down cumulative LAI: for each footprint, find the highest occupied bin and
# sum LAD bins downward over the first `d` metres of canopy
top_bin <- max.col(lad_mat > 1e-6, ties.method = "last")  # index of highest bin
top_bin[rowSums(lad_mat > 1e-6) == 0] <- NA
cum_from_top <- function(d_m) {
  nb <- d_m / 5L
  sapply(seq_len(nrow(lad_mat)), function(i) {
    tb <- top_bin[i]
    if (is.na(tb)) return(NA_real_)
    sum(lad_mat[i, max(1, tb - nb + 1):tb])
  })
}
for (d_m in depths_m) dt[, paste0("lai_top_", d_m) := cum_from_top(d_m)]

dep <- rbindlist(lapply(depths_m, function(d_m) {
  v <- paste0("lai_top_", d_m)
  dt[!is.na(lai_s2) & !is.na(get(v)),
     .(depth_m = d_m,
       cor_gedi = cor(pai_gedi, get(v)),
       cor_s2   = cor(lai_s2, get(v))),
     by = site]
}))
fwrite(dep, file.path(dir_tabs, "Table2_c4_effective_depth.csv"))
print(dep, digits = 3)

dep_l <- melt(dep, id.vars = c("site", "depth_m"),
              variable.name = "sensor", value.name = "r")
dep_l[, sensor := factor(sensor, c("cor_gedi", "cor_s2"),
                         c("GEDI PAI", "Sentinel-2 LAI"))]
p2 <- ggplot(dep_l, aes(depth_m, r, colour = sensor)) +
  annotate("rect", xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "grey40") +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~site) +
  scale_colour_manual(values = c("#7570b3", "#d95f02")) +
  labs(x = "Depth below canopy top (m, cumulative ALS LAI)",
       y = "Pearson r with top-down cumulative ALS LAI",
       title = "How deep into the canopy does each sensor see?",
       subtitle = paste("Grey band: Ch2 optical d_opt (6–8 m).",
                        "40 m = full canopy column."),
       colour = NULL) +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig2_c4_effective_depth.png"), p2,
       width = 9, height = 4.2, dpi = 200)

# ==============================================================================
# 4. MAGNITUDE RECOVERY — R2 of predicting full ALS LAI per site
# ==============================================================================
r2 <- function(m) summary(m)$r.squared
rec <- dt[!is.na(lai_s2),
          .(r2_s2        = r2(lm(lai_als ~ lai_s2)),
            r2_gedi      = r2(lm(lai_als ~ pai_gedi)),
            r2_gedi_s2   = r2(lm(lai_als ~ pai_gedi + lai_s2)),
            r2_gedi_rh   = r2(lm(lai_als ~ pai_gedi + rh98 + cover + fhd_normal)),
            n = .N),
          by = site]
fwrite(rec, file.path(dir_tabs, "Table3_c4_recovery_r2.csv"))
print(rec, digits = 3)

rec_l <- melt(rec, id.vars = c("site", "n"), variable.name = "model",
              value.name = "R2")
rec_l[, model := factor(model,
  c("r2_s2", "r2_gedi", "r2_gedi_s2", "r2_gedi_rh"),
  c("S2 only", "GEDI PAI only", "GEDI + S2", "GEDI PAI+rh98+cover+FHD"))]
p3 <- ggplot(rec_l, aes(model, R2, fill = model)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", R2)), vjust = -0.3, size = 3) +
  facet_wrap(~site) +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  labs(x = NULL, y = expression(R^2 ~ "vs full ALS LAI (footprint level)"),
       title = "Magnitude recovery: does GEDI beat S2 as an anchor?") +
  theme_article() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(dir_figs, "Fig3_c4_recovery_r2.png"), p3,
       width = 9.5, height = 4.4, dpi = 200)

cat("\nDone. Tables in", dir_tabs, "\nFigures in", dir_figs, "\n")
