# ==============================================================================
# Chapter 3 — canonical 3-site S2 saturation stats + FIG1/FIG1b (durable).
# Replaces the lost /tmp figure scripts. Uses the SAME pixel pool convention as
# the SM6a heterogeneity pipeline:
#   lai_als     = sum(Deciduous_Only/ladstack_classic.tif)        [LiDAR LAI]
#   lai_s2_atbd = sm6/{site}/s2lai_summer_atbd_T_res_10_m.tif      [fCover-masked]
#   keep pixels: both non-NA AND lai_als >= 2.0   (deciduous, fCover>=90% baked in)
# Read-only on 03_RESULTS. Writes figures + one stats table.
# Run from NC_Full root:  Rscript scripts/ch3/c3_saturation_3site.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(data.table); library(ggplot2); library(here)
})

sites    <- c("Aigoual", "Blois", "Mormal")
lai_min  <- 2.0
# k convention of the study: PAD precomputed at k_ref=0.5, study uses k_select=0.65
# (steps 06/07/09/11). LAI scales by k_ref/k_select. Must match LAI_ALS_dopt rasters.
k_ref    <- 0.5; k_select <- 0.65
k_scale  <- k_ref / k_select                       # 0.769
als_fn   <- function(s) here("03_RESULTS", s, "Metrics", "Deciduous_Only",
                             "ladstack_classic.tif")
s2_fn    <- function(s) here("output", "intermediate", "sm6", s,
                             "s2lai_summer_atbd_T_res_10_m.tif")

dt <- rbindlist(lapply(sites, function(s) {
  lai_raw <- sum(rast(als_fn(s)), na.rm = TRUE)             # k=0.5 (for the LAI>=2 pool)
  lai_s2  <- rast(s2_fn(s))
  if (!compareGeom(lai_raw, lai_s2, stopOnError = FALSE))
    lai_s2 <- resample(lai_s2, lai_raw, method = "bilinear")
  d <- data.table(site = s,
                  lai_raw = as.numeric(values(lai_raw)),
                  lai_s2  = as.numeric(values(lai_s2)))
  d <- d[is.finite(lai_raw) & is.finite(lai_s2) & lai_raw >= lai_min]  # same pool as SM6a
  d[, lai_als := lai_raw * k_scale]                        # report LAI at k=0.65
  d[, lai_raw := NULL]; d
}))
cat("pixels per site (lai_als>=", lai_min, "):\n", sep = "")
print(dt[, .N, by = site])

# ── Saturation statistics (the numbers the plan must cite) ───────────────────
ov <- dt[, .(sd_als = sd(lai_als), sd_s2 = sd(lai_s2),
             cor_s2_als = cor(lai_s2, lai_als),
             r2_s2_predicts_als = cor(lai_s2, lai_als)^2,
             mean_s2_laile4 = mean(lai_s2[lai_als <= 4]),
             mean_s2_laigt4 = mean(lai_s2[lai_als > 4]),
             n = .N)]
ov[, compression_ratio := sd_als / sd_s2]
persite <- dt[, .(sd_als = sd(lai_als), sd_s2 = sd(lai_s2),
                  slope = coef(lm(lai_s2 ~ lai_als))[2],
                  ratio_mean = mean(lai_als) / mean(lai_s2),
                  n = .N), by = site]

cat("\n=== OVERALL (3-site pooled) ===\n"); print(ov)
cat("\n=== PER SITE ===\n"); print(persite)

stats_out <- rbind(
  data.table(scope = "pooled", site = "ALL", ov),
  data.table(scope = "site", persite,
             cor_s2_als = NA, r2_s2_predicts_als = NA,
             mean_s2_laile4 = NA, mean_s2_laigt4 = NA,
             compression_ratio = persite$sd_als / persite$sd_s2),
  fill = TRUE)
fwrite(stats_out, here("output", "tables", "Table0_saturation_stats_3site.csv"))
fwrite(stats_out, here("manuscripts/ch3", "tables",
                       "Table0_saturation_stats_3site.csv"))

# ── FIG1: scatter LAI_S2 vs LAI_ALS, per-site regression + 1:1 ───────────────
site_cols <- c(Aigoual = "#E8746A", Blois = "#2EA02E", Mormal = "#5B8FF0")
set.seed(1); dplot <- dt[sample(.N, min(.N, 18000))]   # thin for rendering
g1 <- ggplot(dplot, aes(lai_als, lai_s2, colour = site)) +
  geom_point(alpha = 0.06, size = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_colour_manual(values = site_cols) +
  labs(title = "Sentinel-2 saturates in dense canopy",
       subtitle = sprintf("LAI_S2_ATBD vs LiDAR LAI | sd(S2)=%.2f << sd(LiDAR)=%.2f (x%.1f) | dashed = 1:1",
                          ov$sd_s2, ov$sd_als, ov$compression_ratio),
       x = expression(LAI[LiDAR] ~ (m^2/m^2)),
       y = expression(LAI[S2] ~ (PROSAIL)), colour = "Site") +
  coord_equal(xlim = c(0, 8), ylim = c(0, 8)) + theme_minimal(base_size = 12)
ggsave(here("output", "figures", "Fig1_saturation_scatter.png"), g1,
       width = 7, height = 6, dpi = 150)
ggsave(here("manuscripts/ch3", "figures", "Fig1_saturation_scatter.png"), g1,
       width = 7, height = 6, dpi = 150)

# ── FIG1b: distribution LiDAR vs S2 ─────────────────────────────────────────
dlong <- rbind(data.table(source = "LiDAR (ALS)", lai = dt$lai_als),
               data.table(source = "S2 (PROSAIL)", lai = dt$lai_s2))
g2 <- ggplot(dlong, aes(lai, fill = source, colour = source)) +
  geom_density(alpha = 0.35) +
  scale_fill_manual(values = c("LiDAR (ALS)" = "#E8746A", "S2 (PROSAIL)" = "#1FA3A3")) +
  scale_colour_manual(values = c("LiDAR (ALS)" = "#E8746A", "S2 (PROSAIL)" = "#1FA3A3")) +
  labs(title = "S2 saturates: compressed LAI distribution vs LiDAR",
       subtitle = sprintf("sd: LiDAR=%.2f vs S2=%.2f (x%.1f narrower) - S2 cannot resolve high LAI",
                          ov$sd_als, ov$sd_s2, ov$compression_ratio),
       x = expression(LAI ~ (m^2/m^2)), y = "density") +
  theme_minimal(base_size = 12)
ggsave(here("output", "figures", "Fig2_saturation_distribution.png"), g2,
       width = 8, height = 5, dpi = 150)
ggsave(here("manuscripts/ch3", "figures", "Fig2_saturation_distribution.png"), g2,
       width = 8, height = 5, dpi = 150)

cat("\nDONE\n")
