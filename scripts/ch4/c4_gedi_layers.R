# ==============================================================================
# c4_gedi_layers.R — Can GEDI supply the S2-invisible below-d_opt layer?
#
# Chapter 3 showed the canopy splits at d_opt into a top layer (seen by S2) and
# a bottom layer (invisible to S2) that buffers microclimate just as much. Here
# we test, at power-beam footprints, whether the GEDI waveform vertical profile
# (L2B PAI per 5 m height bin) tracks the ALS LAD profile layer by layer:
#  (1) mean vertical profiles GEDI vs ALS per site,
#  (2) per-footprint top/bottom split below the ALS canopy top (5 and 10 m
#      bracketing the optical d_opt ~7 m): r(GEDI_layer, ALS_layer),
#  (3) does adding the GEDI bottom layer to S2 improve full-LAI recovery?
# Read-only. Run from NC_Full root.
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(here)
})
source(here::here("scripts", "_article_style.R"))
k_scale <- 0.5 / 0.65
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

fp <- readRDS(here::here("output", "intermediate", "c4", "gedi_matchup_3site.rds"))
fp <- fp[power == "full" & !is.na(lai_s2)]

hbins    <- sprintf("%d_%d", seq(0, 35, 5), seq(5, 40, 5))
gedi_col <- paste0("pai_", hbins, "m")
als_col  <- paste0("LAD_", hbins, "_cor")
G <- as.matrix(fp[, ..gedi_col]); G[is.na(G)] <- 0
# LAD_* bins are per-metre densities (m2/m3): x5 m bin height -> LAI per bin
A <- as.matrix(fp[, ..als_col]) * k_scale * 5; A[is.na(A)] <- 0

# --- (1) mean vertical profiles ----------------------------------------------
prof <- rbindlist(lapply(seq_along(hbins), function(i)
  fp[, .(h_mid = 5 * i - 2.5,
         GEDI = mean(G[.I, i]), ALS = mean(A[.I, i])), by = site]))
pl <- melt(prof, id.vars = c("site", "h_mid"), variable.name = "sensor",
           value.name = "lai_bin")
p1 <- ggplot(pl, aes(lai_bin, h_mid, colour = sensor)) +
  geom_path(linewidth = 0.9) + geom_point(size = 1.6) +
  facet_wrap(~site) +
  scale_colour_manual(values = c(GEDI = "#7570b3", ALS = "#1b9e77")) +
  labs(x = expression("Mean LAI per 5 m bin (" * m^2 * "/" * m^2 * ")"),
       y = "Height above ground (m)",
       title = "Vertical structure: GEDI waveform PAI vs ALS LAD profile",
       subtitle = "Power-beam leaf-on footprints; ALS at k = 0.65", colour = NULL) +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig5_c4_vertical_profiles.png"), p1,
       width = 9, height = 4.2, dpi = 200)

# --- (2) per-footprint top/bottom split below the ALS canopy top -------------
top_bin <- max.col(A > 1e-6, ties.method = "last")
keep    <- rowSums(A > 1e-6) > 0
layer_split <- function(M, d_m) {
  nb <- d_m / 5L
  t(sapply(seq_len(nrow(M)), function(i) {
    tb <- top_bin[i]
    top <- sum(M[i, max(1, tb - nb + 1):tb])
    c(top = top, bottom = sum(M[i, ]) - top)
  }))
}
tab <- rbindlist(lapply(c(5L, 10L), function(d_m) {
  Gs <- layer_split(G, d_m); As <- layer_split(A, d_m)
  dd <- data.table(site = fp$site, lai_s2 = fp$lai_s2,
                   g_top = Gs[, 1], g_bot = Gs[, 2],
                   a_top = As[, 1], a_bot = As[, 2])[keep]
  dd[, .(depth_m = d_m, n = .N,
         r_top    = cor(g_top, a_top),
         r_bottom = cor(g_bot, a_bot),
         r_s2_top = cor(lai_s2, a_top),
         r_s2_bot = cor(lai_s2, a_bot),
         mean_a_bot = mean(a_bot), mean_g_bot = mean(g_bot)), by = site]
}))
cat("=== layer-by-layer agreement (r), split below ALS canopy top ===\n")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(tab, file.path(dir_tabs, "Table6_c4_layer_agreement.csv"))

# --- (3) does the GEDI bottom layer help recover full ALS LAI? ---------------
r2 <- function(m) summary(m)$r.squared
Gs <- layer_split(G, 10L)
dd <- data.table(site = fp$site, lai_s2 = fp$lai_s2, pai = fp$pai,
                 g_bot = Gs[, 2], lai_als = fp$LAI_lidar_cor * k_scale)[keep]
tab3 <- dd[, .(r2_s2       = r2(lm(lai_als ~ lai_s2)),
               r2_s2_gbot  = r2(lm(lai_als ~ lai_s2 + g_bot)),
               r2_gedi_tot = r2(lm(lai_als ~ pai)),
               n = .N), by = site]
cat("\n=== full-LAI recovery with the GEDI bottom layer as extra feature ===\n")
print(tab3[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(tab3, file.path(dir_tabs, "Table6b_c4_bottom_recovery.csv"))
cat("\nDONE\n")
