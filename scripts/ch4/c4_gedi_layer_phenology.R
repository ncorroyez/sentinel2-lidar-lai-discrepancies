# ==============================================================================
# c4_gedi_layer_phenology.R — Chapter 4: LAYER-RESOLVED GEDI phenology,
# implementing the Oliveira & Zhang (2025, ERL, Amazon) approach on our
# temperate sites, upgraded with per-footprint ALS structure adjustment.
#
# Method transposed: GEDI L2B pai_z (5 m bins); understory = 0-10 m bins,
# upper canopy = bins > 10 m; monthly composites per layer, structure-adjusted
# (lm layer_pai ~ month + lai_als, power beams); phenometrics from a smoothing
# spline (SOS = max of first derivative, EOS = min; + half-amplitude crossings),
# as in Oliveira & Zhang but on adjusted composites instead of 2-degree grids.
#
# Questions: (1) does the understory green up BEFORE the overstory and hold
# LATER in autumn (vernal window / stretched under-storey phenology = measured
# support for the Ch3 s12 hybrid)? (2) EOS_GEDI vs EOS_S2 at Blois, in days.
# Outputs: manuscripts/ch4/tables/Table17*, figures/Fig9_c4_layer_phenology.png
# ==============================================================================
suppressMessages({
  library(data.table); library(sf); library(ggplot2); library(here)
})
source(here::here("scripts", "_article_style.R"))
k_scale  <- 0.5 / 0.65
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

gedi_bins <- sprintf("pai_%d_%dm", seq(0, 35, 5), seq(5, 40, 5))
load_site <- function(s) {
  d <- as.data.table(sf::st_drop_geometry(
    readRDS(here::here("01_DATA", s, "GEDI", paste0(s, "_GEDI_LiDAR_Metrics.rds")))))
  d[, `:=`(site = s, m = as.integer(as.character(month)))]
  d <- d[!is.na(pai) & pai > 0 & sensitivity >= 0.90 & power == "full" &
           !is.na(LAI_lidar_cor) & LAI_lidar_cor >= 2]
  G <- as.matrix(d[, ..gedi_bins]); G[is.na(G)] <- 0
  d[, `:=`(lai_als = LAI_lidar_cor * k_scale,
           pai_under = G[, 1] + G[, 2],          # 0-10 m
           pai_canopy = rowSums(G[, 3:8]))]      # > 10 m
  d
}
dt <- rbindlist(lapply(c("Blois", "Mormal"), load_site))

# --- monthly composites per layer, adjusted to the common mean ALS LAI --------
adj <- function(dd, col) {
  f  <- lm(reformulate(c("factor(m)", "lai_als"), col), data = dd)
  nd <- dd[, .(n = .N), by = m][, lai_als := mean(dd$lai_als)]
  pr <- predict(f, nd, se.fit = TRUE)
  nd[, `:=`(val = pr$fit, se = pr$se.fit, layer = col)][]
}
comp <- dt[, rbind(adj(.SD, "pai_under"), adj(.SD, "pai_canopy")), by = site]
comp[, layer := factor(layer, c("pai_canopy", "pai_under"),
                       c("Upper canopy (>10 m)", "Understory (0-10 m)"))]
fwrite(comp, file.path(dir_tabs, "Table17_c4_layer_monthly.csv"))
print(dcast(comp, site + m ~ layer, value.var = "val")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# --- figure: normalized seasonal course per layer ------------------------------
pd <- copy(comp)
pd[, doy := 30.4 * (m - 0.5)]
pd[, `:=`(vmax = max(val)), by = .(site, layer)]
pd[, val_n := val / vmax]
p <- ggplot(pd, aes(doy, val_n, colour = layer)) +
  geom_line(linewidth = 0.9) +
  geom_pointrange(aes(ymin = (val - 1.96 * se) / vmax,
                      ymax = (val + 1.96 * se) / vmax), size = 0.4) +
  facet_wrap(~site) +
  scale_colour_manual(values = c("#1b9e77", "#8c6d31"), name = NULL) +
  labs(x = "Day of year (2021-2022 pooled by month)",
       y = "PAI, fraction of layer seasonal maximum",
       title = "Layer-resolved GEDI phenology (Oliveira & Zhang method, ALS-adjusted)",
       subtitle = paste("Power beams, adjusted to common ALS LAI.",
                        "Does the understory lead in spring and lag in autumn?")) +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig9_c4_layer_phenology.png"), p,
       width = 9, height = 4.4, dpi = 200)

# --- phenometrics (Mormal, 10 months): spline + derivative extrema ------------
phenometrics <- function(mm, vv) {
  sp <- smooth.spline(30.4 * (mm - 0.5), vv, df = min(6, length(mm) - 1))
  doy <- 15:350
  fit <- predict(sp, doy)$y
  d1  <- predict(sp, doy, deriv = 1)$y
  amp <- max(fit) - min(fit)
  half_up   <- doy[which(fit >= min(fit) + amp / 2)[1]]
  half_down <- rev(doy[which(fit >= min(fit) + amp / 2)])[1]
  list(SOS = doy[which.max(d1)], POS = doy[which.max(fit)],
       EOS = doy[which.min(d1)], mid_green = half_up, mid_fall = half_down)
}
mo <- comp[site == "Mormal"]
tab_ph <- mo[, as.data.table(phenometrics(m, val)), by = layer]
cat("\n=== Mormal phenometrics per layer (DOY) ===\n"); print(tab_ph)

# Blois autumn: half-amplitude fall crossing GEDI vs the S2 series
bl <- comp[site == "Blois" & layer == "Upper canopy (>10 m)"]
s2 <- as.data.table(readRDS(here::here("output", "intermediate",
                                       "blois_s2_lai_ts_2021.rds")))
s2 <- s2[, .(lai = mean(LAI_S2_ATBD, na.rm = TRUE)), by = doy][order(doy)]
half_fall <- function(doy, v) {
  pk <- max(v[doy >= 150 & doy <= 280]); lo <- min(v[doy > 280])
  f  <- approxfun(doy, v)
  dd <- seq(240, max(doy), 1)
  dd[which(f(dd) <= (pk + lo) / 2)[1]]
}
gb <- dt[site == "Blois" & m %in% c(9, 10, 11)]
cat("\nBlois autumn half-amplitude crossing:\n")
cat("  S2 series:", half_fall(s2$doy, s2$lai), "DOY\n")
cat("  GEDI (canopy layer, monthly adjusted):",
    half_fall(30.4 * (bl$m - 0.5), bl$val), "DOY\n")
fwrite(tab_ph, file.path(dir_tabs, "Table17b_c4_phenometrics.csv"))
cat("\nDONE\n")
