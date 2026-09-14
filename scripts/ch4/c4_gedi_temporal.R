# ==============================================================================
# c4_gedi_temporal.R — Chapter 4, the TEMPORAL angle: GEDI as the off-summer
# LiDAR truth that Chapter 3 lacked ("no LiDAR truth outside summer").
#
# GEDI shots span the year (Mormal: 10 months, Blois: 7), so monthly composites
# of power-beam PAI give a waveform-LiDAR seasonal cycle of the canopy:
#  (1) seasonal trajectory per site (median +- IQR), Blois overlaid with the S2
#      annual LAI series (19 dates, ATBD) used by the Ch3 DYN scenarios;
#  (2) winter floor (leaf-off PAI = wood/branch area) and summer plateau by
#      canopy-density stratum (ALS LAI bins): does the GEDI seasonal AMPLITUDE
#      rank density where summer S2 saturates? amplitude = leaf-only LAI proxy;
#  (3) composition control: months sample different footprints, so composites
#      are also reported structure-adjusted (lm pai ~ month + LAI_ALS + FORMS-ish)
#      and on a common ALS-LAI support [2.5, 5.5].
# Quality: sensitivity >= 0.90, PAI > 0, ALS LAI >= 2 (raw k=0.5), power beams.
# Outputs: manuscripts/ch4/{tables/Table12*,figures/Fig6,Fig7}.
# ==============================================================================
suppressMessages({
  library(data.table); library(sf); library(ggplot2); library(patchwork); library(here)
})
source(here::here("scripts", "_article_style.R"))
k_scale  <- 0.5 / 0.65
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

load_site <- function(s) {
  d <- as.data.table(sf::st_drop_geometry(
    readRDS(here::here("01_DATA", s, "GEDI", paste0(s, "_GEDI_LiDAR_Metrics.rds")))))
  d[, `:=`(site = s, m = as.integer(as.character(month)),
           yr = as.integer(as.character(year)))]
  d[!is.na(pai) & pai > 0 & sensitivity >= 0.90 &
      !is.na(LAI_lidar_cor) & LAI_lidar_cor >= 2]
}
dt <- rbindlist(lapply(c("Blois", "Mormal"), load_site))
dt[, lai_als := LAI_lidar_cor * k_scale]
pw <- dt[power == "full"]
cat("power-beam shots by site x month x year:\n")
print(dcast(pw[, .N, by = .(site, m, yr)], site + m ~ yr, value.var = "N", fill = 0))

# --- (1) monthly composites + structure-adjusted month effects ----------------
comp <- pw[, .(n = .N, pai_med = median(pai), pai_q25 = quantile(pai, .25),
               pai_q75 = quantile(pai, .75), lai_als_med = median(lai_als)),
           by = .(site, m)][order(site, m)]
# adjust for the month-to-month composition drift: month effect at the common
# mean structure, from lm(pai ~ month + lai_als), power beams
adj <- pw[, {
  f <- lm(pai ~ factor(m) + lai_als, data = .SD)
  nd <- data.table(m = sort(unique(m)), lai_als = mean(lai_als))
  nd[, pai_adj := predict(f, .SD)][]
}, by = site]
comp <- merge(comp, adj[, .(site, m, pai_adj)], by = c("site", "m"))
fwrite(comp, file.path(dir_tabs, "Table12_c4_gedi_monthly.csv"))
print(comp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# S2 annual series at Blois (site mean of the 19-date ATBD inversion, Ch3)
s2ts <- tryCatch({
  ts <- readRDS(here::here("output", "intermediate", "blois_s2_lai_ts_2021.rds"))
  ts <- as.data.table(ts)
  ts[, .(lai_s2 = mean(LAI_S2_ATBD, na.rm = TRUE)), by = date]
}, error = function(e) NULL)

pd <- copy(comp); pd[, doy_mid := 30.4 * (m - 0.5)]
p1 <- ggplot(pd, aes(doy_mid, pai_med)) +
  geom_ribbon(aes(ymin = pai_q25, ymax = pai_q75), alpha = 0.18, fill = "#7570b3") +
  geom_line(colour = "#7570b3", linewidth = 0.9) +
  geom_point(aes(size = n), colour = "#7570b3") +
  geom_line(aes(y = pai_adj), colour = "#7570b3", linetype = 2) +
  { if (!is.null(s2ts))
      geom_line(data = data.table(site = "Blois",
                                  doy_mid = as.integer(format(s2ts$date, "%j")),
                                  pai_med = s2ts$lai_s2),
                colour = "#d95f02", linewidth = 0.8) } +
  facet_wrap(~site) +
  scale_size_continuous(range = c(1, 3.5), guide = "none") +
  labs(x = "Day of year (2021-2022 pooled by month)",
       y = expression("PAI / LAI (" * m^2 * "/" * m^2 * ")"),
       title = "A waveform-LiDAR seasonal cycle: GEDI monthly composites",
       subtitle = paste("Purple: GEDI power-beam median (IQR ribbon; dashed =",
                        "structure-adjusted). Orange (Blois): S2 ATBD annual",
                        "series. Winter floor = wood/branch area.")) +
  theme_article()
ggsave(file.path(dir_figs, "Fig6_c4_gedi_seasonal.png"), p1,
       width = 10, height = 4.6, dpi = 200)

# --- (2) winter floor / summer plateau / amplitude by density stratum ---------
seas_months <- list(winter = c(1, 2, 3, 12), summer = c(6, 7, 8, 9))
pw[, season := fifelse(m %in% seas_months$winter, "winter",
                fifelse(m %in% seas_months$summer, "summer", NA_character_))]
pw[, dens_bin := cut(lai_als, c(2, 3, 4, 9), labels = c("2-3", "3-4", ">4"),
                     include.lowest = TRUE)]
st <- pw[!is.na(season) & !is.na(dens_bin),
         .(n = .N, pai = median(pai)), by = .(site, dens_bin, season)]
stw <- dcast(st, site + dens_bin ~ season, value.var = c("pai", "n"), fill = NA)
stw[, amplitude := pai_summer - pai_winter]
fwrite(stw, file.path(dir_tabs, "Table13_c4_seasonal_amplitude.csv"))
cat("\nwinter floor / summer plateau / amplitude by ALS-LAI stratum:\n")
print(stw[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

sl <- melt(stw, id.vars = c("site", "dens_bin"),
           measure.vars = c("pai_winter", "pai_summer", "amplitude"),
           variable.name = "what")
sl[, what := factor(what, c("pai_winter", "pai_summer", "amplitude"),
                    c("Winter floor (wood)", "Summer plateau", "Amplitude (leaf)"))]
p2 <- ggplot(sl, aes(dens_bin, value, fill = what)) +
  geom_col(position = "dodge", width = 0.75) +
  facet_wrap(~site) +
  scale_fill_manual(values = c("#8c6d31", "#1b9e77", "#7570b3"), name = NULL) +
  labs(x = expression("ALS LAI stratum (" * m^2 * "/" * m^2 * ", k = 0.65)"),
       y = expression("GEDI PAI (" * m^2 * "/" * m^2 * ")"),
       title = "Does the GEDI seasonal amplitude rank canopy density?",
       subtitle = "Winter = Dec-Mar (Mormal) / leaf-off proxy; Summer = Jun-Sep; power beams") +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig7_c4_seasonal_amplitude.png"), p2,
       width = 9, height = 4.4, dpi = 200)

# --- (3) common-support check (composition drift control) --------------------
cs <- pw[lai_als %between% c(2.5, 5.5),
         .(n = .N, pai_med = median(pai)), by = .(site, m)][order(site, m)]
cat("\ncommon-support [2.5, 5.5] composites (control):\n")
print(cs[, .(site, m, n, pai_med = round(pai_med, 2))])
fwrite(cs, file.path(dir_tabs, "Table12b_c4_monthly_commonsupport.csv"))
cat("\nDONE\n")
