# ==============================================================================
# c4_gedi_blois_summer.R — Chapter 4, Blois only, WITHIN-SUMMER temporal signal.
#
# Chapter 3 assumes a flat summer plateau (all STATIC scenarios; CONST==DYN in
# JJAS). GEDI gives 9 acquisition dates at Blois (2021 + 2022, incl. the 2022
# drought summer): per-date composites, STRUCTURE-ADJUSTED so that dates
# sampling different stands are comparable:
#     lm(pai ~ factor(date) + lai_als), power beams, prediction at the common
#     mean ALS LAI, +-1.96 SE. Robustness: half beams added with a beam term.
# Questions: (1) is the 2021 summer plateau flat (Jun 6 dip real?);
#            (2) Sep 2021 vs Sep 2022 (drought year) at matched structure;
#            (3) agreement with the S2 annual series at the same dates.
# Outputs: manuscripts/ch4/tables/Table14_c4_blois_dates.csv, figures/Fig8_*.png
# ==============================================================================
suppressMessages({
  library(data.table); library(sf); library(ggplot2); library(here)
})
source(here::here("scripts", "_article_style.R"))
k_scale  <- 0.5 / 0.65
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

d <- as.data.table(sf::st_drop_geometry(
  readRDS(here::here("01_DATA", "Blois", "GEDI", "Blois_GEDI_LiDAR_Metrics.rds"))))
d <- d[!is.na(pai) & pai > 0 & sensitivity >= 0.90 &
         !is.na(LAI_lidar_cor) & LAI_lidar_cor >= 2]
d[, `:=`(date = as.Date(sprintf("%s-%s-%s", year, month, day)),
         lai_als = LAI_lidar_cor * k_scale)]
d[, hour := suppressWarnings(as.integer(substr(time_utc, 12, 13)))]
d[, daylight := fifelse(is.na(hour), NA, hour >= 6 & hour <= 20)]

# leaf-on dates only (drop Apr 11 pre-leaf-out and Nov 14 senescent for the
# plateau test; they stay in the printed table for context)
d[, leafon := !format(date, "%m") %in% c("04", "11")]

adj_means <- function(dd, label) {
  f  <- lm(pai ~ factor(date) + lai_als, data = dd)
  nd <- dd[, .(n = .N), by = date][, lai_als := mean(dd$lai_als)]
  pr <- predict(f, nd, se.fit = TRUE)
  nd[, `:=`(pai_adj = pr$fit, lo = pr$fit - 1.96 * pr$se.fit,
            hi = pr$fit + 1.96 * pr$se.fit, pool = label)][]
}
pw  <- d[power == "full"]
tab <- rbind(adj_means(pw, "power"),
             { db <- copy(d); f <- lm(pai ~ factor(date) + lai_als + power, data = db)
               nd <- db[, .(n = .N), by = date][, `:=`(lai_als = mean(db$lai_als),
                                                       power = "full")]
               pr <- predict(f, nd, se.fit = TRUE)
               nd[, `:=`(pai_adj = pr$fit, lo = pr$fit - 1.96 * pr$se.fit,
                         hi = pr$fit + 1.96 * pr$se.fit, pool = "power+half")]
               nd[, power := NULL]; nd })
raw <- d[, .(n_raw = .N, n_power = sum(power == "full"),
             pai_med_power = median(pai[power == "full"]),
             lai_als_med = median(lai_als),
             pct_day = round(100 * mean(daylight, na.rm = TRUE))), by = date]
tab <- merge(tab, raw[, .(date, pai_med_power, lai_als_med, pct_day)], by = "date")
setorder(tab, pool, date)
fwrite(tab, file.path(dir_tabs, "Table14_c4_blois_dates.csv"))
cat("=== Blois per-date composites (structure-adjusted to common LAI_ALS) ===\n")
print(tab[, .(pool, date, n, pct_day, pai_med_power = round(pai_med_power, 2),
              lai_als_med = round(lai_als_med, 2), pai_adj = round(pai_adj, 2),
              ci = sprintf("[%.2f, %.2f]", lo, hi))])

# --- formal contrasts (power beams, leaf-on) ----------------------------------
pw_on <- pw[leafon == TRUE]
f <- lm(pai ~ factor(date) + lai_als, data = pw_on)
ct <- function(d1, d2) {
  b <- coef(f); v <- vcov(f)
  n1 <- paste0("factor(date)", d1); n2 <- paste0("factor(date)", d2)
  g1 <- if (n1 %in% names(b)) b[[n1]] else 0
  g2 <- if (n2 %in% names(b)) b[[n2]] else 0
  vv <- sum(v[intersect(c(n1, n2), rownames(v)), intersect(c(n1, n2), colnames(v))] *
              matrix(c(1, -1)[c(n1 %in% rownames(v), n2 %in% rownames(v))] %o%
                     c(1, -1)[c(n1 %in% colnames(v), n2 %in% colnames(v))],
                     nrow = sum(c(n1, n2) %in% rownames(v))))
  est <- g1 - g2; se <- sqrt(abs(vv))
  sprintf("%s - %s = %+.2f +- %.2f (z = %.1f)", d1, d2, est, 1.96 * se, est / se)
}
cat("\n--- contrasts (adjusted, power beams) ---\n")
cat(ct("2021-06-06", "2021-09-02"), " # June dip vs Sep 2021\n")
cat(ct("2022-09-10", "2021-09-02"), " # drought Sep 2022 vs Sep 2021\n")
cat(ct("2022-10-09", "2022-09-10"), " # early autumn 2022\n")
cat(ct("2022-05-14", "2021-09-02"), " # May 2022 vs Sep 2021\n")

# --- figure: adjusted per-date trajectory + S2 series -------------------------
s2 <- as.data.table(readRDS(here::here("output", "intermediate",
                                       "blois_s2_lai_ts_2021.rds")))
s2 <- s2[, .(lai_s2 = mean(LAI_S2_ATBD, na.rm = TRUE)), by = date]
pt <- tab[pool == "power"][, `:=`(yr = format(date, "%Y"),
                                  doy = as.integer(format(date, "%j")))]
p <- ggplot(pt, aes(doy, pai_adj, colour = yr)) +
  geom_line(data = s2[, .(doy = as.integer(format(date, "%j")), pai_adj = lai_s2,
                          yr = "S2 2021")], linewidth = 0.7) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.5) +
  geom_line(linetype = 3) +
  scale_colour_manual(values = c(`2021` = "#7570b3", `2022` = "#1b9e77",
                                 `S2 2021` = "#d95f02"), name = NULL) +
  labs(x = "Day of year", y = expression("PAI / LAI (" * m^2 * "/" * m^2 * ")"),
       title = "Blois: per-acquisition GEDI composites, structure-adjusted",
       subtitle = paste("Points: GEDI power beams, adjusted to common ALS LAI",
                        "(95% CI). Orange: S2 ATBD annual series (2021).")) +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig8_c4_blois_summer.png"), p,
       width = 8.5, height = 4.8, dpi = 200)
cat("\nDONE\n")
