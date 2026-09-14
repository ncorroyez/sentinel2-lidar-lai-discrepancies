# ==============================================================================
# c4_gedi_crosssite_shape.R — Chapter 4 / Appendix G follow-ups (data only):
#  (1) CROSS-SITE TRANSFERABILITY of the GEDI seasonal shape: build the same
#      leaf-fraction curve at Blois and Mormal (adjusted composites -> weekly
#      interp -> Whittaker -> fraction of amplitude) and quantify their overlap.
#      If the two curves match, a site-agnostic "regional GEDI phenology" is
#      plausible (the perspectives claim).
#  (2) DENSITY-STRATIFIED EOS at Mormal: does denser canopy senesce later?
#      Stand-scale test of the canopy-structure -> autumn-phenology mechanism
#      (Wu et al. 2024). TRS50 on per-stratum adjusted monthly composites.
# Outputs: manuscripts/ch4/tables/Table21*, figures/Fig10_c4_crosssite_shape.png
# ==============================================================================
suppressMessages({
  library(data.table); library(sf); library(ggplot2); library(phenofit); library(here)
})
source(here::here("scripts", "_article_style.R"))
k_scale  <- 0.5 / 0.65
dir_tabs <- here::here("manuscripts/ch4", "tables")
dir_figs <- here::here("manuscripts/ch4", "figures")

load_site <- function(s) {
  d <- as.data.table(sf::st_drop_geometry(
    readRDS(here::here("01_DATA", s, "GEDI", paste0(s, "_GEDI_LiDAR_Metrics.rds")))))
  d[, `:=`(site = s, m = as.integer(as.character(month)),
           date = as.Date(sprintf("%s-%s-%s", year, month, day)))]
  d <- d[!is.na(pai) & pai > 0 & sensitivity >= 0.90 & power == "full" &
           !is.na(LAI_lidar_cor) & LAI_lidar_cor >= 2]
  d[, lai_als := LAI_lidar_cor * k_scale]
  d[site == "Blois" & date == as.Date("2021-05-29"), keep := FALSE]  # n=5 outlier
  d[is.na(keep), keep := TRUE]; d <- d[keep == TRUE]
  d
}
dt <- rbindlist(lapply(c("Blois", "Mormal"), load_site), fill = TRUE)

# --- (1) site-level leaf-fraction curves --------------------------------------
# per-date (Blois) / per-month (Mormal) adjusted composites, then fraction of
# the smoothed amplitude, exactly the construction used for the scenarios
frac_curve <- function(dd, by_date = FALSE) {
  dd <- copy(dd)
  dd[, tvar := if (by_date) as.integer(format(date, "%j")) else round(30.4 * (m - 0.5))]
  f  <- lm(pai ~ factor(tvar) + lai_als, data = dd)
  nd <- dd[, .(n = .N), by = tvar][, lai_als := mean(dd$lai_als)]
  nd[, pai_adj := predict(f, nd)]
  setorder(nd, tvar)
  wk <- seq(min(nd$tvar), max(nd$tvar), by = 7)
  sm <- as.numeric(whit2(approx(nd$tvar, nd$pai_adj, wk)$y, lambda = 10))
  data.table(doy = wk, frac = pmin(pmax((sm - min(sm)) / diff(range(sm)), 0), 1))
}
fb <- frac_curve(dt[site == "Blois"], by_date = TRUE)[, site := "Blois"]
fm <- frac_curve(dt[site == "Mormal"])[, site := "Mormal"]

common <- seq(max(min(fb$doy), min(fm$doy)), min(max(fb$doy), max(fm$doy)), 1)
vb <- approx(fb$doy, fb$frac, common)$y; vm <- approx(fm$doy, fm$frac, common)$y
cat(sprintf("shape overlap (DOY %d-%d): r = %.3f | mean |diff| = %.3f | max |diff| = %.3f\n",
            min(common), max(common), cor(vb, vm), mean(abs(vb - vm)), max(abs(vb - vm))))
fwrite(rbind(fb, fm), file.path(dir_tabs, "Table21_c4_site_shapes.csv"))

p <- ggplot(rbind(fb, fm), aes(doy, frac, colour = site)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(Blois = "#7570b3", Mormal = "#1b9e77"), name = NULL) +
  labs(x = "Day of year", y = "Leaf fraction (of seasonal amplitude)",
       title = "Is the GEDI seasonal shape transferable across sites?",
       subtitle = paste("Structure-adjusted, Whittaker-smoothed leaf-fraction",
                        "curves; Blois (oak) vs Mormal (oak-beech), ~250 km apart")) +
  theme_article() + theme(legend.position = "top")
ggsave(file.path(dir_figs, "Fig10_c4_crosssite_shape.png"), p,
       width = 7.5, height = 4, dpi = 200)

# --- (2) density-stratified EOS at Mormal -------------------------------------
mo <- dt[site == "Mormal"]
mo[, dens_bin := cut(lai_als, c(2, 3, 4, 9), labels = c("2-3", "3-4", ">4"),
                     include.lowest = TRUE)]
trs50_strat <- function(dd) {
  cm <- dd[, .(pai_med = median(pai), n = .N), by = m][order(m)]
  wk <- seq(30.4 * (min(cm$m) - 0.5), 30.4 * (max(cm$m) - 0.5), by = 7)
  sm <- as.numeric(whit2(approx(30.4 * (cm$m - 0.5), cm$pai_med, wk)$y, lambda = 10))
  half <- min(sm) + diff(range(sm)) / 2
  up <- which(sm >= half)
  list(n = nrow(dd), n_months = nrow(cm),
       SOS = round(wk[up[1]]), EOS = round(wk[rev(up)[1]]),
       amp = round(diff(range(sm)), 2))
}
tab_eos <- mo[!is.na(dens_bin), trs50_strat(.SD), by = dens_bin][order(dens_bin)]
cat("\n=== Mormal: TRS50 by density stratum (raw medians; composition within\n",
    "    stratum is narrower so no extra adjustment) ===\n")
print(tab_eos)
fwrite(tab_eos, file.path(dir_tabs, "Table21b_c4_eos_by_density.csv"))
cat("\nDONE\n")
