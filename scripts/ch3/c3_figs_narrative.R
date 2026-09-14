# ==============================================================================
# Chapter 3 — narrative-numbered figures (English), NC_Full side.
#   Fig 3  magnitude recovery (Table6 LOO): mean ok, per-plot pattern not
#   Fig S1 annual S2 LAI series + smoother (secondary, §3.5)
#   Fig S2 per-plot corrected trajectories (secondary, §3.5)
# Fig 1/2 (saturation) and Fig 4 (d_opt transfer) are produced by
# c3_saturation_3site.R and c3_fig5_transfer.R; this script renames them too.
# Run from NC_Full root:  Rscript scripts/ch3/c3_figs_narrative.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
figdir <- here("manuscripts/ch3", "figures")
outfig <- here("output", "figures")
save2 <- function(g, name, w, h) {
  for (d in c(figdir, outfig)) ggsave(file.path(d, name), g, width = w, height = h, dpi = 150)
}

# ── Fig 3 — magnitude recovery (leave-one-plot-out, Table6) ──────────────────
t6 <- fread(here("manuscripts/ch3", "tables", "Table6_LAIrecovery_LOO.csv"))
lab <- c(raw_s2 = "Raw S2", global_rescale = "S2 x mean ratio (ratio-A)",
         rf_s2only = "RF (S2 only)", dopt_physics = "S2 x H/d_opt",
         rf_struct_only = "RF (LiDAR structure)", rf_full = "RF (S2 + LiDAR structure)",
         dopt_lai_rescale = "S2 x ALS/ALS_dopt")
t6 <- t6[method != "dopt_lai_rescale"]                       # degenerate (RMSE 9)
t6[, lab := lab[method]]
t6[, grp := ifelse(applicable_beyond_lidar, "S2-only (works beyond LiDAR)",
                   "needs LiDAR-derived inputs (oracle)")]
t6[, lab := reorder(lab, r2)]
g3 <- ggplot(t6, aes(r2, lab, fill = grp)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", r2)), hjust = -0.15, size = 3.3) +
  scale_fill_manual(values = c("S2-only (works beyond LiDAR)" = "#E8746A",
                               "needs LiDAR-derived inputs (oracle)" = "#9B5DE5"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 0.82), expand = c(0, 0)) +
  labs(title = "Recovering per-plot LiDAR LAI: S2 fixes the mean, not the pattern",
       subtitle = "Leave-one-plot-out R² (Blois, n=47). S2-only methods (incl. the ratio) stay near R²≈0.1;\nonly LiDAR-derived structure recovers the spatial pattern — an ALS-bound oracle.",
       x = expression(R^2~"(predicting LiDAR LAI)"), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
save2(g3, "Fig3_magnitude_recovery.png", 8, 4.3)

# ── Fig S1 — annual S2 LAI series + smoothed trajectory ──────────────────────
ts  <- as.data.table(readRDS(here("output","intermediate","blois_s2_lai_ts_2021.rds")))
cor <- as.data.table(readRDS(here("output","intermediate","blois_s2_lai_corrected_daily_2021.rds")))
ts_summ <- ts[, .(lai = mean(LAI_S2_ATBD, na.rm = TRUE)), by = doy]
gS1 <- ggplot() +
  geom_point(data = ts_summ, aes(doy, lai), colour = "#1FA3A3", size = 2) +
  geom_line(data = cor, aes(doy, LAI_S2_ATBD__raw), colour = "#0B6E6E", linewidth = 1) +
  labs(title = "Annual Sentinel-2 LAI series, Blois 2021 (secondary, §3.5)",
       subtitle = "Points: 19 inverted S2 dates (plot mean). Line: cyclic Whittaker smoothing.",
       x = "Day of year", y = expression(LAI[S2_ATBD]~(m^2/m^2))) +
  theme_minimal(base_size = 11)
save2(gS1, "FigS1_annual_series.png", 8, 4)

# ── Fig S2 — per-plot corrected trajectories (3 magnitude variants) ──────────
pp <- as.data.table(readRDS(here("output","intermediate","blois_perplot_daily_v2_2021.rds")))
ppl <- melt(pp, id.vars = c("id","doy"),
            measure.vars = c("rescaled_ALS","rescaled_dopt","ML_rf"),
            variable.name = "variant", value.name = "lai")
vlab <- c(rescaled_ALS = "x full LiDAR LAI", rescaled_dopt = "x LiDAR LAI (d_opt)", ML_rf = "RF (structure+S2)")
ppl[, variant := vlab[as.character(variant)]]
gS2 <- ggplot(ppl, aes(doy, lai, group = id)) +
  geom_line(alpha = 0.18, colour = "#2EA02E") +
  facet_wrap(~ variant) +
  labs(title = "Per-plot corrected LAI trajectories, Blois 2021 (secondary, §3.5)",
       subtitle = "Greenness-fraction-rescaled S2 shape x three magnitude anchors. One line per plot (n=47).",
       x = "Day of year", y = expression(LAI~(m^2/m^2))) +
  theme_minimal(base_size = 11)
save2(gS2, "FigS2_perplot_corrections.png", 9, 3.6)

cat("DONE Fig3, FigS1, FigS2\n")
