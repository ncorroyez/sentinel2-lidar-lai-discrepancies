# ==============================================================================
# c3_fig_s2_series_archetype.R
#
# Chapter 3 figure: the 2021 Sentinel-2 LAI trajectory read archetype by
# archetype, in the retrieval the chapter forces MuSICA with.
#
# The scenario table scores the whole season in one number per plot. This figure
# shows where that number comes from: whether the four archetypes separate along
# the year at all. It
# is the temporal counterpart of the archetype R2 figure.
#
# Mean +- SD per archetype per date, following the temporal LAI plot shown at
# CSI 2, restricted to the June-September 2021 leaf-on window, which is the one
# the scenarios are scored over.
#
# Inputs : output/intermediate/blois_s2_lai_ts_2021.rds   (19 dates x 60 plots)
#          <rmusica>/out_files/Chapter3/lai_prep/df_plots_real53.rds  (ID -> id_plot)
#          <rmusica>/out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv (archetype)
# Output : manuscripts/ch3/figures/Fig_s2_series_archetype.png
#          Admin/ED_GAIA/CSI3/figures/fig09_ch3_s2_series.png (deck copy)
#   Rscript scripts/ch3/c3_fig_s2_series_archetype.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2) })

RM   <- "/home/corroyez/Documents/z_Example_rmusica_31012025"
ROOT <- "/home/corroyez/Documents/NC_Full"
OUT  <- file.path(ROOT, "manuscripts/ch3/figures/Fig_s2_series_archetype.png")
DECK <- "/home/corroyez/Documents/Admin/ED_GAIA/CSI3/figures/fig09_ch3_s2_series.png"

LV  <- c("P1", "P2", "P3", "P4")
PAL <- c(P1 = "#D7191C", P2 = "#FDAE61", P3 = "#74C476", P4 = "#1A9850")
ALS_DATE <- as.Date("2021-06-17")          # Blois leaf-on flight (Table 1)
SUMMER   <- as.Date(c("2021-06-01", "2021-09-30"))

ts <- as.data.table(readRDS(file.path(ROOT, "output/intermediate/blois_s2_lai_ts_2021.rds")))
ts[, date := as.Date(date)]

# ID (1..60) is the plot index the series carries; df_plots_real53 maps it to the
# plot name, which the Chapter 1 cluster table keys the archetype on.
key <- as.data.table(readRDS(file.path(RM, "out_files/Chapter3/lai_prep/df_plots_real53.rds")))[
        , .(plot = as.integer(ID), id_plot = as.character(id_plot))]
cl  <- fread(file.path(RM, "out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv"))[
        , .(id_plot = as.character(id_plot), P)]
d <- merge(merge(ts, key, by = "plot"), cl, by = "id_plot")
d[, P := factor(P, levels = LV)]
# Only the leaf-on window the chapter scores: 14 of the 19 dates fall in it, and
# the leaf-off tail compresses them against each other.
d <- d[date %between% SUMMER]

n_arch <- d[date == d$date[1], .N, by = P][order(P)]
cat("plots per archetype:\n"); print(n_arch)
stopifnot(nrow(d) > 0, all(LV %in% levels(droplevels(d$P))))

lab <- setNames(sprintf("%s (n = %d)", n_arch$P, n_arch$N), n_arch$P)

# Mean +- SD per archetype per date, following the temporal LAI plot shown at
# CSI 2, restricted to the June-September 2021 leaf-on window that the scenarios
# are scored over.
#
# The acquisitions are irregular: three of them fall within five days in early
# June, then nothing for three weeks. On a calendar axis those three collapse
# into one illegible cluster while the rest of the panel sits empty, so the date
# is drawn as an ordered category, one slot per acquisition. The seasonal shape
# is slightly distorted; the comparison between archetypes, which is what the
# figure is for, becomes readable.
agg <- function(col) d[is.finite(get(col)),
  .(m = mean(get(col)), sd = sd(get(col)), n = .N), by = .(P, date)][order(P, date)]

a <- agg("LAI_S2_ATBD")
a[, slot := factor(format(date, "%d %b"), levels = format(sort(unique(a$date)), "%d %b"))]
dodge <- position_dodge(width = 0.62)

# P1 and P4 carry the comparison: they are the ends of the density gradient, and
# the chapter's claim is about what happens as the canopy closes. If even those
# two do not separate, the intermediate types cannot. P2 and P3 stay in the
# figure, muted, so the reader still sees the gradient is covered.
# L'emphase passe par l'epaisseur et l'opacite, pas par la couleur : P2 et P3
# gardent la leur, sinon le gradient de densite disparait de la figure.
a[, hi := P %in% c("P1", "P4")]

g <- ggplot(a, aes(slot, m, colour = P, group = P)) +
  geom_line(data = a[hi == FALSE], linewidth = 0.5, alpha = 0.5, position = dodge) +
  geom_point(data = a[hi == FALSE], size = 1.3, alpha = 0.6, position = dodge) +
  geom_errorbar(data = a[hi == TRUE], aes(ymin = m - sd, ymax = m + sd),
                width = 0.34, linewidth = 0.5, alpha = 0.8, position = dodge) +
  geom_line(data = a[hi == TRUE], linewidth = 1.0, position = dodge) +
  geom_point(data = a[hi == TRUE], size = 2.4, position = dodge) +
  # sans breaks explicites la legende suit l'ordre des couches (P2, P3, P1, P4)
  scale_colour_manual(values = PAL, labels = lab, breaks = LV, name = NULL) +
  labs(x = NULL, y = expression(LAI[S2] ~ (mean %+-% SD))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
        legend.position = "top", legend.margin = margin(b = 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT, g, width = 9.0, height = 4.2, dpi = 300, bg = "white")
file.copy(OUT, DECK, overwrite = TRUE)

# The summer plateau is what the static scenarios score; report it so the caption
# and the notes can quote a number rather than describe the picture.
pl <- d[date %between% SUMMER, .(ATBD = median(LAI_S2_ATBD, na.rm = TRUE),
                                 opt  = median(LAI_S2_opt,  na.rm = TRUE)), by = P][order(P)]
cat("\nsummer median LAI_S2 by archetype:\n"); print(pl)
message("wrote: ", OUT)
message("wrote: ", DECK)
