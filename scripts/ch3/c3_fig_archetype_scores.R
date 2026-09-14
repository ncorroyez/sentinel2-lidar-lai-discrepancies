# ==============================================================================
# c3_fig_archetype_scores.R
#
# Chapter 3 figure: between-plot R2 of dTmax by canopy archetype, for the five
# scenarios of the chapter, with bootstrap 95 % intervals on all of them.
#
# All five scenarios now cover the 53 plots, the dynamic combination included:
# its annual series was rebuilt on the unmasked Sentinel-2 product and the runs
# redone (c3_rebuild_annual_notmasked.R, c3_rerun_dyn_annual_53.R). The "n.d."
# branch below is kept as a guard, not because a cell is expected to be empty.
#
# The archetype reading replaces the earlier stratification at a near-median
# leaf-area threshold. It shows the one contrast the sample supports without
# qualification: Sentinel-2 collapses in the dense archetype P4, where its
# interval does not overlap the LiDAR interval. Elsewhere the intervals overlap.
#
# Input : manuscripts/ch3/tables/Table_archetype_scores_CHS41.csv (c3_archetype_scores.R)
# Output: manuscripts/ch3/figures/Fig_archetype_scores.png
#         Admin/ED_GAIA/CSI3/figures/fig08_ch3_archetypes.png (report copy)
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

ROOT <- "/home/corroyez/Documents/NC_Full"
TAB  <- file.path(ROOT, "manuscripts/ch3/tables/Table_archetype_scores_CHS41.csv")
OUT  <- file.path(ROOT, "manuscripts/ch3/figures/Fig_archetype_scores.png")
RPT  <- "/home/corroyez/Documents/Admin/ED_GAIA/CSI3/figures/fig08_ch3_archetypes.png"

dt <- fread(TAB)

# Les cinq scenarios du chapitre. Tous portent leur intervalle bootstrap ; les
# deux capteurs seuls, qui portent le contraste, sont traces plus epais et
# etiquetes.
SCEN <- list(
  list(id = "STATIC_ALS",      lab = "LiDAR fixed",         col = "#1F6FA8", hi = TRUE),
  list(id = "STATIC_S2_ATBD",  lab = "Sentinel-2 alone",    col = "#D2691E", hi = TRUE),
  list(id = "STATIC_S2_OPT",   lab = "Sentinel-2 opt",      col = "#E8A87C", hi = FALSE),
  list(id = "DYN_S2_ANNUAL",   lab = "Combination, dynamic", col = "#6FA8DC", hi = FALSE),
  list(id = "RF_RAWBANDS_LOO", lab = "Random forest",       col = "#8FBF7F", hi = FALSE))

groups <- c("P1", "P2", "P3", "P4")
labs   <- c("P1\nopen\nn = 8", "P2\nn = 12", "P3\nn = 13", "P4\ndense\nn = 20")
d <- dt[group %in% groups]

draw <- function() {
  par(mar = c(3.6, 4.2, 2.6, 0.8))
  plot(NA, xlim = c(0.4, 4.6), ylim = c(0, 1),
       xaxt = "n", xlab = "", ylab = "", las = 1, bty = "n")
  mtext(expression("Between-plot " * italic(R)^2 * " on " * Delta * T[max]),
        side = 2, line = 2.7, cex = 0.95)
  abline(h = seq(0, 1, 0.2), col = "grey92", lwd = 0.8)
  rect(3.62, 0, 4.38, 1, col = adjustcolor("#C9A87C", 0.10), border = NA)

  n <- length(SCEN)
  off <- seq(-0.30, 0.30, length.out = n)
  for (i in seq_along(groups)) {
    for (k in seq_len(n)) {
      sc  <- SCEN[[k]]
      row <- d[group == groups[i] & scenario == sc$id]
      x   <- i + off[k]
      # Six des huit placettes P1 n'ont pas tourne sous les scenarios dynamiques :
      # deux placettes ne definissent pas un R2. On marque la case au lieu de
      # laisser un blanc que le lecteur prendrait pour un oubli.
      if (!nrow(row) || !is.finite(row$R2)) {
        text(x, 0.035, "n.d.", cex = 0.6, col = adjustcolor(sc$col, 0.9), font = 3)
        next
      }
      lw  <- if (sc$hi) 2.2 else 1.1
      cap <- if (sc$hi) 0.045 else 0.028
      segments(x, row$R2_lo, x, row$R2_hi,
               col = if (sc$hi) sc$col else adjustcolor(sc$col, 0.75),
               lwd = lw, lend = 1)
      segments(x - cap, row$R2_lo, x + cap, row$R2_lo,
               col = adjustcolor(sc$col, if (sc$hi) 1 else 0.75), lwd = lw * 0.75)
      segments(x - cap, row$R2_hi, x + cap, row$R2_hi,
               col = adjustcolor(sc$col, if (sc$hi) 1 else 0.75), lwd = lw * 0.75)
      points(x, row$R2, pch = 21, bg = sc$col, col = "white",
             cex = if (sc$hi) 1.7 else 1.15, lwd = if (sc$hi) 1.6 else 1.1)
      if (sc$hi)
        text(x, row$R2, sprintf("%.2f", row$R2), pos = if (k == 1) 2 else 4,
             cex = 0.7, col = "grey20", offset = 0.5)
    }
  }
  axis(1, at = 1:4, labels = labs, tick = FALSE, line = 0.4, cex.axis = 0.82,
       col.axis = "grey20", padj = 0.5)
  legend("topright",
         legend = vapply(SCEN, function(z) z$lab, ""),
         pch = 21, pt.bg = vapply(SCEN, function(z) z$col, ""),
         col = "white",
         pt.cex = vapply(SCEN, function(z) if (z$hi) 1.5 else 1.05, 0),
         pt.lwd = 1.3, bty = "n", cex = 0.78, ncol = 3,
         text.font = vapply(SCEN, function(z) if (z$hi) 2L else 1L, 0L))
  mtext("Sentinel-2 fails where the canopy closes; elsewhere the sample cannot separate the sensors",
        side = 3, line = 0.7, adj = 0, cex = 0.86, font = 2, col = "grey15")
}

png(OUT, width = 9.2, height = 4.8, units = "in", res = 300, type = "cairo")
draw(); dev.off()
file.copy(OUT, RPT, overwrite = TRUE)
message("wrote: ", OUT)
message("wrote: ", RPT)
