# ==============================================================================
# c3_fig0_lai_sensitivity_styled_CHS41.R
#
# Figure 1 of the Chapter 3 supervisor document, restyled from the CHS41 run.
#
# The figure in the document dates from June 2026 and was produced under
# in_files/FR-Blo_2021.nc, the forcing Chapter 3 used before the move to the
# CHS41 hybrid. Every other number in the document now comes from CHS41, so the
# figure was the last piece on the old forcing. Under CHS41 the response is
# about three times larger above the Sentinel-2 saturation ceiling, so the
# argument the figure makes is stronger, not weaker.
#
# Style follows the other "_styled" figures of the document: panel tags (a) and
# (b), no plot title, grey band over the saturation range, green for the daytime
# offset and blue for the coupling slope.
#
# Input : manuscripts/ch3/tables/Table0c_lai_sensitivity_CHS41.csv
#           (written by z_Example_rmusica/c3_lai_sensitivity_CHS41.R)
# Output: manuscripts/ch3/figures/Fig0_lai_sensitivity_styled_CHS41.png
#   Rscript scripts/ch3/c3_fig0_lai_sensitivity_styled_CHS41.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork) })

ROOT <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3"
d <- fread(file.path(ROOT, "tables/Table0c_lai_sensitivity_CHS41.csv"))
SAT <- 4                                   # Sentinel-2 saturation ceiling

band <- function() annotate("rect", xmin = SAT, xmax = Inf, ymin = -Inf, ymax = Inf,
                            fill = "grey92", alpha = 0.9)

base <- function(g, ylab, tag) {
  g + band() + geom_line(linewidth = 0.9) + geom_point(size = 1.9) +
    # expression() lirait "one-sided" comme une soustraction : on compose la chaine.
    labs(x = expression("LAI (one-sided, " * m^2 ~ m^-2 * ")"), y = ylab, tag = tag) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
          plot.tag = element_text(face = "bold", size = 15),
          plot.tag.position = c(0.03, 0.97))
}

ga <- base(ggplot(d, aes(LAI, dTmax)) +
             geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60"),
           expression(Delta * italic(T)[max] ~ (degree * C)), "(a)") +
  scale_colour_identity() + aes(colour = "#2E9B57") +
  annotate("text", x = (SAT + max(d$LAI)) / 2, y = max(d$dTmax) * 0.75,
           label = "S2 saturation\n(LAI > 4)", size = 3.6, colour = "grey35", lineheight = 1.05)

gb <- base(ggplot(d, aes(LAI, slope)), "micro–macro coupling slope", "(b)") +
  scale_colour_identity() + aes(colour = "#2C7FB8") +
  annotate("text", x = (SAT + max(d$LAI)) / 2, y = max(d$slope) - 0.02,
           label = "S2 saturation\n(LAI > 4)", size = 3.6, colour = "grey35", lineheight = 1.05)

out <- file.path(ROOT, "figures/Fig0_lai_sensitivity_styled_CHS41.png")
ggsave(out, ga | gb, width = 10.0, height = 4.0, dpi = 300, bg = "white")

span <- d[LAI == SAT]$dTmax - d[LAI == max(LAI)]$dTmax
cat(sprintf("dTmax : %+.2f a LAI 1, pic %+.2f, %+.2f a LAI 9 ; etendue au-dela de LAI %d : %.1f degC\n",
            d[LAI == 1]$dTmax, max(d$dTmax), d[LAI == max(LAI)]$dTmax, SAT, span))
cat(sprintf("pente : %.3f a LAI 1 -> %.3f a LAI 9\n", d[LAI == 1]$slope, d[LAI == max(LAI)]$slope))
message("wrote: ", out)
