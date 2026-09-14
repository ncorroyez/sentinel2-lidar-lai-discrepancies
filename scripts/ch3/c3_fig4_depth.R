# ==============================================================================
# Chapter 3 — Fig 4 (revised): cross-site LiDAR-LAI recovery vs d_opt DEPTH.
# Held-out R² rises monotonically with the depth of LiDAR LAI injected as a
# feature -> the "d_opt benefit" is really "more LiDAR helps", trivial toward the
# full-canopy oracle at deep d. The optical d_opt (~6-8 m, what S2 senses) is the
# shallowest, least-circular point — still ALS-bound. Reads Table3c.
# Run from NC_Full root:  Rscript scripts/ch3/c3_fig4_depth.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
t <- fread(here("manuscripts/ch3","tables","Table3c_dopt_depth_sweep.csv"))
d_persite <- c(Aigoual = 8, Blois = 6, Mormal = 10)

# fixed-depth line
fix <- melt(t, id.vars = c("site","S2_FORMSH"),
            measure.vars = c("dop_5","dop_common7","dop_10","dop_15","dop_20"),
            variable.name = "d", value.name = "r2")
fix[, depth := as.numeric(c(dop_5=5, dop_common7=7, dop_10=10, dop_15=15, dop_20=20)[as.character(d)])]
# per-site optical-depth point
ps <- t[, .(site, depth = d_persite[site], r2 = dop_per_site)]
base <- t[, .(site, S2_FORMSH)]
cols <- c(Aigoual="#E8746A", Blois="#2EA02E", Mormal="#5B8FF0")

g <- ggplot() +
  annotate("rect", xmin = 6, xmax = 8, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = .5) +
  annotate("text", x = 7, y = 0.04, label = "optical d_opt\n(what S2 senses)", size = 2.8, colour = "grey35") +
  geom_hline(data = base, aes(yintercept = S2_FORMSH, colour = site), linetype = "dotted") +
  geom_line(data = fix, aes(depth, r2, colour = site), linewidth = 0.9) +
  geom_point(data = fix, aes(depth, r2, colour = site), size = 2) +
  geom_point(data = ps, aes(depth, r2, colour = site), shape = 17, size = 3.4) +
  scale_colour_manual(values = cols, name = "Held-out site") +
  scale_x_continuous(breaks = c(5,7,10,15,20)) +
  labs(title = "Recovering LiDAR LAI: the d_opt 'benefit' is monotone with depth",
       subtitle = "Held-out R² (S2 + FORMS-H + d_opt depth) climbs to the full-canopy oracle as d_opt deepens.\nDotted = S2 + FORMS-H baseline; triangle = per-site optical d_opt; grey band = optical depth. A non-ALS\npredicted d_opt stays at baseline (Table 3b) — only the TRUE (ALS) depth carries information.",
       x = expression(d[opt]~"depth (m, top of canopy)"),
       y = expression(R^2~"(predicting LiDAR LAI, held-out site)")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 11) + theme(plot.subtitle = element_text(size = 8))

for (dd in c(here("output","figures"), here("manuscripts/ch3","figures")))
  ggsave(file.path(dd, "Fig4_dopt_transfer.png"), g, width = 8.5, height = 5.2, dpi = 150)
cat("DONE\n")
