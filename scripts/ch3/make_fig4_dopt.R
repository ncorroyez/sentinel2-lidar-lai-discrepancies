# Chapter 3 — Fig 4: cross-site LiDAR-LAI recovery vs d_opt depth. ARTICLE style.
# Target: Remote Sensing of Environment (single column 85 mm; vector PDF). Reads Table3c.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here("scripts","_article_style.R"))
t <- fread(here("manuscripts/ch3","tables","Table3c_dopt_depth_sweep.csv"))
d_persite <- c(Aigoual=8, Blois=6, Mormal=10)
fix <- melt(t, id.vars=c("site","S2_FORMSH"),
            measure.vars=c("dop_5","dop_common7","dop_10","dop_15","dop_20"),
            variable.name="d", value.name="r2")
fix[, depth := as.numeric(c(dop_5=5,dop_common7=7,dop_10=10,dop_15=15,dop_20=20)[as.character(d)])]
ps <- t[, .(site, depth=d_persite[site], r2=dop_per_site)]
base <- t[, .(site, S2_FORMSH)]
g <- ggplot() +
  annotate("rect", xmin=6, xmax=8, ymin=-Inf, ymax=Inf, fill="grey88", alpha=.6) +
  annotate("text", x=7, y=0.06, label="optical~italic(d)[opt]", parse=TRUE, size=2.7, colour="grey35") +
  geom_hline(data=base, aes(yintercept=S2_FORMSH, colour=site), linetype="dotted", linewidth=0.6) +
  geom_line(data=fix, aes(depth, r2, colour=site), linewidth=0.9) +
  geom_point(data=fix, aes(depth, r2, colour=site), size=1.8) +
  geom_point(data=ps, aes(depth, r2, colour=site, shape="per-site optical"), size=3) +
  scale_colour_manual(values=PAL_SITE, name="Held-out site") +
  scale_shape_manual(values=c("per-site optical"=17), name=NULL) +
  scale_x_continuous(breaks=c(5,7,10,15,20)) +
  labs(x=expression(italic(d)[opt]~"depth (m, top of canopy)"),
       y=expression(R^2~"(predicting LiDAR LAI, held-out site)")) +
  coord_cartesian(ylim=c(0,1)) + theme_article(11) + legend_corner(0.98,0.02)
outdir <- here("outputs","figures_article_ch3")
saveRDS(list(fix=fix, ps=ps, base=base), file.path(outdir,"fig4_dopt_data.rds"))
ggsave_article(file.path(outdir,"Fig4_dopt_transfer"), g, 6.7, 4.2)
ggsave_article(here("manuscripts/ch3","figures","Fig4_dopt_transfer_styled"), g, 6.7, 4.2)
cat("DONE Fig4\n")
