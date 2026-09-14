# Chapter 3 — summer 53-HOBO validation: R² by scenario, 2 metrics × 2 periods.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })
source(here("scripts","_article_style.R"))
t <- fread(here("manuscripts/ch3","tables","Table4f_summer_2metrics_2periods.csv"))
lab <- c(STATIC_ALS="LiDAR full (static)", STATIC_ALS_DOPT="LiDAR d_opt (static)",
  DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
  DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full (dyn)", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt (dyn)")
t[, lbl := lab[scenario]]
m <- melt(t, id.vars=c("lbl","period"), measure.vars=c("dT_r2","sl_r2"), variable.name="metric", value.name="R2")
m[, metric := ifelse(metric=="dT_r2","ΔT_max","slope-and-equilibrium")]
m[, period := ifelse(period=="all","all summer days","hottest days (Tmacro ≥ p90)")]
ord <- t[period=="all"][order(dT_r2), lbl]
m[, lbl := factor(lbl, levels=ord)]
m[, grp := fifelse(grepl("LiDAR", lbl), "LiDAR", "S2")]
g <- ggplot(m, aes(R2, lbl, fill=grp)) +
  geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%.2f", R2)), hjust=-0.15, size=2.6) +
  facet_grid(metric ~ period, switch="y") +
  scale_fill_manual(values=c(LiDAR="#1A9850", S2="#FDAE61"), name=NULL) +
  scale_x_continuous(limits=c(0,1.0), expand=c(0,0)) +
  labs(x=expression(R^2~"(simulated vs observed, n = 53 HOBO)"), y=NULL) +
  theme_article(11) + theme(legend.position="bottom", panel.spacing=unit(0.6,"lines"))
outdir <- here("outputs","figures_article_ch3")
ggsave_article(file.path(outdir,"Fig_validation_2metrics_2periods"), g, 9, 6.2)
ggsave_article(here("manuscripts/ch3","figures","Fig_validation_2metrics_2periods"), g, 9, 6.2)
cat("DONE validation figure\n")
