# Chapter 3 — Fig 0 (LAI sensitivity), Fig 4e (summer scenarios), Fig Sh (Shapley) — ARTICLE style.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here); library(patchwork) })
source(here("scripts","_article_style.R"))
T <- function(n) fread(here("manuscripts/ch3","tables",n))
outdir <- here("outputs","figures_article_ch3")
sv <- function(g,nm,w,h){ ggsave_article(file.path(outdir,nm),g,w,h); ggsave_article(here("manuscripts/ch3","figures",paste0(nm,"_styled")),g,w,h) }

# Fig 0 — LAI sensitivity sweep (one-sided LAI on x). Panel a: ΔTmax ; b: coupling slope.
t0 <- T("Table0c_lai_sensitivity.csv")
band <- annotate("rect", xmin=4, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey88", alpha=.6)
a <- ggplot(t0, aes(LAI, dTmax)) + band +
  annotate("text", x=6.6, y=1.2, label="S2 saturation\n(LAI > 4)", size=2.7, colour="grey35") +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60") + geom_line(linewidth=0.9, colour="#1A9850") + geom_point(size=1.4, colour="#1A9850") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.3, vjust=1.4, label="(a)", fontface="bold", size=4) +
  labs(x="LAI (one-sided, m² m⁻²)", y=expression(Delta*italic(T)[max]~(degree*C))) + theme_article(11)
b <- ggplot(t0, aes(LAI, slope)) + band + geom_line(linewidth=0.9, colour="#3182bd") + geom_point(size=1.4, colour="#3182bd") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.3, vjust=1.4, label="(b)", fontface="bold", size=4) +
  labs(x="LAI (one-sided, m² m⁻²)", y="micro–macro coupling slope") + theme_article(11)
sv(a+b, "Fig0_lai_sensitivity", 9, 3.6)

# Fig 4e — summer MuSICA scenarios: RMSE bars coloured by R²
t4 <- T("Table4e_summer_lai_scenarios.csv")
lab <- c(SUM_S2_ATBD="S2 ATBD (raw)", SUM_S2_opt="S2 opt (raw)", SUM_ATBD_corr="S2 ATBD ×ratio",
  SUM_opt_corr="S2 opt ×ratio", SUM_ALS="LiDAR full", SUM_ALS_dopt="LiDAR d_opt",
  SUM_S2_rescaled="S2 rescaled", SUM_hybrid="hybrid (layered)", SUM_RF="RF", SUM_ratioA="S2 × ratio-A")
t4[, lbl:=reorder(lab[scenario], -rmse)]
g4 <- ggplot(t4, aes(lbl, rmse, fill=r2)) + geom_col(width=0.7) +
  geom_text(aes(label=sprintf("R²=%.2f", r2)), vjust=-0.3, size=2.7) +
  scale_fill_viridis_c(option="D", name=expression(italic(R)^2)) +
  labs(x=NULL, y=expression(RMSE~Delta*italic(T)[max]~(degree*C))) +
  theme_article(11) + theme(axis.text.x=element_text(angle=25, hjust=1)) + coord_cartesian(ylim=c(0, max(t4$rmse)*1.08))
sv(g4, "Fig4e_summer_lai_scenarios", 8, 4.2)

# Fig Sh — per-point cLHS Shapley: factor contributions by cluster (boxplots)
sh <- T("Table6c_shapley_clhs_2x.csv")
ml <- melt(sh, id.vars=c("pid","Cluster"), measure.vars=c("LAI","Hmax","fCover","LAD"),
           variable.name="factor", value.name="phi")
ml[, Cluster:=paste0("Cluster ", Cluster)]
gS <- ggplot(ml, aes(factor, phi, fill=factor)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60") +
  geom_boxplot(outlier.size=0.3, alpha=0.85, linewidth=0.4) + facet_wrap(~Cluster, nrow=1) +
  scale_fill_brewer(palette="Set2", guide="none") +
  labs(x=NULL, y=expression(Shapley~phi~(degree*C~on~Delta*italic(T)[max]))) + theme_article(11)
sv(gS, "FigSh_shapley_clhs_clusters", 9, 3.4)
cat("DONE Fig0 / Fig4e / FigSh\n")
