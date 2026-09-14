# ==============================================================================
# Fig K — regime-switching hybrid (LiDAR dense + S2 open, gated at LAI 3.86, near the
# median LiDAR LAI of the 53 plots; this is a chosen split, not a Ch2-derived crossover).
# Shows (i) the hybrid beats either sensor alone on ΔTmax RMSE / bias / amplitude,
# (ii) correcting S2 per Ch2 (d_opt / rescaled) in the OPEN regime HURTS: raw ATBD
# is the right S2 to use where S2 is used. Data: figK_hybrid_ch2.csv (c3_hybrid_ch2_wc.R).
#   Rscript scripts/ch3/make_ch3_figK.R   (from NC_Full root)
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("scripts/_article_style.R")
source("scripts/_figure_label_check.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
K<-fread(file.path(TAB,"figK_hybrid_ch2.csv"))
# keep the LAI-crossover oracle switch + the two single-sensor baselines
keep<-c("LiDAR full (pooled)","S2 ATBD (pooled)",
        "Hybrid LAI-abs: LiDAR+S2 ATBD","Hybrid LAI-abs: LiDAR+S2 d_opt","Hybrid LAI-abs: LiDAR+S2 rescaled")
lab<-c("LiDAR full (pooled)"="LiDAR alone","S2 ATBD (pooled)"="Sentinel-2 alone",
       "Hybrid LAI-abs: LiDAR+S2 ATBD"="Hybrid: +S2 raw (ATBD)",
       "Hybrid LAI-abs: LiDAR+S2 d_opt"="Hybrid: +S2 d_opt-corr.",
       "Hybrid LAI-abs: LiDAR+S2 rescaled"="Hybrid: +S2 rescaled")
grp<-c("LiDAR alone"="single","Sentinel-2 alone"="single","Hybrid: +S2 raw (ATBD)"="best",
       "Hybrid: +S2 d_opt-corr."="corrected","Hybrid: +S2 rescaled"="corrected")
K<-K[product%in%keep]; K[,prod:=factor(lab[product],levels=rev(unname(lab)))]; K[,grp:=grp[as.character(prod)]]
Km<-melt(K,id.vars=c("prod","grp"),measure.vars=c("dT_R2","dT_RMSE","dT_SDrec"),variable.name="metric",value.name="v")
Km[,metric:=factor(metric,levels=c("dT_R2","dT_RMSE","dT_SDrec"),
   labels=c("Ranking R² (higher better)","RMSE °C (lower better)","Amplitude SD-rec. (higher better)"))]
pal<-c("single"="grey55","best"=PAL_SENSOR[["Fusion"]],"corrected"="#E6A33E")
check_labels(Km$v, "manuscripts/ch3/Chapter3_article_standalone_EN.md",
             what="Fig K — fusion par hauteur")

pK<-ggplot(Km,aes(v,prod,fill=grp))+geom_col(width=0.68)+
  geom_text(aes(label=fmt2(v)),hjust=-0.15,size=2.8,colour="grey25")+
  facet_wrap(~metric,scales="free_x")+scale_fill_manual(values=pal,guide="none")+
  scale_x_continuous(expand=expansion(mult=c(0,0.18)))+
  labs(x=NULL,y=NULL,
       title="Best microclimate product: LiDAR in dense + raw S2 in open (gated at LAI 3.86)",
       subtitle="Correcting S2 per Ch2 in the OPEN regime hurts — S2 is used exactly where it needs no correction")+
  theme_article(11)+theme(legend.position="none")
ggsave_article(file.path(FIG,"FigK_hybrid_ch2"),pK,9.4,3.8)
cat("DONE: FigK_hybrid_ch2 -> manuscripts/ch3/figures/\n")
