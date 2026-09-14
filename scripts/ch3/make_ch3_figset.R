# ==============================================================================
# Chapter 3 figure set (v3.2.3, genuine-53). Honest case that LiDAR (static ALS)
# is the reference microclimate driver.
#   Fig A  LAI-space: which LAI product predicts observed buffering (LiDAR wins).
#   Fig B  MuSICA HOBO validation per scenario: R2 / RMSE / bias (ΔTmax).
#   Fig C  Magnitude axis: warm bias vs mean LAI fed (truncating to d_opt under-buffers).
#   Fig D  Why R2 misleads: R2 vs SD-recovery (amplitude) per scenario.
# Data from c3_figset_data.R (figA/figB/figC csv). English, article style.
#   Rscript scripts/ch3/make_ch3_figset.R   (run from NC_Full root)
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
source("scripts/_article_style.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
FAM<-c("LiDAR"="#0072B2","Sentinel-2"="#D55E00","Fusion"="#CC79A7","d_opt truncated"="#666666")

## ---- Fig A : LAI-space prediction of observed microclimate --------------------
A<-fread(file.path(TAB,"figA_laispace_cor.csv"))
Al<-melt(A,id.vars="product",measure.vars=c("cor_dTmax","cor_slope"),variable.name="metric",value.name="r")
Al[,metric:=factor(metric,levels=c("cor_dTmax","cor_slope"),labels=c("Daytime offset (ΔT[max])","Coupling slope"))]
Al[,product:=factor(product,levels=A[order(cor_dTmax)]$product)]
Al[,fam:=fifelse(grepl("LiDAR",product),"LiDAR",fifelse(grepl("FORMS",product),"Fusion",fifelse(grepl("Corrected",product),"Fusion","Sentinel-2")))]
pA<-ggplot(Al,aes(r,product,colour=fam))+
  geom_vline(xintercept=0,linetype="dotted",colour="grey50")+
  geom_segment(aes(x=0,xend=r,yend=product),linewidth=0.9)+geom_point(size=3)+
  facet_wrap(~metric)+scale_colour_manual(values=FAM,guide="none")+
  scale_x_continuous(limits=c(-1,0.45))+
  labs(x="Correlation with observed buffering\n(r < 0 = correct: more leaf → cooler)",y=NULL,
       title="LiDAR LAI best predicts observed buffering",
       subtitle="Optical Sentinel-2 LAI is weak and inverts in dense canopy (r > 0)")+
  theme_article(12)
ggsave_article(file.path(FIG,"FigA_laispace_prediction"),pA,7.4,4.4)

## ---- Fig B : MuSICA HOBO validation per scenario (ΔTmax) ----------------------
B<-fread(file.path(TAB,"figB_scenario_metrics.csv")); B[,scenario:=factor(scenario,levels=B[order(dT_bias)]$scenario)]
Bl<-melt(B,id.vars=c("scenario","family"),measure.vars=c("dT_bias","dT_RMSE","dT_R2"),variable.name="metric",value.name="v")
Bl[,metric:=factor(metric,levels=c("dT_bias","dT_RMSE","dT_R2"),labels=c("Warm bias (°C)  ↓ better","RMSE (°C)  ↓ better","R²  (ranking, not amplitude)"))]
pB<-ggplot(Bl,aes(v,scenario,fill=family))+geom_col(width=0.72)+
  facet_wrap(~metric,scales="free_x")+scale_fill_manual(values=FAM,name=NULL)+
  labs(x=NULL,y=NULL,title="LiDAR full has the lowest warm bias vs 53 HOBO loggers",
       subtitle="R² measures ranking; the pooled S2 lead is an open-canopy effect, not superior skill (Fig D, J)")+
  theme_article(11)+legend_corner(0.99,0.02)
ggsave_article(file.path(FIG,"FigB_hobo_validation"),pB,8.6,4.6)

## ---- Fig C : magnitude axis --------------------------------------------------
C<-fread(file.path(TAB,"figC_magnitude_axis.csv"))
C[,fam:=fifelse(grepl("trunc",scenario),"d_opt truncated",fifelse(grepl("LiDAR",scenario),"LiDAR","Sentinel-2"))]
pC<-ggplot(C,aes(meanLAI,dT_bias,colour=fam))+
  geom_smooth(aes(group=1),method="lm",se=FALSE,colour="grey60",linewidth=0.6,linetype="dashed")+
  geom_point(size=4)+ggrepel::geom_text_repel(aes(label=scenario),size=3.3,seg.color="grey70",box.padding=0.6)+
  scale_colour_manual(values=FAM,guide="none")+
  labs(x=expression("Mean LAI fed to MuSICA"~(m^2~m^-2)),y="Warm bias ΔTmax (°C)",
       title="Truncating the leaf column to d_opt under-buffers",
       subtitle="Feeding only the S2-sensed depth drops the LAI → +1 °C warm bias")+
  theme_article(12)
ggsave_article(file.path(FIG,"FigC_magnitude_axis"),pC,7.0,4.4)

## ---- Fig D : R2 vs amplitude (why R2 misleads) -------------------------------
pD<-ggplot(B,aes(dT_SDrec,dT_R2,colour=family))+
  geom_point(size=3.2)+ggrepel::geom_text_repel(aes(label=scenario),size=3,seg.color="grey70",max.overlaps=20)+
  scale_colour_manual(values=FAM,name=NULL)+
  labs(x="Amplitude recovered  SD(sim) / SD(obs)   (1 = realistic)",
       y="Between-plot R² (ΔTmax)",
       title="R² measures ranking, not amplitude",
       subtitle="Every scenario recovers ≤ ¼ of the observed spread (ABL binary); high R² with low amplitude is not skill")+
  theme_article(12)+legend_corner(0.02,0.99)
ggsave_article(file.path(FIG,"FigD_r2_vs_amplitude"),pD,7.2,4.6)

cat("DONE: FigA/B/C/D -> manuscripts/ch3/figures/ (png+pdf)\n")
