# ==============================================================================
# Fig L — mean vertical PAD profile, dense vs open plots, with the effective
# optical depth d_opt = 7 m marked. Quantifies the mechanism behind the truncation
# warm bias: 63% of the dense-canopy leaf area lies below what Sentinel-2 senses.
# Data: manuscripts/ch3/tables/figL_pad_profile.csv (c3_pad_profile.R).
#   Rscript scripts/ch3/make_ch3_figL.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("scripts/_article_style.R")
TAB<-"manuscripts/ch3/tables"; FIG<-"manuscripts/ch3/figures"
P<-fread(file.path(TAB,"figL_pad_profile.csv")); S<-fread(file.path(TAB,"figL_pad_summary.csv"))
P[,stratum:=factor(stratum,levels=grep("Dense",unique(P$stratum),value=TRUE) |>
   c(grep("Open",unique(P$stratum),value=TRUE)))]
dtop<-S$dense_top; dopt<-S$dopt
pL<-ggplot(P,aes(pad,height,colour=stratum))+
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=dtop-dopt,ymax=dtop,alpha=0.10,fill="#0072B2")+
  annotate("text",x=Inf,y=dtop-dopt/2,label=sprintf("d_opt window\n(top %.0f m S2 senses)",dopt),
           hjust=1.05,vjust=0.5,size=3,colour="#0072B2")+
  geom_path(linewidth=1)+
  scale_colour_manual(values=c("#1A9850","#D7191C"),name=NULL)+
  labs(x=expression("Plant area density"~(m^2~m^-3)),y="Height above ground (m)",
       title="Most dense-canopy leaf area lies below the Sentinel-2 optical depth",
       subtitle=sprintf("%.0f%% of the leaf area in dense plots is below the d_opt window (open plots: %.0f%%)",
                        100*S$frac_below_dense,100*S$frac_below_open))+
  theme_article(12)+legend_corner(0.99,0.02)
ggsave_article(file.path(FIG,"FigL_pad_profile"),pL,7.0,5.2)
cat("DONE: FigL_pad_profile -> manuscripts/ch3/figures/\n")
