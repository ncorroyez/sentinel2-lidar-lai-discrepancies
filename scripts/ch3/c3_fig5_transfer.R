# ==============================================================================
# Chapter 3 — FIG5 cross-site LAI transfer, regenerated from Table3 so it shows
# the ORACLE COLLAPSE: true d_opt lifts held-out R², predicted (non-ALS) d_opt
# does not. Durable; reads Table3_FORMSH_operationnalisation_R2.csv.
# Run from NC_Full root:  Rscript scripts/ch3/c3_fig5_transfer.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(here) })

t3 <- fread(here("manuscripts/ch3", "tables",
                 "Table3_FORMSH_operationnalisation_R2.csv"))
m <- melt(t3, id.vars = "site", variable.name = "feat", value.name = "r2")
lab <- c(S2_seul = "S2 only",
         S2_FORMSH = "S2 + FORMS-H",
         S2_FORMSH_doptVRAI = "S2 + FORMS-H + d_opt (TRUE, ALS)",
         S2_FORMSH_doptPREDIT = "S2 + FORMS-H + d_opt (predicted, non-ALS)")
m[, feat := factor(lab[as.character(feat)], levels = rev(lab))]
m[, site := factor(site, levels = c("Aigoual", "Blois", "Mormal"))]

cols <- c("S2 only" = "#F08C7A",
          "S2 + FORMS-H" = "#1FA3A3",
          "S2 + FORMS-H + d_opt (TRUE, ALS)" = "#9B5DE5",
          "S2 + FORMS-H + d_opt (predicted, non-ALS)" = "#C9B8E8")

g <- ggplot(m, aes(r2, feat, fill = feat)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey60") +
  facet_wrap(~ site) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(title = "Cross-site LAI transfer (leave-one-site-out)",
       subtitle = "R² predicting LiDAR LAI on the held-out site. TRUE (ALS) d_opt lifts it; a satellite-predicted d_opt collapses back to baseline.",
       x = expression(R^2), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8.5))

ggsave(here("output", "figures", "Fig4_dopt_transfer.png"), g,
       width = 9.5, height = 4, dpi = 150)
ggsave(here("manuscripts/ch3", "figures", "Fig4_dopt_transfer.png"), g,
       width = 9.5, height = 4, dpi = 150)
cat("DONE\n")
