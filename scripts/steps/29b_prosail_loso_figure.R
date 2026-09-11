# ---
# title:  29b_prosail_loso_figure.R
# desc:   Figure A.19 — dispersion of the inter-sensor RMSE over the
#         225-configuration PROSAIL ensemble at each site's Pareto d_opt, with
#         the ATBD baseline, the site-specific Pareto optimum, the best
#         LiDAR-free configuration and the leave-one-site-out transferred
#         configuration marked on the distribution.
#
#         Prerequisite: scripts/steps/29_prosail_loso_transfer.R
#
#         Writes: output/figures/Fig_A19_loso_rmse_spread.{png,pdf}
# ---

suppressMessages({
  library(here); library(data.table); library(ggplot2)
})
source(here::here("R", "paths.R"))

norm_select <- "DSM_keepTrees"; h_min_select <- 10L
lai_scenario_select <- "common"; sites <- c("Aigoual", "Blois", "Mormal")
ATBD_COL <- "LIDFa=1_lai=1_LMA=1_BROWN=1"

sm5_dir <- file.path(paths$output, "intermediate", "sm5")
tab_dir <- file.path(paths$output, "tables")
fig_dir <- file.path(paths$output, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

metrics <- fread(file.path(sm5_dir, "all_results_combined_LIDFa_lai_LMA_BROWN.csv"))
metrics <- metrics[Norm == norm_select & h_min == h_min_select &
                     lai_scenario == lai_scenario_select & Site %in% sites]
metrics[, lai_lvl := as.integer(sub(".*_lai=([0-9]+)_.*", "\\1", Column))]

dopt_ref <- fread(file.path(sm5_dir, "dopt_reference.csv"))
dopt_ref <- dopt_ref[Norm == norm_select & method_dopt == "pareto" & Site %in% sites]
d_opt_site <- setNames(as.integer(dopt_ref$d_opt), dopt_ref$Site)

pro_opt <- fread(file.path(sm5_dir, "prosail_opt.csv"))
pro_opt <- pro_opt[d_opt_source == "per_site" & Norm == norm_select]
col_pareto <- setNames(pro_opt$Column_opt, pro_opt$Site)

a16  <- fread(file.path(tab_dir, "Table_A16_partial_attribution_check.csv"))
loso <- fread(file.path(tab_dir, "Table_A17_loso_transfer.csv"))

cloud <- rbindlist(lapply(sites, function(s) {
  x <- metrics[Site == s & Depth == d_opt_site[[s]]]
  data.table(Site = s, RMSE = round(x$RMSE, 2),
             Prior = fifelse(x$lai_lvl %in% 4:5,
                             "LAI prior derived from LiDAR",
                             "LAI prior not derived from LiDAR"))
}))

marks <- rbindlist(lapply(sites, function(s) {
  d <- d_opt_site[[s]]
  get_rmse <- function(cc, dd = d)
    round(metrics[Site == s & Depth == dd & Column == cc, RMSE], 2)
  lo <- loso[Held_out == s & grepl("ALS-free priors", Mode)]
  data.table(
    Site = s,
    What = factor(c("ATBD baseline", "Site-specific Pareto optimum",
                    "Best LiDAR-free configuration",
                    "Transferred from the two other sites"),
                  levels = c("ATBD baseline", "Site-specific Pareto optimum",
                             "Best LiDAR-free configuration",
                             "Transferred from the two other sites")),
    RMSE = c(get_rmse(ATBD_COL), get_rmse(col_pareto[[s]]),
             a16[Site == s & Configuration == "Best ALS-free", RMSE],
             lo$RMSE))
}))

lab_of <- function(s) sprintf("%s~(italic(d)[opt] == %d~m)", s, d_opt_site[s])
cloud[, Site_lab := lab_of(Site)]
marks[, Site_lab := lab_of(Site)]
marks[, x := 1.15 + 0.13 * (as.integer(What) - 1L)]

set.seed(123)
p <- ggplot(cloud, aes(x = 0.85, y = RMSE)) +
  geom_jitter(aes(colour = Prior), width = 0.22, height = 0, size = 0.9,
              alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = NA, colour = "grey30",
               linewidth = 0.35) +
  geom_point(data = marks, aes(x = x, y = RMSE, shape = What, fill = What),
             size = 3.2, colour = "black", stroke = 0.5) +
  scale_x_continuous(limits = c(0.58, 1.65), expand = c(0, 0)) +
  facet_wrap(~ Site_lab, labeller = label_parsed) +
  scale_colour_manual(values = c("LAI prior not derived from LiDAR" = "#4D4D4D",
                                 "LAI prior derived from LiDAR"     = "#D55E00"),
                      name = NULL) +
  scale_shape_manual(values = c(22, 21, 24, 23), name = NULL) +
  scale_fill_manual(values = c("ATBD baseline" = "white",
                               "Site-specific Pareto optimum" = "#0072B2",
                               "Best LiDAR-free configuration" = "#009E73",
                               "Transferred from the two other sites" = "#CC79A7"),
                    name = NULL) +
  labs(x = NULL, y = expression("Inter-sensor RMSE in LAI (m"^2~"m"^-2*")")) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 2, alpha = 0.6)),
         shape  = guide_legend(order = 2), fill = guide_legend(order = 2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "pt"))

ggsave(file.path(fig_dir, "Fig_A19_loso_rmse_spread.png"), p,
       width = 8.5, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Fig_A19_loso_rmse_spread.pdf"), p,
       width = 8.5, height = 5.4)
cat("written:", file.path(fig_dir, "Fig_A19_loso_rmse_spread.png"), "\n")
