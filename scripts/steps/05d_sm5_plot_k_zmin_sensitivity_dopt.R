# ---
# title:  05d_sm5_plot_k_zmin_sensitivity_dopt.R
# desc:   Figures — sensitivity of d_opt (Pareto) to:
#           (A) extinction coefficient k and scan angle theta
#               Tile heatmap: x = k, y = theta, fill = d_opt (discrete),
#               faceted by site (DSM_keepTrees only).
#           (B) minimum height z_min (2, 3, 5 m)
#               Dot-plot: x = z_min, y = d_opt, one coloured line per site.
#
#         Reads (output of 05c):
#           output/intermediate/sm5/dopt_k_theta_sensitivity.csv
#           output/intermediate/sm5/dopt_zmin_sensitivity.csv
#
#         Outputs:
#           output/figures/sm_dopt_k_theta_heatmap.pdf/.png
#           output/figures/sm_dopt_zmin.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("scripts/05d_sm5_plot_k_zmin_sensitivity_dopt.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(cli)

source(here::here("R", "paths.R"))

# ── I/O ────────────────────────────────────────────────────────────────────────

kt_csv   <- file.path(paths$output, "intermediate", "sm5",
                        "dopt_k_theta_sensitivity.csv")
zmin_csv <- file.path(paths$output, "intermediate", "sm5",
                        "dopt_zmin_sensitivity.csv")
out_dir  <- file.path(paths$output, "figures")

if (!file.exists(kt_csv))
  stop("Missing: ", kt_csv, "\nRun scripts/05c first.")
if (!file.exists(zmin_csv))
  stop("Missing: ", zmin_csv, "\nRun scripts/05c first.")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dopt_kt   <- data.table::fread(kt_csv)
dopt_zmin <- data.table::fread(zmin_csv)
cli::cli_alert_success("Loaded kt: {nrow(dopt_kt)} rows, zmin: {nrow(dopt_zmin)} rows")

# ── Parameters ─────────────────────────────────────────────────────────────────

sites        <- c("Aigoual", "Blois", "Mormal")
# Combined-site row label written by script 05c (matches combined_mode
# from script 06). Both underscore and space variants are accepted.
combined_label_candidates <- c("Sites_averaged", "Sites averaged",
                               "Sites_combined", "Sites combined")
site_colours <- c(Aigoual = "#1b7837", Blois = "#762a83", Mormal = "#d6604d",
                  `Sites averaged` = "grey25")
norm_sel     <- "DSM_keepTrees"

theme_pub <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor  = ggplot2::element_blank(),
    strip.background  = ggplot2::element_rect(fill = "grey92"),
    strip.text        = ggplot2::element_text(face = "bold", size = 11),
    legend.position   = "right"
  )

# ── Figure A: k × theta → d_opt heatmap ───────────────────────────────────────

# Include the combined-site row in the heatmap (4th panel). Normalise its
# label to "Sites averaged" for consistency across plots.
kt_present_combined <- intersect(combined_label_candidates, unique(dopt_kt$Site))
if (length(kt_present_combined) > 0L) {
  combined_in_csv <- kt_present_combined[[1L]]
  dopt_kt[Site == combined_in_csv, Site := "Sites averaged"]
}
sites_plot <- c(sites, "Sites averaged")
kt_sub <- dopt_kt[Site %in% sites_plot & Norm == norm_sel]
kt_sub[, site  := factor(Site,  levels = sites_plot)]
kt_sub[, theta := factor(theta, levels = sort(unique(theta)))]

# d_opt as factor for discrete colour scale
dopt_range <- sort(unique(kt_sub$d_opt))
kt_sub[, dopt_f := factor(d_opt, levels = dopt_range)]

# Colour scale: sequential blue-orange palette for low-to-high d_opt
n_levels  <- length(dopt_range)
tile_cols <- scales::hue_pal()(n_levels)  # replaced below with viridis-like

p_heat <- ggplot2::ggplot(
  kt_sub,
  ggplot2::aes(x = k, y = theta, fill = dopt_f)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::facet_wrap(~ site, nrow = 1) +
  ggplot2::scale_fill_viridis_d(
    name   = expression(italic(d)[opt] ~ (m)),
    option = "D",
    direction = 1
  ) +
  ggplot2::scale_x_continuous(
    name   = expression(italic(k) ~ (extinction ~ coefficient)),
    breaks = sort(unique(kt_sub$k)),
    labels = function(x) sub("^0", "", sprintf("%.2f", x))
  ) +
  ggplot2::scale_y_discrete(
    name = expression(theta ~ (scan ~ angle ~ deg))
  ) +
  theme_pub +
  ggplot2::theme(
    legend.title = ggplot2::element_text(size = 10),
    axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
  )

w_heat <- 26; h_heat <- 10

ggplot2::ggsave(
  file.path(out_dir, "sm_dopt_k_theta_heatmap.pdf"),
  p_heat, width = w_heat, height = h_heat, units = "cm", device = cairo_pdf
)
ggplot2::ggsave(
  file.path(out_dir, "sm_dopt_k_theta_heatmap.png"),
  p_heat, width = w_heat, height = h_heat, units = "cm",
  device = "png", dpi = 600, bg = "white"
)
cli::cli_alert_success("Written: sm_dopt_k_theta_heatmap.pdf / .png")

# ── Figure B: z_min → d_opt dot-line plot ──────────────────────────────────────

zmin_present_combined <- intersect(combined_label_candidates,
                                    unique(dopt_zmin$Site))
if (length(zmin_present_combined) > 0L) {
  dopt_zmin[Site == zmin_present_combined[[1L]], Site := "Sites averaged"]
}
zmin_sub <- dopt_zmin[Site %in% sites_plot & Norm == norm_sel]
zmin_sub[, site  := factor(Site, levels = sites_plot)]
zmin_sub[, z_min := as.numeric(z_min)]

p_zmin <- ggplot2::ggplot(
  zmin_sub,
  ggplot2::aes(x = z_min, y = d_opt, colour = site, group = site)
) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 3.5, shape = 16) +
  ggplot2::scale_colour_manual(
    values = site_colours,
    name   = NULL
  ) +
  ggplot2::scale_x_continuous(
    name   = expression(italic(z)[min] ~ (m)),
    breaks = c(2, 3, 5)
  ) +
  ggplot2::scale_y_continuous(
    name   = expression(italic(d)[opt] ~ (m)),
    breaks = scales::breaks_pretty(n = 6)
  ) +
  ggplot2::labs(
    caption = paste0(
      "k = 0.6 (value used in this study), theta = 0 deg. ",
      "DSM_keepTrees normalisation."
    )
  ) +
  theme_pub +
  ggplot2::theme(
    legend.position  = "bottom",
    legend.key.width = ggplot2::unit(1.2, "cm"),
    plot.caption     = ggplot2::element_text(size = 8, colour = "grey40",
                                              hjust = 0)
  )

w_zmin <- 12; h_zmin <- 10

ggplot2::ggsave(
  file.path(out_dir, "sm_dopt_zmin.pdf"),
  p_zmin, width = w_zmin, height = h_zmin, units = "cm", device = cairo_pdf
)
ggplot2::ggsave(
  file.path(out_dir, "sm_dopt_zmin.png"),
  p_zmin, width = w_zmin, height = h_zmin, units = "cm",
  device = "png", dpi = 600, bg = "white"
)
cli::cli_alert_success("Written: sm_dopt_zmin.pdf / .png")

cli::cli_h1("Done — d_opt sensitivity figures")
