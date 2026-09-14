# ---
# title:  21_sm6_density_scatter.R
# desc:   Density scatter plots — evolution of three ALS-S2 LAI relationships
#         with horizontal heterogeneity (DSM_SD). One figure per site:
#           rows    = 3 LAI combinations (ATBD_vs_ALS, ATBD_vs_ALS_dopt,
#                     opt_vs_ALS_dopt)
#           columns = 4 heterogeneity classes (Low, Medium, High, Total)
#
#         Each panel: 2D histogram (geom_bin_2d, log10 colour scale),
#         dashed 1:1 line, OLS regression line, n count in top-left corner.
#
#         Thresholds read from DSM_sd_thresholds.csv (script 15 output).
#         Rasters reloaded from SM6a + external results directories.
#
#         Outputs (one pair per site):
#           output/figures/density_scatter_Aigoual.pdf / .png
#           output/figures/density_scatter_Blois.pdf   / .png
#           output/figures/density_scatter_Mormal.pdf  / .png
#
# Run from the project root (NC_Full/):
#   source("scripts/21_sm6_density_scatter.R")
# ---

library(here)
library(terra)
library(data.table)
library(ggplot2)
library(scales)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_load_rasters.R"))
source(here::here("R", "sm6_dataframe.R"))
source(here::here("R", "sm6_classify.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites               <- c("Aigoual", "Blois", "Mormal")
norm_ref            <- "DSM_keepTrees"
d_opt_source_select <- "all_sites"
s2_opt_fn           <- "s2lai_summer_opt_common_res_10_m.tif"

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")
ext_dir  <- paths$ext_results
out_dir  <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── d_opt per site ─────────────────────────────────────────────────────────────

prosail_opt_csv <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv"
)
if (!file.exists(prosail_opt_csv))
  stop("prosail_opt.csv not found — run step 12_sm5_select_prosail_opt.R first:\n  ",
       prosail_opt_csv, call. = FALSE)

prosail_opt  <- data.table::fread(prosail_opt_csv)
opt_ref      <- prosail_opt[Norm == norm_ref & Site %in% sites &
                             d_opt_source == d_opt_source_select]
if (nrow(opt_ref) != length(sites))
  stop("prosail_opt.csv missing rows for some sites (norm=", norm_ref,
       ", source=", d_opt_source_select, "). Got: ",
       paste(opt_ref$Site, collapse = ", "), call. = FALSE)

dopt_by_site <- setNames(opt_ref$d_opt, opt_ref$Site)[sites]
cli::cli_alert_info(
  "d_opt: {paste(names(dopt_by_site), dopt_by_site, sep='=', collapse=', ')}"
)

# ── DSM thresholds (Low / High boundary) ──────────────────────────────────────

thr_csv <- file.path(paths$output, "intermediate", "sm6",
  paste0("DSM_sd_thresholds_", d_opt_source_select, ".csv")
)
if (!file.exists(thr_csv))
  stop("Missing: ", thr_csv,
       "\nRun script 15 first (it now writes mode-tagged thresholds).")
thr_dt   <- data.table::fread(thr_csv)
low_thr  <- thr_dt$low_threshold[1L]
high_thr <- thr_dt$high_threshold[1L]
cli::cli_alert_info("DSM thresholds: Low < {low_thr} m <= Medium < {high_thr} m <= High")

# ── Load rasters ───────────────────────────────────────────────────────────────

cli::cli_h2("Loading rasters (may take ~30 s)")

lai_s2_opt_paths <- setNames(
  file.path(sm6a_dir, sites, s2_opt_fn), sites
)
lai_als_dopt_paths <- setNames(
  file.path(paths$output, "intermediate", "lai_als_dopt", sites,
             "LAI_ALS_dopt_per_site.tif"),
  sites
)
dt_full <- build_multisite_dt(
  sites              = sites,
  dopt_by_site       = dopt_by_site,
  sm6a_dir           = sm6a_dir,
  ext_dir            = ext_dir,
  lai_s2_opt_paths   = lai_s2_opt_paths,
  lai_als_dopt_paths = lai_als_dopt_paths
)
cli::cli_alert_success("Total pixels loaded: {nrow(dt_full)}")

# ── Classify pixels into Low / Medium / High ───────────────────────────────────

classify_heterogeneity(dt_full, "dsm_sd",
                       low_threshold  = low_thr,
                       high_threshold = high_thr)

# ── Build long-format scatter data ─────────────────────────────────────────────

combo_defs <- list(
  list(name = "ATBD_vs_ALS",      x_col = "lai_als",      y_col = "lai_s2_atbd"),
  list(name = "ATBD_vs_ALS_dopt", x_col = "lai_als_dopt", y_col = "lai_s2_atbd"),
  list(name = "opt_vs_ALS_dopt",  x_col = "lai_als_dopt", y_col = "lai_s2_opt")
)

het_levels   <- c("Low", "Medium", "High", "Total")
combo_levels <- c("ATBD_vs_ALS", "ATBD_vs_ALS_dopt", "opt_vs_ALS_dopt")

parts <- vector("list", length(combo_defs))

for (ci in seq_along(combo_defs)) {
  cd <- combo_defs[[ci]]

  # Classed pixels (Low / Medium / High)
  sub_cls <- dt_full[!is.na(het_class),
                      .(site, het_class,
                        x_val       = get(cd$x_col),
                        y_val       = get(cd$y_col),
                        combination = cd$name)]

  # Total: all pixels regardless of class
  sub_tot <- dt_full[, .(site, het_class = "Total",
                          x_val       = get(cd$x_col),
                          y_val       = get(cd$y_col),
                          combination = cd$name)]

  parts[[ci]] <- data.table::rbindlist(list(sub_cls, sub_tot))
}

scatter_dt <- data.table::rbindlist(parts)
scatter_dt <- scatter_dt[!is.na(x_val) & !is.na(y_val)]
scatter_dt[, het_class   := factor(het_class,   levels = het_levels)]
scatter_dt[, combination := factor(combination, levels = combo_levels)]

# ── Per-panel n count ──────────────────────────────────────────────────────────

n_dt <- scatter_dt[, .(n = .N), by = .(site, combination, het_class)]
n_dt[, lbl := paste0("n = ", formatC(n, format = "d", big.mark = ","))]

# ── Combination strip labels (plotmath) ────────────────────────────────────────

combo_strip_raw <- c(
  ATBD_vs_ALS      = 'LAI[S2_ATBD] ~ "vs" ~ LAI[ALS]',
  ATBD_vs_ALS_dopt = 'LAI[S2_ATBD] ~ "vs" ~ LAI["ALS_dopt"]',
  opt_vs_ALS_dopt  = 'LAI[S2_opt] ~ "vs" ~ LAI["ALS_dopt"]'
)

combo_labeller <- ggplot2::as_labeller(
  combo_strip_raw,
  default = ggplot2::label_parsed
)

# ── Per-site figures ───────────────────────────────────────────────────────────

for (s in sites) {
  cli::cli_h2("Plotting site: {s}")

  dt_s <- scatter_dt[site == s]
  n_s  <- n_dt[site == s]

  # Caption: n pixels per het_class (use ATBD_vs_ALS as representative count)
  cap_parts <- sapply(het_levels, function(hc) {
    n_val <- n_s[het_class == hc & combination == "ATBD_vs_ALS", n]
    paste0(hc, ": n = ", formatC(n_val, format = "d", big.mark = ","))
  })
  caption_text <- paste(cap_parts, collapse = "     ")

  # Common axis range (all combinations pooled) so 1:1 line is interpretable
  xy_lo <- floor(min(dt_s$x_val, dt_s$y_val, na.rm = TRUE))
  xy_hi <- ceiling(max(dt_s$x_val, dt_s$y_val, na.rm = TRUE))

  p <- ggplot2::ggplot(dt_s, ggplot2::aes(x = x_val, y = y_val)) +
    ggplot2::geom_bin_2d(bins = 60) +
    ggplot2::scale_fill_viridis_c(
      name   = "pixels",
      trans  = "log10",
      labels = scales::label_comma(),
      option = "D"
    ) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "red", linewidth = 0.7
    ) +
    ggplot2::geom_smooth(
      method = "lm", formula = y ~ x,
      se = FALSE, colour = "dodgerblue", linewidth = 0.75, na.rm = TRUE
    ) +
    ggplot2::facet_grid(
      combination ~ het_class,
      labeller = ggplot2::labeller(
        combination = combo_labeller,
        het_class   = ggplot2::label_value
      )
    ) +
    ggplot2::scale_x_continuous(
      name   = expression(LAI[ALS] ~ (m^2 ~ m^{-2})),
      limits = c(xy_lo, xy_hi),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::scale_y_continuous(
      name   = expression(LAI[S2] ~ (m^2 ~ m^{-2})),
      limits = c(xy_lo, xy_hi),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = NULL, caption = caption_text) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background  = ggplot2::element_rect(fill = "grey92"),
      strip.text.x      = ggplot2::element_text(face = "bold", size = 9),
      strip.text.y      = ggplot2::element_text(face = "bold", size = 8),
      legend.position   = "right",
      panel.grid.minor  = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(
        size = 10, face = "bold", hjust = 0.5
      ),
      plot.caption      = ggplot2::element_text(
        hjust = 0.5, size = 10, colour = "grey20"
      )
    )

  out_base <- file.path(out_dir, paste0("density_scatter_", s))
  ggplot2::ggsave(paste0(out_base, ".pdf"), p,
                  width = 30, height = 22, units = "cm", device = cairo_pdf)
  ggplot2::ggsave(paste0(out_base, ".png"), p,
                  width = 30, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: density_scatter_{s}.pdf / .png")
}

cli::cli_h1("Done — SM6 density scatter plots (3 sites)")
