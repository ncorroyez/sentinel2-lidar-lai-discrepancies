# ---
# title:  20f_sm6a_heter_continuous.R
# desc:   Companion to 20d_sm6_heter_continuous.R: continuous heterogeneity
#         view (sliding window over DSM_sd) but computed on SM6a pixels
#         (rasters from output/intermediate/sm6/{site}/...) instead of the
#         SM5 stratified-uniform sample points.
#
#         For comparability with 20d, an SM6a-pixel uniform-LAI sample of
#         5,000 points per site is drawn before windowing.
#
#         Robustness check: does the qualitative continuous story hold on
#         the SM6a pixel pool (Total ≈ 50k masked pixels/site, deciduous-only
#         + low-veg-majority + artefacts masked at SM6a stage)?
#
#         Reads:
#           output/intermediate/sm6/{site}/lai_*.tif, dsm_sd_*.tif,
#             chm_sd_*.tif, max_*.tif
#           output/intermediate/sm5/prosail_opt.csv  (Column_opt + d_opt)
#
#         Writes (to output/figures/):
#           het_metrics_continuous_SM6a_DSM_per_site.{pdf,png}
#           het_metrics_continuous_SM6a_DSM_all_sites.{pdf,png}
#
# Run from project root (NC_Full/):
#   source("scripts/steps/20f_sm6a_heter_continuous.R")
# ---

library(here)
library(terra)
library(data.table)
library(ggplot2)
library(patchwork)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm6_load_rasters.R"))
source(here::here("R", "sm6_dataframe.R"))
source(here::here("R", "sm6_metrics.R"))         # sample_uniform_lai_dt()
source(here::here("R", "sm6_heterogeneity.R"))   # window_metrics()

# ── Parameters ─────────────────────────────────────────────────────────────────

sites              <- c("Aigoual", "Blois", "Mormal")
norm_ref           <- "DSM_keepTrees"
dopt_modes_to_plot <- c("per_site", "all_sites")
dsm_sd_filename    <- "dsm_sd_res_10_m.tif"
chm_sd_filename    <- "chm_sd_res_10_m.tif"
s2_opt_fn          <- "s2lai_summer_opt_res_10_m.tif"

# Sampling on the SM6a pixels.
#   uniform_lai_sampling = FALSE → use all pixels (post LAI ≥ lai_min filter)
#   uniform_lai_sampling = TRUE  → uniform-LAI sample of n_per_site per site
#                                  (comparable to SM5 n=5000/site)
uniform_lai_sampling <- FALSE
n_per_site           <- 5000L
lai_min              <- 2.0
lai_max_quantile     <- 0.98
bin_width            <- 1.0
sampling_seed        <- 42L

# Sliding window. With all pixels (~50k+ per site), the same window_n=500
# corresponds to a much narrower fraction of the DSM_sd range than with the
# 5k-point sample, so consider scaling window_n/window_step up if the curves
# look noisy at the chosen value.
window_n           <- 500L   # number of points per window (main figure)
window_step        <- 50L    # step between window centres (in rank)
min_n_per_window   <- 200L   # window dropped if fewer valid pairs

# Sensitivity check on window width (SI figure).
window_n_sensitivity <- c(300L, 500L, 1000L)

# ── Combinations ──────────────────────────────────────────────────────────────

combinations <- list(
  list(name = "ATBD_vs_ALS",      lidar_col = "lai_als",      s2_col = "lai_s2_atbd"),
  list(name = "ATBD_vs_ALS_dopt", lidar_col = "lai_als_dopt", s2_col = "lai_s2_atbd"),
  list(name = "opt_vs_ALS_dopt",  lidar_col = "lai_als_dopt", s2_col = "lai_s2_opt")
)

combo_colours <- c(
  ATBD_vs_ALS      = "#004d00",
  ATBD_vs_ALS_dopt = "#008000",
  opt_vs_ALS_dopt  = "#4CBB17"
)
combo_labels <- c(
  ATBD_vs_ALS      = expression(LAI[S2_ATBD] ~ "vs" ~ LAI[ALS]),
  ATBD_vs_ALS_dopt = expression(LAI[S2_ATBD] ~ "vs" ~ LAI["ALS_dopt"]),
  opt_vs_ALS_dopt  = expression(LAI[S2_opt]  ~ "vs" ~ LAI["ALS_dopt"])
)

metric_parsed_levels <- c("italic(r)", "RMSE~(m^2/m^2)",
                          "Bias~(m^2/m^2)", "Slope")

# ── Load prosail_opt.csv to get d_opt + Column_opt per site/mode ──────────────

opt_csv <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv")
if (!file.exists(opt_csv))
  stop("Missing: ", opt_csv, "\nRun step 12 first.")
opt_dt <- data.table::fread(opt_csv)

# ── Output directory ──────────────────────────────────────────────────────────

out_dir  <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sm6a_dir <- file.path(paths$output, "intermediate", "sm6")
ext_dir  <- paths$ext_results

# ── One figure per dopt mode ──────────────────────────────────────────────────

for (current_mode in dopt_modes_to_plot) {

  cli::cli_h1("Mode: {current_mode}")

  opt_ref <- opt_dt[Norm == norm_ref & Site %in% sites &
                     d_opt_source == current_mode]
  if (nrow(opt_ref) != length(sites)) {
    cli::cli_warn("prosail_opt rows missing for mode '{current_mode}' — skipping.")
    next
  }

  dopt_by_site <- setNames(opt_ref$d_opt, opt_ref$Site)
  cli::cli_alert_info(
    "d_opt per site: ",
    paste(names(dopt_by_site), dopt_by_site, sep = "=", collapse = ", ")
  )

  # ── Build the multi-site pixel data.table from SM6a rasters ───────────────
  lai_s2_opt_paths   <- setNames(file.path(sm6a_dir, sites, s2_opt_fn), sites)
  lai_als_dopt_paths <- setNames(
    file.path(paths$output, "intermediate", "lai_als_dopt", sites,
              paste0("LAI_ALS_dopt_", current_mode, ".tif")),
    sites
  )

  dt_full <- build_multisite_dt(
    sites              = sites,
    dopt_by_site       = dopt_by_site,
    sm6a_dir           = sm6a_dir,
    ext_dir            = ext_dir,
    lai_s2_opt_paths   = lai_s2_opt_paths,
    lai_als_dopt_paths = lai_als_dopt_paths,
    dsm_sd_filename    = dsm_sd_filename,
    chm_sd_filename    = chm_sd_filename
  )

  # ── Apply LAI ≥ lai_min and drop NA DSM_sd ─────────────────────────────────
  dt_full <- dt_full[lai_als >= lai_min & !is.na(dsm_sd)]
  cli::cli_alert_info(
    "n per site (after LAI ≥ {lai_min}): ",
    paste(table(dt_full$site), collapse = " / ")
  )

  # ── Sampling: all-pixels (default) OR uniform-LAI sample ──────────────────
  if (isTRUE(uniform_lai_sampling)) {
    dt_sampled <- sample_uniform_lai_dt(
      dt_full, n_target_per_site = n_per_site,
      lai_min = lai_min, lai_max_quantile = lai_max_quantile,
      bin_width = bin_width, seed = sampling_seed
    )
    cli::cli_alert_info(
      "n per site after uniform-LAI sampling: ",
      paste(table(dt_sampled$site), collapse = " / ")
    )
  } else {
    dt_sampled <- dt_full
    cli::cli_alert_info("Using ALL pixels (no LAI sampling).")
  }

  # ── Sliding-window metrics: per (site × combination), looped over window_n
  # for the SI sensitivity overlay. Main figure uses window_n only.
  compute_long_dt <- function(win_n) {
    rows <- list()
    for (s in sites) {
      sub_s <- dt_sampled[site == s]
      for (combo in combinations) {
        wm <- window_metrics(
          x_var            = sub_s[[combo$lidar_col]],
          y_var            = sub_s[[combo$s2_col]],
          dsm_sd           = sub_s$dsm_sd,
          window_n         = win_n,
          window_step      = window_step,
          min_n_per_window = min_n_per_window
        )
        if (nrow(wm) == 0L) next
        wm[, site        := s]
        wm[, combination := combo$name]
        wm[, window_n    := win_n]
        rows[[length(rows) + 1L]] <- wm
      }
    }
    if (length(rows) == 0L) return(data.table::data.table())
    data.table::rbindlist(rows)
  }

  long_dt <- compute_long_dt(window_n)
  if (nrow(long_dt) == 0L) {
    cli::cli_warn("No windowed metrics — skipping mode '{current_mode}'.")
    next
  }

  sensi_dt <- data.table::rbindlist(
    lapply(window_n_sensitivity, compute_long_dt),
    use.names = TRUE, fill = TRUE
  )

  # ── Pivot to long format for facet_grid(metric × site) ────────────────────
  plot_dt <- data.table::melt(
    long_dt,
    id.vars       = c("site", "combination", "dsm_sd_center", "n_window"),
    measure.vars  = c("R", "RMSE", "Bias", "Slope"),
    variable.name = "metric",
    value.name    = "value"
  )
  plot_dt[, metric := factor(
    as.character(metric),
    levels = c("R", "RMSE", "Bias", "Slope"),
    labels = metric_parsed_levels
  )]
  plot_dt[, site        := factor(site,        levels = sites)]
  plot_dt[, combination := factor(combination, levels = names(combo_colours))]

  ref_lines <- data.table::data.table(
    metric     = factor(c("Bias~(m^2/m^2)", "Slope"),
                        levels = metric_parsed_levels),
    yintercept = c(0, 1)
  )

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = dsm_sd_center, y = value, colour = combination,
                 group = combination)
  ) +
    ggplot2::geom_hline(
      data        = ref_lines,
      ggplot2::aes(yintercept = yintercept),
      linetype    = "dashed",
      colour      = "grey40",
      linewidth   = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::facet_grid(
      metric ~ site,
      scales   = "free_y",
      labeller = ggplot2::labeller(
        metric = ggplot2::label_parsed,
        site   = ggplot2::label_value
      )
    ) +
    ggplot2::scale_colour_manual(
      values = combo_colours,
      labels = combo_labels,
      name   = "Relationship",
      guide  = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 3.5, linewidth = 1.1)
      )
    ) +
    ggplot2::labs(
      x = bquote("Horizontal Heterogeneity (" * DSM[SD] * ")"),
      y = NULL,
      subtitle = if (isTRUE(uniform_lai_sampling))
        paste0("SM6a pixels -- uniform-LAI sample (n=", n_per_site,
               " per site)")
      else
        "SM6a pixels -- all pixels (LAI >= lai_min)"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30")
    )

  out_pdf <- file.path(out_dir,
    paste0("het_metrics_continuous_SM6a_DSM_", current_mode, ".pdf"))
  out_png <- file.path(out_dir,
    paste0("het_metrics_continuous_SM6a_DSM_", current_mode, ".png"))

  ggplot2::ggsave(out_pdf, p, width = 26, height = 22, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_png, p, width = 26, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_pdf}")
  cli::cli_alert_success("Written: {out_png}")

  # ── SI figure: sensitivity to window_n ────────────────────────────────────
  sensi_plot_dt <- data.table::melt(
    sensi_dt,
    id.vars       = c("site", "combination", "window_n",
                      "dsm_sd_center", "n_window"),
    measure.vars  = c("R", "RMSE", "Bias", "Slope"),
    variable.name = "metric",
    value.name    = "value"
  )
  sensi_plot_dt[, metric := factor(
    as.character(metric),
    levels = c("R", "RMSE", "Bias", "Slope"),
    labels = metric_parsed_levels
  )]
  sensi_plot_dt[, site        := factor(site,        levels = sites)]
  sensi_plot_dt[, combination := factor(combination, levels = names(combo_colours))]
  sensi_plot_dt[, window_n    := factor(window_n,    levels = window_n_sensitivity,
                                        labels = paste0("n = ", window_n_sensitivity))]

  p_sensi <- ggplot2::ggplot(
    sensi_plot_dt,
    ggplot2::aes(x = dsm_sd_center, y = value, colour = combination,
                 linetype = window_n, group = interaction(combination, window_n))
  ) +
    ggplot2::geom_hline(
      data        = ref_lines,
      ggplot2::aes(yintercept = yintercept),
      linetype    = "dashed",
      colour      = "grey40",
      linewidth   = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::facet_grid(
      metric ~ site,
      scales   = "free_y",
      labeller = ggplot2::labeller(
        metric = ggplot2::label_parsed,
        site   = ggplot2::label_value
      )
    ) +
    ggplot2::scale_colour_manual(
      values = combo_colours,
      labels = combo_labels,
      name   = "Relationship",
      guide  = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 3.5, linewidth = 1.1)
      )
    ) +
    ggplot2::scale_linetype_manual(
      values = c("dotted", "solid", "longdash"),
      name   = "Window width"
    ) +
    ggplot2::labs(
      x = bquote("Horizontal Heterogeneity (" * DSM[SD] * ")"),
      y = NULL,
      subtitle = if (isTRUE(uniform_lai_sampling))
        paste0("SM6a pixels (n=", n_per_site, "/site) -- ",
               "sensitivity to window width")
      else
        "SM6a pixels (all) -- sensitivity to window width"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "vertical",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30")
    )

  out_sensi_pdf <- file.path(out_dir,
    paste0("het_metrics_continuous_SM6a_DSM_winsensi_", current_mode, ".pdf"))
  out_sensi_png <- file.path(out_dir,
    paste0("het_metrics_continuous_SM6a_DSM_winsensi_", current_mode, ".png"))
  ggplot2::ggsave(out_sensi_pdf, p_sensi, width = 26, height = 22, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_sensi_png, p_sensi, width = 26, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_sensi_pdf}")
  cli::cli_alert_success("Written: {out_sensi_png}")
}

cli::cli_h1("Done -- SM6a-based continuous heterogeneity figures")
