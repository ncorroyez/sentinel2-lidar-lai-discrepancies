# ---
# title:  20d_sm6_heter_continuous.R
# desc:   Alternative to 20_sm6_plot_heterogeneity.R — continuous view:
#         goodness-of-fit metrics (r, RMSE, Bias, Slope) plotted against
#         DSM_sd as a continuous covariate (no Low/Medium/High classes).
#
#         Metrics are computed inside overlapping sliding windows over the
#         DSM_sd-ordered SM5 sample points (5,000 per site, h_min=10,
#         LAI-uniform, fCover > 0.9).
#
#         Three LAI comparison combinations:
#           ATBD_vs_ALS       LAI[S2_ATBD] vs LAI[ALS]_full
#           ATBD_vs_ALS_dopt  LAI[S2_ATBD] vs LAI[ALS_dopt]
#           opt_vs_ALS_dopt   LAI[S2_opt]  vs LAI[ALS_dopt]
#
#         Reads (per site):
#           output/intermediate/sampling/{site}/
#             Sampling_stratified_uniform_hmin10_nbSamples_5000.GPKG
#             PAD_DSM_Depth_{d_opt}_Samples_*hmin10*.csv
#             PAD_DSM_Depth_38_Samples_*hmin10*.csv   # full LAI_ALS
#           output/intermediate/PROSAIL_Models/{site}/.../atbd/
#             LAI_estimated_atbd_stratified_uniform_hmin10_nbSamples_5000.csv
#           output/intermediate/PROSAIL_Models/{site}/.../common/
#             LAI_estimated_common_stratified_uniform_hmin10_nbSamples_5000.csv
#           output/intermediate/sm6/{site}/dsm_sd_*.tif
#           output/intermediate/sm5/prosail_opt.csv   (Column_opt + d_opt)
#
#         Figures written to output/figures/:
#           het_metrics_continuous_DSM_per_site.{pdf,png}
#           het_metrics_continuous_DSM_all_sites.{pdf,png}
#
# Run from project root (NC_Full/):
#   source("scripts/steps/20d_sm6_heter_continuous.R")
# ---

library(here)
library(terra)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(cli)

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_sample_loader.R"))
source(here::here("R", "sm6_heterogeneity.R"))   # window_metrics()

# ── Parameters ─────────────────────────────────────────────────────────────────

sites              <- c("Aigoual", "Blois", "Mormal")
norm_ref           <- "DSM_keepTrees"
dopt_modes_to_plot <- c("per_site", "all_sites")
dsm_sd_filename    <- "dsm_sd_res_10_m.tif"

# k rescaling (PADs precomputed at k_ref, study uses k_select)
k_ref              <- 0.5
k_select           <- 0.65
k_scale            <- k_ref / k_select   # ≈ 0.833

# Sliding window on DSM_sd-ordered sample points
window_n           <- 500L   # number of points per window (main figure)
window_step        <- 50L    # step between window centres (in rank)
min_n_per_window   <- 200L   # window dropped if fewer valid pairs

# Sensitivity check on window width (SI figure).
# The main figure uses window_n above; the SI overlay shows that the
# qualitative shape is invariant to the window-width choice in the range
# tested below.
window_n_sensitivity <- c(300L, 500L, 1000L)

# Show "Total" pooled metric stars (one per site x combination), drawn to
# the right of each panel at x = max(dsm_sd_center) + offset.
show_total_stars <- TRUE

# Sampling file naming
h_min              <- 10L
nb_samples         <- 5000L
fn_suffix          <- sprintf("stratified_uniform_hmin%d_nbSamples_%d",
                              h_min, nb_samples)

# Per-site sample loader: build_sm5_site_dt() from R/sm5_sample_loader.R

# Sliding-window helper: window_metrics() from R/sm6_heterogeneity.R

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

# ── Load prosail_opt.csv to get Column_opt + d_opt per site/mode ──────────────

opt_csv <- file.path(paths$output, "intermediate", "sm5", "prosail_opt.csv")
if (!file.exists(opt_csv))
  stop("Missing: ", opt_csv, "\nRun step 12 first.")
opt_dt <- data.table::fread(opt_csv)

# ── Output directory ──────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── One figure per dopt mode ──────────────────────────────────────────────────

for (current_mode in dopt_modes_to_plot) {

  cli::cli_h1("Mode: {current_mode}")

  opt_ref <- opt_dt[Norm == norm_ref & Site %in% sites &
                     d_opt_source == current_mode]
  if (nrow(opt_ref) != length(sites)) {
    cli::cli_warn("prosail_opt rows missing for mode '{current_mode}' — skipping.")
    next
  }

  per_site <- vector("list", length(sites))
  for (i in seq_along(sites)) {
    s   <- sites[i]
    row <- opt_ref[Site == s]
    cli::cli_alert_info(
      "{s}: d_opt = {row$d_opt}  |  Column_opt = {row$Column_opt}"
    )
    per_site[[i]] <- build_sm5_site_dt(
      site            = s,
      d_opt_per_site  = row$d_opt,
      column_opt      = row$Column_opt,
      k_select        = k_select,
      k_ref           = k_ref,
      h_min           = h_min,
      nb_samples      = nb_samples,
      paths_obj       = paths,
      ext_results     = paths$ext_results,
      dsm_sd_filename = dsm_sd_filename
    )
  }
  dt_full <- data.table::rbindlist(per_site)

  # Apply LAI ≥ 2 filter (same as compute_metrics_by_class)
  # dt_full <- dt_full[lai_als >= 2 & !is.na(dsm_sd)]

  # ── DSM_sd histogram per site (only on first mode — DSM_sd does not depend
  # on d_opt / Column_opt, so the second mode would write an identical figure)
  if (current_mode == dopt_modes_to_plot[1L]) {
    hist_dt <- data.table::copy(dt_full)
    hist_dt[, site := factor(site, levels = sites)]

    med_dt <- hist_dt[, .(med = stats::median(dsm_sd, na.rm = TRUE),
                          n   = .N), by = site]

    p_hist <- ggplot2::ggplot(
      hist_dt, ggplot2::aes(x = dsm_sd)
    ) +
      ggplot2::geom_histogram(binwidth = 0.25, fill = "grey60",
                              colour = "white", linewidth = 0.2,
                              na.rm = TRUE) +
      ggplot2::geom_vline(
        data        = med_dt,
        ggplot2::aes(xintercept = med),
        colour      = "firebrick", linetype = "dashed", linewidth = 0.6,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data        = med_dt,
        ggplot2::aes(x = Inf, y = Inf,
                     label = sprintf("median = %.2f m\nn = %d", med, n)),
        hjust       = 1.05, vjust = 1.4, size = 3.8, colour = "grey20",
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(~ site, ncol = 3L, scales = "free_y") +
      ggplot2::labs(
        x = bquote(DSM[SD] ~ (m)),
        y = "Count"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey90"),
        strip.text       = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )

    out_h_pdf <- file.path(out_dir, "het_DSM_sd_histogram_sm5_sample.pdf")
    out_h_png <- file.path(out_dir, "het_DSM_sd_histogram_sm5_sample.png")
    ggplot2::ggsave(out_h_pdf, p_hist, width = 26, height = 10, units = "cm",
                    device = cairo_pdf)
    ggplot2::ggsave(out_h_png, p_hist, width = 26, height = 10, units = "cm",
                    device = "png", dpi = 600, bg = "white")
    cli::cli_alert_success("Written: {out_h_pdf}")
    cli::cli_alert_success("Written: {out_h_png}")
  }
  cli::cli_alert_info(
    "n per site (after LAI ≥ 2): {paste(table(dt_full$site), collapse=' / ')}"
  )

  # Sliding-window metrics: per (site, combination), for each window_n tested.
  # The main figure uses only window_n; the SI overlay uses all window_n values.
  compute_long_dt <- function(win_n) {
    rows <- list()
    for (s in sites) {
      sub_s <- dt_full[site == s]
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

  # Main figure data: window_n only
  long_dt <- compute_long_dt(window_n)
  if (nrow(long_dt) == 0L) {
    cli::cli_warn("No windowed metrics — skipping mode '{current_mode}'.")
    next
  }

  # Sensitivity data: all window_n_sensitivity values stacked
  sensi_dt <- data.table::rbindlist(
    lapply(window_n_sensitivity, compute_long_dt),
    use.names = TRUE, fill = TRUE
  )

  # Pivot to long format for facet_grid(metric × site)
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

  # ── Total (pooled) metrics per (site, combination), drawn as stars at the
  # right of each panel. Computed on the full per-site sample (no windowing),
  # matching the "Total" column of the discrete-class plot (script 20).
  if (isTRUE(show_total_stars)) {
    total_rows <- list()
    for (s in sites) {
      sub_s <- dt_full[site == s]
      for (combo in combinations) {
        x <- sub_s[[combo$lidar_col]]
        y <- sub_s[[combo$s2_col]]
        ok <- !is.na(x) & !is.na(y)
        if (sum(ok) < 10L) next
        xv <- x[ok]; yv <- y[ok]
        diff <- yv - xv
        r_val    <- stats::cor(xv, yv)
        rmse_val <- sqrt(mean(diff^2))
        bias_val <- mean(diff)
        slope_val <- stats::coef(stats::lm(yv ~ xv))[[2L]]
        total_rows[[length(total_rows) + 1L]] <- data.table::data.table(
          site = s, combination = combo$name,
          R = r_val, RMSE = rmse_val, Bias = bias_val, Slope = slope_val
        )
      }
    }
    total_dt <- data.table::rbindlist(total_rows)

    # Star x-position per site: max(dsm_sd_center for that site) + offset
    x_total_dt <- plot_dt[, .(x_total = max(dsm_sd_center, na.rm = TRUE) *
                                          1.06), by = site]
    total_long <- data.table::melt(
      total_dt,
      id.vars       = c("site", "combination"),
      measure.vars  = c("R", "RMSE", "Bias", "Slope"),
      variable.name = "metric",
      value.name    = "value"
    )
    total_long[, metric := factor(
      as.character(metric),
      levels = c("R", "RMSE", "Bias", "Slope"),
      labels = metric_parsed_levels
    )]
    total_long[, site        := factor(site,        levels = sites)]
    total_long[, combination := factor(combination, levels = names(combo_colours))]
    total_long <- merge(total_long, x_total_dt, by = "site")
  } else {
    total_long <- NULL
  }

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
    {if (isTRUE(show_total_stars) && !is.null(total_long) && nrow(total_long) > 0L)
      ggplot2::geom_point(
        data = total_long,
        ggplot2::aes(x = x_total, y = value, colour = combination),
        shape       = 8,
        size        = 3.0,
        stroke      = 1.3,
        inherit.aes = FALSE,
        na.rm       = TRUE
      )
    } +
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
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  out_pdf <- file.path(out_dir,
    paste0("het_metrics_continuous_DSM_", current_mode, ".pdf"))
  out_png <- file.path(out_dir,
    paste0("het_metrics_continuous_DSM_", current_mode, ".png"))
  
  print(out_png)

  ggplot2::ggsave(out_pdf, p, width = 26, height = 22, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_png, p, width = 26, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_pdf}")
  cli::cli_alert_success("Written: {out_png}")

  # ── SI figure: sensitivity to window_n ─────────────────────────────────────
  # Overlay of curves at each window_n; linetype encodes window width.
  # Shows that the qualitative shape of metric = f(DSM_sd) does not depend on
  # the choice of window width in the range tested.

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
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.box       = "vertical",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  out_sensi_pdf <- file.path(out_dir,
    paste0("het_metrics_continuous_DSM_winsensi_", current_mode, ".pdf"))
  out_sensi_png <- file.path(out_dir,
    paste0("het_metrics_continuous_DSM_winsensi_", current_mode, ".png"))
  ggplot2::ggsave(out_sensi_pdf, p_sensi, width = 26, height = 22, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_sensi_png, p_sensi, width = 26, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_sensi_pdf}")
  cli::cli_alert_success("Written: {out_sensi_png}")

  # ── Build the long-format paired data for the hex figure below ────────────
  # (Previously also used by a scatter version; scatter dropped — hex is the
  # retained 9-panel view.)

  scatter_long <- data.table::rbindlist(lapply(combinations, function(combo) {
    data.table::data.table(
      site        = dt_full$site,
      dsm_sd      = dt_full$dsm_sd,
      lidar       = dt_full[[combo$lidar_col]],
      s2          = dt_full[[combo$s2_col]],
      combination = combo$name
    )
  }))
  scatter_long <- scatter_long[!is.na(lidar) & !is.na(s2) & !is.na(dsm_sd)]
  scatter_long[, site        := factor(site,        levels = sites)]
  scatter_long[, combination := factor(combination, levels = names(combo_colours),
                                        labels = c(
    "LAI[S2_ATBD]~vs~LAI[ALS]",
    "LAI[S2_ATBD]~vs~LAI['ALS_dopt']",
    "LAI[S2_opt]~vs~LAI['ALS_dopt']"
  ))]

  # Shared LAI axis range across all panels for visual comparison
  lai_max_axis <- max(scatter_long$lidar, scatter_long$s2,
                       na.rm = TRUE) * 1.02
  lai_min_axis <- 0

  # DSM_sd colour scale: cap at the maximum window centre reached by the
  # continuous-curve figure (max of dsm_sd_center across all sites /
  # combinations). This makes the hex palette and the continuous-curve x-axis
  # share the same upper bound, so the hex stops "showing" extreme DSM_sd
  # exactly where the continuous curves stop (e.g. ~5.5 m for Mormal Slope).
  # Values above the cap are squished to the top colour.
  dsm_cap <- max(long_dt$dsm_sd_center, na.rm = TRUE)

  # ── Per-panel stats annotations (site × combination) ──────────────────────
  # Three text lines drawn in the top-left of each panel:
  #   1) y = ax + b
  #   2) R = .., R^2 = ..
  #   3) RMSE = .., Bias = ..
  # Convention lm(s2 ~ lidar), matching compute_metrics_by_class.
  hex_stats <- scatter_long[
    , {
      ok <- !is.na(lidar) & !is.na(s2)
      x  <- lidar[ok]; y <- s2[ok]
      if (length(x) < 10L) {
        list(a = NA_real_, b = NA_real_, R = NA_real_, R2 = NA_real_,
             RMSE = NA_real_, Bias = NA_real_)
      } else {
        fit  <- stats::lm(y ~ x)
        co   <- stats::coef(fit)
        diff <- y - x
        r    <- stats::cor(x, y)
        list(a    = unname(co[[2L]]),
             b    = unname(co[[1L]]),
             R    = r,
             R2   = r^2,
             RMSE = sqrt(mean(diff^2)),
             Bias = mean(diff))
      }
    },
    by = .(site, combination)
  ]
  hex_stats[, b_str := ifelse(b >= 0,
                              sprintf(" + %.2f", b),
                              sprintf(" - %.2f", abs(b)))]
  hex_stats[, label1 := sprintf("y = %.2fx%s",            a, b_str)]
  hex_stats[, label2 := sprintf("R = %.2f,  R² = %.2f", R, R2)]
  hex_stats[, label3 := sprintf("RMSE = %.2f,  Bias = %.2f", RMSE, Bias)]

  x_ann  <- lai_min_axis + (lai_max_axis - lai_min_axis) * 0.03
  y_top  <- lai_min_axis + (lai_max_axis - lai_min_axis) * 0.96
  dy_ann <- (lai_max_axis - lai_min_axis) * 0.060

  # ── 9-panel layout, hex binning coloured by median DSM_sd ─────────────────
  # Each hexagonal cell aggregates points falling inside it; cell colour =
  # MEDIAN DSM_sd in the cell. Eliminates overplotting at the cost of point-
  # level resolution. Bins=25 → larger hexes, fewer empty cells; no min-
  # population threshold (all non-empty cells shown).

  p_hex <- ggplot2::ggplot(
    scatter_long,
    ggplot2::aes(x = lidar, y = s2, z = dsm_sd)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey15", linewidth = 0.7) +
    ggplot2::stat_summary_hex(
      bins = 25L,
      fun  = function(x) stats::median(x, na.rm = TRUE)
    ) +
    ggplot2::geom_text(
      data        = hex_stats,
      ggplot2::aes(label = label1),
      x = x_ann, y = y_top,
      hjust = 0, vjust = 1, size = 3.0, colour = "grey10",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data        = hex_stats,
      ggplot2::aes(label = label2),
      x = x_ann, y = y_top - dy_ann,
      hjust = 0, vjust = 1, size = 3.0, colour = "grey10",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data        = hex_stats,
      ggplot2::aes(label = label3),
      x = x_ann, y = y_top - 2 * dy_ann,
      hjust = 0, vjust = 1, size = 3.0, colour = "grey10",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(
      combination ~ site,
      labeller = ggplot2::labeller(
        combination = ggplot2::label_parsed,
        site        = ggplot2::label_value
      )
    ) +
    ggplot2::scale_fill_viridis_c(
      option    = "viridis",
      name      = bquote(median(DSM[SD]) ~ (m)),
      limits    = c(0, dsm_cap),
      oob       = scales::squish,
      na.value  = "transparent"
    ) +
    ggplot2::coord_equal(xlim = c(lai_min_axis, lai_max_axis),
                         ylim = c(lai_min_axis, lai_max_axis)) +
    ggplot2::labs(
      x = expression(LAI[ALS] ~ (m^2/m^2)),
      y = expression(LAI[S2]  ~ (m^2/m^2))
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position  = "right",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(barheight = grid::unit(8, "cm"))
    )

  out_hex_pdf <- file.path(out_dir,
    paste0("het_hex_LAI_S2_vs_ALS_by_DSM_sd_", current_mode, ".pdf"))
  out_hex_png <- file.path(out_dir,
    paste0("het_hex_LAI_S2_vs_ALS_by_DSM_sd_", current_mode, ".png"))
  ggplot2::ggsave(out_hex_pdf, p_hex, width = 26, height = 24, units = "cm",
                  device = cairo_pdf)
  ggplot2::ggsave(out_hex_png, p_hex, width = 26, height = 24, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cli::cli_alert_success("Written: {out_hex_pdf}")
  cli::cli_alert_success("Written: {out_hex_png}")
}

cli::cli_h1("Done -- continuous-heterogeneity figures")
