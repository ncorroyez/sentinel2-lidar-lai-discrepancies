# ---
# title:  05b_sm5_plot_dopt_metrics.R
# desc:   Visualisation — R (Pearson r) vs canopy integration depth, faceted
#         by h_min × Site, matching the style of Fig. 5 from the submitted
#         paper (Main_01_analysis_optim_depth.R).
#
#         Per panel:
#           - Curves coloured by Norm (DSM_keepTrees blue, DTM_keepTrees red)
#           - Gray shaded rectangle for depths > h_min (redundancy zone)
#           - Downward arrows at Pareto d_opt per norm, coloured by norm
#           - Dashed black horizontal: r_all = R at Depth=38 (total LAI ref)
#           - Dotted black horizontal: r_uniform = R at Depth=h_min (sampled
#             LAI ref — boundary of redundancy zone)
#
#         Reads:
#           revision/output/intermediate/sm5/all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#
#         Figures written to revision/output/figures/sm5/.
#
# Run from the project root (NC_Full/):
#   source("revision/scripts/05b_sm5_plot_dopt_metrics.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("terra")
library("cli")

source(here::here("revision", "R", "paths.R"))
source(here::here("revision", "R", "sm5_dopt.R"))
source(here::here("revision", "R", "sm5_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites_individual <- c("Aigoual", "Blois", "Mormal")
norms_select  <- c("DSM_keepTrees", "DTM_keepTrees")
h_min_values  <- c(10L, 15L, 20L)
h_min_10      <- 10L

norm_colours   <- c(DSM_keepTrees = "#d6604d", DTM_keepTrees = "#2166ac")
norm_labels    <- c(DSM_keepTrees = "CHM", DTM_keepTrees = "DTM")
norm_linetypes <- c(DSM_keepTrees = "solid",  DTM_keepTrees = "dashed")

# ── Sensitivity parameters ─────────────────────────────────────────────────────
# k_select = 0.5 and theta_select = 0 → use 05a's CSV (h_min 10/15/20 available)
# any other (k, theta)               → use 05c's kt CSV (h_min = 10 only)
k_ref         <- 0.5
k_select      <- 0.5   # extinction coefficient to display
theta_select  <- 0     # scan angle (degrees) to display

# "pool"    — use pooled-pixel combined row ("Sites combined")
# "average" — compute site-averaged metrics post-hoc ("Sites averaged")
combined_mode <- "average"

# ── Load and prepare data ──────────────────────────────────────────────────────

use_reference <- isTRUE(k_select == k_ref && theta_select == 0)

if (use_reference) {
  atbd_csv <- here::here(
    "revision", "output", "intermediate", "sm5",
    "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(atbd_csv))
    stop("Run 05a_sm5_compute_metrics_atbd.R first:\n  ", atbd_csv)
  dt <- data.table::fread(atbd_csv)
  dt <- dt[ATBD == TRUE & Norm %in% norms_select & h_min %in% h_min_values]
} else {
  kt_csv <- here::here(
    "revision", "output", "intermediate", "sm5",
    "kt_sensitivity_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(kt_csv))
    stop("Run 05c_sm5_k_zmin_sensitivity_dopt.R first:\n  ", kt_csv)
  dt_kt <- data.table::fread(kt_csv)
  dt <- dt_kt[k == k_select & theta == theta_select & Norm %in% norms_select]
  if (nrow(dt) == 0L)
    stop("No rows found for k=", k_select, ", theta=", theta_select,
         " in:\n  ", kt_csv)
  dt[, h_min := h_min_10]
  h_min_values <- h_min_10   # only h_min=10 from 05c
  cli::cli_alert_info(
    "Using 05c CSV — k={k_select}, theta={theta_select}°, h_min=10 only"
  )
}

# ── Apply combined_mode ────────────────────────────────────────────────────────
# Always derive the combined row consistently regardless of CSV source.

combined_row_names <- c("All_sites", "Sites_combined", "Sites_averaged",
                        "Sites combined", "Sites averaged")
dt_individual <- dt[!Site %in% combined_row_names]

if (combined_mode == "average") {
  combined_label <- "Sites averaged"
  avg_dt <- dt_individual[Site %in% sites_individual,
    .(R     = mean(R,     na.rm = TRUE),
      R2    = mean(R2,    na.rm = TRUE),
      RMSE  = mean(RMSE,  na.rm = TRUE),
      Bias  = mean(Bias,  na.rm = TRUE),
      Slope = mean(Slope, na.rm = TRUE),
      ATBD  = TRUE, Column = "LAI_atbd"),
    by = .(Norm, Depth, h_min)
  ]
  avg_dt[, Site := combined_label]
  dt <- data.table::rbindlist(list(dt_individual, avg_dt),
                               use.names = TRUE, fill = TRUE)
} else {
  combined_label <- "Sites combined"
  pool_dt <- dt[Site %in% combined_row_names][, Site := combined_label]
  dt <- data.table::rbindlist(list(dt_individual, pool_dt),
                               use.names = TRUE, fill = TRUE)
}

sites       <- c(sites_individual, combined_label)
hmin_levels <- paste0("h_min = ", h_min_values, " m")

dt[, Site := factor(Site, levels = sites)]
dt[, h_min_label := factor(paste0("h_min = ", h_min, " m"), levels = hmin_levels)]

# ── Derived columns used in multi-metric figures ───────────────────────────────

dt[, abs_Bias   := abs(Bias)]
dt[, abs_Slope1 := abs(Slope - 1)]

# Pareto distance to utopia (one value per Site × Norm × h_min × Depth)
compute_dist_utopia <- function(sub) {
  m1 <- 1 - sub$R2; m2 <- sub$RMSE
  m3 <- abs(sub$Bias); m4 <- abs(sub$Slope - 1)
  norm01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / diff(r)
  }
  data.table::data.table(
    Depth = sub$Depth,
    dist_utopia = sqrt(norm01(m1)^2 + norm01(m2)^2 + norm01(m3)^2 + norm01(m4)^2)
  )
}
dist_dt <- dt[, compute_dist_utopia(.SD), by = .(Site, Norm, h_min)]
dt      <- merge(dt, dist_dt, by = c("Site", "Norm", "Depth", "h_min"))

# ── Normalization bounds (per Site × Norm × h_min) — reused for external refs ─
extract_norm_bounds <- function(sub) {
  m1 <- 1 - sub$R2; m2 <- sub$RMSE
  m3 <- abs(sub$Bias); m4 <- abs(sub$Slope - 1)
  data.table::data.table(
    m1_min = min(m1, na.rm = TRUE), m1_max = max(m1, na.rm = TRUE),
    m2_min = min(m2, na.rm = TRUE), m2_max = max(m2, na.rm = TRUE),
    m3_min = min(m3, na.rm = TRUE), m3_max = max(m3, na.rm = TRUE),
    m4_min = min(m4, na.rm = TRUE), m4_max = max(m4, na.rm = TRUE)
  )
}
bounds_dt <- dt[, extract_norm_bounds(.SD), by = .(Site, Norm, h_min)]

project_metrics <- function(R2, RMSE, Bias, Slope, b) {
  norm01 <- function(x, mn, mx) if (mx > mn) (x - mn) / (mx - mn) else 0
  sqrt(
    norm01(1 - R2,        b$m1_min, b$m1_max)^2 +
    norm01(RMSE,          b$m2_min, b$m2_max)^2 +
    norm01(abs(Bias),     b$m3_min, b$m3_max)^2 +
    norm01(abs(Slope - 1),b$m4_min, b$m4_max)^2
  )
}

compute_ref_metrics <- function(als_v, s2_v) {
  if (length(als_v) < 5L) return(NULL)
  m     <- compute_metrics_s2_lidar(s2_v, als_v)
  slope <- stats::coef(stats::lm(als_v ~ s2_v))[[2L]]
  list(R = m$R, RMSE = m$RMSE, Bias = m$Bias, R2 = m$R2,
       abs_Bias = abs(m$Bias), Slope = slope, abs_Slope1 = abs(slope - 1))
}

# ── Pareto d_opt per (Site × Norm × h_min) ────────────────────────────────────
# Computed here so arrows can be drawn at the chosen d_opt.

cli::cli_alert_info("Computing Pareto d_opt for each h_min...")
dopt_by_hmin <- data.table::rbindlist(lapply(h_min_values, function(hm) {
  sub <- dt[h_min == hm]
  d   <- select_dopt(sub, methods = "pareto", max_depth = hm,
                     prosail_filter = "ATBD")
  d[method_dopt == "pareto"][, h_min := hm]
}))

dopt_by_hmin[, h_min_label := factor(
  paste0("h_min = ", h_min, " m"), levels = hmin_levels
)]

# Arrow data: join d_opt with metric values at that depth
arrow_dt <- merge(
  dopt_by_hmin[, .(Site, Norm, h_min, h_min_label, d_opt, R, RMSE, Bias, Slope)],
  dt[, .(Site, Norm, Depth, h_min, dist_utopia)],
  by.x = c("Site", "Norm", "d_opt", "h_min"),
  by.y = c("Site", "Norm", "Depth", "h_min"),
  all.x = TRUE
)
arrow_dt[, abs_Bias   := abs(Bias)]
arrow_dt[, abs_Slope1 := abs(Slope - 1)]
arrow_dt[, Site := factor(Site, levels = sites)]

# ── Reference horizontal lines for the R figure ───────────────────────────────
# r_all     (dashed): mean R at Depth=38 across norms — "total LAI" reference
# r_uniform (dotted): mean R at Depth=h_min across norms — "sampled LAI" ref,
#            i.e. the last valid depth before the redundancy zone.

rref_all <- dt[Depth == 38L & Norm %in% norms_select,
               .(r_all = mean(R, na.rm = TRUE)),
               by = .(Site, h_min, h_min_label)]

rref_uniform <- data.table::rbindlist(lapply(h_min_values, function(hm) {
  sub <- dt[Depth == hm & h_min == hm & Norm %in% norms_select,
            .(r_uniform = mean(R, na.rm = TRUE)),
            by = .(Site, h_min)]
  sub[, h_min_label := factor(paste0("h_min = ", hm, " m"), levels = hmin_levels)]
  sub
}))

hline_dt <- merge(rref_all, rref_uniform, by = c("Site", "h_min", "h_min_label"),
                  all = TRUE)
hline_dt[, Site := factor(Site, levels = sites)]

# ── Plot helper ───────────────────────────────────────────────────────────────

make_metric_figure <- function(plot_dt, yvar, ylab,
                                ylim = NULL, show_legend = FALSE,
                                arrow_dt = NULL,
                                hline_dt = NULL,
                                ref_line_labels = NULL,
                                base_size = 13,
                                row_facet = "h_min_label") {

  shade_df <- unique(plot_dt[, .(h_min_label, h_min_val = h_min)])

  has_random <- !is.null(hline_dt) && "r_random" %in% names(hline_dt)

  # Expand colour/linetype scales to include reference hlines when given
  if (!is.null(hline_dt)) {
    def_rl <- if (has_random) expression(r[all], r[uniform], r[random]) else expression(r[all], r[uniform])
    rl <- if (!is.null(ref_line_labels)) ref_line_labels else def_rl
    clr_values <- c(norm_colours, r_all = "grey20", r_uniform = "grey20")
    lty_values <- c(norm_linetypes, r_all = "dashed", r_uniform = "dotted")
    cl_breaks  <- c("DSM_keepTrees", "DTM_keepTrees", "r_all", "r_uniform")
    if (has_random) {
      clr_values <- c(clr_values, r_random = "grey50")
      lty_values <- c(lty_values, r_random = "dotdash")
      cl_breaks  <- c(cl_breaks, "r_random")
    }
    cl_labels  <- c(expression("CHM", "DTM"), rl)
  } else {
    clr_values <- norm_colours
    lty_values <- norm_linetypes
    cl_breaks  <- c("DSM_keepTrees", "DTM_keepTrees")
    cl_labels  <- c("CHM", "DTM")
  }

  p <- ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = Depth, y = .data[[yvar]],
                 colour = Norm, linetype = Norm)
  ) +
    ggplot2::geom_rect(
      data = shade_df,
      ggplot2::aes(xmin = h_min_val, xmax = Inf, ymin = -Inf, ymax = Inf,
                   fill = "Redundancy zone"),
      inherit.aes = FALSE, alpha = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = clr_values, breaks = cl_breaks, labels = cl_labels, name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = lty_values, breaks = cl_breaks, labels = cl_labels, name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Redundancy zone" = "grey85"), name = NULL
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 10, 20, 30), limits = c(0, NA)) +
    ggplot2::labs(x = "Depth (m)", y = ylab) +
    ggplot2::facet_grid(stats::as.formula(paste(row_facet, "~ Site"))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = if (show_legend) "bottom" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92"),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      strip.text.y     = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::guides(
      colour   = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1),
      fill     = ggplot2::guide_legend(order = 2,
                                        override.aes = list(alpha = 0.5))
    )

  if (!is.null(ylim))
    p <- p + ggplot2::coord_cartesian(ylim = ylim)

  # Reference hlines with mapped aesthetics so they appear in the legend
  if (!is.null(hline_dt)) {
    hline_parts <- list(
      hline_dt[, .(Site, h_min, h_min_label, yintercept = r_all,
                   colour_key = "r_all",     linetype_key = "r_all")],
      hline_dt[, .(Site, h_min, h_min_label, yintercept = r_uniform,
                   colour_key = "r_uniform", linetype_key = "r_uniform")]
    )
    if (has_random)
      hline_parts <- c(hline_parts, list(
        hline_dt[, .(Site, h_min, h_min_label, yintercept = r_random,
                     colour_key = "r_random", linetype_key = "r_random")]
      ))
    hline_long <- data.table::rbindlist(hline_parts)
    hline_long[, Site := factor(Site, levels = levels(plot_dt$Site))]
    p <- p + ggplot2::geom_hline(
      data      = hline_long,
      ggplot2::aes(yintercept = yintercept,
                   colour     = colour_key,
                   linetype   = linetype_key),
      linewidth = 0.5
    )
  }

  # Downward arrows at d_opt per norm — same colour as the respective curve
  if (!is.null(arrow_dt) && nrow(arrow_dt) > 0 && yvar %in% names(arrow_dt)) {
    y_range <- range(plot_dt[[yvar]], na.rm = TRUE)
    offset  <- diff(y_range) * 0.20

    for (norm_i in norms_select) {
      sub_a <- data.table::copy(arrow_dt[Norm == norm_i])
      if (nrow(sub_a) == 0L) next
      sub_a[, y_tip   := get(yvar)]
      sub_a[, y_base  := y_tip - offset]
      p <- p + ggplot2::geom_segment(
        data = sub_a,
        ggplot2::aes(x = d_opt, xend = d_opt,
                     y = y_base, yend = y_tip),
        colour      = norm_colours[[norm_i]],
        linewidth   = 1.1,
        arrow       = ggplot2::arrow(length = ggplot2::unit(0.50, "cm"),
                                     type = "closed"),
        inherit.aes = FALSE
      )
    }
  }

  p
}

# ── Output directory ──────────────────────────────────────────────────────────

out_dir <- here::here("revision", "output", "figures", "sm5")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Figure 5 — R only, with reference lines and arrows ───────────────────────

cat("Plotting R figure (all h_min)...\n")

p_R_fig <- make_metric_figure(
  plot_dt     = dt,
  yvar        = "R",
  ylab        = "Pearson r",
  show_legend = TRUE,
  arrow_dt    = arrow_dt,
  hline_dt    = hline_dt
)

out_r <- file.path(out_dir, "dopt_R_ATBD.pdf")
ggplot2::ggsave(out_r, p_R_fig, width = 28, height = 18, units = "cm")
cat("  Written:", out_r, "\n")

# ── Figure 5 — h_min = 10 m only ──────────────────────────────────────────────

cat("Plotting Pareto figure (h_min = 10 m only)...\n")

show_ref_lines <- FALSE   # set TRUE to add d_all / d_uniform hlines

dt10       <- dt[h_min == h_min_10]
arrow_dt10 <- arrow_dt[h_min == h_min_10]

# ── Raster-based reference metrics for h_min=10 figures ─────────────────────
# Three spatial samples per site: all pixels, 5000 stratified uniform, 5000 random.
# Computed only when k_select == k_ref and theta_select == 0 (rasters are at
# the reference k/theta). Skipped otherwise — hlines are omitted from figures.

ref_met     <- list()
d_all_rows  <- list()
site_pixels <- list()

if (use_reference) {

  site_meta_raster <- list(
    Aigoual = list(s2_date = "2021-07-11"),
    Blois   = list(s2_date = "2021-06-14"),
    Mormal  = list(s2_date = "2021-06-14")
  )
  n_random <- 5000L
  set.seed(42L)

  for (s in names(site_meta_raster)) {
    meta     <- site_meta_raster[[s]]
    base     <- file.path(paths$ext_results, s, "Metrics", "Deciduous_Only")
    als_path <- file.path(base, "lidarlai_res_10_m.tif")
    s2_path  <- file.path(base, paste0("s2lai_", meta$s2_date, "_atbd_res_10_m.tif"))
    if (!file.exists(als_path) || !file.exists(s2_path)) {
      cli::cli_warn("Raster missing for {s}, skipping reference metrics.")
      next
    }
    r_als <- terra::rast(als_path)
    r_s2  <- terra::rast(s2_path)
    if (!terra::compareGeom(r_als, r_s2, stopOnError = FALSE))
      r_s2 <- terra::resample(r_s2, r_als, method = "bilinear")
    als_v <- terra::values(r_als, mat = FALSE)
    s2_v  <- terra::values(r_s2,  mat = FALSE)
    valid <- !is.na(als_v) & !is.na(s2_v) & als_v >= 0 & s2_v >= 0
    als_v <- als_v[valid]; s2_v <- s2_v[valid]
    site_pixels[[s]] <- list(als = als_v, s2 = s2_v)

    rm_all <- compute_ref_metrics(als_v, s2_v)
    idx_r  <- sample(length(als_v), min(n_random, length(als_v)))
    rm_rnd <- compute_ref_metrics(als_v[idx_r], s2_v[idx_r])

    gpkg_path <- here::here(
      "03_RESULTS", s, "PROSAIL_Optimization", "sampling",
      "Sampling_stratified_uniform_hmin10_nbSamples_5000.GPKG"
    )
    rm_uni <- NULL
    if (file.exists(gpkg_path)) {
      pts   <- terra::vect(gpkg_path)
      e_als <- terra::extract(r_als, pts, ID = FALSE)[[1L]]
      e_s2  <- terra::extract(r_s2,  pts, ID = FALSE)[[1L]]
      vu    <- !is.na(e_als) & !is.na(e_s2) & e_als >= 0 & e_s2 >= 0
      rm_uni <- compute_ref_metrics(e_als[vu], e_s2[vu])
    }

    ref_met[[s]] <- list(all = rm_all, uniform = rm_uni, random = rm_rnd)

    slope  <- rm_all$Slope
    d_vals <- vapply(norms_select, function(norm) {
      b <- bounds_dt[Site == s & Norm == norm & h_min == h_min_10]
      if (nrow(b) == 0L) return(NA_real_)
      project_metrics(rm_all$R2, rm_all$RMSE, rm_all$Bias, slope, b)
    }, numeric(1L))
    d_all_rows <- c(d_all_rows, list(data.table::data.table(
      Site        = s,
      h_min       = h_min_10,
      h_min_label = factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels),
      r_all       = mean(d_vals, na.rm = TRUE)
    )))
  }

  # Combined reference: pool or average depending on combined_mode
  if (combined_mode == "pool" && length(site_pixels) > 0L) {
    comb_als  <- unlist(lapply(site_pixels, `[[`, "als"))
    comb_s2   <- unlist(lapply(site_pixels, `[[`, "s2"))
    rm_comb   <- compute_ref_metrics(comb_als, comb_s2)
    idx_rc    <- sample(length(comb_als), min(n_random, length(comb_als)))
    rm_comb_r <- compute_ref_metrics(comb_als[idx_rc], comb_s2[idx_rc])
    ref_met[[combined_label]] <- list(all = rm_comb, uniform = NULL, random = rm_comb_r)

    d_vals_comb <- vapply(norms_select, function(norm) {
      b <- bounds_dt[Site == combined_label & Norm == norm & h_min == h_min_10]
      if (nrow(b) == 0L) return(NA_real_)
      project_metrics(rm_comb$R2, rm_comb$RMSE, rm_comb$Bias, rm_comb$Slope, b)
    }, numeric(1L))
    d_all_rows <- c(d_all_rows, list(data.table::data.table(
      Site        = combined_label,
      h_min       = h_min_10,
      h_min_label = factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels),
      r_all       = mean(d_vals_comb, na.rm = TRUE)
    )))
  } else if (combined_mode == "average" && length(ref_met) > 0L) {
    # Average per-site reference metrics
    avg_metric <- function(field, subfield) {
      vals <- lapply(sites_individual, function(s) ref_met[[s]][[field]][[subfield]])
      vals <- Filter(Negate(is.null), vals)
      if (length(vals) == 0L) return(NULL)
      mean(unlist(vals), na.rm = TRUE)
    }
    rm_avg <- lapply(c("R", "RMSE", "Bias", "R2", "abs_Bias", "Slope", "abs_Slope1"),
                     function(fld) avg_metric("all", fld))
    names(rm_avg) <- c("R", "RMSE", "Bias", "R2", "abs_Bias", "Slope", "abs_Slope1")
    ref_met[[combined_label]] <- list(all = rm_avg, uniform = NULL, random = NULL)
  }

} else {
  cli::cli_alert_info(
    "Raster-based reference metrics skipped (k={k_select}, theta={theta_select}° ≠ reference)."
  )
}

d_all_site    <- data.table::rbindlist(d_all_rows)
d_uniform_site <- dt10[Depth == 38L,
                        .(r_uniform = mean(dist_utopia, na.rm = TRUE)),
                        by = .(Site, h_min, h_min_label)]

pareto_hline_dt10 <- merge(d_all_site, d_uniform_site,
                            by = c("Site", "h_min", "h_min_label"))
pareto_hline_dt10[, Site := factor(Site, levels = sites)]

# Per-metric hline table builder (for individual metric h_min=10 figures)
# r_uniform comes from dt10[Depth==38] (same source as the curves) so the line
# lands exactly at the curve endpoint averaged across norms.
build_metric_hline <- function(metric_col, hmin_lbl) {
  unif_vals <- dt10[Depth == 38L,
                    .(r_uniform = mean(get(metric_col), na.rm = TRUE)),
                    by = .(Site)]
  rows <- lapply(names(ref_met), function(s) {
    rm  <- ref_met[[s]]
    r_u <- unif_vals[Site == s, r_uniform]
    data.table::data.table(
      Site        = s,
      h_min       = h_min_10,
      h_min_label = hmin_lbl,
      r_all       = if (!is.null(rm$all))    rm$all[[metric_col]]    else NA_real_,
      r_uniform   = if (length(r_u) > 0L)    r_u                     else NA_real_,
      r_random    = if (!is.null(rm$random)) rm$random[[metric_col]] else NA_real_
    )
  })
  d <- data.table::rbindlist(rows)
  d[, Site := factor(Site, levels = sites)]
  d
}

p_pareto_10 <- make_metric_figure(
  plot_dt         = dt10,
  yvar            = "dist_utopia",
  ylab            = "Pareto distance to utopia",
  show_legend     = TRUE,
  arrow_dt        = arrow_dt10,
  hline_dt        = if (show_ref_lines) pareto_hline_dt10 else NULL,
  ref_line_labels = expression(d[all], d[uniform]),
  base_size       = 14,
  row_facet       = "."
)

out_pareto10 <- file.path(out_dir, "dopt_Pareto_ATBD_hmin10.png")
ggplot2::ggsave(out_pareto10, p_pareto_10, width = 24, height = 12, units = "cm",
                device = "png", dpi = 300)
cat("  Written:", out_pareto10, "\n")

# ── Per-metric figures — h_min = 10 m, with all/uniform/random reference lines ─

cat("Plotting per-metric figures (h_min = 10 m, with reference lines)...\n")

hmin_lbl10 <- factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels)

metrics_ref <- list(
  list(var = "R",          col = "R",          ylab = "Pearson r",
       labels = expression(r[all], r[uniform], r[random]),          file = "R"),
  list(var = "RMSE",       col = "RMSE",       ylab = "RMSE (m²/m²)",
       labels = expression(RMSE[all], RMSE[uniform], RMSE[random]), file = "RMSE"),
  list(var = "abs_Bias",   col = "abs_Bias",   ylab = "|Bias| (m²/m²)",
       labels = expression(Bias[all], Bias[uniform], Bias[random]), file = "absBias"),
  list(var = "abs_Slope1", col = "abs_Slope1", ylab = "|Slope - 1|",
       labels = expression(Slope[all], Slope[uniform], Slope[random]), file = "absSlope1")
)

for (spec in metrics_ref) {
  hl <- build_metric_hline(spec$col, hmin_lbl10)
  fig10 <- make_metric_figure(
    plot_dt         = dt10,
    yvar            = spec$var,
    ylab            = spec$ylab,
    show_legend     = TRUE,
    arrow_dt        = arrow_dt10,
    hline_dt        = hl,
    ref_line_labels = spec$labels,
    base_size       = 14,
    row_facet       = "."
  )
  out10 <- file.path(out_dir, paste0("dopt_", spec$file, "_ATBD_hmin10.png"))
  ggplot2::ggsave(out10, fig10, width = 24, height = 12, units = "cm",
                  device = "png", dpi = 300)
  cat("  Written:", out10, "\n")
}

# ── Multi-metric exploratory figures ─────────────────────────────────────────

metrics_spec <- list(
  list(var = "R",           ylab = "Pearson R",             file = "R"),
  list(var = "RMSE",        ylab = "RMSE (m²/m²)",          file = "RMSE"),
  list(var = "abs_Bias",    ylab = "|Bias| (m²/m²)",        file = "absBias"),
  list(var = "abs_Slope1",  ylab = "|Slope - 1|",           file = "absSlope1"),
  list(var = "dist_utopia", ylab = "Pareto dist. to utopia", file = "Pareto")
)

for (spec in metrics_spec) {
  cat("Plotting metric:", spec$var, "\n")
  fig <- make_metric_figure(
    plot_dt     = dt,
    yvar        = spec$var,
    ylab        = spec$ylab,
    show_legend = (spec$var == "dist_utopia"),
    arrow_dt    = arrow_dt
  )
  out_file <- file.path(out_dir, paste0("dopt_metric_", spec$file, "_ATBD.pdf"))
  ggplot2::ggsave(out_file, fig, width = 28, height = 18, units = "cm")
  cat("  Written:", out_file, "\n")
}

cat("Done.\n")
