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
#           output/intermediate/sm5/all_results_atbd_LIDFa_lai_LMA_BROWN.csv
#
#         Figures written to output/figures/.
#
# Run from the project root (NC_Full/):
#   source("scripts/05b_sm5_plot_dopt_metrics.R")
# ---

library("here")
library("data.table")
library("ggplot2")
library("patchwork")
library("terra")
library("cli")

source(here::here("R", "paths.R"))
source(here::here("R", "sm5_dopt.R"))
source(here::here("R", "sm5_metrics.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites_individual <- c("Aigoual", "Blois", "Mormal")
norms_select  <- c("DSM_keepTrees", "DTM_keepTrees")
h_min_values  <- c(10L, 15L, 20L)
h_min_10      <- 10L

# Pareto selection criteria — MUST match scripts/steps/06 and 12.
# Subset of c("R","RMSE","Bias","Slope").
pareto_criteria <- c("R", "RMSE", "Bias", "Slope")
# pareto_criteria <- c("R", "RMSE", "Bias")
pareto_norm_method <- "minmax"
# pareto_norm_method <- "max"

norm_colours       <- c(DSM_keepTrees = "#d6604d", DTM_keepTrees = "#2166ac")
norm_arrow_colours <- c(DSM_keepTrees = "#8b1a0e", DTM_keepTrees = "#0d2d5e")
norm_labels        <- c(DSM_keepTrees = "CHM", DTM_keepTrees = "DTM")
norm_linetypes     <- c(DSM_keepTrees = "solid",  DTM_keepTrees = "dashed")

# ── Sensitivity parameters ─────────────────────────────────────────────────────
# k_select = 0.5 and theta_select = 0 → use 05a's CSV (h_min 10/15/20 available)
# any other (k, theta)               → use 05c's kt CSV (h_min = 10 only)
k_ref         <- 0.5
k_select      <- 0.65  # extinction coefficient to display
theta_select  <- 0     # scan angle (degrees) to display

# "pool"    — use pooled-pixel combined row ("Sites combined")
# "average" — compute site-averaged metrics post-hoc ("Sites averaged")
combined_mode <- "average"

# ── Load and prepare data ──────────────────────────────────────────────────────

use_reference <- isTRUE(k_select == k_ref && theta_select == 0)

if (use_reference) {
  atbd_csv <- here::here(
    "output", "intermediate", "sm5",
    "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(atbd_csv))
    stop("Run 05a_sm5_compute_metrics_atbd.R first:\n  ", atbd_csv)
  dt <- data.table::fread(atbd_csv)
  dt <- dt[ATBD == TRUE & Norm %in% norms_select & h_min %in% h_min_values]
} else {
  kt_csv <- here::here(
    "output", "intermediate", "sm5",
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

# Avg score (continuous): mean of normalised minimisable criteria (subset
# given by `pareto_criteria`), expressed in "higher = better" convention
# (1 = best, 0 = worst). Normalisation bounds are derived from depths <= hm
# (the valid selection zone), then applied to ALL depths so the curve is
# continuous over the full depth range. NOTE: this differs from
# compute_pareto_front() in R/sm5_dopt.R, which renormalises on Front-1
# only and is used for selection; here we need a continuous score for plotting.
compute_avg_score <- function(sub, hm) {
  raw_cols <- list(
    R     = 1 - sub$R,
    RMSE  = sub$RMSE,
    Bias  = abs(sub$Bias),
    Slope = abs(sub$Slope - 1)
  )
  use <- raw_cols[pareto_criteria]

  sel <- sub$Depth <= hm
  norm_proj <- function(x, xref) {
    hi <- max(xref, na.rm = TRUE)
    if (pareto_norm_method == "max") {
      if (hi == 0) return(rep(0, length(x)))
      return(x / hi)
    }
    lo <- min(xref, na.rm = TRUE)
    if (hi == lo) return(rep(0, length(x)))
    (x - lo) / (hi - lo)
  }
  m_norm <- vapply(use, function(v) norm_proj(v, v[sel]),
                    numeric(length(sub$Depth)))
  data.table::data.table(
    Depth     = sub$Depth,
    Avg_score = 1 - rowMeans(m_norm)
  )
}
dist_dt <- dt[, compute_avg_score(.SD, hm = h_min[1L]), by = .(Site, Norm, h_min)]
dt      <- merge(dt, dist_dt, by = c("Site", "Norm", "Depth", "h_min"))

# ── Normalization bounds (per Site × Norm × h_min) — derived from depths <= hm
extract_norm_bounds <- function(sub, hm) {
  sel <- sub$Depth <= hm
  m1 <- 1 - sub$R[sel]; m2 <- sub$RMSE[sel]
  m3 <- abs(sub$Bias[sel]); m4 <- abs(sub$Slope[sel] - 1)
  data.table::data.table(
    m1_min = min(m1, na.rm = TRUE), m1_max = max(m1, na.rm = TRUE),
    m2_min = min(m2, na.rm = TRUE), m2_max = max(m2, na.rm = TRUE),
    m3_min = min(m3, na.rm = TRUE), m3_max = max(m3, na.rm = TRUE),
    m4_min = min(m4, na.rm = TRUE), m4_max = max(m4, na.rm = TRUE)
  )
}
bounds_dt <- dt[, extract_norm_bounds(.SD, hm = h_min[1L]), by = .(Site, Norm, h_min)]

project_metrics <- function(R, RMSE, Bias, Slope, b) {
  norm01 <- function(x, mn, mx) if (mx > mn) (x - mn) / (mx - mn) else 0
  sqrt(
    norm01(1 - R,         b$m1_min, b$m1_max)^2 +
    norm01(RMSE,          b$m2_min, b$m2_max)^2 +
    norm01(abs(Bias),     b$m3_min, b$m3_max)^2 +
    norm01(abs(Slope - 1),b$m4_min, b$m4_max)^2
  )
}

compute_ref_metrics <- function(als_v, s2_v) {
  if (length(als_v) < 5L) return(NULL)
  m <- compute_metrics_s2_lidar(s2_v, als_v)   # lm(s2 ~ lidar) convention
  list(R = m$R, RMSE = m$RMSE, Bias = m$Bias, R2 = m$R2,
       abs_Bias = abs(m$Bias), Slope = m$Slope, abs_Slope1 = abs(m$Slope - 1))
}

# ── Pareto d_opt per (Site × Norm × h_min) ────────────────────────────────────
# Computed here so arrows can be drawn at the chosen d_opt.

cli::cli_alert_info("Computing Pareto d_opt for each h_min...")
dopt_by_hmin <- data.table::rbindlist(lapply(h_min_values, function(hm) {
  sub <- dt[h_min == hm]
  d   <- select_dopt(sub, methods = "pareto", max_depth = hm,
                     prosail_filter = "ATBD",
                     pareto_criteria = pareto_criteria,
                     pareto_norm_method = pareto_norm_method)
  d[method_dopt == "pareto"][, h_min := hm]
}))

dopt_by_hmin[, h_min_label := factor(
  paste0("h_min = ", h_min, " m"), levels = hmin_levels
)]

# Arrow data: join d_opt with metric values at that depth
arrow_dt <- merge(
  dopt_by_hmin[, .(Site, Norm, h_min, h_min_label, d_opt, R, RMSE, Bias, Slope)],
  dt[, .(Site, Norm, Depth, h_min, Avg_score)],
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
                                row_facet = "h_min_label",
                                fade_redundancy = FALSE) {

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

  # When fade_redundancy: valid zone is solid for all norms; redundancy zone
  # uses a fixed longdash (not in aes), so norm linetypes become redundant there.
  if (fade_redundancy) {
    for (nk in names(norm_linetypes)) {
      if (nk %in% names(lty_values)) lty_values[[nk]] <- "solid"
    }
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
    ggplot2::scale_colour_manual(
      values = clr_values, breaks = cl_breaks, labels = cl_labels,
      name = "Normalization"
    ) +
    ggplot2::scale_linetype_manual(
      values = lty_values, breaks = cl_breaks, labels = cl_labels,
      name = "Normalization"
    ) +
    { if (fade_redundancy) list(
        ggplot2::geom_line(
          data = plot_dt[Depth <= h_min],
          linewidth = 0.8, na.rm = TRUE
        ),
        ggplot2::geom_line(
          data      = plot_dt[Depth > h_min],
          linewidth = 0.5, alpha = 0.3, linetype = "longdash", na.rm = TRUE
        )
      ) else
        ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE)
    } +
    ggplot2::scale_fill_manual(
      values = c("Redundancy zone" = "grey85"), name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks       = c(0, 5, 10, 20, 30),
      minor_breaks = c(1:4, 6:9),
      limits       = c(0, NA)
    ) +
    ggplot2::labs(x = "Depth (m)", y = ylab) +
    ggplot2::facet_grid(stats::as.formula(paste(row_facet, "~ Site"))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position    = if (show_legend) "bottom" else "none",
      panel.grid.minor.x = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
      panel.grid.minor.y = ggplot2::element_blank(),
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
      sub_a[, y_base  := y_tip + offset]
      p <- p + ggplot2::geom_segment(
        data = sub_a,
        ggplot2::aes(x = d_opt, xend = d_opt,
                     y = y_base, yend = y_tip),
        colour      = norm_arrow_colours[[norm_i]],
        linewidth   = 1.1,
        arrow       = ggplot2::arrow(length = ggplot2::unit(0.30, "cm"),
                                     type = "closed"),
        inherit.aes = FALSE
      )
    }
  }

  p
}

# ── Output directory ──────────────────────────────────────────────────────────

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Figure 5 — R only, with reference lines and arrows ───────────────────────

cat("Plotting R figure (all h_min)...\n")

p_R_fig <- make_metric_figure(
  plot_dt         = dt,
  yvar            = "R",
  ylab            = expression(Pearson ~ italic(r)),
  show_legend     = TRUE,
  arrow_dt        = arrow_dt,
  hline_dt        = hline_dt,
  fade_redundancy = TRUE
)

out_r <- file.path(out_dir, "dopt_R_ATBD.pdf")
ggplot2::ggsave(out_r, p_R_fig, width = 28, height = 18, units = "cm", device = cairo_pdf)
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

# ── Raster reference metrics — always computed at ATBD k=0.5 ──────────────────
# These baselines (all pixels / stratified uniform / 5000 random) are independent
# of k_select and theta_select: the rasters represent the full-pixel LAI at the
# reference parametrisation and are valid reference lines for any k/theta curve.

ref_met     <- list()
site_pixels <- list()
n_random    <- 5000L
set.seed(42L)

site_meta_raster <- list(
  Aigoual = list(s2_date = "2021-07-11"),
  Blois   = list(s2_date = "2021-06-14"),
  Mormal  = list(s2_date = "2021-06-14")
)

for (s in names(site_meta_raster)) {
  meta     <- site_meta_raster[[s]]
  base     <- file.path(paths$ext_results, s, "Metrics", "Deciduous_Only")
  als_path <- file.path(base, "lidarlai_res_10_m.tif")
  s2_path  <- file.path(paths$output, "intermediate", "sm6",
                        s, "s2lai_summer_atbd_T_res_10_m.tif")
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

  gpkg_path <- file.path(
    paths$sampling, s,
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
}

# Combined label: average per-site reference metrics (always "average" approach)
if (length(ref_met) > 0L) {
  flds <- c("R", "RMSE", "Bias", "R2", "abs_Bias", "Slope", "abs_Slope1")
  avg_ref_field <- function(fld) {
    vals <- lapply(sites_individual, function(s) ref_met[[s]][["all"]][[fld]])
    vals <- Filter(Negate(is.null), vals)
    if (length(vals) == 0L) return(NA_real_)
    mean(unlist(vals), na.rm = TRUE)
  }
  ref_met[[combined_label]] <- list(
    all     = setNames(lapply(flds, avg_ref_field), flds),
    uniform = NULL,
    random  = NULL
  )
}

# ── Pareto-distance d_all reference — only meaningful at reference k/theta ─────
# (requires bounds_dt from the current dt, so cannot be pre-computed above)

d_all_rows <- list()

if (use_reference) {
  for (s in names(site_meta_raster)) {
    rm_all <- ref_met[[s]]$all
    if (is.null(rm_all)) next
    d_vals <- vapply(norms_select, function(norm) {
      b <- bounds_dt[Site == s & Norm == norm & h_min == h_min_10]
      if (nrow(b) == 0L) return(NA_real_)
      project_metrics(rm_all$R, rm_all$RMSE, rm_all$Bias, rm_all$Slope, b)
    }, numeric(1L))
    d_all_rows <- c(d_all_rows, list(data.table::data.table(
      Site        = s,
      h_min       = h_min_10,
      h_min_label = factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels),
      r_all       = mean(d_vals, na.rm = TRUE)
    )))
  }

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
      project_metrics(rm_comb$R, rm_comb$RMSE, rm_comb$Bias, rm_comb$Slope, b)
    }, numeric(1L))
    d_all_rows <- c(d_all_rows, list(data.table::data.table(
      Site        = combined_label,
      h_min       = h_min_10,
      h_min_label = factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels),
      r_all       = mean(d_vals_comb, na.rm = TRUE)
    )))
  }
} else {
  cli::cli_alert_info(
    "Pareto distance reference lines skipped (k={k_select}, theta={theta_select} != reference k=0.5)."
  )
}

d_all_site    <- data.table::rbindlist(d_all_rows)
d_uniform_site <- dt10[Depth == 38L,
                        .(r_uniform = mean(Avg_score, na.rm = TRUE)),
                        by = .(Site, h_min, h_min_label)]

pareto_hline_dt10 <- if (nrow(d_all_site) > 0L) {
  out <- merge(d_all_site, d_uniform_site, by = c("Site", "h_min", "h_min_label"))
  out[, Site := factor(Site, levels = sites)]
  out
} else {
  NULL
}

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
  if (nrow(d) == 0L) return(NULL)
  d[, Site := factor(Site, levels = sites)]
  d
}

p_pareto_10 <- make_metric_figure(
  plot_dt          = dt10,
  yvar             = "Avg_score",
  ylab             = "Avg score (mean of normalised criteria)",
  show_legend      = TRUE,
  arrow_dt         = arrow_dt10,
  hline_dt         = if (show_ref_lines) pareto_hline_dt10 else NULL,
  ref_line_labels  = expression(d[all], d[uniform]),
  base_size        = 14,
  row_facet        = ".",
  fade_redundancy  = TRUE
)

# NOTE: ylim removed after Avg_score convention flip (higher = better);
# the old c(0.85, 1.85) zoom no longer applies. Re-tune visually if needed.

out_pareto10 <- file.path(out_dir, "dopt_Pareto_ATBD_hmin10.png")
ggplot2::ggsave(out_pareto10, p_pareto_10, width = 24, height = 12, units = "cm",
                device = "png", dpi = 600, bg = "white")
cat("  Written:", out_pareto10, "\n")

# ── Per-metric figures — h_min = 10 m, with all/uniform/random reference lines ─

cat("Plotting per-metric figures (h_min = 10 m, with reference lines)...\n")

hmin_lbl10 <- factor(paste0("h_min = ", h_min_10, " m"), levels = hmin_levels)

metrics_ref <- list(
  list(var = "R",          col = "R",          ylab = expression(Pearson ~ italic(r)),
       labels = expression(r[all], r[uniform], r[random]),          file = "R"),
  list(var = "RMSE",       col = "RMSE",       ylab = "RMSE (m²/m²)",
       labels = expression(RMSE[all], RMSE[uniform], RMSE[random]), file = "RMSE"),
  list(var = "abs_Bias",   col = "abs_Bias",   ylab = "|Bias| (m²/m²)",
       labels = expression(Bias[all], Bias[uniform], Bias[random]), file = "absBias"),
  list(var = "abs_Slope1", col = "abs_Slope1", ylab = "|Slope - 1|",
       labels = expression(Slope[all], Slope[uniform], Slope[random]), file = "absSlope1")
)

for (spec in metrics_ref) {
  fig10 <- make_metric_figure(
    plot_dt         = dt10,
    yvar            = spec$var,
    ylab            = spec$ylab,
    show_legend     = TRUE,
    arrow_dt        = arrow_dt10,
    hline_dt        = NULL,
    ref_line_labels = NULL,
    base_size       = 14,
    row_facet       = ".",
    fade_redundancy = TRUE
  )
  out10 <- file.path(out_dir, paste0("dopt_", spec$file, "_ATBD_hmin10.png"))
  ggplot2::ggsave(out10, fig10, width = 24, height = 12, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cat("  Written:", out10, "\n")
}

# ── Multi-metric exploratory figures ─────────────────────────────────────────

metrics_spec <- list(
  list(var = "R",           ylab = "Pearson R",             file = "R"),
  list(var = "RMSE",        ylab = "RMSE (m²/m²)",          file = "RMSE"),
  list(var = "abs_Bias",    ylab = "|Bias| (m²/m²)",        file = "absBias"),
  list(var = "abs_Slope1",  ylab = "|Slope - 1|",           file = "absSlope1"),
  list(var = "Avg_score",   ylab = "Avg score (mean of normalised criteria)", file = "Pareto")
)

for (spec in metrics_spec) {
  cat("Plotting metric:", spec$var, "\n")
  fig <- make_metric_figure(
    plot_dt         = dt,
    yvar            = spec$var,
    ylab            = spec$ylab,
    show_legend     = (spec$var == "Avg_score"),
    arrow_dt        = arrow_dt,
    fade_redundancy = TRUE
  )
  out_file <- file.path(out_dir, paste0("dopt_metric_", spec$file, "_ATBD.pdf"))
  ggplot2::ggsave(out_file, fig, width = 28, height = 18, units = "cm", device = cairo_pdf)
  cat("  Written:", out_file, "\n")
}

# ── Variant: facet_grid(Site x h_min) — rows = sites (incl. averaged) ────────
# 4 rows (Aigoual, Blois, Mormal, Sites averaged) x 3 cols (h_min = 10/15/20).
# One curve per Norm (CHM solid / DTM dashed). Redundancy zone fade preserved.
# Arrows at d_opt (per Site x Norm x h_min) drawn from arrow_dt.

make_metric_figure_site_hmin <- function(plot_dt, yvar, ylab,
                                          arrow_dt    = NULL,
                                          show_legend = TRUE,
                                          base_size   = 13) {
  shade_df <- unique(plot_dt[, .(h_min_label, h_min_val = h_min)])

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
    ggplot2::geom_line(
      data      = plot_dt[Depth <= h_min],
      linewidth = 0.85, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data      = plot_dt[Depth > h_min],
      linewidth = 0.5, alpha = 0.3, linetype = "longdash", na.rm = TRUE
    ) +
    ggplot2::facet_grid(Site ~ h_min_label, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = norm_colours,
      breaks = c("DSM_keepTrees", "DTM_keepTrees"),
      labels = c("CHM", "DTM"),
      name   = "Normalization"
    ) +
    ggplot2::scale_linetype_manual(
      values = c(DSM_keepTrees = "solid", DTM_keepTrees = "dashed"),
      breaks = c("DSM_keepTrees", "DTM_keepTrees"),
      labels = c("CHM", "DTM"),
      name   = "Normalization"
    ) +
    ggplot2::scale_fill_manual(
      values = c("Redundancy zone" = "grey85"), name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks       = c(0, 5, 10, 20, 30),
      minor_breaks = c(1:4, 6:9),
      limits       = c(0, NA)
    ) +
    ggplot2::labs(x = "Depth (m)", y = ylab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position    = if (show_legend) "bottom" else "none",
      legend.box         = "horizontal",
      panel.grid.minor.x = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
      panel.grid.minor.y = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "grey92"),
      strip.text         = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::guides(
      colour   = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
      linetype = ggplot2::guide_legend(order = 1),
      fill     = ggplot2::guide_legend(order = 2, override.aes = list(alpha = 0.5))
    )

  # Optional d_opt arrows (one per Site x Norm x h_min)
  if (!is.null(arrow_dt) && nrow(arrow_dt) > 0L && yvar %in% names(arrow_dt)) {
    y_range  <- range(plot_dt[[yvar]], na.rm = TRUE)
    offset   <- diff(y_range) * 0.18
    arrow_pl <- arrow_dt[, .(Site, Norm, h_min, h_min_label, d_opt,
                              y_tip = get(yvar))]
    arrow_pl[, y_base := y_tip + offset]
    p <- p + ggplot2::geom_segment(
      data        = arrow_pl,
      ggplot2::aes(x = d_opt, xend = d_opt, y = y_base, yend = y_tip,
                   colour = Norm),
      arrow       = ggplot2::arrow(length = ggplot2::unit(0.20, "cm"),
                                    type = "closed"),
      linewidth   = 0.8,
      inherit.aes = FALSE,
      na.rm       = TRUE
    )
  }
  p
}

# For this variant we need the three h_min values (10, 15, 20). The current
# `dt` may be limited to h_min = 10 when k_select != k_ref (05c CSV). We
# therefore reload 05a's all-h_min CSV here. R and R² are invariant to the
# k rescaling; Slope/RMSE/Bias shown below are at k_ref = 0.5 (noted in the
# filename). Set use_kref_for_site_hmin <- FALSE to skip the reload and
# only plot the h_min values present in `dt`.
use_kref_for_site_hmin <- TRUE
dt_sh <- dt
arrow_dt_sh <- arrow_dt

if (use_kref_for_site_hmin) {
  atbd_csv <- here::here(
    "output", "intermediate", "sm5",
    "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (file.exists(atbd_csv)) {
    dt_kref <- data.table::fread(atbd_csv)
    dt_kref <- dt_kref[ATBD == TRUE & Norm %in% norms_select &
                         h_min %in% c(10L, 15L, 20L)]
    dt_ind  <- dt_kref[!Site %in% combined_row_names]
    avg_kref <- dt_ind[Site %in% sites_individual,
      .(R = mean(R, na.rm = TRUE), R2 = mean(R2, na.rm = TRUE),
        RMSE = mean(RMSE, na.rm = TRUE), Bias = mean(Bias, na.rm = TRUE),
        Slope = mean(Slope, na.rm = TRUE),
        ATBD = TRUE, Column = "LAI_atbd"),
      by = .(Norm, Depth, h_min)]
    avg_kref[, Site := combined_label]
    dt_sh <- data.table::rbindlist(list(dt_ind, avg_kref),
                                    use.names = TRUE, fill = TRUE)
    sites_sh <- c(sites_individual, combined_label)
    hmin_lvl <- paste0("h_min = ", c(10L, 15L, 20L), " m")
    dt_sh[, Site := factor(Site, levels = sites_sh)]
    dt_sh[, h_min_label := factor(paste0("h_min = ", h_min, " m"),
                                  levels = hmin_lvl)]
    dt_sh[, abs_Bias   := abs(Bias)]
    dt_sh[, abs_Slope1 := abs(Slope - 1)]
    # Avg score: continuous (same logic as compute_avg_score above),
    # bounds derived from depths <= h_min per (Site, Norm, h_min).
    dist_kref <- dt_sh[, compute_avg_score(.SD, hm = h_min[1L]),
                       by = .(Site, Norm, h_min)]
    dt_sh <- merge(dt_sh, dist_kref,
                   by = c("Site", "Norm", "Depth", "h_min"))
    # Arrows: recompute d_opt per (Site x Norm x h_min) at k_ref
    arrow_kref <- data.table::rbindlist(lapply(c(10L, 15L, 20L), function(hm) {
      sub <- dt_sh[h_min == hm]
      d   <- select_dopt(sub, methods = "pareto", max_depth = hm,
                         prosail_filter = "ATBD",
                         pareto_criteria    = pareto_criteria,
                         pareto_norm_method = pareto_norm_method)
      d[method_dopt == "pareto"][, h_min := hm]
    }))
    arrow_kref[, h_min_label := factor(paste0("h_min = ", h_min, " m"),
                                        levels = hmin_lvl)]
    arrow_dt_sh <- merge(
      arrow_kref[, .(Site, Norm, h_min, h_min_label, d_opt, R, RMSE,
                     Bias, Slope)],
      dt_sh[, .(Site, Norm, Depth, h_min, Avg_score)],
      by.x = c("Site", "Norm", "d_opt", "h_min"),
      by.y = c("Site", "Norm", "Depth", "h_min"),
      all.x = TRUE
    )
    arrow_dt_sh[, abs_Bias   := abs(Bias)]
    arrow_dt_sh[, abs_Slope1 := abs(Slope - 1)]
    arrow_dt_sh[, Site := factor(Site, levels = sites_sh)]
    cli::cli_alert_info(
      "Site x h_min figures: using 05a CSV (k = {k_ref}, all h_min 10/15/20)"
    )
  } else {
    cli::cli_warn("05a CSV not found — falling back to current dt (likely h_min = 10 only)")
  }
}

for (spec in metrics_spec) {
  cat("Plotting metric (Site x h_min facet):", spec$var, "\n")
  fig_sh <- make_metric_figure_site_hmin(
    plot_dt     = dt_sh,
    yvar        = spec$var,
    ylab        = spec$ylab,
    arrow_dt    = arrow_dt_sh,
    show_legend = TRUE
  )
  out_base <- file.path(out_dir, paste0("dopt_metric_", spec$file,
                                          "_site_hmin"))
  ggplot2::ggsave(paste0(out_base, ".pdf"), fig_sh,
                  width = 28, height = 22, units = "cm", device = cairo_pdf)
  ggplot2::ggsave(paste0(out_base, ".png"), fig_sh,
                  width = 28, height = 22, units = "cm",
                  device = "png", dpi = 600, bg = "white")
  cat("  Written:", out_base, ".pdf / .png\n")
}

# ── SM Appendix Figure 1: Pareto, all h_min (double facet h_min × Site) ───────
# Always uses reference data (k=0.5). If k_select ≠ k_ref, reload from CSV.

cat("Plotting SM Figure 1 — Pareto all h_min...\n")

if (use_reference) {
  dt_sm1       <- dt
  arrow_dt_sm1 <- arrow_dt
} else {
  atbd_csv_ref <- here::here(
    "output", "intermediate", "sm5",
    "all_results_atbd_LIDFa_lai_LMA_BROWN.csv"
  )
  if (!file.exists(atbd_csv_ref))
    stop("Missing reference CSV for SM Figure 1:\n  ", atbd_csv_ref)

  hm_ref <- c(10L, 15L, 20L)
  hl_ref <- paste0("h_min = ", hm_ref, " m")

  dt_sm1 <- data.table::fread(atbd_csv_ref)
  dt_sm1 <- dt_sm1[ATBD == TRUE & Norm %in% norms_select & h_min %in% hm_ref]

  comb_names <- c("All_sites", "Sites_combined", "Sites_averaged",
                  "Sites combined", "Sites averaged")
  dt_sm1_ind <- dt_sm1[!Site %in% comb_names]
  avg_sm1 <- dt_sm1_ind[Site %in% sites_individual,
    .(R=mean(R,na.rm=TRUE), R2=mean(R2,na.rm=TRUE), RMSE=mean(RMSE,na.rm=TRUE),
      Bias=mean(Bias,na.rm=TRUE), Slope=mean(Slope,na.rm=TRUE),
      ATBD=TRUE, Column="LAI_atbd"),
    by=.(Norm,Depth,h_min)]
  avg_sm1[, Site := combined_label]
  dt_sm1 <- data.table::rbindlist(list(dt_sm1_ind, avg_sm1), use.names=TRUE, fill=TRUE)
  dt_sm1[, Site := factor(Site, levels=sites)]
  dt_sm1[, h_min_label := factor(paste0("h_min = ",h_min," m"), levels=hl_ref)]
  dt_sm1[, abs_Bias   := abs(Bias)]
  dt_sm1[, abs_Slope1 := abs(Slope-1)]

  dist_sm1 <- dt_sm1[, compute_avg_score(.SD, hm = h_min[1L]), by=.(Site,Norm,h_min)]
  dt_sm1   <- merge(dt_sm1, dist_sm1, by=c("Site","Norm","Depth","h_min"))

  dopt_sm1 <- data.table::rbindlist(lapply(hm_ref, function(hm) {
    d <- select_dopt(dt_sm1[h_min==hm], methods="pareto",
                     max_depth=hm, prosail_filter="ATBD",
                     pareto_criteria = pareto_criteria,
                     pareto_norm_method = pareto_norm_method)
    d[method_dopt=="pareto"][, h_min := hm]
  }))
  dopt_sm1[, h_min_label := factor(paste0("h_min = ",h_min," m"), levels=hl_ref)]

  arrow_dt_sm1 <- merge(
    dopt_sm1[, .(Site,Norm,h_min,h_min_label,d_opt,R,RMSE,Bias,Slope)],
    dt_sm1[, .(Site,Norm,Depth,h_min,Avg_score)],
    by.x=c("Site","Norm","d_opt","h_min"), by.y=c("Site","Norm","Depth","h_min"),
    all.x=TRUE
  )
  arrow_dt_sm1[, abs_Bias   := abs(Bias)]
  arrow_dt_sm1[, abs_Slope1 := abs(Slope-1)]
  arrow_dt_sm1[, Site := factor(Site, levels=sites)]
}

p_sm1 <- make_metric_figure(
  plot_dt         = dt_sm1,
  yvar            = "Avg_score",
  ylab            = "Avg score (mean of normalised criteria)",
  show_legend     = TRUE,
  arrow_dt        = arrow_dt_sm1,
  base_size       = 13,
  row_facet       = "h_min_label",
  fade_redundancy = TRUE
)

out_sm1 <- file.path(out_dir, "sm_dopt_pareto_all_hmin.png")
ggplot2::ggsave(out_sm1, p_sm1, width = 32, height = 22, units = "cm",
                device = "png", dpi = 600, bg = "white")
cat("  Written:", out_sm1, "\n")

# ── SM Appendix Figure 2: facet_grid(metric × site), légende unifiée ──────────
#
# Layout: rows = metrics (r, RMSE, |Bias|, |Slope-1|), cols = sites.
# Shared y scale per row (free_y), shared x scale.
# Single legend: line style encodes Norm (CHM/DTM) and reference type
# (all pixels / stratified uniform / random 5 000), independent of metric.

cat("Plotting SM Figure 2 — facet_grid(metric x site), h_min=10...\n")

metrics_sm2_def <- list(
  list(var = "R",          col = "R",          label = "Pearson~italic(r)"),
  list(var = "RMSE",       col = "RMSE",       label = "RMSE~(m^2/m^2)"),
  list(var = "abs_Bias",   col = "abs_Bias",   label = "'|Bias|'~(m^2/m^2)"),
  list(var = "abs_Slope1", col = "abs_Slope1", label = "'|Slope - 1|'")
)
metric_levels_sm2 <- sapply(metrics_sm2_def, `[[`, "label")

# ── Long-format curves ─────────────────────────────────────────────────────────

dt10_long <- data.table::rbindlist(lapply(metrics_sm2_def, function(m) {
  sub <- dt10[, .(Site, Norm, Depth, h_min, Value = get(m$var))]
  sub[, metric_label := m$label]
  sub
}))
dt10_long[, metric_label := factor(metric_label, levels = metric_levels_sm2)]
dt10_long[, Site         := factor(Site,         levels = sites)]

# ── Long-format reference hlines ───────────────────────────────────────────────
# Each metric provides r_all / r_uniform / r_random values per site.
# Column names kept generic ("r_all" etc.) — labels handled by scale_*.

hline_long_sm2 <- data.table::rbindlist(lapply(metrics_sm2_def, function(m) {
  hl <- build_metric_hline(m$col, hmin_lbl10)
  if (is.null(hl) || nrow(hl) == 0L) return(NULL)
  present <- intersect(c("r_all", "r_uniform", "r_random"), names(hl))
  melted  <- data.table::melt(
    hl[, c("Site", present), with = FALSE],
    id.vars      = "Site",
    measure.vars = present,
    variable.name = "line_key",
    value.name    = "yintercept"
  )
  melted[, metric_label := m$label]
  melted
}), fill = TRUE)

if (nrow(hline_long_sm2) > 0L) {
  hline_long_sm2[, metric_label := factor(metric_label, levels = metric_levels_sm2)]
  hline_long_sm2[, Site         := factor(Site,         levels = sites)]
  hline_long_sm2[, line_key     := as.character(line_key)]
}

# ── Long-format arrows ─────────────────────────────────────────────────────────

arrow_long_sm2 <- data.table::rbindlist(lapply(metrics_sm2_def, function(m) {
  if (!m$var %in% names(arrow_dt10)) return(NULL)
  sub <- arrow_dt10[, .(Site, Norm, d_opt, Value = get(m$var))]
  sub[, metric_label := m$label]
  sub
}))
if (nrow(arrow_long_sm2) > 0L) {
  arrow_long_sm2[, metric_label := factor(metric_label, levels = metric_levels_sm2)]
  arrow_long_sm2[, Site         := factor(Site,         levels = sites)]
  y_ranges_sm2  <- dt10_long[, .(ymin = min(Value, na.rm = TRUE),
                                  ymax = max(Value, na.rm = TRUE)),
                               by = metric_label]
  arrow_long_sm2 <- merge(arrow_long_sm2, y_ranges_sm2, by = "metric_label")
  arrow_long_sm2[, y_offset := (ymax - ymin) * 0.20]
  arrow_long_sm2[, y_tip    := Value]
  arrow_long_sm2[, y_base   := y_tip + y_offset]
}

# ── Per-metric individual d_opt (drawn as dashed arrows on top of Pareto) ────
# For each metric (R: max; RMSE/|Bias|/|Slope-1|: min), the depth that would
# have been selected if THAT metric alone had been optimised (restricted to
# Depth <= h_min, so outside the redundancy zone).
indiv_dir_main <- c(
  "Pearson~italic(r)"   = "max",
  "RMSE~(m^2/m^2)"      = "min",
  "'|Bias|'~(m^2/m^2)"  = "min",
  "'|Slope - 1|'"       = "min"
)
arrow_indiv_long_main <- dt10_long[Depth <= h_min_10, {
  lbl <- as.character(unique(metric_label))
  dir <- indiv_dir_main[[lbl]]
  idx <- if (dir == "max") which.max(Value) else which.min(Value)
  .(d_opt_indiv = Depth[idx], Value = Value[idx])
}, by = .(Site, Norm, metric_label)]
if (nrow(arrow_indiv_long_main) > 0L) {
  arrow_indiv_long_main <- merge(arrow_indiv_long_main, y_ranges_sm2,
                                  by = "metric_label")
  arrow_indiv_long_main[, y_offset := (ymax - ymin) * 0.20]
  arrow_indiv_long_main[, y_tip    := Value]
  # Arrow points UP from below the data line — keeps Pareto (top, down-pointing)
  # and individual (bottom, up-pointing) visually distinct on the same panel.
  arrow_indiv_long_main[, y_base   := y_tip - y_offset]
  arrow_indiv_long_main[, metric_label :=
                           factor(metric_label, levels = metric_levels_sm2)]
  arrow_indiv_long_main[, Site := factor(Site, levels = sites)]
}

# ── Unified colour + linetype scale ───────────────────────────────────────────
# Five line types, two groups: Norm (CHM/DTM) + reference baselines.

sm2_line_keys <- c("DSM_keepTrees", "DTM_keepTrees")
sm2_colours <- c(
  DSM_keepTrees = "#d6604d", DTM_keepTrees = "#2166ac"
)
sm2_linetypes <- c(
  DSM_keepTrees = "solid",   DTM_keepTrees = "solid"
)
sm2_line_labels <- c(
  DSM_keepTrees = "CHM",
  DTM_keepTrees = "DTM"
)
sm2_linewidths <- c(
  DSM_keepTrees = 0.9, DTM_keepTrees = 0.9
)

# ── Shade rectangles (redundancy zone, per metric row) ────────────────────────

shade_sm2 <- unique(dt10_long[, .(metric_label, xmin = h_min)])

# ── Build figure ──────────────────────────────────────────────────────────────

fig_sm2 <- ggplot2::ggplot(
  dt10_long,
  ggplot2::aes(x = Depth, y = Value, colour = Norm, linetype = Norm)
) +
  ggplot2::geom_rect(
    data        = shade_sm2,
    ggplot2::aes(xmin = xmin, xmax = Inf, ymin = -Inf, ymax = Inf,
                 fill = "Redundancy zone"),
    inherit.aes = FALSE, alpha = 0.45
  ) +
  ggplot2::geom_line(
    data = dt10_long[Depth <= h_min_10],
    linewidth = 0.9, na.rm = TRUE
  ) +
  ggplot2::geom_line(
    data      = dt10_long[Depth > h_min_10],
    linewidth = 0.5, alpha = 0.3, linetype = "longdash", na.rm = TRUE
  ) +
  ggplot2::facet_grid(metric_label ~ Site, scales = "free_y",
                      switch = "y",
                      labeller = ggplot2::labeller(metric_label = ggplot2::label_parsed)) +
  ggplot2::scale_colour_manual(
    values = sm2_colours,   labels = sm2_line_labels,
    breaks = sm2_line_keys, name   = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = sm2_linetypes, labels = sm2_line_labels,
    breaks = sm2_line_keys, name   = NULL
  ) +
  ggplot2::scale_fill_manual(
    values = c("Redundancy zone" = "grey85"), name = NULL
  ) +
  ggplot2::scale_x_continuous(
    breaks       = c(0, 5, 10, 20, 30),
    minor_breaks = c(1:4, 6:9),
    limits       = c(0, NA)
  ) +
  ggplot2::labs(x = "Integration depth (m)", y = NULL) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position    = "bottom",
    legend.key.width   = ggplot2::unit(1.8, "cm"),
    panel.grid.minor.x = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background   = ggplot2::element_rect(fill = "grey92"),
    strip.text.x       = ggplot2::element_text(face = "bold", size = 11),
    strip.text.y.left  = ggplot2::element_text(face = "bold", size = 9, angle = 0),
    strip.placement    = "outside"
  ) +
  ggplot2::guides(
    colour   = ggplot2::guide_legend(order = 1, nrow = 2,
                                      override.aes = list(linewidth = 1.2)),
    linetype = ggplot2::guide_legend(order = 1, nrow = 2),
    fill     = ggplot2::guide_legend(order = 2,
                                      override.aes = list(alpha = 0.5))
  )

# Arrows at d_opt
#   Solid closed arrowhead → Pareto d_opt (joint minimisation, all 4 criteria)
#   Dashed open arrowhead  → per-metric individual d_opt (single-criterion
#                             optimum, restricted to Depth <= h_min)
if (nrow(arrow_long_sm2) > 0L) {
  for (norm_i in norms_select) {
    sub_a <- arrow_long_sm2[Norm == norm_i]
    if (nrow(sub_a) > 0L) {
      fig_sm2 <- fig_sm2 + ggplot2::geom_segment(
        data        = sub_a,
        ggplot2::aes(x = d_opt, xend = d_opt, y = y_base, yend = y_tip),
        colour      = norm_arrow_colours[[norm_i]],
        linewidth   = 1.0,
        linetype    = "solid",
        arrow       = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"),
                                     type   = "closed"),
        inherit.aes = FALSE
      )
    }

    sub_i <- arrow_indiv_long_main[Norm == norm_i]
    if (nrow(sub_i) > 0L) {
      fig_sm2 <- fig_sm2 + ggplot2::geom_segment(
        data        = sub_i,
        ggplot2::aes(x = d_opt_indiv, xend = d_opt_indiv,
                     y = y_base,      yend = y_tip),
        colour      = norm_arrow_colours[[norm_i]],
        linewidth   = 1.0,
        linetype    = "dashed",
        arrow       = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"),
                                     type   = "open"),
        inherit.aes = FALSE
      )
    }
  }
}

out_sm2     <- file.path(out_dir, "sm_dopt_4metrics_hmin10.png")
out_sm2_pdf <- file.path(out_dir, "sm_dopt_4metrics_hmin10.pdf")
ggplot2::ggsave(out_sm2,     fig_sm2, width = 30, height = 30, units = "cm",
                device = "png", dpi = 600, bg = "white")
ggplot2::ggsave(out_sm2_pdf, fig_sm2, width = 30, height = 30, units = "cm", device = cairo_pdf)
cat("  Written:", out_sm2, "\n")
cat("  Written:", out_sm2_pdf, "\n")

# ── SM Figure 2b: combined — Pareto arrows (solid) + individual optima (dashed) ─
# Solid   closed arrow → Pareto d_opt (joint minimisation of 4 criteria)
# Dashed  open   arrow → individual per-metric d_opt (single-criterion optimum,
#                        restricted to Depth ≤ h_min to stay outside redundancy zone)

cat("Plotting SM Figure 2b — combined Pareto + individual arrows...\n")

indiv_dir <- c(
  "Pearson~italic(r)"   = "max",
  "RMSE~(m^2/m^2)"      = "min",
  "'|Bias|'~(m^2/m^2)"  = "min",
  "'|Slope - 1|'"       = "min"
)

arrow_indiv_long <- dt10_long[Depth <= h_min_10, {
  lbl <- as.character(unique(metric_label))
  dir <- indiv_dir[[lbl]]
  idx <- if (dir == "max") which.max(Value) else which.min(Value)
  .(d_opt_indiv = Depth[idx], Value = Value[idx])
}, by = .(Site, Norm, metric_label)]

arrow_indiv_long <- merge(arrow_indiv_long, y_ranges_sm2, by = "metric_label")
arrow_indiv_long[, y_offset := (ymax - ymin) * 0.20]
arrow_indiv_long[, y_tip    := Value]
arrow_indiv_long[, y_base   := y_tip - y_offset]   # below → arrow points UP ↑
arrow_indiv_long[, metric_label := factor(metric_label, levels = metric_levels_sm2)]
arrow_indiv_long[, Site         := factor(Site, levels = sites)]

# Base figure — same as fig_sm2 before arrows were added
fig_sm2_comb <- ggplot2::ggplot(
  dt10_long,
  ggplot2::aes(x = Depth, y = Value, colour = Norm, linetype = Norm)
) +
  ggplot2::geom_rect(
    data        = shade_sm2,
    ggplot2::aes(xmin = xmin, xmax = Inf, ymin = -Inf, ymax = Inf,
                 fill = "Redundancy zone"),
    inherit.aes = FALSE, alpha = 0.45
  ) +
  ggplot2::geom_line(
    data = dt10_long[Depth <= h_min_10],
    linewidth = 0.9, na.rm = TRUE
  ) +
  ggplot2::geom_line(
    data      = dt10_long[Depth > h_min_10],
    linewidth = 0.5, alpha = 0.3, linetype = "longdash", na.rm = TRUE
  ) +
  ggplot2::facet_grid(metric_label ~ Site, scales = "free_y", switch = "y",
                      labeller = ggplot2::labeller(metric_label = ggplot2::label_parsed)) +
  ggplot2::scale_colour_manual(
    values = sm2_colours,   labels = sm2_line_labels,
    breaks = sm2_line_keys, name   = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = sm2_linetypes, labels = sm2_line_labels,
    breaks = sm2_line_keys, name   = NULL
  ) +
  ggplot2::scale_fill_manual(
    values = c("Redundancy zone" = "grey85"), name = NULL
  ) +
  ggplot2::scale_x_continuous(
    breaks       = c(0, 5, 10, 20, 30),
    minor_breaks = c(1:4, 6:9),
    limits       = c(0, NA)
  ) +
  ggplot2::labs(x = "Integration depth (m)", y = NULL) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    legend.position    = "bottom",
    legend.key.width   = ggplot2::unit(1.8, "cm"),
    panel.grid.minor.x = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
    panel.grid.minor.y = ggplot2::element_blank(),
    strip.background   = ggplot2::element_rect(fill = "grey92"),
    strip.text.x       = ggplot2::element_text(face = "bold", size = 11),
    strip.text.y.left  = ggplot2::element_text(face = "bold", size = 9, angle = 0),
    strip.placement    = "outside"
  ) +
  ggplot2::guides(
    colour   = ggplot2::guide_legend(order = 1, nrow = 2,
                                      override.aes = list(linewidth = 1.2)),
    linetype = ggplot2::guide_legend(order = 1, nrow = 2),
    fill     = ggplot2::guide_legend(order = 2,
                                      override.aes = list(alpha = 0.5))
  )

for (norm_i in norms_select) {
  # Pareto d_opt — solid segment, closed filled arrowhead
  sub_p <- arrow_long_sm2[Norm == norm_i]
  if (nrow(sub_p) > 0L)
    fig_sm2_comb <- fig_sm2_comb + ggplot2::geom_segment(
      data        = sub_p,
      ggplot2::aes(x = d_opt, xend = d_opt, y = y_base, yend = y_tip),
      colour      = norm_arrow_colours[[norm_i]],
      linewidth   = 1.1,
      linetype    = "solid",
      arrow       = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"),
                                    type   = "closed"),
      inherit.aes = FALSE
    )

  # Individual d_opt — dashed segment, open arrowhead
  sub_i <- arrow_indiv_long[Norm == norm_i]
  if (nrow(sub_i) > 0L)
    fig_sm2_comb <- fig_sm2_comb + ggplot2::geom_segment(
      data        = sub_i,
      ggplot2::aes(x = d_opt_indiv, xend = d_opt_indiv,
                   y = y_base,      yend = y_tip),
      colour      = norm_arrow_colours[[norm_i]],
      linewidth   = 1.1,
      linetype    = "dashed",
      arrow       = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"),
                                    type   = "open"),
      inherit.aes = FALSE
    )
}

out_sm2b     <- file.path(out_dir, "sm_dopt_4metrics_hmin10_combined.png")
out_sm2b_pdf <- file.path(out_dir, "sm_dopt_4metrics_hmin10_combined.pdf")
ggplot2::ggsave(out_sm2b,     fig_sm2_comb, width = 30, height = 30, units = "cm",
                device = "png", dpi = 600, bg = "white")
ggplot2::ggsave(out_sm2b_pdf, fig_sm2_comb, width = 30, height = 30, units = "cm", device = cairo_pdf)
cat("  Written:", out_sm2b, "\n")
cat("  Written:", out_sm2b_pdf, "\n")

cat("Done.\n")
