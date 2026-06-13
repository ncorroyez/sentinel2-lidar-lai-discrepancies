# ---
# title:  24_plot_s2_atbd_validation.R
# desc:   Heatmap comparing r and RMSE between three 10-m LAI products:
#           LAIS2_ATBD_T (codistribution_lai = TRUE),
#           LAIS2_ATBD_F (codistribution_lai = FALSE),
#           LAI_SNAP     (ESA SNAP Biophysical Toolbox).
#         One row per site (Aigoual, Blois, Mormal), two facets (r / RMSE).
#         Pixels masked by the artifacts + low-vegetation mask.
#
#         Prerequisite: run 25_train_apply_prosail_atbd_rasters.R first to
#         generate ATBD_T and ATBD_F rasters with prosail 3.0.0.
#
#         Reads:
#           output/intermediate/sm6/{site}/
#             s2lai_summer_atbd_T_res_10_m.tif   (ATBD_T, from script 25)
#             s2lai_summer_atbd_F_res_10_m.tif   (ATBD_F, from script 25)
#           03_RESULTS/{site}/LiDAR/Heterogeneity_Masks/
#             artifacts_low_vegetation_majority_90_p_res_10_m.envi
#           ~/Documents/Softwares/SNAP_Toolbox_Results/{site}/
#             {snap_image}_resampled_biophysical10m.data/lai.img
#
#         Outputs:
#           output/figures/s2_atbd_validation_heatmap.pdf/.png
#
# Run from the project root (NC_Full/):
#   source("scripts/24_plot_s2_atbd_validation.R")
# ---

library(here)
library(data.table)
library(ggplot2)
library(terra)
library(cli)

source(here::here("R", "paths.R"))

# ── Parameters ─────────────────────────────────────────────────────────────────

sites <- c("Aigoual", "Blois", "Mormal")

snap_ids <- c(
  Aigoual = "S2A_MSIL2A_20210711T104031_N0301_R008_T31TEJ_20210711T134956",
  Blois   = "S2A_MSIL2A_20210614T105031_N0300_R051_T31TCN_20210614T140120",
  Mormal  = "S2A_MSIL2A_20210614T105031_N0300_R051_T31UER_20210614T140120"
)

snap_root <- paths$snap_lai

out_dir <- file.path(paths$output, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Compute r and RMSE per site ────────────────────────────────────────────────

cli::cli_h1("Computing ATBD validation metrics...")

rows <- vector("list", length(sites) * 3L)
row_i <- 1L

sm6_dir <- file.path(paths$output, "intermediate", "sm6")

for (site in sites) {

  atbd_path   <- file.path(sm6_dir, site, "s2lai_summer_atbd_T_res_10_m.tif")
  codist_path <- file.path(sm6_dir, site, "s2lai_summer_atbd_F_res_10_m.tif")
  mask_path   <- file.path(paths$ext_results, site, "LiDAR", "Heterogeneity_Masks",
                            "artifacts_low_vegetation_majority_90_p_res_10_m.envi")
  snap_id   <- snap_ids[[site]]
  snap_path  <- file.path(snap_root, site,
                           paste0(snap_id, "_resampled_biophysical10m.data"),
                           "lai.img")

  for (p in c(atbd_path, codist_path, mask_path, snap_path)) {
    if (!file.exists(p)) stop("Missing file:\n  ", p)
  }

  cli::cli_alert_info("{site}: loading rasters...")

  mask_r   <- terra::rast(mask_path)
  atbd_r   <- terra::rast(atbd_path)  |> terra::project(mask_r) |> terra::mask(mask_r)
  codist_r <- terra::rast(codist_path)|> terra::project(mask_r) |> terra::mask(mask_r)
  snap_r   <- terra::rast(snap_path)  |> terra::project(mask_r) |> terra::mask(mask_r)

  v_atbd   <- terra::values(atbd_r,   na.rm = FALSE)[, 1L]
  v_codist <- terra::values(codist_r, na.rm = FALSE)[, 1L]
  v_snap   <- terra::values(snap_r,   na.rm = FALSE)[, 1L]

  good <- complete.cases(cbind(v_atbd, v_codist, v_snap))
  v_atbd   <- v_atbd[good]
  v_codist <- v_codist[good]
  v_snap   <- v_snap[good]

  cli::cli_alert_success("{site}: {sum(good)} valid pixels")

  rmse_fn <- function(x, y) sqrt(mean((x - y)^2))

  comps <- list(
    list(
      label = "LAI[S2_ATBD_T]*\" vs \"*LAI[S2_ATBD_F]",
      r     = cor(v_atbd, v_codist),
      rmse  = rmse_fn(v_atbd, v_codist)
    ),
    list(
      label = "LAI[S2_ATBD_T]*\" vs \"*LAI[SNAP]",
      r     = cor(v_atbd, v_snap),
      rmse  = rmse_fn(v_atbd, v_snap)
    ),
    list(
      label = "LAI[S2_ATBD_F]*\" vs \"*LAI[SNAP]",
      r     = cor(v_codist, v_snap),
      rmse  = rmse_fn(v_codist, v_snap)
    )
  )

  for (comp in comps) {
    rows[[row_i]] <- data.table(
      Site       = site,
      Comparison = comp$label,
      r          = comp$r,
      RMSE       = comp$rmse
    )
    row_i <- row_i + 1L
  }
}

results <- data.table::rbindlist(rows)
cli::cli_alert_success("Metrics table ({nrow(results)} rows):")
print(results)

# ── Pivot to long format ────────────────────────────────────────────────────────

results_long <- data.table::melt(
  results,
  id.vars       = c("Site", "Comparison"),
  measure.vars  = c("r", "RMSE"),
  variable.name = "Metric",
  value.name    = "Value"
)
results_long[, Site := factor(Site, levels = rev(c("Aigoual", "Blois", "Mormal")))]

metric_labeller <- ggplot2::as_labeller(
  c("r" = "italic(r)", "RMSE" = "RMSE"),
  default = ggplot2::label_parsed
)

# ── Plot ───────────────────────────────────────────────────────────────────────

p <- ggplot2::ggplot(
  results_long,
  ggplot2::aes(x = Comparison, y = Site, fill = Value)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = round(Value, 2)),
    color = "black", size = 4
  ) +
  ggplot2::facet_wrap(
    ~ Metric, nrow = 1, scales = "free",
    labeller = metric_labeller
  ) +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  ) +
  ggplot2::labs(
    x    = "Compared S2 ATBD Distributions",
    y    = "Site",
    fill = "Value"
  ) +
  ggplot2::scale_x_discrete(labels = function(x) parse(text = x))

# ── Save ───────────────────────────────────────────────────────────────────────

out_pdf <- file.path(out_dir, "s2_atbd_validation_heatmap.pdf")
out_png <- file.path(out_dir, "s2_atbd_validation_heatmap.png")

ggplot2::ggsave(out_pdf, p, width = 10, height = 10, units = "in", device = cairo_pdf)
ggplot2::ggsave(out_png, p, width = 10, height = 10, units = "in",
                device = "png", dpi = 600, bg = "white")

cli::cli_alert_success("Written: {out_pdf}")
cli::cli_alert_success("Written: {out_png}")
cli::cli_h1("Done — S2 ATBD validation heatmap")
