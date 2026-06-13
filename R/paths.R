# ---
# title:  paths.R
# desc:   Reads config.yml at the project root and exports a `paths` list
#         with absolute paths to all external data directories AND the output
#         root.  Source this file at the top of any script that reads from
#         data/ or external sources, or writes to output/.
#
#         config.yml supports two profiles (switch with `profile: local|smb`):
#           local  — data in data/, outputs in output/
#           smb    — everything on the shared INRAE server (SMB must be mounted)
#
#         Usage:
#           source(here::here("R", "paths.R"))
#           # read inputs:  paths$ext_results, paths$raw_data, …
#           # write outputs: paths$output  (root of output/)
#
#         Fallback: if config.yml is absent, defaults to local NC_Full layout
#         (data_root = project root, output_root = output/).
# ---

local({
  cfg_path <- here::here("config.yml")

  if (file.exists(cfg_path)) {
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("Package 'yaml' is required to read config.yml.\n",
           "Install with: install.packages('yaml')")

    cfg     <- yaml::read_yaml(cfg_path)
    profile <- cfg$profile %||% "local"

    if (!profile %in% names(cfg$profiles))
      stop("Profile '", profile, "' not found in config.yml.\n",
           "Available profiles: ", paste(names(cfg$profiles), collapse = ", "))

    p           <- cfg$profiles[[profile]]
    data_root   <- p$data_root
    output_root <- p$output_root
    dirs        <- p$dirs

    message("[paths.R] profile: ", profile,
            " | data_root: ", data_root,
            " | output_root: ", output_root)

  } else {
    # No config.yml — assume data lives at the project root (default layout).
    data_root   <- here::here()
    output_root <- here::here("output")
    dirs <- list(
      raw_data      = "data/raw",
      ext_results   = "data/results",
      prosail_lidar = "data/prosail_lidar",
      snap_lai      = "data/snap_lai"
    )
    message("[paths.R] no config.yml — using project root as data_root")
  }

  profile_name <- if (file.exists(cfg_path)) profile else "default"

  bootstrap_hint <- function(reason) {
    smb_hint <- if (identical(profile_name, "smb"))
      "Profile 'smb' is selected — make sure the INRAE share is mounted.\n"
    else ""
    paste0(
      reason, "\n",
      smb_hint,
      "Bootstrap (profile '", profile_name, "'):\n",
      "  1. Download data_bundle.zip from the release / Zenodo page\n",
      "     listed in README.md.\n",
      "  2. Extract it so that data/raw, data/results and ",
      "data/prosail_lidar exist.\n",
      "  3. Copy config.yml.template to config.yml and adjust ",
      "data_root / output_root."
    )
  }

  if (!dir.exists(data_root))
    stop(bootstrap_hint(paste0("data_root does not exist: ", data_root)))

  # snap_lai is optional — fall back to a default RELATIVE to data_root.
  # When data_root already points to a "data/" folder (local profile),
  # the default below resolves to data_root/snap_lai.
  snap_lai_dir <- dirs$snap_lai %||% "snap_lai"

  expected <- file.path(
    data_root,
    c(dirs$raw_data, dirs$ext_results, dirs$prosail_lidar)
  )
  if (!any(dir.exists(expected)))
    stop(bootstrap_hint(paste0(
      "data_root exists but is empty: ", data_root, "\n",
      "None of ", paste(basename(expected), collapse = " / "), " were found."
    )))

  paths <<- list(
    data_root     = data_root,
    output        = output_root,
    raw_data      = file.path(data_root, dirs$raw_data),
    ext_results   = file.path(data_root, dirs$ext_results),
    prosail_lidar = file.path(data_root, dirs$prosail_lidar),
    snap_lai      = file.path(data_root, snap_lai_dir),
    # Always local — SQLite/GPKG cannot be written to SMB shares
    sampling      = file.path(output_root, "intermediate", "sampling")
  )

  invisible(paths)
})

# Minimal null-coalescing operator (avoids rlang dependency)
`%||%` <- function(a, b) if (!is.null(a)) a else b
