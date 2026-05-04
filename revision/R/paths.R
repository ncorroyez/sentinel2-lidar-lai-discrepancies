# ---
# title:  paths.R
# desc:   Reads config.yml at the project root and exports a `paths` list
#         with absolute paths to all external data directories AND the output
#         root.  Source this file at the top of any script that reads from
#         01_DATA / 03_RESULTS / PROSAIL-Optimization, or writes to
#         revision/output.
#
#         config.yml supports two profiles (switch with `profile: local|smb`):
#           local  — data in revision/data/, outputs in revision/output/
#           smb    — everything on the shared INRAE server (SMB must be mounted)
#
#         Usage:
#           source(here::here("revision", "R", "paths.R"))
#           # read inputs:  paths$ext_results, paths$raw_data, …
#           # write outputs: paths$output  (root of revision/output/)
#
#         Fallback: if config.yml is absent, defaults to local NC_Full layout
#         (data_root = project root, output_root = revision/output/).
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
    output_root <- here::here("revision", "output")
    dirs <- list(
      raw_data      = "01_DATA",
      ext_results   = "03_RESULTS",
      prosail_lidar = "PROSAIL-Optimization/01_DATA",
      prosail_codes = "PROSAIL-Optimization/02_CODES"
    )
    message("[paths.R] no config.yml — using project root as data_root")
  }

  if (!dir.exists(data_root))
    stop(
      "data_root does not exist: ", data_root, "\n",
      "Check config.yml (or copy config.yml.template to config.yml ",
      "and set the correct data_root).\n",
      "For SMB: make sure the share is mounted before running the pipeline."
    )

  paths <<- list(
    data_root     = data_root,
    output        = output_root,
    raw_data      = file.path(data_root, dirs$raw_data),
    ext_results   = file.path(data_root, dirs$ext_results),
    prosail_lidar = file.path(data_root, dirs$prosail_lidar),
    prosail_codes = file.path(data_root, dirs$prosail_codes)
  )

  invisible(paths)
})

# Minimal null-coalescing operator (avoids rlang dependency)
`%||%` <- function(a, b) if (!is.null(a)) a else b
