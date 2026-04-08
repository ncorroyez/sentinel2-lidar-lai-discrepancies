# ==============================================================================
# Download HR-VPP (LAI, FCOVER) from Copernicus DataSpace Ecosystem
# Description: Uses STAC API and OAuth2 to fetch High Resolution Vegetation 
#              Phenology and Productivity products based on a shapefile extent.
# ==============================================================================

# ----------------------------- Clear environment ------------------------------
rm(list=ls(all=TRUE))
gc()

# ----------------------------- Working Directory ------------------------------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  cat("Working directory set to:", getwd(), "\n")
}

# ----------------------------- Libraries --------------------------------------
library(httr2)
library(rstac)
library(yaml)
library(sf)

# ==============================================================================
# FUNCTIONS (Kept from original with minor English translations)
# ==============================================================================

get_copdataspace_token <- function(credential_file){
  creds <- yaml::read_yaml(credential_file)
  cop_user <- creds$username
  cop_pass <- creds$password
  cop_token_url <- "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
  
  client <- httr2::oauth_client(id = "cdse-public", token_url = cop_token_url)
  token <- httr2::oauth_flow_password(client, username = cop_user, password = cop_pass)
  return(list(client = client, token = token))
}

download_copdataspace_products <- function(x, outdir, credential_file, dry_run=TRUE, overwrite=FALSE){
  if(!inherits(x, "doc_items")){
    stop("x must be of class 'rstac::doc_items'")
  }
  if(!dir.exists(outdir)){
    stop("Directory not found: ", outdir)
  }
  
  res <- get_copdataspace_token(credential_file)
  client <- res$client
  token <- res$token
  
  if(dry_run){
    message("Running as dry run: printing only the destination paths, without downloading files.")
  }
  
  zipfiles <- c()
  for(f in x$features){
    # HR-VPP products are usually delivered as .tar or .zip depending on the exact asset
    outfile <- file.path(outdir, f$id)
    zipfile <- paste0(outfile, ".tar") 
    
    if(file.exists(zipfile) && !overwrite){
      message(sprintf("Skipping (already exists): %s", f$id))
      next
    }
    message(sprintf("Downloading: %s --> %s", f$id, basename(zipfile)))
    
    product_url <- f$assets[["PRODUCT"]]$href
    
    if(Sys.time() >= token$expires_at){
      token <- httr2::oauth_flow_refresh(client, token$refresh_token)
    }
    
    req <- httr2::request(product_url) |>
      httr2::req_auth_bearer_token(token = token$access_token) |>
      req_options(unrestricted_auth = 1) 
    
    if(!dry_run){
      resp <- httr2::req_perform(req, zipfile)
    }
    zipfiles <- c(zipfiles, zipfile)
  }
  return(zipfiles)
}

# ==============================================================================
# MAIN EXECUTION: STAC SEARCH & DOWNLOAD
# ==============================================================================

credential_file <- "../../01_DATA/copernicus-credentials.yml"

# --- 1. Load AOI and extract WGS84 Bounding Box ---
# STAC queries require WGS84 (EPSG:4326) coordinates
v_aoi <- st_read("../../01_DATA//utm_init.shp", quiet = TRUE)
v_aoi_wgs84 <- st_transform(v_aoi, 4326)
bbox <- st_bbox(v_aoi_wgs84)

cat("AOI Bounding Box (WGS84):\n")
print(bbox)

# --- 2. Build and execute STAC query ---
s_obj <- stac("https://catalogue.dataspace.copernicus.eu/stac")

# collection "hrvpp-vi" contains Vegetation Indices (LAI, FCOVER, FAPAR, etc.)
products <- s_obj |> 
  stac_search(
    collections = "hrvpp-vi",
    bbox = c(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]),
    datetime = "2025-06-01T00:00:00Z/2024-09-30T23:59:59Z"
  ) |> 
  post_request() |>
  items_fetch() |> 
  assets_select(asset_names = "PRODUCT")

cat(sprintf("Found %d HR-VPP products intersecting the AOI.\n", length(products$features)))

# --- 3. Download the products ---
outdir <- file.path(getwd(), "out_files", "HRVPP_Raw")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Set dry_run = FALSE to actually download the files
files <- download_copdataspace_products(
  x = products, 
  outdir = outdir, 
  credential_file = credential_file, 
  dry_run = TRUE, 
  overwrite = FALSE
)