# ----------------------------- (Optional) Clear the environment and free memory -------------------------------------

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory

# --------------------------------------------------------------------------------------------------------------------

# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path));getwd()
}


library(httr2)
library(rstac)
library(yaml)

get_copdataspace_token <- function(credential_file){
  creds = yaml::read_yaml(credential_file)
  cop_user = creds$username
  cop_pass = creds$password
  cop_token_url <- "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
  
  client <- httr2::oauth_client(id = "cdse-public", token_url = cop_token_url)
  token <- httr2::oauth_flow_password(client, username = cop_user, password = cop_pass)
  return(list(client=client, token=token))
}

# Download Copernicus DataSpace products listed from STAC search.
#
# @param x rstac::doc_items
# @param outdir output dir
# @param credential_file yaml file with Copernicus DataSpace credentials.
# The expected Copernicus DataSpace credential file content:
#   username: "xxxx@email.org" 
#   password: "xxxxx"
# @param dry_run logical. If TRUE, the authetication is tried but the download is not performed.
# @param overwrite logical. If TRUE, existing files are overwritten.
#
# @return list of downloaded files
download_copdataspace_products <- function(x, outdir, credential_file, dry_run=T, overwrite=F){
  if(!inherits(x, "doc_items")){
    stop("x must be of class 'rstac::doc_items'")
  }
  
  if(!dir.exists(outdir)){
    stop("Directory not found: ", outdir)
  }
  
  res = get_copdataspace_token(credential_file)
  client = res$client
  token = res$token
  
  if(dry_run){
    message("Running as dry run: printing only the destination paths, without downloading files.")
  }
  zipfiles = c()
  for(f in x$features){
    outfile = file.path(outdir, f$id)
    zipfile = paste0(outfile, ".zip")
    if(file.exists(zipfile) && !overwrite){
      message(sprintf("Skipping: %s", f$id))
      next
    }
    message(sprintf("Downloading: %s --> %s", f$id, zipfile))
    
    product_url = f$assets[["PRODUCT"]]$href
    # print(product_url)
    # if expired:
    if(Sys.time() >= token$expires_at)
      token = httr2::oauth_flow_refresh(client, token$refresh_token)
    
    req <- httr2::request(product_url) |>
      httr2::req_auth_bearer_token(token = token$access_token) |>
      # With r-curl < 5.2.1, unrestricted_auth=1 was the default.
      # Changed since r-curl 5.2.1, see https://github.com/jeroen/curl/blob/master/NEWS
      # 
      # It seems necessary to add unrestricted_auth
      # as Copernicus Dataspace would systematically redirect from
      # https://catalogue.dataspace.copernicus.eu/.... to
      # https://download.dataspace.copernicus.eu/....
      # set verbosity=1 at req_perform to see that.
      # See https://github.com/r-lib/httr2/issues/475 for security details.
      # See also https://superuser.com/questions/936042/how-can-i-instruct-curl-to-reuse-credentials-after-it-followed-a-redirect
      req_options(unrestricted_auth = 1) 
    
    if(!dry_run){
      resp <- httr2::req_perform(req, zipfile) # add verbosity=1 to see the details
    }
    zipfiles = c(zipfiles, zipfile)
  }
  return(zipfiles)
}

####### Exemple #########
credential_file = "../../01_DATA/copernicus-credentials.yml"
# expected Copernicus DataSpace credential file content:
# username: "xxxx@email.org" 
# password: "xxxxx"

s_obj <- stac("https://catalogue.dataspace.copernicus.eu/stac")
products <- res <- s_obj |> stac_search(
  collections="SENTINEL-2",
  datetime = "2021-11-17/2021-11-20") |> # "2021-11-10/2021-11-20"
  ext_filter(tileId == "31TEJ" && 
               productType == "S2MSI2A" && 
               id %like% "%\\_N0301\\_%" # % is the wildcard, _ is the single character wildcard, \\ is the escape
  ) |> 
  post_request() |>
  items_fetch() |> 
  assets_select(asset_names="PRODUCT")

outdir = getwd()
# unlink(outdir, recursive=T)
# dir.create(outdir)

files <- download_copdataspace_products(products, outdir, credential_file, dry_run=F)