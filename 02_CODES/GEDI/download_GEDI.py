import os
import geopandas as gpd
import earthaccess

# PYTHONPATH="" /home/corroyez/miniconda3/envs/gedi_env/bin/python /home/corroyez/Documents/NC_Full/02_CODES/GEDI/download_GEDI.py

# 1. Earthdata login (uses ~/.netrc or direct credentials)
auth = earthaccess.login(strategy="netrc")  # assumes ~/.netrc exists with login info

# Or if you prefer to define them explicitly:
# auth = earthaccess.login(strategy="environment")

# Environment variables must be set:
# os.environ["EARTHDATA_USER"] = "your_username"
# os.environ["EARTHDATA_PASS"] = "your_password"

# 2. Your main data directory
data_dir = "/home/corroyez/Documents/NC_Full/01_DATA"
hard_dir = "/media/corroyez/MyPassport/01_DATA"
# sites = ["Aigoual", "Blois", "Mormal"]
sites = ["Aigoual", "Blois", "Hayes", "Mormal", "Reine"]

# 3. Loop over sites
for site in sites:
    print(f"Processing {site}...")

    # --- Read AOI ---
    aoi_path = os.path.join(hard_dir, site, "Geo_Files", "utm_from_lidar.gpkg")
    aoi = gpd.read_file(aoi_path)
    aoi_wgs = aoi.to_crs(epsg=4326)
    bbox = aoi_wgs.total_bounds  # [minx, miny, maxx, maxy]

    # Extract bounding box values
    ul_lon, lr_lon = bbox[0], bbox[2]
    lr_lat, ul_lat = bbox[1], bbox[3]

    # --- Output directory ---
    # outdir = os.path.join(data_dir, site, "GEDI")
    outdir = os.path.join(hard_dir, site, "GEDI")
    os.makedirs(outdir, exist_ok=True)

    # --- Date range ---
    daterange = ("2019-04-18", "2023-03-16")

    # --- Products to query ---
    # products = ["GEDI02_B"]
    products = ["GEDI02_A", "GEDI02_B"]

    for product in products:
        print(f"  Searching {product} data...")

        # Query data from LP DAAC
        results = earthaccess.search_data(
            short_name=product,
            version="002",
            temporal=daterange,
            bounding_box=(ul_lon, lr_lat, lr_lon, ul_lat)
        )

        if not results:
            print(f"  ⚠️  No {product} data found for {site}.")
            continue

        print(f"  Found {len(results)} {product} files — downloading...")

        # Download data
        earthaccess.download(results, outdir)

    print(f"✅ Finished downloading GEDI data for {site}\n")
