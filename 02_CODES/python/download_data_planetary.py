import os
import requests
import pandas as pd

# Define your API key and endpoint
API_KEY = '9e1a5f8a49df4346986d94f82fb13fcc'
ENDPOINT = 'https://planetarycomputer.microsoft.com/api/stac/v1/'

# Define function to download image
def download_image(url, filename):
    response = requests.get(url)
    if response.status_code == 200:
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {filename}")
    else:
        print(f"Failed to download {filename}")

# Define function to search for items
def search_items(query_params):
    response = requests.get(f"{ENDPOINT}/collections/sentinel-2-l2a-boa/items", params=query_params, headers={"Authorization": f"Bearer {API_KEY}"})
    return response.json()

# Define function to download images for given items
def download_images(data):
    for index, row in data.iterrows():
        query_params = {
            'bbox': f"{row['Lon_IS']},{row['Lat_IS']},{row['Lon_IS']},{row['Lat_IS']}",  
            'datetime': row['TIME_IS'],
            'limit': 1  
        }
        print(query_params)
        items = search_items(query_params)
        if 'features' in items and len(items['features']) > 0:
            item = items['features'][0]
            asset_url = item['assets']['data']['href']
            download_image(asset_url, f"{row['GBOV_ID']}.tif")
        else:
            print(f"No image found for {row['GBOV_ID']}")

def main():
    # Load CSV data
    df = pd.read_csv("/home/corroyez/Documents/NC_Full/01_DATA/COPERNICUS_GBOV_RM7_20240227101352/Deciduous_Forests/june_unique_plot_data.csv")

    # Download images
    download_images(df)

if __name__ == "__main__":
    main()
