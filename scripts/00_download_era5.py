import cdsapi
import numpy as np
import os
import zipfile

# -----------------------------
# USER SETTINGS
# -----------------------------

# Starting point of the grid
lat_start = 54.2
lon_start = -3.5

# Grid size: 1 degrees × 1 degrees
grid_size = 1

# Resolution: 0.1 degree
step = 0.1

# Variables to download
variables = [
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_solar_radiation_downwards",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind"
]

# Output folder
output_root = "data/era5_grid/"
os.makedirs(output_root, exist_ok=True)

# -----------------------------
# GENERATE GRID
# -----------------------------

lats = np.arange(lat_start, lat_start + grid_size, step)
lons = np.arange(lon_start, lon_start + grid_size, step)

# -----------------------------
# CDS API CLIENT
# -----------------------------
client = cdsapi.Client()

# -----------------------------
# LOOP THROUGH GRID POINTS
# -----------------------------

for lat in lats:
    for lon in lons:

        # Create folder for this point
        point_folder = f"{output_root}/{lat:.2f}_{lon:.2f}/"
        os.makedirs(point_folder, exist_ok=True)

        # Output ZIP file
        zip_path = f"{point_folder}/download.zip"
        extract_folder = f"{point_folder}/nc/"
        os.makedirs(extract_folder, exist_ok=True)
        
        # --- SKIP IF DATA ALREADY DOWNLOADED ---
        nc_files = [f for f in os.listdir(extract_folder) if f.endswith(".nc")]
        if len(nc_files) > 0:
            print(f"Skipping {lat:.2f}, {lon:.2f} — already downloaded.")
            continue
        
        print(f"Downloading ERA5-Land for point {lat:.2f}, {lon:.2f}")


        dataset = "reanalysis-era5-land-timeseries"
        request = {
            "variable": variables,
            "location": {"longitude": lon, "latitude": lat},
            "date": ["1995-01-01/2026-03-13"],
            "data_format": "netcdf"
        }

        client.retrieve(dataset, request).download(zip_path)


        # -----------------------------
        # UNZIP AND ORGANIZE FILES
        # -----------------------------
        
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(extract_folder)

        print(f"Extracted files for {lat:.2f}, {lon:.2f}")

        # Optional: delete ZIP to save space
        os.remove(zip_path)
