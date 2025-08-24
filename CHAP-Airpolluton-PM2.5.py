import xarray as xr
import pandas as pd
import numpy as np
import os
from datetime import datetime

def extract_pm25_for_cities(file_path, cities_coords):
    """
    Extract PM2.5 values for specific cities from a NetCDF file.
    
    Args:
        file_path (str): Path to the NetCDF file.
        cities_coords (dict): Dict of city names and (lat, lon) tuples, e.g., {'city': (lat, lon)}
        
    Returns:
        dict: PM2.5 values and metadata for each city.
    """
    try:
        ds = xr.open_dataset(file_path)
        print("Variables in dataset:", list(ds.data_vars.keys()))
        
        # Extract date from filename
        filename = os.path.basename(file_path)
        date_str = filename.split('_')[3]
        date = datetime.strptime(date_str, '%Y%m%d').strftime('%Y-%m-%d')
        
        # Try to find the PM2.5 variable name
        pm25_var_name = None
        possible_names = ['PM25', 'pm25', 'PM2_5', 'PM2.5', 'pm2_5', 'pm2.5']
        for var in ds.data_vars.keys():
            if var in possible_names or 'pm25' in var.lower() or 'pm2' in var.lower():
                pm25_var_name = var
                break
        if pm25_var_name is None:
            print(f"PM2.5 variable not found in {filename}. Available variables: {list(ds.data_vars.keys())}")
            ds.close()
            return None
        print(f"Using PM2.5 variable: {pm25_var_name}")
        
        results = {}
        for city, (target_lat, target_lon) in cities_coords.items():
            try:
                # Find nearest grid point
                lat_idx = np.abs(ds.lat.values - target_lat).argmin()
                lon_idx = np.abs(ds.lon.values - target_lon).argmin()
                pm25_value = ds[pm25_var_name].isel(lat=lat_idx, lon=lon_idx).values
                if pm25_value.ndim > 0:
                    pm25_value = pm25_value.item() if pm25_value.size == 1 else pm25_value[0]
                actual_lat = ds.lat.isel(lat=lat_idx).values
                actual_lon = ds.lon.isel(lon=lon_idx).values
                results[city] = {
                    'date': date,
                    'pm25': float(pm25_value) if not np.isnan(pm25_value) else None,
                    'target_lat': target_lat,
                    'target_lon': target_lon,
                    'actual_lat': float(actual_lat),
                    'actual_lon': float(actual_lon),
                    'distance_km': calculate_distance(target_lat, target_lon, float(actual_lat), float(actual_lon))
                }
            except Exception as city_e:
                print(f"Error extracting data for city {city}: {city_e}")
                results[city] = {
                    'date': date,
                    'pm25': None,
                    'target_lat': target_lat,
                    'target_lon': target_lon,
                    'actual_lat': None,
                    'actual_lon': None,
                    'distance_km': None
                }
        ds.close()
        return results
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None

def calculate_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the distance between two (lat, lon) points in kilometers.
    """
    from math import radians, cos, sin, asin, sqrt
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    distance = 2 * asin(sqrt(a)) * 6371
    return round(distance, 2)

def batch_process_pm25_files(folder_path, cities_coords, output_file):
    """
    Batch process NetCDF files in a folder for PM2.5 extraction.
    
    Args:
        folder_path (str): Folder containing NetCDF files.
        cities_coords (dict): Dict of city names and (lat, lon) tuples.
        output_file (str): Path to save the output CSV.
        
    Returns:
        DataFrame or None
    """
    all_results = []
    nc_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.nc')])
    print(f"Found {len(nc_files)} NetCDF files.")
    for i, filename in enumerate(nc_files, 1):
        file_path = os.path.join(folder_path, filename)
        print(f"Processing ({i}/{len(nc_files)}): {filename}")
        file_results = extract_pm25_for_cities(file_path, cities_coords)
        if file_results:
            for city, data in file_results.items():
                row = {
                    'date': data['date'],
                    'city': city,
                    'pm25': data['pm25'],
                    'target_lat': data['target_lat'],
                    'target_lon': data['target_lon'],
                    'actual_lat': data['actual_lat'],
                    'actual_lon': data['actual_lon'],
                    'distance_km': data['distance_km']
                }
                all_results.append(row)
    if all_results:
        df = pd.DataFrame(all_results)
        df = df.sort_values(['city', 'date'])
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\nData saved to {output_file}")
        print(f"Processed {len(nc_files)} files, {len(cities_coords)} cities")
        print(f"Data shape: {df.shape}")
        print("\nPreview of first 5 rows:")
        print(df.head())
        print("\nPM2.5 statistics:")
        print(df['pm25'].describe())
        return df
    else:
        print("No data extracted.")
        return None

# Example usage
if __name__ == "__main__":
    # Define city coordinates
    cities_coords = {
        'Jiangnanhospotal-Observation': (31.484705, 120.27716),
        # Add more cities as needed
    }
    folder_path = r"C:\Users\xieru\Desktop\PM2.5ASD\2023"
    output_file = r"C:\Users\xieru\Desktop\pm25_data_2023.csv"
    df = batch_process_pm25_files(folder_path, cities_coords, output_file)
    if df is not None:
        print("\nAverage PM2.5 by city:")
        print(df.groupby('city')['pm25'].mean().sort_values(ascending=False))
        print("\nMissing data count by city:")
        print(df.groupby('city')['pm25'].apply(lambda x: x.isnull().sum()))
