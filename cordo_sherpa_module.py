import os
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import unicodedata


# Helper to normalize input string
def normalize_string(value):
    return unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('utf-8').lower()

def coordo_sherpa(sc="", pol="", year=""):
    # Normalize and validate the scenario input
    sc = normalize_string(sc)
    pol = normalize_string(pol)

    valid_scenarios = ["s1", "s2", "s3", "s4"]
    if sc not in valid_scenarios:
        raise ValueError(f"The specified scenario ('sc') must be one of the following: {valid_scenarios}")

    # Standardize pollutant names
    pol = {
        "ug_pm25_rh50_high": "ug_pm25_rh50",
        "ug_pm25_rh50_low": "ug_pm25_rh50",
        "ug_no2_high": "ug_no2",
        "ug_no2_low": "ug_no2"
    }.get(pol, pol)

    # Paths to data (ensure UTF-8 compatibility)
    path_SHERPA_delta = rf"data/1-processed-data/SHERPA/scenarios/{sc}/DCconc_{sc}_{year}_SURF_{pol}.nc"
    path_SHERPA_2019 = rf"data/1-processed-data/SHERPA/conc-2019/BCconc_emepV4_45_cams61_withCond_01005_2019_SURF_{pol}.nc"
    path_AREA = r"data/1-processed-data/SHERPA/emiRedOn_01005_France_new.nc"

    # Load files using xarray
    ds_delta = xr.open_dataset(path_SHERPA_delta)
    ds_2019 = xr.open_dataset(path_SHERPA_2019)
    ds_area = xr.open_dataset(path_AREA)

    # Extract and process data
    delta_conc = ds_delta["delta_conc"].values
    conc_2019 = ds_2019["conc"].values
    area = ds_area["AREA"].values
    latitude = ds_delta["latitude"].values
    longitude = ds_delta["longitude"].values

    # Filter and reshape the data
    mask = (area > 0) & ~np.isnan(delta_conc)
    filtered_conc = np.where(mask, conc_2019, 0)
    filtered_delta_conc = np.where(mask, delta_conc, 0)

    longitude, latitude = np.meshgrid(longitude, latitude)
    lon_lat_pairs = np.vstack([longitude[mask], latitude[mask]]).T
    values = np.vstack([filtered_conc[mask], filtered_delta_conc[mask]]).T

    # Create GeoDataFrame (ensure encoding if French names are added)
    conc_points = gpd.GeoDataFrame(values, columns=["conc", "delta_conc"],
                                   geometry=[Point(xy) for xy in lon_lat_pairs])
    conc_points.set_crs("EPSG:4326", inplace=True)

    print(f"Concentrations in 2019 and {year} are calculated for the pollutant '{pol}' ({sc}).")

    return conc_points

#Cordo_sherpa based on new grid of EMEP and SRR 2021
import xarray as xr
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

def coordo_sherpa(sc="", pol="", year=""):
    sc = normalize_string(sc)
    pol = normalize_string(pol)
    valid_scenarios = ["s1", "s2", "s3", "s4"]
    if sc not in valid_scenarios:
        raise ValueError(f"The specified scenario ('sc') must be one of the following: {valid_scenarios}")
    pol = {
        "ug_pm25_rh50_high": "ug_pm25_rh50",
        "ug_pm25_rh50_low": "ug_pm25_rh50",
        "ug_no2_high": "ug_no2",
        "ug_no2_low": "ug_no2"
    }.get(pol, pol)
    path_SHERPA_delta = rf"data/1-processed-data/SHERPA/scenarios/{sc}/DCconc_{sc}_{year}_SURF_{pol}.nc"
    path_SHERPA_2019 = rf"data/1-processed-data/SHERPA/conc-2019/BCconc_emepV45_cams80_SURF_{pol}.nc"
    path_AREA = r"data/1-processed-data/SHERPA/emiRedOn_01005_France.nc"
    ds_delta = xr.open_dataset(path_SHERPA_delta)
    ds_2019 = xr.open_dataset(path_SHERPA_2019)
    ds_area = xr.open_dataset(path_AREA)

    # Interpolate base concentration and area to delta grid
    conc_2019 = ds_2019['conc'].interp_like(ds_delta).values
    area = ds_area['AREA'].interp_like(ds_delta).values
    delta_conc = ds_delta["delta_conc"].values
    latitude = ds_delta["latitude"].values
    longitude = ds_delta["longitude"].values

    # 2D mask
    mask = (area > 0) & ~np.isnan(delta_conc)

    # Find indices where mask is True
    idx_y, idx_x = np.where(mask)
    # Use meshgrid to get correct coordinates for each cell
    lon2d, lat2d = np.meshgrid(longitude, latitude)

    # Get values at valid indices
    lon_valid = lon2d[mask]
    lat_valid = lat2d[mask]
    conc_valid = conc_2019[mask]
    delta_valid = delta_conc[mask]

    # Stack arrays for GeoDataFrame
    values = np.vstack([conc_valid, delta_valid]).T
    lon_lat_pairs = np.vstack([lon_valid, lat_valid]).T

    conc_points = gpd.GeoDataFrame(values, columns=["conc", "delta_conc"],
                                   geometry=[Point(xy) for xy in lon_lat_pairs])
    conc_points.set_crs("EPSG:4326", inplace=True)
    print(f"Concentrations in 2019 and {year} are calculated for the pollutant '{pol}' ({sc}).")
    return conc_points

#For GNFR analysis, loop over sectors (g) in cordo_sherpa
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point

def coordo_sherpa_new(sc="s1", pol="ug_PM25_RH50", year="2019", g="agri"):
    # Validate the scenario input
    valid_scenarios = ["s1", "s2", "s3", "s4"]
    if sc not in valid_scenarios:
        raise ValueError("The specified scenario ('sc') must be one of the following: 's1', 's2', 's3', 's4'.")

    # Standardize pollutant names
    if pol in ["ug_PM25_RH50_high", "ug_PM25_RH50_low"]:
        pol = "ug_PM25_RH50"
    elif pol in ["ug_NO2_high", "ug_NO2_low"]:
        pol = "ug_NO2"

    # Paths to data
    path_SHERPA_delta = f"data/1-processed-data/SHERPA/scenarios/{sc}/DCconc_{sc}_{year}_{g}_SURF_{pol}.nc"
    path_SHERPA_2019 = f"data/1-processed-data/SHERPA/conc-2019/BCconc_emepV4_45_cams61_withCond_01005_2019_SURF_{pol}.nc"
    path_AREA = "data/1-processed-data/SHERPA/emiRedOn_01005_France.nc"

    try:
        # Load files using xarray
        ds_delta = xr.open_dataset(path_SHERPA_delta)
        ds_2019 = xr.open_dataset(path_SHERPA_2019)
        ds_area = xr.open_dataset(path_AREA)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return None

    delta_conc = ds_delta["delta_conc"].values
    conc_2019 = ds_2019["conc"].values
    area = ds_area["AREA"].values
    latitude = ds_delta["latitude"].values
    longitude = ds_delta["longitude"].values

    # Check if data frames have the expected dimensions
    if delta_conc.shape != area.shape:
        raise ValueError("Dimension mismatch between area and delta_conc.")

    # Filter and reshape the data
    mask = (area > 0) & ~np.isnan(delta_conc)
    filtered_conc = np.where(mask, conc_2019, 0)
    filtered_delta_conc = np.where(mask, delta_conc, 0)

    longitude, latitude = np.meshgrid(longitude, latitude)
    lon_lat_pairs = np.vstack([longitude[mask], latitude[mask]]).T
    values = np.vstack([filtered_conc[mask], filtered_delta_conc[mask]]).T

    # Create GeoDataFrame
    conc_points = gpd.GeoDataFrame(values, columns=["conc", "delta_conc"],
                                   geometry=[Point(xy) for xy in lon_lat_pairs])
    conc_points.set_crs(epsg=4326, inplace=True)
    print(f"Concentrations in 2019 and {year} are calculated for the pollutant {pol} according to scenario {sc}.")

    # Save the conc_points as a shapefile if needed
    # shapefile_path = f"data/output/conc_points_{sc}_{pol}_{year}.shp"
    # conc_points.to_file(shapefile_path, driver='ESRI Shapefile')

    # List the files in the working directory
    print("Files in the current working directory:")
    print(os.listdir())

    return conc_points


def coordo_sherpa_gnfr(sc="s1", pol="ug_PM25_RH50", year="2019", g="GNFR1"):
    # Validate the scenario input
    valid_scenarios = ["s1", "s2", "s3", "s4"]
    if sc not in valid_scenarios:
        raise ValueError("The specified scenario ('sc') must be one of the following: 's1', 's2', 's3', 's4'.")

    # Standardize pollutant names
    if pol in ["ug_PM25_RH50_high", "ug_PM25_RH50_low"]:
        pol = "ug_PM25_RH50"
    elif pol in ["ug_NO2_high", "ug_NO2_low"]:
        pol = "ug_NO2"

    # Paths to data
    path_SHERPA_delta = f"data/1-processed-data/SHERPA/scenarios/sectorbysector/decrease_{g}_SURF_{pol}.nc"
    path_SHERPA_2019 = f"data/1-processed-data/SHERPA/conc-2019/BCconc_emepV4_45_cams61_withCond_01005_2019_SURF_{pol}.nc"
    path_AREA = "data/1-processed-data/SHERPA/emiRedOn_01005_France.nc"

    try:
        # Load files using xarray
        ds_delta = xr.open_dataset(path_SHERPA_delta)
        ds_2019 = xr.open_dataset(path_SHERPA_2019)
        ds_area = xr.open_dataset(path_AREA)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return None

    delta_conc = ds_delta["delta_conc"].values
    conc_2019 = ds_2019["conc"].values
    area = ds_area["AREA"].values
    latitude = ds_delta["latitude"].values
    longitude = ds_delta["longitude"].values

    # Check if data frames have the expected dimensions
    if delta_conc.shape != area.shape:
        raise ValueError("Dimension mismatch between area and delta_conc.")

    # Filter and reshape the data
    mask = (area > 0) & ~np.isnan(delta_conc)
    filtered_conc = np.where(mask, conc_2019, 0)
    filtered_delta_conc = np.where(mask, delta_conc, 0)

    longitude, latitude = np.meshgrid(longitude, latitude)
    lon_lat_pairs = np.vstack([longitude[mask], latitude[mask]]).T
    values = np.vstack([filtered_conc[mask], filtered_delta_conc[mask]]).T

    # Create GeoDataFrame
    conc_points = gpd.GeoDataFrame(values, columns=["conc", "delta_conc"],
                                   geometry=[Point(xy) for xy in lon_lat_pairs])
    conc_points.set_crs(epsg=4326, inplace=True)
    print(f"Concentrations in 2019 and {year} are calculated for the pollutant {pol} according to scenario {sc}.")

    # Save the conc_points as a shapefile if needed
    # shapefile_path = f"data/output/conc_points_{sc}_{pol}_{year}.shp"
    # conc_points.to_file(shapefile_path, driver='ESRI Shapefile')

    # List the files in the working directory
    print("Files in the current working directory:")
    print(os.listdir())

    return conc_points

def coordo_sherpa_EU(sc="s1", pol="ug_PM25_RH50", year="2019", g="agri"):
    # Validate the scenario input
    valid_scenarios = ["s1", "s2", "s3", "s4"]
    if sc not in valid_scenarios:
        raise ValueError("The specified scenario ('sc') must be one of the following: 's1', 's2', 's3', 's4'.")

    # Standardize pollutant names
    if pol in ["ug_PM25_RH50_high", "ug_PM25_RH50_low"]:
        pol = "ug_PM25_RH50"
    elif pol in ["ug_NO2_high", "ug_NO2_low"]:
        pol = "ug_NO2"

    # Paths to data
    path_SHERPA_delta = f"data/1-processed-data/SHERPA/scenarios/{sc}/DCconc_EU_{sc}_{year}_{g}_SURF_{pol}.nc"
    path_SHERPA_2019 = f"data/1-processed-data/SHERPA/conc-2019/BCconc_emepV4_45_cams61_withCond_01005_2019_SURF_{pol}.nc"
    path_AREA = "data/1-processed-data/SHERPA/emiRedOn_01005_France.nc"

    try:
        # Load files using xarray
        ds_delta = xr.open_dataset(path_SHERPA_delta)
        ds_2019 = xr.open_dataset(path_SHERPA_2019)
        ds_area = xr.open_dataset(path_AREA)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return None

    delta_conc = ds_delta["delta_conc"].values
    conc_2019 = ds_2019["conc"].values
    area = ds_area["AREA"].values
    latitude = ds_delta["latitude"].values
    longitude = ds_delta["longitude"].values

    # Check if data frames have the expected dimensions
    if delta_conc.shape != area.shape:
        raise ValueError("Dimension mismatch between area and delta_conc.")

    # Filter and reshape the data
    mask = (area > 0) & ~np.isnan(delta_conc)
    filtered_conc = np.where(mask, conc_2019, 0)
    filtered_delta_conc = np.where(mask, delta_conc, 0)

    longitude, latitude = np.meshgrid(longitude, latitude)
    lon_lat_pairs = np.vstack([longitude[mask], latitude[mask]]).T
    values = np.vstack([filtered_conc[mask], filtered_delta_conc[mask]]).T

    # Create GeoDataFrame
    conc_points = gpd.GeoDataFrame(values, columns=["conc", "delta_conc"],
                                   geometry=[Point(xy) for xy in lon_lat_pairs])
    conc_points.set_crs(epsg=4326, inplace=True)
    print(f"Concentrations in 2019 and {year} are calculated for the pollutant {pol} according to scenario {sc}.")

    # Save the conc_points as a shapefile if needed
    # shapefile_path = f"data/output/conc_points_{sc}_{pol}_{year}.shp"
    # conc_points.to_file(shapefile_path, driver='ESRI Shapefile')

    # List the files in the working directory
    print("Files in the current working directory:")
    print(os.listdir())

    return conc_points


