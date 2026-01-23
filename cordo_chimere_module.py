#CHIMERE code to extract data points
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import os

#Cordo ACRA codes to extract data points
def build_acra_filename(year, SC, pol):
    """
    Builds the filename for ACRA data based on year, scenario (SC), and pollutant.
    SC can be None for base year.
    """
    # Standardize pollutant code
    pol_code = {
        "ug_PM25_RH50_high": "PM25",
        "ug_PM25_RH50_low": "PM25",
        "ug_NO2_high": "NO2",
        "ug_NO2_low": "NO2",
        "PM25": "PM25",
        "NO2": "NO2"
    }.get(pol, pol)

    # FRA01 for NO2, FRA02 for PM25
    fra_code = "FRA01" if pol_code == "NO2" else "FRA02"

    # Build filename based on year and scenario
    if year == "2018":
        filename = f"outl.2018_{fra_code}_{pol_code}_analysis_yravg.nc"
    elif year == "2030" and SC in ["AME", "AMS"]:
        filename = f"outl.2030{SC}_{fra_code}_{pol_code}_analysis_yravg.nc"
    else:
        raise ValueError(f"No file defined for year={year}, scenario={SC}, pol={pol_code}")

    return filename, pol_code

def coordo_acra(pol="NO2", year="2018", SC=None):
    print("Starting coordo_acra function")

    # Build filename and standardized pollutant code
    filename, pol_code = build_acra_filename(year, SC, pol)
    path_ACRA = os.path.join("data/1-processed-data/SHERPA/ACRA", filename)

    var = pol_code

    # Check if the file exists
    if not os.path.exists(path_ACRA):
        raise FileNotFoundError(f"File not found: {path_ACRA}")

    # Load the data file
    print(f"Loading data from {path_ACRA}")
    with xr.open_dataset(path_ACRA) as ds:
        # Read coordinates
        latitude = ds["lat"].values
        longitude = ds["lon"].values

        # Read the concentration variable
        if var in ds.data_vars:
            conc = ds[var].squeeze().values
        else:
            # Fallback: use the first data variable if names differ
            var_name = [v for v in ds.data_vars if v not in ['Times_bnds']]
            if not var_name:
                raise ValueError("No suitable data variable found in dataset.")
            conc = ds[var_name[0]].squeeze().values

    # Create a DataFrame with longitude, latitude, and concentration values
    lon, lat = np.meshgrid(longitude, latitude)
    df = pd.DataFrame({
        'x': lon.ravel(),
        'y': lat.ravel(),
        'conc': conc.ravel()
    }).dropna()

    # Convert to a GeoDataFrame
    conc_acra = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.x, df.y), crs="EPSG:2154")

    print("Finished processing coordo_acra function")
    return conc_acra

# ---------------------------
# CHIMERE code to extract data points (All codes)
# ---------------------------

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from netCDF4 import Dataset

# ---------------------------
# Helper: missing value handler
# ---------------------------
def process_conc_array(arr):
    """
    Replace -999 and invalid values with NaN.
    Do NOT fill NaNs with zero.
    """
    arr = np.array(arr, dtype=float)
    arr[arr <= -999] = np.nan
    return arr


# ---------------------------
# Filename builder (unchanged logic)
# ---------------------------
def build_chimere_filename(year, SC, pol):
    if pol.startswith("ug_PM25_RH50"):
        pol_code = "PM25"
    elif pol.startswith("ug_NO2"):
        pol_code = "NO2"
    elif pol in ["PM25", "NO2"]:
        pol_code = pol
    else:
        pol_code = pol

    sc_part = "" if str(year) == "2019" else f"_scen{SC}"
    fra_code = "FRA01" if pol_code == "NO2" else "FRA02"

    filename = f"outl.{year}{sc_part}_{fra_code}_{pol_code}_analysis_yravg.nc"
    return filename, pol_code


# ---------------------------
# FIXED: CHIMERE loader
# ---------------------------
def coordo_chimere(pol="NO2", year="2019", SC="S3"):
    print("Starting coordo_chimere function")

    filename, pol_code = build_chimere_filename(year, SC, pol)
    path_CHIMERE = os.path.join("data", "1-processed-data", "SHERPA", "CHIMERE", filename)

    if not os.path.exists(path_CHIMERE):
        raise FileNotFoundError(f"File not found: {path_CHIMERE}")

    print(f"Loading data from {path_CHIMERE}")

    with xr.open_dataset(path_CHIMERE) as ds:
        # Safer variable selection
        if pol_code in ds.data_vars:
            var = ds[pol_code].squeeze()
        else:
            data_vars = [v for v in ds.data_vars if v.lower() not in ['times_bnds']]
            if not data_vars:
                raise ValueError("No suitable data variable found in dataset.")
            var = ds[data_vars[0]].squeeze()

        # Convert to DataFrame safely (preserves lat/lon mapping)
        df = var.to_dataframe(name="conc").reset_index()

    # Standardize coordinate column names
    if "lat" in df.columns:
        df = df.rename(columns={"lat": "y"})
    if "lon" in df.columns:
        df = df.rename(columns={"lon": "x"})

    # Handle missing values correctly
    df["conc"] = process_conc_array(df["conc"].values)

    # Drop missing rows
    df = df.dropna(subset=["conc", "x", "y"])

    # Convert to GeoDataFrame
    conc_chimere = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.x, df.y),
        crs="EPSG:4326"
    )

    print("Finished processing coordo_chimere function")
    return conc_chimere


# ---------------------------
# FIXED: INERIS (baseline) loader
# ---------------------------
def coordo_ineris_chimere(pol, year=2019):
    print("Starting coordo_ineris function")

    pol = {
        "ug_PM25_RH50_high": "ug_PM25_RH50",
        "ug_PM25_RH50_low": "ug_PM25_RH50",
        "ug_NO2_high": "ug_NO2",
        "ug_NO2_low": "ug_NO2"
    }.get(pol, pol)

    path_INERIS_2019 = {
        "ug_PM25_RH50": "data/1-processed-data/SHERPA/CHIMERE/outl.2019_FRA02_PM25_analysis_yravg.nc",
        "ug_NO2": "data/1-processed-data/SHERPA/CHIMERE/outl.2019_FRA01_NO2_analysis_yravg.nc"
    }.get(pol)

    if path_INERIS_2019 is None:
        raise ValueError(f"Unsupported pollutant: {pol}")

    var_name = "PM25" if pol == "ug_PM25_RH50" else "NO2"

    if not os.path.exists(path_INERIS_2019):
        raise FileNotFoundError(f"File not found: {path_INERIS_2019}")

    print(f"Loading data from {path_INERIS_2019}")

    # Use xarray for consistency and correct mapping
    with xr.open_dataset(path_INERIS_2019) as ds:
        if var_name not in ds.data_vars:
            raise ValueError(f"Variable '{var_name}' not found in dataset.")

        var = ds[var_name].squeeze()
        df = var.to_dataframe(name="conc19").reset_index()

    # Standardize coordinate column names
    if "lat" in df.columns:
        df = df.rename(columns={"lat": "y"})
    if "lon" in df.columns:
        df = df.rename(columns={"lon": "x"})

    # Handle missing values safely
    df["conc19"] = process_conc_array(df["conc19"].values)

    # Remove missing rows
    df = df.dropna(subset=["conc19", "x", "y"])

    # Convert to GeoDataFrame
    conc_ineris = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.x, df.y),
        crs="EPSG:4326"
    )

    print("Finished processing coordo_ineris function")
    return conc_ineris

# ---------------------------
# FIXED: Correction loader
# ---------------------------
from scipy.spatial import cKDTree
def correction_chimere(conc_points, conc_ineris):
    """
    Computes delta as: scenario - baseline (2019)
    """

    # Ensure same CRS
    if conc_points.crs != conc_ineris.crs:
        conc_ineris = conc_ineris.to_crs(conc_points.crs)

    # Extract coordinates
    points_coords = np.column_stack((conc_points.geometry.x, conc_points.geometry.y))
    ineris_coords = np.column_stack((conc_ineris.geometry.x, conc_ineris.geometry.y))

    # KDTree match
    tree = cKDTree(ineris_coords)
    _, idx = tree.query(points_coords, k=1)

    # Values
    scen_conc = conc_points['conc'].values
    base_conc = conc_ineris.iloc[idx]['conc19'].values

    # Correct delta direction
    #delta_conc = scen_conc - base_conc
    delta_conc = base_conc - scen_conc

    return gpd.GeoDataFrame({
        'conc': scen_conc,
        'delta_conc': delta_conc,
        'geometry': conc_points.geometry
    }, crs=conc_points.crs)
