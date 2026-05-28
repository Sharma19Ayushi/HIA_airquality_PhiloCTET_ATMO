import geopandas as gpd
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import box

import geopandas as gpd
import numpy as np
from scipy.spatial import cKDTree

def correction(conc_points, conc_ineris):
    """
    Correction function equivalent to the R version using nearest-neighbor matching.
    Adjusts delta_conc and conc based on INERIS conc19.
    """

    # Ensure both GeoDataFrames use the same CRS
    if conc_points.crs != conc_ineris.crs:
        conc_ineris = conc_ineris.to_crs(conc_points.crs)

    # --- Nearest Neighbor Search (equivalent to st_nn in R) ---
    # Extract coordinates
    points_coords = np.column_stack((conc_points.geometry.x, conc_points.geometry.y))
    ineris_coords = np.column_stack((conc_ineris.geometry.x, conc_ineris.geometry.y))

    # Build KDTree for conc_ineris
    tree = cKDTree(ineris_coords)

    # Query nearest INERIS point for each point in conc_points
    _, idx = tree.query(points_coords, k=1)

    # Extract nearest conc19 values
    nearest_conc19 = conc_ineris.iloc[idx]['conc19'].values

    # Original values
    conc = conc_points['conc'].values
    delta_conc = conc_points['delta_conc'].values

    # --- Apply exact same correction formula as in R ---
    # new_delta_conc = delta_conc * nearest_conc19 / conc
    # new_conc = nearest_conc19 * (1 - delta_conc / conc)

    with np.errstate(divide='ignore', invalid='ignore'):
        new_delta_conc = np.where(conc != 0, delta_conc * nearest_conc19 / conc, 0)
        new_conc = np.where(conc != 0, nearest_conc19 * (1 - delta_conc / conc), nearest_conc19)
        # Final safety clamp
        new_conc = np.maximum(new_conc, 0.0)
        new_delta_conc = np.maximum(new_delta_conc, 0.0)

    # --- Create corrected GeoDataFrame (st_sf equivalent) ---
    corrected = gpd.GeoDataFrame({
        'conc': new_conc,
        'delta_conc': new_delta_conc,
        'geometry': conc_points.geometry
    }, crs=conc_points.crs)
    # --- Diagnostic prints (like R interactive diagnostics) ---
    print("Correction diagnostics:")
    print(f"Mean conc before correction: {np.nanmean(conc):.4f}")
    print(f"Mean delta_conc before correction: {np.nanmean(delta_conc):.4f}")
    print(f"Mean conc after correction: {np.nanmean(new_conc):.4f}")
    print(f"Mean delta_conc after correction: {np.nanmean(new_delta_conc):.4f}")

    return corrected