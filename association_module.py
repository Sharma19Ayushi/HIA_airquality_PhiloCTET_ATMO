import geopandas as gpd
from shapely.geometry import Polygon
import pandas as pd
from geopandas import GeoDataFrame
from shapely.ops import unary_union
import multiprocessing
from functools import partial
import os


def generate_points(polygone, conc_points):
    """
    Generate grid points for a given polygon and associate them with concentration points.
    """
    if polygone.geometry.is_empty or not polygone.geometry.is_valid:
        return gpd.GeoDataFrame()  # Return empty GeoDataFrame if invalid

    if conc_points is None or conc_points.empty:
        return gpd.GeoDataFrame()  # Return empty GeoDataFrame if conc_points is empty

    # Filter concentration points within the polygon
    points_within = conc_points[conc_points.geometry.within(polygone.geometry)].copy()
    points_within['iriscod'] = polygone['iriscod']

    return points_within


def calculate_intersection_percentages(grille_combinee, donnees_exportees_transformed):
    """
    Iteratively calculate the intersection percentages for each grid point's rectangle and associated polygon.
    Returns:
        - GeoDataFrame: grille_combinee with calculated 'perc' values.
    """
    # Initialize 'perc' column with zeros, explicitly as float
    grille_combinee['perc'] = 0.0

    for idx, row in grille_combinee.iterrows():
        x, y = row.geometry.x, row.geometry.y

        # Create a rectangle around the point
        rectangle = Polygon([
            (x - 0.05, y - 0.025),
            (x + 0.05, y - 0.025),
            (x + 0.05, y + 0.025),
            (x - 0.05, y + 0.025),
            (x - 0.05, y - 0.025)
        ])

        # Ensure geometry is consistent
        rectangle = gpd.GeoSeries([rectangle], crs=grille_combinee.crs)

        # Filter polygon corresponding to the current grid point
        iriscod = row['iriscod']
        polygon_row = donnees_exportees_transformed[
            donnees_exportees_transformed['iriscod'] == iriscod
            ]

        if polygon_row.empty:
            continue  # Skip if no corresponding polygon exists

        polygon_geom = polygon_row.geometry.iloc[0]

        # Calculate intersection
        intersection = rectangle.iloc[0].intersection(polygon_geom)

        # Check if the intersection is valid and not empty
        if not intersection.is_empty:  # Check for single geometry
            area_intersection = intersection.area
            area_polygon = polygon_geom.area

            if area_polygon > 0:  # Avoid division by zero
                perc = round(area_intersection / area_polygon, 3)
                grille_combinee.at[idx, 'perc'] = perc

    # Filter out rows where intersection percentage is zero
    return grille_combinee[grille_combinee['perc'] > 0]


def process_subset(subset, conc_points):
    """
    Process a subset of polygons with concentration points.
    """
    subset_results = []
    for idx, row in subset.iterrows():
        points = generate_points(row, conc_points)
        subset_results.append(points)

    return GeoDataFrame(pd.concat(subset_results, ignore_index=True) if subset_results else [])


def worker_function(donnees_subset, conc_points):
    """
    Wrapper to process subsets of polygons using multiprocessing.
    """
    try:
        return process_subset(donnees_subset, conc_points)
    except Exception as e:
        print(f"Worker function error: {e}")
        return gpd.GeoDataFrame()


def association(donnees_exportees_transformed, conc_points):
    """
    Perform the association operation and return the resulting GeoDataFrame.
    """
    # Ensure input GeoDataFrames are valid and CRS is assigned
    if not isinstance(donnees_exportees_transformed, GeoDataFrame):
        raise TypeError("donnees_exportees_transformed must be a GeoDataFrame")
    if not isinstance(conc_points, GeoDataFrame):
        raise TypeError("conc_points must be a GeoDataFrame")

    if donnees_exportees_transformed.crs is None or conc_points.crs is None:
        # Assign default CRS if missing
        print("Assigning default CRS 'EPSG:4326' to inputs...")
        donnees_exportees_transformed.set_crs("EPSG:4326", inplace=True)
        conc_points.set_crs("EPSG:4326", inplace=True)

    # Split the data into smaller subsets for multiprocessing
    num_cores = multiprocessing.cpu_count()
    subset_size = max(1, len(donnees_exportees_transformed) // num_cores)
    subsets = [
        donnees_exportees_transformed.iloc[i:i + subset_size]
        for i in range(0, len(donnees_exportees_transformed), subset_size)
    ]

    if not subsets:
        print("No subsets created. Returning empty GeoDataFrame.")
        return GeoDataFrame()

    # Use multiprocessing to process the subsets
    with multiprocessing.Pool(num_cores) as pool:
        worker_func = partial(worker_function, conc_points=conc_points)
        grids = pool.map(worker_func, subsets)

    grids = [grid for grid in grids if not grid.empty]

    # If grids are empty after filtering, return an empty GeoDataFrame
    if not grids:
        print("No grids generated. Returning empty GeoDataFrame.")
        return GeoDataFrame()

    # Combine results into a single GeoDataFrame
    result = GeoDataFrame(pd.concat(grids, ignore_index=True))

    # Calculate intersection percentages finally
    result = calculate_intersection_percentages(result, donnees_exportees_transformed)

    return result
