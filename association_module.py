import geopandas as gpd
from shapely.geometry import Polygon
import pandas as pd
from geopandas import GeoDataFrame
from shapely.ops import unary_union
import multiprocessing
from functools import partial
import os
from shapely.geometry import box

def generate_points(polygone, conc_points):
    """
    Generate grid points for a given polygon and associate them with concentration points
    based on rectangle-polygon intersection instead of centroid inclusion.
    """
    if polygone.geometry.is_empty or not polygone.geometry.is_valid:
        return gpd.GeoDataFrame()
    if conc_points is None or conc_points.empty:
        return gpd.GeoDataFrame()

    polygon_geom = polygone.geometry
    # Build rectangles around all points
    rect_geoms = conc_points.geometry.apply(
        lambda pt: Polygon([
            (pt.x - 0.05, pt.y - 0.025),
            (pt.x + 0.05, pt.y - 0.025),
            (pt.x + 0.05, pt.y + 0.025),
            (pt.x - 0.05, pt.y + 0.025),
            (pt.x - 0.05, pt.y - 0.025)
        ])
    )
    rect_gdf = gpd.GeoDataFrame(
        conc_points.drop(columns="geometry").copy(),
        geometry=rect_geoms,
        crs=conc_points.crs
    )
    # Fast spatial filtering: keep rectangles that intersect the polygon (don't work on centroid)
    intersects_mask = rect_gdf.intersects(polygon_geom)
    points_within = conc_points.loc[intersects_mask].copy()
    points_within["iriscod"] = polygone["iriscod"]
    return points_within

import geopandas as gpd
from shapely.geometry import box
def calculate_intersection_percentages(
    grille_points: gpd.GeoDataFrame,
    iris_polygons: gpd.GeoDataFrame,
    cell_size_m: float = 1000.0,
    denominator: str = "cell",   # "cell" or "iris"
):
    if grille_points.crs != iris_polygons.crs:
        raise ValueError(f"CRS mismatch: grille={grille_points.crs}, iris={iris_polygons.crs}")
    if "iriscod" not in grille_points.columns:
        raise ValueError("grille_points must contain 'iriscod'")

    half = cell_size_m / 2.0

    g = grille_points.copy()
    iris = iris_polygons[["iriscod", "geometry"]].copy()

    # Build cell polygons around each point
    xs = g.geometry.x.to_numpy()
    ys = g.geometry.y.to_numpy()
    cell_geom = [box(x-half, y-half, x+half, y+half) for x, y in zip(xs, ys)]
    g["cell_geom"] = cell_geom

    # attach IRIS geometry
    g = g.merge(iris, on="iriscod", how="left", suffixes=("", "_iris"))

    # IMPORTANT: make them GeoSeries so .intersection works
    cell_gs = gpd.GeoSeries(g["cell_geom"], crs=g.crs)
    iris_gs = gpd.GeoSeries(g["geometry_iris"], crs=g.crs)

    inter_area = cell_gs.intersection(iris_gs).area

    if denominator == "cell":
        denom = float(cell_size_m * cell_size_m)
    elif denominator == "iris":
        denom = iris_gs.area
    else:
        raise ValueError("denominator must be 'cell' or 'iris'")

    g["perc"] = (inter_area / denom).astype(float)

    g = g[g["perc"] > 0].copy()
    g = g.drop(columns=["cell_geom", "geometry_iris"])
    return g

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
    result = result.to_crs(2154)
    donnees_exportees_transformed = donnees_exportees_transformed.to_crs(2154)

    # Calculate intersection percentages finally
    result = calculate_intersection_percentages(result, donnees_exportees_transformed, 1000, denominator="iris")
    return result