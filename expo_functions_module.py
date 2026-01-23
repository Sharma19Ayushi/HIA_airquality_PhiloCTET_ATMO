import pandas as pd
import geopandas as gpd
import numpy as np
from multiprocessing import Pool, cpu_count
from shapely.geometry import Point
from scipy.spatial import cKDTree
import logging


def expo(donnees_exportees_transformed, conc_corrigee, grille_combinee):
    """
    Parallelized and optimized function for processing exposure metrics.
    """
    logging.info("Starting optimized expo function")

    try:
        # Validate input data
        required_columns = ['geometry', 'conc', 'delta_conc']
        for col in required_columns:
            if col not in conc_corrigee.columns:
                raise ValueError(f"Missing column '{col}' in conc_corrigee GeoDataFrame")

        # Build spatial index (cKDTree) once for efficiency
        coords = np.array(list(zip(conc_corrigee.geometry.x, conc_corrigee.geometry.y)))
        tree = cKDTree(coords)

        # Limit number of cores for multiprocessing
        num_cores = min(cpu_count() - 1, 4)  # Max 4 cores for stability
        subsets = np.array_split(donnees_exportees_transformed, num_cores)

        # Prepare arguments for parallel processing
        args = [(subset, conc_corrigee, grille_combinee, tree) for subset in subsets]

        # Use multiprocessing pool to process subsets
        with Pool(num_cores) as pool:
            results = pool.map(process_expo_subset, args)

        # Combine processed results into a single DataFrame
        processed_results = pd.concat(results, ignore_index=True)

        logging.info("Expo processing completed successfully")
        return processed_results

    except Exception as e:
        logging.error(f"Error in expo function: {e}")
        return pd.DataFrame()  # Ensure a DataFrame, even if empty, is returned


def process_expo_subset(args):
    """
    Optimized subset processing for computing meanconc and meandelta using spatial data.
    """
    donnees, conc_corrigee, grille_combinee, tree = args
    logging.info(f"Processing {len(donnees)} rows in subset")

    # Work on a copy of the data
    donnees = donnees.copy()

    # Group grille_combinee by 'iriscod' (optimize filtering)
    grille_grouped = grille_combinee[grille_combinee['perc'] > 0].groupby('iriscod')

    # Initialize results storage
    meanconc_results = []
    meandelta_results = []
    donnees_indices = []

    # Default values if no valid neighbors are found
    default_meanconc = 0  # Default mean concentration
    default_meandelta = 0  # Default mean delta concentration

    # Process in bulk # Go through all rows in donnees
    for idx, row in donnees.iterrows():
        iriscode = row['iriscod']
        if iriscode not in grille_grouped.groups:
            # No weights for this code; assign default
            donnees.at[idx, 'meanconc'] = default_meanconc
            donnees.at[idx, 'meandelta'] = default_meandelta
            continue

        points = grille_grouped.get_group(iriscode)
        point_coords = np.array(list(zip(points.geometry.x, points.geometry.y)))
        distances, indices = tree.query(point_coords, k=1)

        distances_test, _ = tree.query(np.array(list(zip(grille_combinee.geometry.x, grille_combinee.geometry.y))), k=1)
        max_distance_threshold = min(np.percentile(distances_test, 95), 1000)
        valid_mask = distances <= max_distance_threshold

        if not valid_mask.any():
            donnees.at[idx, 'meanconc'] = default_meanconc
            donnees.at[idx, 'meandelta'] = default_meandelta
            continue

        indices = indices[valid_mask]
        nearest_points = conc_corrigee.iloc[indices]
        conc_values = nearest_points['conc'].values
        delta_values = nearest_points['delta_conc'].values
        weights = points.iloc[valid_mask]['perc'].values
        weights = np.clip(weights, 0, 1)
        weight_sum = np.sum(weights)
        if weight_sum > 0:
            weights = weights / weight_sum
            meanconc = np.dot(conc_values, weights)
            meandelta = np.dot(delta_values, weights)
        else:
            meanconc = default_meanconc
            meandelta = default_meandelta
        donnees.at[idx, 'meanconc'] = meanconc
        donnees.at[idx, 'meandelta'] = meandelta
        donnees_indices.append(idx)
    # Return the processed subset with computed meanconc and meandelta
    return donnees[['iriscod', 'meanconc', 'meandelta']].loc[donnees_indices].reset_index(drop=True)

