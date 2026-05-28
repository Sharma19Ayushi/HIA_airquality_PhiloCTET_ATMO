import pandas as pd
import geopandas as gpd
import numpy as np
from multiprocessing import Pool, cpu_count
from shapely.geometry import Point
from scipy.spatial import cKDTree
import logging

TARGET_CRS = "EPSG:2154"  # Lambert-93, meters — standard for France
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree
import logging

def expo_fast(
    donnees_exportees_transformed: gpd.GeoDataFrame,
    conc_corrigee: gpd.GeoDataFrame,
    grille_combinee: gpd.GeoDataFrame,
    max_distance_m: float | None = None,   # e.g. None (no filter) or 5000
    fallback_for_missing: bool = True,     # fill IRIS without grid using representative_point
):
    """
    Fast exposure computation:
      meanconc(iriscod)  = sum_i( conc(nearest(grid_i))  * w_i ) / sum_i(w_i)
      meandelta(iriscod) = sum_i( delta(nearest(grid_i)) * w_i ) / sum_i(w_i)
    where grid_i are points in grille_combinee (perc>0), w_i=perc.
    Returns DataFrame with columns: iriscod, meanconc, meandelta
    """
    logging.info("Starting expo_fast")

    # --- validation ---
    for col in ("geometry", "conc", "delta_conc"):
        if col not in conc_corrigee.columns:
            raise ValueError(f"Missing '{col}' in conc_corrigee")

    if not (donnees_exportees_transformed.crs == conc_corrigee.crs == grille_combinee.crs):
        raise ValueError(
            f"CRS mismatch: donnees={donnees_exportees_transformed.crs}, "
            f"conc={conc_corrigee.crs}, grille={grille_combinee.crs}"
        )

    # --- Build KDTree for conc points ---
    conc_xy = np.column_stack([
        conc_corrigee.geometry.x.to_numpy(),
        conc_corrigee.geometry.y.to_numpy()
    ])
    tree = cKDTree(conc_xy)

    # --- Use only grid points with positive weights ---
    g = grille_combinee.loc[grille_combinee["perc"] > 0, ["iriscod", "perc", "geometry"]].copy()
    if g.empty:
        logging.warning("No grid points with perc>0. Returning empty result.")
        return pd.DataFrame(columns=["iriscod", "meanconc", "meandelta"])

    grid_xy = np.column_stack([g.geometry.x.to_numpy(), g.geometry.y.to_numpy()])

    # nearest conc for all grid points
    distances, idx = tree.query(grid_xy, k=1)

    if max_distance_m is not None:
        valid = distances <= max_distance_m
        g = g.loc[valid].copy()
        idx = idx[valid]
        if g.empty:
            logging.warning("All grid points filtered out by max_distance_m.")
            return pd.DataFrame(columns=["iriscod", "meanconc", "meandelta"])

    # attach conc/delta to each grid point
    g["conc"] = conc_corrigee["conc"].to_numpy()[idx]
    g["delta_conc"] = conc_corrigee["delta_conc"].to_numpy()[idx]
    w = g["perc"].to_numpy(dtype=float)
    w = np.clip(w, 0, None)
    g["w"] = w

    # weighted sums by IRIS (vectorized)
    grouped = g.groupby("iriscod", as_index=False).agg(
        wsum=("w", "sum"),
        conc_w=("conc", lambda s: 0.0),  # placeholder; we'll compute below faster
        delta_w=("delta_conc", lambda s: 0.0),
    )

    # Faster than lambda: compute weighted sums using temporary columns
    g["conc_w"] = g["conc"] * g["w"]
    g["delta_w"] = g["delta_conc"] * g["w"]

    sums = g.groupby("iriscod", as_index=False)[["w", "conc_w", "delta_w"]].sum()
    sums = sums.rename(columns={"w": "wsum"})

    sums["meanconc"] = sums["conc_w"] / sums["wsum"]
    sums["meandelta"] = sums["delta_w"] / sums["wsum"]

    out = sums[["iriscod", "meanconc", "meandelta"]].copy()

    # --- ensure ALL IRIS are present (left join) ---
    all_iris = donnees_exportees_transformed[["iriscod", "geometry"]].drop_duplicates().copy()
    out = all_iris.merge(out, on="iriscod", how="left").drop(columns="geometry")

    # --- fallback for missing IRIS (no grid points) ---
    if fallback_for_missing:
        missing = out["meanconc"].isna()
        if missing.any():
            miss_iris = donnees_exportees_transformed.loc[
                donnees_exportees_transformed["iriscod"].isin(out.loc[missing, "iriscod"]),
                ["iriscod", "geometry"]
            ].drop_duplicates()

            rp = miss_iris.geometry.representative_point()
            miss_xy = np.column_stack([rp.x.to_numpy(), rp.y.to_numpy()])
            _, miss_idx = tree.query(miss_xy, k=1)

            fill_df = pd.DataFrame({
                "iriscod": miss_iris["iriscod"].to_numpy(),
                "meanconc": conc_corrigee["conc"].to_numpy()[miss_idx],
                "meandelta": conc_corrigee["delta_conc"].to_numpy()[miss_idx],
            })

            out = out.drop(columns=["meanconc", "meandelta"]).merge(fill_df, on="iriscod", how="left")

    logging.info("expo_fast completed")
    return out

def process_expo_subset(args):
    """
    Optimized subset processing for computing meanconc and meandelta using spatial data.
    Fixed: no rows dropped, threshold computed once outside loop.
    """
    donnees, conc_corrigee, grille_combinee, tree = args
    logging.info(f"Processing {len(donnees)} rows in subset")

    donnees = donnees.copy()
    donnees["meanconc"]  = 0.0   # initialise ALL rows to default upfront
    donnees["meandelta"] = 0.0   # → no rows will be lost regardless of branch taken

    # Group grille_combinee by iriscod once
    grille_grouped = grille_combinee[grille_combinee["perc"] > 0].groupby("iriscod")

    # ── Compute distance threshold ONCE outside the loop ─────────────────────
    grille_coords = np.array(
        list(zip(grille_combinee.geometry.x, grille_combinee.geometry.y))
    )
    distances_global, _ = tree.query(grille_coords, k=1)
    max_distance_threshold = min(np.percentile(distances_global, 95), 1000)

    # ── Main loop ─────────────────────────────────────────────────────────────
    for idx, row in donnees.iterrows():
        iriscode = row["iriscod"]

        if iriscode not in grille_grouped.groups:
            # default 0 already set above — just continue
            continue

        points       = grille_grouped.get_group(iriscode)
        point_coords = np.array(list(zip(points.geometry.x, points.geometry.y)))
        distances, indices = tree.query(point_coords, k=1)

        valid_mask = distances <= max_distance_threshold

        if not valid_mask.any():
            continue   # default 0 already set

        indices        = indices[valid_mask]
        nearest_points = conc_corrigee.iloc[indices]
        conc_values    = nearest_points["conc"].values
        delta_values   = nearest_points["delta_conc"].values
        weights        = points.iloc[valid_mask]["perc"].values
        weights        = np.clip(weights, 0, 1)
        weight_sum     = np.sum(weights)

        if weight_sum > 0:
            weights    = weights / weight_sum
            donnees.at[idx, "meanconc"]  = np.dot(conc_values,  weights)
            donnees.at[idx, "meandelta"] = np.dot(delta_values, weights)
        # else: weight_sum == 0 → default 0 already set, nothing to do

    # Return ALL rows (no donnees_indices filtering needed anymore)
    return donnees[["iriscod", "meanconc", "meandelta"]].reset_index(drop=True)

