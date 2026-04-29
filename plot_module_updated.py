import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt

#Helper functions to plot maps for PM2.5 and NO2 concentrations
def pretty_pol(pol):
    return "PM2.5" if pol == "ug_PM25_RH50" else "NO2"

def _square_extent_from_bounds(bounds, pad_frac=0.02):
    xmin, ymin, xmax, ymax = bounds
    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    half = max(dx, dy) / 2.0
    half = half * (1.0 + pad_frac)
    return (cx - half, cx + half, cy - half, cy + half)

def _get_iris_code_col(gdf):
    for c in ["CODE_IRIS", "iriscod", "code_iris", "IRIS_CODE", "IRIS", "iris"]:
        if c in gdf.columns:
            return c
    return None

def _points_to_iris_polygons(points_gdf, value_col, iris_polys, iris_code_col, fill_missing=True):
    pts = points_gdf[[value_col, "geometry"]].copy()
    if pts.crs != iris_polys.crs:
        pts = pts.to_crs(iris_polys.crs)

    joined = gpd.sjoin(pts, iris_polys[[iris_code_col, "geometry"]], how="left", predicate="intersects")

    missing_pts = joined[iris_code_col].isna()
    if bool(missing_pts.any()):
        try:
            nearest = gpd.sjoin_nearest(
                pts.loc[missing_pts],
                iris_polys[[iris_code_col, "geometry"]],
                how="left",
                distance_col="__dist",
            )
            joined.loc[missing_pts, iris_code_col] = nearest[iris_code_col].values
        except Exception:
            pass

    agg = (
        joined.dropna(subset=[iris_code_col])
        .groupby(iris_code_col)[value_col]
        .mean()
        .reset_index()
    )

    out = iris_polys.merge(agg, on=iris_code_col, how="left")

    if fill_missing and out[value_col].isna().any():
        s = pd.to_numeric(out[value_col], errors="coerce").astype(float)
        if s.notna().any():
            out[value_col] = s.fillna(s.mean()).values

    return out

def get_pw_mean_commune(iris_vals_gdf, val_col, pop_commune):
    """
    Population-weighted mean by commune:
    - mean(val_col) across IRIS within each comcod
    - weight commune means by commune population
    """
    commune_vals = iris_vals_gdf.groupby("comcod")[val_col].mean().reset_index()
    merged = pd.merge(commune_vals, pop_commune, on="comcod", how="inner")
    total_pop = merged["pop_val"].sum()
    if total_pop == 0:
        return 0.0
    return float((merged[val_col] * merged["pop_val"]).sum() / total_pop)

# -------------------------
# map plot function (polygons)
# -------------------------
def fast_plot(ax, poly_gdf, col, title, cmap, norm, extent, pw_mean=None, unit_label="µg/m³"):
    poly_gdf.plot(
        column=col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        linewidth=0.0,
        edgecolor="none",
        missing_kwds={"color": "lightgrey", "edgecolor": "none"},
        rasterized=True,
    )

    title_text = f"{title}"
    if pw_mean is not None:
        title_text += f"\n$\\mathbf{{Population\\ Weighted\\ Mean\\ =\\ {pw_mean:.2f}\\ {unit_label}}}$"
    ax.set_title(title_text, pad=10)
    xmin, xmax, ymin, ymax = extent
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.tick_params(axis="both", which="major", length=4, width=0.8, direction="out")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sm

#Plotting map codes for PM2.5 and NO2 concentrations
import logging
from cordo_chimere_module import *
plt.rcParams.update({
    "figure.dpi": 300,
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11
})

scenario_colors = {
    "S1": "#6C8EBF",
    "S2": "#B47CC7",
    "S3": "#6ABF69",
    "S4": "#E88D67",
}

scenario_dark_colors = {
    "S1": "#355C9A",
    "S2": "#8E5EA2",
    "S3": "#3F8F44",
    "S4": "#C96A3D",
}
baseline_year = "2019"
baseline_label = "Baseline (2019)"
scenarios = ["s1", "s2", "s3", "s4"]
pollutants = ["ug_PM25_RH50", "ug_NO2"]
years = ["2030", "2050"]
def _safe_numeric_series(series):
    vals = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return vals.astype(float).to_numpy()
def load_chimere_absolute_data(pol):
    """
    Returns a tidy dataframe with absolute concentrations for:
    - Baseline 2019 (from coordo_ineris_chimere, column conc19)
    - Scenarios 2030/2050 (from coordo_chimere, column conc)
    """
    records = []
    pol_lbl = pretty_pol(pol)

    logging.info(f"Loading baseline data for {pol_lbl}")
    try:
        base_pts = coordo_ineris_chimere(pol, year=baseline_year)
        if "conc19" not in base_pts.columns:
            logging.warning(f"'conc19' column not found in baseline data for {pol}")
        else:
            base_vals = _safe_numeric_series(base_pts["conc19"])
            for val in base_vals:
                records.append({
                    "Pollutant": pol_lbl,
                    "Group": baseline_label,
                    "Scenario": "Baseline",
                    "Year": baseline_year,
                    "Concentration": float(val),
                })
    except Exception as e:
        logging.warning(f"Could not load baseline data for {pol}: {e}")

    for year in years:
        for sc in scenarios:
            label = f"{sc.upper()} - {year}"
            logging.info(f"Loading {label} for {pol_lbl}")
            try:
                scen_pts = coordo_chimere(pol, year=year, SC=sc.upper())
                if "conc" not in scen_pts.columns:
                    logging.warning(f"'conc' column not found for {label} / {pol}")
                    continue

                scen_vals = _safe_numeric_series(scen_pts["conc"])
                for val in scen_vals:
                    records.append({
                        "Pollutant": pol_lbl,
                        "Group": label,
                        "Scenario": sc.upper(),
                        "Year": year,
                        "Concentration": float(val),
                    })
            except Exception as e:
                logging.warning(f"Could not load {label} for {pol}: {e}")

    return pd.DataFrame.from_records(records)
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import logging

# Helper: pretty pollutant label
def pretty_pol(pol):
    return "PM2.5" if pol == "ug_PM25_RH50" else "NO2"

def _square_extent_from_bounds(bounds, pad_frac=0.02):
    xmin, ymin, xmax, ymax = bounds
    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    half = max(dx, dy) / 2.0 * (1.0 + pad_frac)
    return (cx - half, cx + half, cy - half, cy + half)

def fast_plot(ax, poly_gdf, col, cmap, norm, extent):
    poly_gdf.plot(
        column=col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        linewidth=0,
        edgecolor="none",
        missing_kwds={"color": "lightgrey", "edgecolor": "none"},
        rasterized=True,
    )
    xmin, xmax, ymin, ymax = extent
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

def plot_maps_matrix(
    data,     # dict keyed as (year, sc, pol) → dict with keys "scen" (gdf), "pw_mean_conc" (float)
    years,    # list of years (row order)
    scenarios, pollutant_labels, # e.g. ["s2", "s3"], ["ug_PM25_RH50", "ug_NO2"]
    extent,   # map extent
    vmin, vmax,
    output_path
):
    nrows = len(years)
    ncols = len(scenarios) * len(pollutant_labels)
    pol_lbls = [pretty_pol(pol) for pol in pollutant_labels]
    combo_order = [(sc, pol) for pol in pollutant_labels for sc in scenarios]

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(2.5*ncols, 2.9*nrows),
        subplot_kw=dict(aspect="equal")
    )
    if nrows == 1:
        axes = axes[None, :]
    if ncols == 1:
        axes = axes[:, None]

    plt.subplots_adjust(
        left=0.06, right=0.91, bottom=0.07, top=0.92,
        wspace=0.07, hspace=0.05
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = "viridis"

    # Plot each panel; add a text for pop-w mean
    for i, year in enumerate(years):
        for j, (pol, sc) in enumerate(combo_order):
            ax = axes[i, j]
            entry = data[(year, sc, pol)]
            scen = entry["scen"]
            fast_plot(ax, scen, "conc", cmap, norm, extent)
            popw = entry["pw_mean_conc"]
            # (Add pop-weighted mean conc above map, bold/centered/small)
            ax.text(
                0.5, 1.05,
                f"Pop-weighted mean: {popw:.2f} µg/m³",
                ha="center", va="bottom",
                fontsize=10, fontweight="bold",
                transform=ax.transAxes,
            )

    # Set only unique year labels (row) and unique column labels (top)
    for i, year in enumerate(years):
        axes[i, 0].set_ylabel(str(year), fontsize=13, fontweight="bold", labelpad=18)

    for j, (pol, sc) in enumerate(combo_order):
        axes[0, j].set_title(f"{pretty_pol(pol)}\n{sc.upper()}", fontsize=12, fontweight="bold", pad=12)

    # Unified colorbar (right)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(
        sm, ax=axes, orientation="vertical",
        fraction=0.03, pad=0.02, aspect=34
    )
    cbar.set_label("Mean Concentration (µg/m³)", fontsize=13)
    cbar.ax.tick_params(labelsize=10)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

