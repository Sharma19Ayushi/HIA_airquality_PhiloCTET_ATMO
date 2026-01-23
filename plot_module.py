
import os
import seaborn as sns
from shapely.geometry import Point
import numpy as np
from shapely.ops import transform
from functools import partial
import pyproj
import logging
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import geopandas as gpd
import contextily as ctx
from joypy import joyplot

# Function to read shapefiles
def read_shapefile(path, title):
    return gpd.read_file(os.path.join(path, f"{title}.shp"))

# Function to align CRS of spatial objects
def align_crs(data, target_crs):
    if data.crs != target_crs:
        data = data.to_crs(target_crs)
    return data

# Function to export data to shapefile
def export_data_shp(data, path, title):
    shp_path = os.path.join(path, f"{title}.shp")
    data.to_file(shp_path, driver='ESRI Shapefile')

# Function to plot the number of IRIS polygons based on the characteristic distance
def plot_distance(donnees_exportees):
    donnees_exportees['area'] = donnees_exportees.geometry.area
    donnees_exportees['sqrt_area'] = donnees_exportees['area'].apply(np.sqrt)
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    
    sns.histplot(donnees_exportees['sqrt_area'], bins=30, ax=ax[0], color='blue')
    ax[0].axvline(1200, color='red', linestyle='--', linewidth=1)
    ax[0].set_title('Number of IRIS polygons based on the square root of the polygon area')
    ax[0].set_xlabel('Square root of the polygon area (m)')
    ax[0].set_ylabel('Number of IRIS polygons')
    
    donnees_exportees.plot(column='sqrt_area', ax=ax[1], legend=True, cmap='viridis')
    ctx.add_basemap(ax[1], crs=donnees_exportees.crs.to_string())
    ax[1].set_title('Map of IRIS polygons colored by the square root of the polygon area')
    
    plt.tight_layout()
    plt.show()

# Function to plot population map of over 30s for a year
def plot_carte_iris(donnees_mixtes, annee):
    pop = f"pop{annee}"
    donnees_mixtes.plot(column=pop, legend=True)
    plt.title(f"Population map of over 30s for {annee}")
    plt.show()

# Function to save population map of over 30s for a year
def save_carte_iris(donnees_mixtes, annee, path):
    pop = f"pop{annee}"
    donnees_mixtes.plot(column=pop, legend=True)
    plt.title(f"Population map of over 30s for {annee}")
    plt.savefig(path)
    plt.close()

# Function to plot exposure maps using concentration data
def plot_carte_expo(result, col, n):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    # Plot exposure data from the GeoDataFrame
    result.plot(column=col, cmap="viridis", linewidth=0, ax=ax, legend=True,
                norm=Normalize(vmin=0, vmax=n))
    # Add gridlines and labels for lat/lon
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color="gray")  # Enable grid
    ax.set_axisbelow(True)  # Place grid below the polygons
    ax.tick_params(labelsize=12)  # Adjust tick label size
    ax.set_xlabel("Longitude", fontsize=12)  # X-axis label
    ax.set_ylabel("Latitude", fontsize=12)  # Y-axis label
    # Set the title
    ax.set_title("Exposure in μg/m³", fontsize=14)
    # Display the plot
    plt.tight_layout()
    plt.show()

# Function to save exposure map
def save_carte_expo(result, path, col, n):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    # Plot exposure data from the GeoDataFrame
    result.plot(column=col, cmap="viridis", linewidth=0, ax=ax, legend=True,
                norm=Normalize(vmin=0, vmax=n))
    # Add gridlines and labels for lat/lon
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, color="gray")  # Enable grid
    ax.set_axisbelow(True)  # Place grid below the polygons
    ax.tick_params(labelsize=12)  # Adjust tick label size
    ax.set_xlabel("Longitude", fontsize=12)  # X-axis label
    ax.set_ylabel("Latitude", fontsize=12)  # Y-axis label
    # Set the title
    ax.set_title("Exposure in μg/m³", fontsize=14)
    # Save the plot to the specified path
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close(fig)


# Function to get scale value n1 based on pollutant
def echelle_n1(pol):
    if pol in ["ug_PM25_RH50_high", "ug_PM25_RH50_low"]:
        pol = "ug_PM25_RH50"
    elif pol in ["ug_NO2_high", "ug_NO2_low"]:
        pol = "ug_NO2"
    return 13 if pol == "ug_PM25_RH50" else 34
#return 13 if pol == "ug_PM25_RH50" else 34

# Function to get scale value n2 based on pollutant
def echelle_n2(pol):
    if pol in ["ug_PM25_RH50_high", "ug_PM25_RH50_low"]:
        pol = "ug_PM25_RH50"
    elif pol in ["ug_NO2_high", "ug_NO2_low"]:
        pol = "ug_NO2"
    return 9 if pol == "ug_PM25_RH50" else 25
#return 9 if pol == "ug_PM25_RH50" else 25

# Function to create map
def create_map(data, scenario, variable, n, titles_colors):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    data.plot(column=variable, ax=ax, legend=True, cmap='viridis')
    ax.set_title(titles_colors[scenario]["title"])
    plt.savefig(f"{scenario}_{variable}.png")  # Save map instead of displaying

print('Successfully loaded plotting command')

#Plot the number of IRIS mortality maps based on mortality results
def plot_multiple_iris_maps(
        iris_level_results, iris_shapefile, columns_to_plot, titles=None,
        cmap="RdYlBu_r", base_path="data/2-output-data"
):
    """
    Plot multiple choropleth maps of IRIS-level mortality results with improved aesthetics.
    Parameters:
    - iris_level_results: DataFrame with mortality analysis results
    - iris_shapefile: GeoDataFrame with 'iriscod' and 'geometry'
    - columns_to_plot: List of columns to visualize
    - titles: List of titles for the plots
    - cmap: Colormap for all plots
    - base_path: Path where the plot will be saved
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import os

    # Merge result with shapefile
    plot_data = iris_shapefile.merge(iris_level_results, on="iriscod", how="inner")
    plot_data = plot_data.fillna(0)

    # Set up the figure with adjusted dimensions
    n_maps = len(columns_to_plot)
    fig, axes = plt.subplots(1, n_maps, figsize=(20, 8))

    # Ensure axes is always an array
    if n_maps == 1:
        axes = [axes]

    if titles is None:
        titles = columns_to_plot

    # Define consistent style parameters
    TITLE_PARAMS = {'fontsize': 12, 'fontweight': 'bold', 'pad': 15}
    ANNOTATION_PARAMS = {'fontsize': 11, 'fontweight': 'bold', 'color': 'darkblue'}
    LEGEND_PARAMS = {'bbox_to_anchor': (1.05, 0.5), 'loc': 'center left'}
    MAP_PARAMS = {
        'linewidth': 0.2,
        'edgecolor': '0.8',
        'missing_kwds': {'color': 'lightgrey'},
    }

    for i, column in enumerate(columns_to_plot):
        if pd.api.types.is_numeric_dtype(plot_data[column]):
            # Calculate value range for consistent color scaling
            vmin, vmax = plot_data[column].quantile([0.01, 0.99])

            # Plot the map
            im = plot_data.plot(
                column=column,
                cmap=cmap,
                legend=True,
                ax=axes[i],
                vmin=vmin,
                vmax=vmax,
                legend_kwds={
                    'label': column,
                    'orientation': 'vertical',
                    'shrink': 0.8,
                    'fraction': 0.046,
                    'aspect': 20
                },
                **MAP_PARAMS
            )

            # Adjust legend
            legend = axes[i].get_legend()
            if legend:
                legend.set_bbox_to_anchor((1.05, 0.5))
                legend.set_frame_on(True)
                legend.set_framealpha(0.8)
                legend.set_edgecolor('0.8')

            # Calculate and format annotation
            is_total = column.startswith(("total", "dlg", "ci"))
            value = plot_data[column].sum() if is_total else plot_data[column].sum()
            label = "Total" if is_total else "Sum"
            format_str = "{:,.0f}" if is_total else "{:,.2f}"
            annotation_text = f"{label}: {format_str.format(value)}"

            # Add title and annotation
            axes[i].set_title(titles[i], **TITLE_PARAMS)
            axes[i].text(
                0.5, -0.1,
                annotation_text,
                transform=axes[i].transAxes,
                ha='center',
                va='center',
                **ANNOTATION_PARAMS
            )
        else:
            axes[i].set_title(f"{titles[i]}\n(Non-numeric data)", **TITLE_PARAMS)

        axes[i].axis('off')

    # Adjust layout
    plt.tight_layout(w_pad=3)

    # Save with high resolution
    output_path = os.path.join(base_path, "iris_mortality_maps_continuous.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.5)

    return fig, axes


# Main execution
if __name__ == "__main__":
    donnees_exportees = gpd.read_file("path_to_donnees_exportees_shapefile.shp")
    grille_combinee = gpd.read_file("path_to_grille_combinee_shapefile.shp")
    donnees_merged = pd.read_csv("path_to_donnees_merged.csv")

    # Ensure proper CRS alignment
    target_crs = "EPSG:2154"
    donnees_exportees = align_crs(donnees_exportees, target_crs)
    grille_combinee = align_crs(grille_combinee, target_crs)

    # Example call to plot_distance function
    plot_distance(donnees_exportees)

# Function to plot boxplots for comparison across all 4 scenarios for 2030, 2050 and Baseline
def plot_boxplot_comparison(data, title, output_path):
    plt.figure(figsize=(6, 12))  # Adjusted to fit all scenarios and years

    # Filter data for Baseline, 2030, and 2050
    data = data[(data["Year"].isin(["2019", "2030", "2050"]))].copy()

    # Replace Baseline year for better readability
    data["Year"] = data["Year"].replace("2019", "Baseline")

    # Combine legend/scenario information instead of x-axis labels
    data["Scenario_Legend"] = data["Scenario"].replace({
        "Baseline": "Baseline",
        "s1": "s1",
        "s2": "s2",
        "s3": "s3",
        "s4": "s4"
    })

    # Define consistent colors for the scenarios, including baseline
    base_colors = ["grey", "#5C9AD6", "#B47CC7", "#6ABF69", "#E48A4C"]
    palette = dict(zip(["Baseline", "s1", "s2", "s3", "s4"], base_colors))

    # Plot boxplot
    sns.boxplot(
        x="Year",
        y="Average_Concentration",
        hue="Scenario_Legend",
        data=data,
        palette=palette,
        width=0.8,
        showfliers=False  # Hide outliers beyond whiskers
    )

    # Add title and labels
    plt.title(title, pad=20)  # Add some padding for the title
    plt.xlabel("Year")
    plt.ylabel("Average Concentration (µg/m³)")

    # Move the legend outside the plot
    plt.legend(title="Scenario", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()


def plot_combined_ridgeplot(data, pollutant, title, output_path):
    plt.rcParams.update({'font.size': 14})

    # Convert "Year" and "Scenario" columns to strings before concatenation
    data["Year"] = data["Year"].astype(str)
    data["Scenario"] = data["Scenario"].astype(str)

    # Combine "Year" and "Scenario" into a string grouping category
    data["Year_Scenario"] = data["Year"] + " - " + data["Scenario"]

    # Create ridge plot (distributions grouped by Year and Scenario)
    fig, axes = joyplot(
        data=data,
        column="Average_Concentration",
        by="Year_Scenario",
        colormap=base_colors,  # General colormap
        kind="kde",
        xlabels=True,
        ylabels=True,
        legend=False,  # Hide legend as categories are grouped in labels
        linewidth=1.0,
    )

    # Add title and axis labels
    fig.suptitle(title, fontsize=14)
    plt.xlabel("Average Concentration (µg/m³)")
    plt.ylabel("Year - Scenario")
    # Adjust figure size to be tall and narrow
    fig.set_size_inches(6, 10)  # Width = 4 inches, Height = 10 inches (adjust as needed)
    # Save the plot and display
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    
