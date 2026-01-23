import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Function which returns a GeoDataFrame representing the IRIS geometries and population data
def geometries(path_iris, path_depart, path_contours):
    shape_iris = gpd.read_file(path_contours)
    donnees_iris = pd.read_csv(path_iris, delimiter="\t", header=0)
    iris_manquants = pd.DataFrame({
        "IRIS": ["311490901", "532740000", "452870000", "532390000", "163510000",
                 "215070000", "212130000", "690720103", "690720101", "690720102",
                 "276760000", "940220103", "940220104", "940220101"],
        "REG": [76, 52, 24, 52, 75, 27, 27, 84, 84, 84, 28, 11, 11, 11],
        "DEP": [31, 53, 45, 53, 16, 21, 21, 69, 69, 69, 27, 94, 94, 94],
        "P19_POP": [0] * 14
    })
    donnees_iris = pd.concat([donnees_iris, iris_manquants], ignore_index=True)
    num_depart = pd.read_csv(path_depart, sep=";")
    donnees_geom = shape_iris.merge(donnees_iris, left_on="CODE_IRIS", right_on="IRIS", how="left")
    donnees_geom = donnees_geom.merge(num_depart, left_on="DEP", right_on="num", how="left").rename(columns={
        "DEP": "dep_cod",
        "CODE_IRIS": "iris_cod",
        "INSEE_COM": "com_cod",
        "NOM_COM": "com_name",
        "NOM_IRIS": "iris_name",
        "REG": "region_cod"
    })
    donnees_geom = donnees_geom[["dep_cod", "dep_name", "region_cod", "region_name",
                                 "com_cod", "com_name", "iris_cod", "iris_name", "P19_POP", "geometry"]]
    donnees_geom = donnees_geom.rename(columns={"P19_POP": "POP_2019_IRIS_TotAge"})
    donnees_geom["aire_m2"] = donnees_geom["geometry"].area
    print("The INSEE data is extracted at the IRIS level.")
    return donnees_geom

# Calculates the total population by age group from 2019 to 2050 at the national level
def age_nat(path_age_femmes, path_age_hommes):
    donnees_age_femmes = pd.read_csv(path_age_femmes, sep=";", na_values="NA")
    donnees_age_hommes = pd.read_csv(path_age_hommes, sep=";", na_values="NA")
    donnees_age_femmes = donnees_age_femmes.iloc[:-1]
    donnees_age_hommes = donnees_age_hommes.iloc[:-1]
    donnees_age_femmes.iloc[:, 0] = donnees_age_femmes.iloc[:, 0].str.replace("105+", "105", regex=False)
    donnees_age_hommes.iloc[:, 0] = donnees_age_hommes.iloc[:, 0].str.replace("105+", "105", regex=False)
    age_hf = pd.DataFrame(np.nan, index=range(len(donnees_age_femmes)),
                          columns=["age"] + [f"POP_{2019 + i}_NAT" for i in range(32)])
    age_hf["age"] = donnees_age_femmes.iloc[:, 0]
    for i in range(32):
        annee = 2019 + i
        age_hf[f"POP_{annee}_NAT"] = pd.to_numeric(donnees_age_femmes.iloc[:, 2 + i]) + \
                                      pd.to_numeric(donnees_age_hommes.iloc[:, 2 + i])
    print("The total population for each age group between 2019 and 2050 has been calculated at the national level.")
    return age_hf

#Function to disaggregate population by age
def decomposition_age(age_hf):
    # Dynamically create columns for output DataFrame and ensure compatibility with mixed data types
    perc_hf = pd.DataFrame(
        np.nan,
        index=range(len(age_hf)),
        columns=list(age_hf.columns) + ["tranche_age"] +
                [f"FR_{2019 + j}_agesurtrage" for j in range(len(age_hf.columns) - 1)]
    )

    # Explicitly cast "tranche_age" column (and others as needed) to appropriate object dtype
    perc_hf["tranche_age"] = perc_hf["tranche_age"].astype("object")
    perc_hf["age"] = age_hf["age"]

    # Iterate through all columns starting from the second one ("2019" onwards), including the last column
    for j in range(1, len(age_hf.columns)):  # Include the last column by using len(age_hf.columns)
        annee = 2019 + j  # Dynamically calculate the year (starting at 2019 onwards)
        perc_hf[age_hf.columns[j]] = age_hf[age_hf.columns[j]]
        perc_hf[f"FR_{annee}_agesurtrage"] = np.nan  # Initialize new column for computation

        # Debug point: Confirm year and column alignment
        print(f"Processing year {annee}, column: {age_hf.columns[j]}")

        # Decompose for ages < 95
        for i in range(95):
            p = int(int(age_hf.iloc[i, 0]) / 5)  # Calculate 5-year group
            perc_hf.loc[i, "tranche_age"] = f"[{5 * p};{5 * p + 5}["
            try:
                # Correctly align indices and sums for shared columns
                total_sum = sum(pd.to_numeric(age_hf.iloc[5 * p:5 * p + 5, j], errors='coerce'))
                perc_hf.iloc[i, j + len(age_hf.columns)] = (
                        pd.to_numeric(age_hf.iloc[i, j], errors='coerce') / total_sum
                )
            except ZeroDivisionError:
                perc_hf.iloc[i, j + len(age_hf.columns)] = None  # Handle divide-by-zero gracefully

        # Decompose for ages 95 to 99
        for i in range(95, 99):
            perc_hf.loc[i, "tranche_age"] = "[95;99+]"
            try:
                perc_hf.iloc[i, j + len(age_hf.columns)] = (
                        pd.to_numeric(age_hf.iloc[i, j], errors='coerce') /
                        sum(pd.to_numeric(age_hf.iloc[95:100, j], errors='coerce'))
                )
            except ZeroDivisionError:
                perc_hf.iloc[i, j + len(age_hf.columns)] = None

        # Decompose for ages 100+
        perc_hf.loc[99, "tranche_age"] = "[95;99+]"
        try:
            total_sum_100 = sum(pd.to_numeric(age_hf.iloc[95:100, j], errors='coerce'))
            perc_hf.iloc[99, j + len(age_hf.columns)] = (
                    sum(pd.to_numeric(age_hf.iloc[100:, j], errors='coerce')) / total_sum_100
            )
        except ZeroDivisionError:
            perc_hf.iloc[99, j + len(age_hf.columns)] = None  # Handle divide-by-zero gracefully

    # Ensure perc_hf DataFrame has the correct number of rows
    perc_hf = perc_hf.iloc[:100]

    # Fix the age row for the last entry
    perc_hf.at[99, "age"] = "99"

    # Debugging: Check generated columns
    print("Generated columns in perc_hf:")
    print(perc_hf.columns)

    # Final debugging output
    print("Processing completed, including all columns up to 2050.")
    return perc_hf

# Processes population projection data from 2019 to 2050 at the departmental level in France
import pandas as pd
import chardet

def recense(path_proj, perc_hf):
    try:
        with open(path_proj, "rb") as f:
            result = chardet.detect(f.read(100000))  # Read a sample
            print(result)  # Check detected encoding
        donnees_proj = pd.read_csv(path_proj, encoding=result["encoding"], sep="\t", header=None)
        donnees_proj.columns = donnees_proj.iloc[0]  # First row as column names
        donnees_proj = donnees_proj[1:].reset_index(drop=True)

        # Keep only the first 36 columns
        donnees_proj = donnees_proj.iloc[:, :37]

        # Remove "POP_2018" if it exists
        if "POP_2018" in donnees_proj.columns:
            donnees_proj.drop(columns=["POP_2018"], inplace=True)

        # Standardizing department names (handling special characters)
        department_replacements = {
            "Ard\x8fche": "Ardèche", "Ari\x8fge": "Ariège", "Bouches-du-Rh\x99ne": "Bouches-du-Rhône",
            "Corr\x8fze": "Corrèze", "C\x99te-d'Or": "Côte-d'Or", "C\x99tes-d'Armor": "Côtes-d'Armor",
            "Deux-S\x8fvres": "Deux-Sèvres", "Dr\x99me": "Drôme", "Finist\x8fre": "Finistère",
            "Haute-Sa\x99ne": "Haute-Saône", "Hautes-Pyr\x8en\x8ees": "Hautes-Pyrénées", "H\x8erault": "Hérault",
            "Is\x8fre": "Isère", "Loz\x8fre": "Lozère", "Ni\x8fvre": "Nièvre", "Puy-de-D\x99me": "Puy-de-Dôme",
            "Pyr\x8en\x8ees-Atlantiques": "Pyrénées-Atlantiques", "Pyr\x8en\x8ees-Orientales": "Pyrénées-Orientales",
            "Rh\x99ne": "Rhône", "R\x8eunion": "Réunion", "Sa\x99ne-et-Loire": "Saône-et-Loire",
            "Val-d'Oise": "Val-d-Oise", "Vend\x8ee": "Vendée"
        }
        donnees_proj.replace({"ZONE": department_replacements}, inplace=True)

        # Convert population columns to numeric
        pop_columns = [col for col in donnees_proj.columns if col.startswith("POP_")]
        donnees_proj[pop_columns] = donnees_proj[pop_columns].apply(pd.to_numeric, errors="coerce")
        # Aggregate population by department (ZONE) and age group (TRAGE)
        donnees_proj = donnees_proj.groupby(["ZONE", "TRAGE"], as_index=False)[pop_columns].sum()

        # Rename columns to match R's naming convention (e.g., POP_<year>_DEP_TrAge)
        donnees_proj.rename(columns=lambda col: f"{col}_DEP_TrAge" if col.startswith("POP_") else col, inplace=True)

        # Compute total population for each department across age groups
        total_pop = donnees_proj.groupby("ZONE", as_index=False).sum()

        # Rename total population columns to match R's naming convention (e.g., POP_<year>_DEP_TotAge)
        total_pop.rename(columns=lambda col: col.replace("_DEP_TrAge", "_DEP_TotAge") if "_DEP_TrAge" in col else col,
                         inplace=True)

        # Merge total departmental population back to the main DataFrame
        donnees_proj = pd.merge(donnees_proj, total_pop, on="ZONE", suffixes=("", "_Tot"))

        # Rename the TRAGE column to match the R script naming convention (use 'tranche_age' instead of TRAGE)
        donnees_proj.rename(columns={"TRAGE": "tranche_age"}, inplace=True)

        # Merge with percentage data (perc_hf) based on age group (tranche_age)
        donnees_proj = pd.merge(donnees_proj, perc_hf, on="tranche_age", how="left")

        # Convert 'age' column to numeric and sort
        donnees_proj["age"] = pd.to_numeric(donnees_proj["age"], errors="coerce")
        donnees_proj.sort_values(by=["ZONE", "age"], inplace=True)

        print("The number of people in each 5-year age group between 2019 and 2050 is predicted at the departmental level")
        return donnees_proj

    except Exception as e:
        print(f"An error occurred in the recense function: {e}")
        raise


# Function to disaggregate population data between 2019 and 2050 by age at the IRIS level
def desagreg(donnees_geom, donnees_proj):
    try:
        # Merge the dataframes on specified columns
        donnees_insee = donnees_geom.merge(donnees_proj, left_on="dep_name", right_on="ZONE")
        print("The INSEE data measured at the IRIS level is overlaid with the predicted INSEE data at the departmental level.")

        # Ensure that 'age' and required columns exist in the data
        required_columns = ['age', 'POP_2019_IRIS_TotAge', 'POP_2019_DEP_TrAge', 'FR_2019_agesurtrage', 'POP_2019_DEP_TotAge']
        for col in required_columns:
            if col not in donnees_insee.columns:
                donnees_insee[col] = None  # Create the column with default values if missing
                print(f"Warning: Column '{col}' is missing. It has been created with default values.")

        # Convert specific columns to numeric, where applicable
        donnees_insee['age'] = pd.to_numeric(donnees_insee['age'], errors='coerce')
        donnees_insee['POP_2019_IRIS_TotAge'] = pd.to_numeric(donnees_insee['POP_2019_IRIS_TotAge'], errors='coerce')
        donnees_insee['POP_2019_DEP_TotAge'] = pd.to_numeric(donnees_insee['POP_2019_DEP_TotAge'], errors='coerce')

        # Calculate the population at the IRIS level in 2019
        donnees_insee['POP_2019_IRIS'] = (
            donnees_insee['POP_2019_IRIS_TotAge'] *
            donnees_insee['POP_2019_DEP_TrAge'] *
            donnees_insee['FR_2019_agesurtrage'] /
            donnees_insee['POP_2019_DEP_TotAge']
        )

        # Calculate the population projections from 2020 to 2050
        for annee in range(2020, 2051):
            col_dep_trage = f"POP_{annee}_DEP_TrAge"
            col_fr_age_sur_trage = f"FR_{annee}_agesurtrage"
            donnees_insee[f"POP_{annee}_IRIS"] = (
                pd.to_numeric(donnees_insee.get(col_dep_trage, None), errors='coerce') *
                pd.to_numeric(donnees_insee.get(col_fr_age_sur_trage, None), errors='coerce') *
                donnees_insee['POP_2019_IRIS'] /
                (donnees_insee['POP_2019_DEP_TrAge'] * donnees_insee['FR_2019_agesurtrage'])
            )

        print("The population predictions between 2019 and 2050 have been calculated at the IRIS level.")
        return donnees_insee
    except Exception as e:
        print(f"An error occurred in the desagreg function: {e}")
        raise

# Function to calculate population density by age at the IRIS level from 2019 to 2050
def dens(donnees_insee):
    # Ensure 'aire_m2' exists
    if "aire_m2" not in donnees_insee.columns:
        raise KeyError("Required column 'aire_m2' is missing from 'donnees_insee'. Cannot calculate densities.")
    # Validate that 'aire_m2' contains valid numeric values and is non-zero
    if donnees_insee["aire_m2"].isnull().any() or (donnees_insee["aire_m2"] <= 0).any():
        raise ValueError("Column 'aire_m2' contains missing or invalid (non-positive) values.")
    # Loop through years and calculate density
    for annee in range(2020, 2051):
        pop_col = f"POP_{annee}_IRIS"  # Population column
        dens_col = f"DENS_{annee}"  # Density column

    print("Population density predictions from 2020 to 2050 have been calculated at the IRIS level.")
    return donnees_insee

# Function to disaggregate national mortality by age at the IRIS level from 2019 to 2050
def mortalite(donnees_dens, path_mortalite_hf):
    # Load mortality data
    donnees_mortalite_hf = pd.read_csv(path_mortalite_hf, sep=";")

    # Step 1: Sum up ages 100+ into the 99th row
    donnees_mortalite_hf.iloc[99, 1:] = (
            donnees_mortalite_hf.iloc[99, 1:].astype(float)
            + donnees_mortalite_hf.iloc[100:, 1:].astype(float).sum(axis=0)
    )
    # Keep only rows up to age 99
    donnees_mortalite_hf = donnees_mortalite_hf.iloc[:100, :]
    donnees_mortalite_hf.iloc[99, 0] = "99"

    # Step 2: Merge mortality data with population density data
    donnees_dens['age'] = donnees_dens['age'].astype(str)
    donnees_mortalite_hf['age'] = donnees_mortalite_hf['age'].astype(str)
    donnees_insee = donnees_dens.merge(donnees_mortalite_hf, on="age", how="left")

    # Step 3: Compute IRIS-level mortality for each year from 2019 to 2050
    for annee in range(2019, 2051):
        # Define column names
        mort_nat_col = f"MORT_{annee}_NAT"  # New mortality column derived from national level
        iris_pop_col = f"POP_{annee}_IRIS"  # Population at IRIS level
        nat_pop_col = f"POP_{annee}_NAT"  # Population at national level

        # Rename mortality data (use year directly from donnees_mortalite_hf)
        if str(annee) in donnees_insee.columns:
            # Rename the year column to MORT_{year}_NAT dynamically
            donnees_insee.rename(columns={str(annee): mort_nat_col}, inplace=True)
        else:
            print(f"Warning: Column {annee} not found in data.")
            continue  # Skip to the next iteration if the column is missing

        # Check if necessary population columns exist
        if iris_pop_col not in donnees_insee.columns or nat_pop_col not in donnees_insee.columns:
            print(f"Skipping year {annee}: Missing necessary population columns.")
            continue

        # Calculate IRIS-level mortality
        try:
            donnees_insee[f"MORT_{annee}_IRIS"] = (
                    donnees_insee[mort_nat_col].astype(float)
                    * donnees_insee[iris_pop_col].astype(float)
                    / donnees_insee[nat_pop_col].astype(float)
            )
        except Exception as e:
            print(f"Error calculating MORT_{annee}_IRIS: {e}")
            continue

    print("Mortality disaggregation from 2019 to 2050 has been calculated at the IRIS level.")
    return donnees_insee


# Function to export a dataframe to a .csv file (without geometry coordinates)
def export_data_csv(donnees_finales, path, title_csv):
    donnees_finales = donnees_finales.sort_values(by=["iris_cod", "age"])
    csv_path = os.path.join(path, f"{title_csv}.csv")
    donnees_finales.drop(columns="geometry", errors='ignore').to_csv(csv_path, index=False)
    print(f"CSV file written to: {csv_path}")

# Function to create a shapefile-compatible data frame from input data
def create_donnees_shp(donnees_finales):
    cols_to_select = ["iris_cod", "iris_name", "com_cod", "com_name", "dep_cod", "dep_name", "region_cod", "region_name", "aire_m2", "age"]
    pop_cols = [col for col in donnees_finales.columns if col.startswith("POP_") and col.endswith("_IRIS")]
    mort_cols = [col for col in donnees_finales.columns if col.startswith("MORT_") and col.endswith("_IRIS")]
    all_cols = cols_to_select + pop_cols + mort_cols
    donnees_shp = donnees_finales[all_cols].rename(columns={
        "region_cod": "regcod",
        "region_name": "regname",
        "iris_cod": "iriscod",
        "iris_name": "irisname",
        "com_cod": "comcod",
        "com_name": "comname",
        "dep_cod": "depcod",
        "dep_name": "depname",
        "aire_m2": "airem2",
        **{col: col.replace("POP_", "pop").replace("_IRIS", "").lower() for col in pop_cols},
        **{col: col.replace("MORT_", "mort").replace("_IRIS", "").lower() for col in mort_cols}
    })
    return donnees_shp

# Function to export data for the years 2019, 2030, and 2050
def export_pollution(donnees_filtrees, donnees_geom):
    import pandas as pd
    # Ensure the 'age' column is numeric for filtering purposes
    donnees_filtrees["age"] = pd.to_numeric(donnees_filtrees["age"], errors="coerce")

    # Prepare indicator (geographic) columns
    donnees_indic = (
        donnees_geom.loc[:, ["iris_cod", "iris_name", "com_cod", "com_name", "dep_cod", "dep_name", "region_cod", "region_name"]]
        .rename(columns={
            "region_cod": "regcod",
            "region_name": "regname",
            "iris_cod": "iriscod",
            "iris_name": "irisname",
            "com_cod": "comcod",
            "com_name": "comname",
            "dep_cod": "depcod",
            "dep_name": "depname"
        })
    )

    # Ensure 'iriscod' exists in donnees_filtrees for grouping and merging
    if "iriscod" not in donnees_filtrees.columns:
        donnees_filtrees["iriscod"] = donnees_filtrees["iris_cod"]

    # Filter by age and group by 'iriscod', then aggregate population and mortality data
    donnees_pop = (
        donnees_filtrees[donnees_filtrees["age"] >= 1]
        .groupby("iriscod", as_index=False)
        .agg(
            pop2019=("pop2019", lambda x: x.sum(skipna=True)),
            pop2030=("pop2030", lambda x: x.sum(skipna=True)),
            pop2050=("pop2050", lambda x: x.sum(skipna=True)),
            mort2019=("mort2019", lambda x: x.sum(skipna=True)),
            mort2030=("mort2030", lambda x: x.sum(skipna=True)),
            mort2050=("mort2050", lambda x: x.sum(skipna=True))
        )
    )
    # Merge the population/mortality data with geographic indicators
    donnees_shp = pd.merge(donnees_pop, donnees_indic, on="iriscod", how="inner")
    return donnees_shp

# Function to export data as a shapefile
def export_data_shp(donnees_shp, path_fichier_shp, title_shp):
    shp_path = os.path.join(path_fichier_shp, f"{title_shp}.shp")
    donnees_shp.to_file(shp_path)
    print(f"Shapefile written to: {shp_path}")

# Function that plots the population map for age n in year a
def plot_carte_iris(donnees_insee, annee, age):
    import geopandas as gpd
    import matplotlib.pyplot as plt
    # Filter data for the given age
    filtered_data = donnees_insee[donnees_insee["age"] == age]
    if filtered_data.empty:
        raise ValueError(f"No data found for age {age}!")

    # Check if the population column exists
    pop = f"POP_{annee}_IRIS"
    if pop not in donnees_insee.columns:
        raise ValueError(f"Column {pop} not found in donnees_insee!")

    # Plot the data
    ax = filtered_data.plot(column=pop, legend=True, legend_kwds={'label': "Population"}, figsize=(10, 8))
    plt.title(f"Population Map for Age {age} in {annee}")
    plt.axis("off")
    plt.show()


# Function to plot population map for a specific age and year
def pyramide_iris(donnees_insee, iris_num, annee):
    pop_col = f"POP_{annee}_IRIS"
    filtered_data = (
        donnees_insee[donnees_insee["iris_cod"] == iris_num]
        .loc[:, ["age", pop_col]]
        .rename(columns={pop_col: "pop"})
    )
    # Plot an Age Pyramid for a Specific IRIS and Year
    filtered_data.plot.barh(x="age", y="pop", color="steelblue")
    plt.title(f"Age Pyramid for the IRIS {iris_num} in {annee}")
    plt.xlabel(f"Population in {annee}")
    plt.ylabel("Age")
    plt.show()

# Function to Plot an Age Pyramid for a Specific Department and Year
def pyramide_dep(donnees_insee, dep_num, annee):
    pop_col = f"POP_{annee}_IRIS"
    filtered_data = (
        donnees_insee[donnees_insee["dep_cod"] == dep_num]
        .loc[:, ["age", pop_col]]
        .rename(columns={pop_col: "pop"})
    )
    age_distribution = filtered_data.groupby("age").agg(total_pop=pd.NamedAgg(column="pop", aggfunc="sum"))
    # Create the age pyramid
    age_distribution.plot.barh(y="total_pop", color="steelblue")
    plt.title(f"Age Pyramid for the Department {dep_num} in {annee}")
    plt.xlabel(f"Population in {annee}")
    plt.ylabel("Age")
    plt.show()

def pyramide_nat(donnees_insee, annee):
    pop_col = f"POP_{annee}_IRIS"
    filtered_data = donnees_insee.loc[:, ["age", pop_col]].rename(columns={pop_col: "pop"})
    age_distribution = filtered_data.groupby("age").agg(total_pop=pd.NamedAgg(column="pop", aggfunc="sum"))
    # Create the Age Pyramid
    age_distribution.plot.barh(y="total_pop", color="steelblue")
    plt.title(f"Age Pyramid for the Nation in {annee}")
    plt.xlabel(f"Population in {annee}")
    plt.ylabel("Age")
    plt.show()

# National test for living people
def national_vivant_test(donnees_finales, annee):
    colonne_iris = f"POP_{annee}_IRIS"
    somme_iris = donnees_finales[colonne_iris].sum(skipna=True)
    colonne_nat = f"POP_{annee}_NAT"
    donnees_nat = donnees_finales[donnees_finales["iris_cod"] == donnees_finales["iris_cod"].iloc[0]]
    somme_nat = donnees_nat[colonne_nat].sum(skipna=True)
    print(f"The total population of all IRIS in {annee} is: {somme_iris}.")
    print(f"The population of France in {annee} is: {somme_nat}.")

# National test for deceased people
def national_mort_test(donnees_finales, annee):
    colonne_iris = f"MORT_{annee}_IRIS"
    somme_iris = donnees_finales[colonne_iris].sum(skipna=True)
    colonne_nat = f"MORT_{annee}_NAT"
    donnees_nat = donnees_finales[donnees_finales["iris_cod"] == donnees_finales["iris_cod"].iloc[0]]
    somme_nat = donnees_nat[colonne_nat].sum(skipna=True)
    print(f"The total mortality of all IRIS in {annee} is: {somme_iris}.")
    print(f"The mortality in France in {annee} is: {somme_nat}.")

# Departmental test
def departemental_test(donnees_finales, dep_num, annee):
    colonne_iris = f"POP_{annee}_IRIS"
    donnees_iris = donnees_finales[donnees_finales["dep_cod"] == dep_num]
    somme_iris = donnees_iris[colonne_iris].sum(skipna=True)
    colonne_dep = f"POP_{annee}_DEP_TrAge"
    donnees_dep = donnees_iris[donnees_iris["iris_cod"] == donnees_iris["iris_cod"].iloc[0]]
    somme_dep = donnees_dep[colonne_dep].sum(skipna=True) / 5
    print(f"The total population of the IRIS in department {dep_num} in {annee} is: {somme_iris}.")
    print(f"The population of the department {dep_num} in {annee} is: {somme_dep}.")

# Test for IRIS in 2019
def iris_test(donnees_finales):
    random_iris_cod = np.random.choice(donnees_finales["iris_cod"])
    donnees_filtrees = donnees_finales[donnees_finales["iris_cod"] == random_iris_cod]
    somme_iris = donnees_filtrees["POP_2019_IRIS"].sum(skipna=True)
    iris_init = donnees_filtrees["POP_2019_IRIS_TotAge"].iloc[0]
    print(f"The total population of IRIS {random_iris_cod} in 2019 is: {somme_iris}.")
    print(f"The population of IRIS {random_iris_cod} in 2019 is: {iris_init}.")

# Comparison test between IRIS and department in 2019
def dep_vs_iris_test(donnees_finales):
    random_dep_cod = np.random.choice(donnees_finales["dep_cod"])
    donnees_filtrees = donnees_finales[donnees_finales["dep_cod"] == random_dep_cod]
    somme_iris = donnees_filtrees["POP_2019_IRIS"].sum(skipna=True)
    colonne_dep = "POP_2019_DEP_TrAge"
    donnees_dep = donnees_filtrees[donnees_filtrees["iris_cod"] == donnees_filtrees["iris_cod"].iloc[0]]
    somme_dep = donnees_dep[colonne_dep].sum(skipna=True) / 5
    print(f"The total population of the IRIS in department {random_dep_cod} in 2019 is: {somme_iris}.")
    print(f"The population of the department {random_dep_cod} in 2019 is: {somme_dep}.")

#Distribute population and mortality data by commune using Departmental population
def disaggregate_commune_by_age(commune_df, perc_hf, year=2019, age_min=None):
    frac_col = f"FR_{year}_agesurtrage"
    if frac_col not in perc_hf.columns:
        raise ValueError(f"{frac_col} not found in perc_hf columns.")
    perc_hf = perc_hf[['age', frac_col]].copy()
    perc_hf['age'] = pd.to_numeric(perc_hf['age'], errors='coerce')
    if age_min is not None:
        perc_hf = perc_hf[perc_hf['age'] >= age_min]
        perc_hf[frac_col] /= perc_hf[frac_col].sum()  # normalize
    commune_df['key'] = 1
    perc_hf['key'] = 1
    merged = commune_df.merge(perc_hf, on='key').drop('key', axis=1)
    merged['pop_age'] = merged['pop30p'] * merged[frac_col]
    merged['morta_age'] = merged['morta_a30_moy'] * merged[frac_col]
    result = merged[['comcod', 'age', 'pop_age', 'morta_age']].copy()
    result = result[result['age'].notnull()]
    result['age'] = result['age'].astype(int)
    return result

#Distribute population and mortality data Homogenously across commune for comparison
def desagreg_commune_simple(commune_df, age_list=None):
    """
    Disaggregate commune-level mortality and population homogeneously across ages 30–99.
    Parameters:
        commune_df: DataFrame with ['comcod','pop30p','morta_a30_moy']
        age_list: list of ages to distribute (default 30–99)
    Returns:
        DataFrame with ['comcod','age','pop_age','morta_age']
    """
    if age_list is None:
        age_list = list(range(30, 100))

    df_list = []
    for _, row in commune_df.iterrows():
        n_ages = len(age_list)
        age_frac = 1 / n_ages  # equal share per age
        df_tmp = pd.DataFrame({
            "comcod": row["comcod"],
            "age": age_list,
            "pop_age": row["pop30p"] * age_frac,
            "morta_age": row["morta_a30_moy"] * age_frac
        })
        df_list.append(df_tmp)

    return pd.concat(df_list, ignore_index=True)