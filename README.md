# HIA_airquality_PhiloCTET_ATMO
This project assesses the health impacts of air quality in France in 2050 under alternative decarbonization strategies. These strategies correspond to the four carbon neutrality scenarios published by ADEME (2021) in its foresight project Transition(s): https://www.ademe.fr/les-futurs-en-transition/.
Using projected changes in national emissions of major air pollutants for each scenario, emission-reduction targets for 2030 and 2050 (relative to 2019) are summarized and used as inputs to atmospheric modeling. Air pollutant concentrations are then estimated using the CHIMERE chemistry-transport model provided by INERIS. CHIMERE provides spatially resolved concentration fields that reflect scenario-driven changes in emissions and atmospheric processes. These modeled concentrations are subsequently combined with INSEE spatial population data (IRIS, municipality, or departmental level) to produce exposure indicators.
Finally, the project quantifies health impacts using a Health Impact Assessment (HIA) approach by applying concentration–response functions for specific pollutants. Health benefits especially avoided mortality (and associated metrics such as life years gained) are computed using established epidemiological relationships.
Importance: The project supports integrating quantified health co-benefits into public decision-making for climate and air-quality policies, strengthening the case for rapid action toward carbon neutrality in France.

This repository supports comprehensive Health Impact Assessment (HIA) studies on air quality changes, leveraging CHIMERE model outputs. The workflow is modular, ensuring clear data provenance from raw pollutant concentration to final health outcomes and spatial visualizations.
**Python Modules Overview**
**Population Module**
Prepares demographic, geospatial, and mortality datasets by cleaning, aggregating, and merging input sources. Outputs structured, geo-enabled datasets essential for exposure and HIA calculations.__ file name: pop_functions_module.py
**CHIMERE Concentration Module**
Loads and processes CHIMERE concentration outputs (e.g., netCDF) for a given pollutant, year, and scenario. It reshapes gridded concentration fields into geospatial point datasets (with coordinates and concentration metrics) that can be intersected with administrative geometries. file name: cordo_chimere_module.py. Processes baseline annual pollutant concentration data (like 2019 reference data). Standardizes pollutant names, validates inputs, and extracts georeferenced data from netCDF files. Provides GeoDataFrames for spatial harmonization and comparisons.
**ASSOCIATION Module**
generate_points(): Creates grid points within polygons and links concentration points to administrative units via an identifier (e.g., IRIS code).
calculate_perc(): Computes polygon–grid overlap percentages to weight grid points by intersection area, improving spatial allocation accuracy.
**CORRECTION Module**
Adjusts CHIMERE concentration points using proximity-based correction against reference concentrations (e.g., baseline INERIS/CHIMERE reference). Uses nearest-neighbor logic and parallelization for efficiency, updating concentrations and deltas while ensuring robust error handling.
**Pourcent_file_export.ipynb**
Process and export percentage data (Area weighted % of intersection) for associating concentrations to IRIS using area-weighing approach.
**EXPO Module**
Computes exposure metrics by spatially combining corrected concentrations with administrative polygons and population data. Uses optimized spatial indexing (e.g., KD-trees) and parallel processing to compute mean concentrations and mean deltas by geographic identifier, including weighted averaging where required.
**MORTALITY Analysis Module**
Quantifies mortality impacts attributable to concentration changes using concentration–response functions and relative risks by pollutant. Produces outcomes such as avoided deaths and life years gained, including age-stratified outputs and uncertainty propagation via mc simulations (1000 iterations).
**MORBIDITY Analysis Module**
Quantifies morbidity impacts of selected diseases attributable to concentration changes using concentration–response functions and relative risks by pollutant. Produces outcomes such as avoided cases and years lived with disability, including age-stratified outputs and uncertainty propagation via mc simulations (1000 iterations).**helper_morbidity_analysis.ipynb**: A support notebook for generating Disability Weights, and Durations from GBD 2019 files for morbidity analysis.
**eqis_chimere_cluster.ipynb & eqis_chimere_new.ipynb**
Workflows related to HIAs using Chimere outputs, for overall scenario and cluster analysis.
**Visualisation of HIA**
Figures_Delayed impact.ipynb: Generates figures and analysis regarding the delayed impact of air quality changes on outcomes, using pre-computed csv files. 
Figures_Immediate Impact.ipynb: Focuses on visualizing and analyzing the immediate impacts of air quality interventions.
**Sensitivity_analysis.ipynb**
Carries out sensitivity analyses on model parameters and assumptions to test result robustness.
**National Life Expectancy Gains.ipynb**
Estimates or visualizes changes in national life expectancy due to modifications in air quality.
**cessation_lag_morbidity.ipynb & cessation_lag_mortality.ipynb**
Analyze how time lags after cessation of pollution affect morbidity and mortality outcomes.
**marginal_effect_checks.ipynb**
Checks or visualizes marginal effects—such as the additional effect of changing one factor while holding others constant.
**sensitivity_cessation_lags.ipynb**
Sensitivity analyses focused on how different assumptions about cessation lag affect results.
**temporal_lag_progression.ipynb**
Studies how effects (e.g., health outcomes) evolve over time after exposure changes.
