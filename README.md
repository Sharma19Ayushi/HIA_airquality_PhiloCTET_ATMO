# HIA_airquality_PhiloCTET_ATMO

This project assesses the health impacts of air quality in France in 2050 under alternative decarbonization strategies. These strategies correspond to the four carbon neutrality scenarios published by ADEME (2021) in its foresight project Transition(s): https://www.ademe.fr/les-futurs-en-transition/.
Using projected changes in national emissions of major air pollutants for each scenario, emission-reduction targets for 2030 and 2050 (relative to 2019) are summarized and used as inputs to atmospheric modeling. Air pollutant concentrations are then estimated using the CHIMERE chemistry-transport model provided by INERIS. CHIMERE provides spatially resolved concentration fields that reflect scenario-driven changes in emissions and atmospheric processes. These modeled concentrations are subsequently combined with INSEE spatial population data (IRIS, municipality, or departmental level) to produce exposure indicators.
Finally, the project quantifies health impacts using a Health Impact Assessment (HIA) approach by applying concentration–response functions for specific pollutants. Health benefits especially avoided mortality (and associated metrics such as life years gained) are computed using established epidemiological relationships.
Importance: The project supports integrating quantified health co-benefits into public decision-making for climate and air-quality policies, strengthening the case for rapid action toward carbon neutrality in France.
 
Description of Python Modules (CHIMERE)
**Population Module**
Prepares demographic, geospatial, and mortality datasets for analysis and mapping. It cleans, aggregates, and merges multiple input sources and outputs structured (geo-enabled) datasets for downstream exposure and HIA calculations.
**CHIMERE Concentration Module**
Loads and processes CHIMERE concentration outputs (e.g., netCDF) for a given pollutant, year, and scenario. It reshapes gridded concentration fields into geospatial point datasets (with coordinates and concentration metrics) that can be intersected with administrative geometries.
**INERIS (Baseline Concentrations) Module**
Processes baseline annual pollutant concentration data (e.g., 2019) used as a reference for harmonization/correction. It standardizes pollutant naming, validates inputs, extracts lat/lon/concentration arrays from netCDF, and returns a GeoDataFrame suitable for spatial comparisons.
**ASSOCIATION Module**
generate_points(): Creates grid points within polygons and links concentration points to administrative units via an identifier (e.g., IRIS code).
calculate_perc(): Computes polygon–grid overlap percentages to weight grid points by intersection area, improving spatial allocation accuracy.
**CORRECTION Module**
Adjusts CHIMERE concentration points using proximity-based correction against reference concentrations (e.g., baseline INERIS/CHIMERE reference). Uses nearest-neighbor logic and parallelization for efficiency, updating concentrations and deltas while ensuring robust error handling.
**EXPO Module**
Computes exposure metrics by spatially combining corrected concentrations with administrative polygons and population data. Uses optimized spatial indexing (e.g., KD-trees) and parallel processing to compute mean concentrations and mean deltas by geographic identifier, including weighted averaging where required.
**MORTALITY Analysis Module**
Quantifies mortality impacts attributable to concentration changes using concentration–response functions and relative risks by pollutant. Produces outcomes such as avoided deaths and life years gained, including age-stratified outputs and uncertainty propagation via mc simulations (1000 iterations).
**MORBIDITY Analysis Module**
Quantifies morbidity impacts of selected diseases attributable to concentration changes using concentration–response functions and relative risks by pollutant. Produces outcomes such as avoided cases and years lived with disability, including age-stratified outputs and uncertainty propagation via mc simulations (1000 iterations).
**PLOT Module**
Supports geospatial visualization and export: reading shapefiles, CRS alignment checks, exporting results, and producing maps/histograms for population, exposure, and polygon characteristics with configurable styling and validation. **Note: ** .ipynb files are execution codes to use above modules and visualise the generated outputs
