import os
import unicodedata
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Morbidity filters & RR dictionary ----------------------------------
# Morbidity config (medical_cost_per_person = avg_annual_cost * duration)
morb_config = [
    {
        "short": "hypertension",
        "disease": "Hypertension",
        "pm25_age_min": 18,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.17, 1.05, 1.30),
        "rr_ug_NO2": None,
        "disability_weight": (0.08, 0.07, 0.09),
        "duration": 15.6,  # From SPF report
        "discounted_mean_duration": 12.91772502,
        "medical_cost_per_person": 621 * 15.6,
        "lag": 5
    },
    {
        "short": "lung_cancer",
        "disease": "Lung Cancer",
        "pm25_age_min": 35,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.16, 1.10, 1.23),
        "rr_ug_NO2": None,
        "disability_weight": (0.13, 0.10, 0.15),
        "duration": 2.12,
        "discounted_mean_duration": 2.064799501,
        "medical_cost_per_person": 15458 * 2.12,
        "lag": 20
    },
    {
        "short": "asthma_child",
        "disease": "Asthma in children",
        "pm25_age_min": 0,
        "pm25_age_max": 18,
        "rr_ug_PM25_RH50": (1.34, 1.10, 1.63),
        "no2_age_min": 0,
        "no2_age_max": 18,
        "rr_ug_NO2": (1.10, 1.05, 1.18),
        "disability_weight": (0.04, 0.03, 0.05),
        "duration": 16.7,
        "discounted_mean_duration": 13.65234033,
        "medical_cost_per_person": 328 * 16.7,
        "lag": 0
    },
    {
        "short": "asthma_adult",
        "disease": "Asthma in adult",
        "pm25_age_min": 18,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": None,
        "no2_age_min": 18,
        "no2_age_max": 99,
        "rr_ug_NO2": (1.10, 1.01, 1.21),
        "disability_weight": (0.04, 0.03, 0.05),
        "duration": 16.7,
        "discounted_mean_duration": 13.65234033,
        "medical_cost_per_person": 152 * 16.7,
        "lag": 0
    },
    {
        "short": "copd",
        "disease": "COPD",
        "pm25_age_min": 40,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.18, 1.13, 1.23),
        "rr_ug_NO2": None,
        "disability_weight": (0.031, 0.027, 0.035),
        "duration": 14.4,
        "discounted_mean_duration": 12.09294696,
        "medical_cost_per_person": 549 * 14.4,
        "lag": 10
    },
    {
        "short": "alri",
        "disease": "ALRI in children",
        "pm25_age_min": 0,
        "pm25_age_max": 12,
        "rr_ug_PM25_RH50": None,
        "no2_age_min": 0,
        "no2_age_max": 12,
        "rr_ug_NO2": (1.09, 1.03, 1.16),
        "disability_weight": (0.06, 0.04, 0.08),
        "duration": 0.02,
        "discounted_mean_duration": 0.019995001,
        "medical_cost_per_person": 0,
        "lag": 0
    },
    {
        "short": "stroke",
        "disease": "Stroke",
        "pm25_age_min": 35,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.16, 1.12, 1.20),
        "rr_ug_NO2": None,
        "disability_weight": (0.15, 0.11, 0.18),
        "duration": 9.7,
        "discounted_mean_duration": 8.613450097,
        "medical_cost_per_person": 4046 * 9.7,
        "lag": 5
    },
    {
        "short": "ihd",
        "disease": "IHD events",
        "pm25_age_min": 30,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.13, 1.05, 1.22),
        "rr_ug_NO2": None,
        "disability_weight": (0.03, 0.02, 0.04),
        "duration": 6.9,
        "discounted_mean_duration": 6.337668448,
        "medical_cost_per_person": 2157 * 6.9,
        "lag": 5
    },
    {
        "short": "diabetes_2",
        "disease": "Diabetes T2",
        "pm25_age_min": 45,
        "pm25_age_max": 99,
        "rr_ug_PM25_RH50": (1.10, 1.03, 1.18),
        "rr_ug_NO2": None,
        "disability_weight": (0.07, 0.05, 0.09),
        "duration": 24.5,
        "discounted_mean_duration": 18.32023246,
        "medical_cost_per_person": 2237 * 24.5,
        "lag": 10
    }
]


#Helper functions for morbidity analysis
def get_rr_key(pol: str) -> str:
    """
    Get the relative risk key corresponding to a pollutant string.
    """
    if pol.startswith("ug_PM25_RH50"):
        return "rr_ug_PM25_RH50"
    elif pol.startswith("ug_NO2"):
        return "rr_ug_NO2"
    raise ValueError(f"Unknown pollutant format: {pol}")

def get_pollutant_base(pol: str) -> str:
    """
    Normalize any pollutant string that may include suffixes (e.g., '..._mean')
    to the base key used in dictionaries and thresholds.
    Returns one of: 'ug_PM25_RH50' or 'ug_NO2'.
    """
    if pol.startswith("ug_PM25_RH50"):
        return "ug_PM25_RH50"
    if pol.startswith("ug_NO2"):
        return "ug_NO2"
    raise ValueError(f"Unknown pollutant format: {pol}")

def normalize_string(value) -> str:
    """
    Normalize a string to lowercase ASCII, replacing spaces/hyphens with underscores.
    Returns an empty string if input is None or not a string.
    """
    if not isinstance(value, str):
        return ""
    return (
        unicodedata.normalize("NFKD", value)
        .encode("ascii", "ignore")
        .decode("utf-8")
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )

def find_matching_morbidity_config(disease: str, pollutant: str, config: list):
    """
    Find the morbidity configuration for a given disease and pollutant.
    Returns:
        (config_dict, rr_key) if found and RR exists,
        (None, rr_key) if disease exists but RR for pollutant is missing,
        (None, None) if no matching disease found at all.
    """
    disease_norm = normalize_string(disease)
    rr_key = get_rr_key(pollutant)

    for cfg in config:
        cfg_disease = normalize_string(cfg.get("disease", ""))
        if cfg_disease == disease_norm:
            if rr_key not in cfg or cfg[rr_key] is None:
                return None, rr_key
            return cfg, rr_key

    return None, None
def get_rr_key(pol: str) -> str:
    """
    Get the relative risk key corresponding to a pollutant string.
    """
    if pol.startswith("ug_PM25_RH50"):
        return "rr_ug_PM25_RH50"
    elif pol.startswith("ug_NO2"):
        return "rr_ug_NO2"
    raise ValueError(f"Unknown pollutant format: {pol}")

def get_pollutant_base(pol: str) -> str:
    """
    Normalize any pollutant string that may include suffixes (e.g., '..._mean')
    to the base key used in dictionaries and thresholds.
    Returns one of: 'ug_PM25_RH50' or 'ug_NO2'.
    """
    if pol.startswith("ug_PM25_RH50"):
        return "ug_PM25_RH50"
    if pol.startswith("ug_NO2"):
        return "ug_NO2"
    raise ValueError(f"Unknown pollutant format: {pol}")

def normalize_string(value) -> str:
    """
    Normalize a string to lowercase ASCII, replacing spaces/hyphens with underscores.
    Returns an empty string if input is None or not a string.
    """
    if not isinstance(value, str):
        return ""
    return (
        unicodedata.normalize("NFKD", value)
        .encode("ascii", "ignore")
        .decode("utf-8")
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )

def find_matching_morbidity_config(disease: str, pollutant: str, config: list):
    """
    Find the morbidity configuration for a given disease and pollutant.
    Returns:
        (config_dict, rr_key) if found and RR exists,
        (None, rr_key) if disease exists but RR for pollutant is missing,
        (None, None) if no matching disease found at all.
    """
    disease_norm = normalize_string(disease)
    rr_key = get_rr_key(pollutant)

    for cfg in config:
        cfg_disease = normalize_string(cfg.get("disease", ""))
        if cfg_disease == disease_norm:
            if rr_key not in cfg or cfg[rr_key] is None:
                return None, rr_key
            return cfg, rr_key

    return None, None

# Match SpF thresholds and life expectancy
LIFE_EXPECTANCY = {2019: 80, 2030: 84, 2050: 86}
#Main function with mc simulation for propagation of uncertainties
#SPF_THRESHOLDS = {"ug_PM25_RH50": 5, "ug_NO2": 10}
SPF_THRESHOLDS = {"ug_PM25_RH50": 0, "ug_NO2": 0}
#Main function with mc simulation for propagation of uncertainties
def morbidity_mortality_mc_by_age_spf_comcod(
        donnees_merged_new,
        donnees_expo,
        annee,
        pol,
        disease,
        morb_config,
        n_mc=1000,
        seed=42,
        chunk_size=10000
):
    """
    Monte Carlo morbidity + mortality HIA
    - Optimized for memory: Uses vectorized operations and minimizes intermediate large object creation.
    """

    import numpy as np
    import pandas as pd
    import logging
    import unicodedata

    logger = logging.getLogger(__name__)
    rng = np.random.default_rng(seed)
    annee = int(annee)

    RR_MORT = {
        "ug_PM25_RH50": (1.10, 1.06, 1.13),
        "ug_NO2": (1.05, 1.03, 1.07)
    }

    if "comcod" not in donnees_merged_new.columns:
        return pd.DataFrame({"error": ["Missing comcod"]})
    #AGE IS DYNAMIC FOR ALL POLLUTANTS
    pop_col = f"pop{annee}" if f"pop{annee}" in donnees_merged_new.columns else "pop2019"
    mort_col = f"mort{annee}" if f"mort{annee}" in donnees_merged_new.columns else "mort2019"

    # Optimization: Filter by disease immediately to reduce data size
    def norm(s):
        if not isinstance(s, str): return ""
        return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower().replace(" ", "_").replace(
            "-", "_")

    target_disease_norm = norm(disease)

    # Efficiently filter the population dataframe first
    df_pop = donnees_merged_new.loc[
        donnees_merged_new["disease"].apply(norm) == target_disease_norm,
        ["comcod", "iriscod", "disease", "age", pop_col, mort_col, "absolute_incidence_iris"]
    ].copy()

    if df_pop.empty:
        return pd.DataFrame()

    df_pop["comcod"] = df_pop["comcod"].astype(str).str.upper()
    df_pop[pop_col] = pd.to_numeric(df_pop[pop_col], errors="coerce").fillna(0.0)
    df_pop[mort_col] = pd.to_numeric(df_pop[mort_col], errors="coerce").fillna(0.0)
    df_pop["absolute_incidence_iris"] = pd.to_numeric(df_pop["absolute_incidence_iris"], errors="coerce").fillna(0.0)

    # Exposure aggregation
    if "meandelta" not in donnees_expo.columns:
        return pd.DataFrame({"error": ["Missing meandelta"]})

    if "comcod" not in donnees_expo.columns and "iriscod" in donnees_expo.columns:
        tmp = donnees_expo[["iriscod", "meandelta"]].merge(
            df_pop[["iriscod", "comcod", pop_col]], on="iriscod", how="inner"
        )
        tmp["comcod"] = tmp["comcod"].astype(str).str.upper()
        tmp["weighted_delta"] = tmp["meandelta"] * tmp[pop_col]
        grp = tmp.groupby("comcod")
        df_expo_comcod = (grp["weighted_delta"].sum() / grp[pop_col].sum().replace(0, np.nan)).fillna(0).reset_index()
        df_expo_comcod.columns = ["comcod", "meandelta"]
    else:
        df_expo_comcod = (
            donnees_expo[["comcod", "meandelta"]]
            .assign(comcod=lambda d: d["comcod"].astype(str).str.upper())
            .groupby("comcod", as_index=False)["meandelta"]
            .mean()
        )

    merged = df_pop.merge(df_expo_comcod, on="comcod", how="left")
    merged["meandelta"] = pd.to_numeric(merged["meandelta"], errors="coerce").fillna(0.0)

    cfg, rr_key = find_matching_morbidity_config(disease, pol, morb_config)
    if cfg is None:
        return pd.DataFrame()

    RRm, RRm_lo, RRm_hi = cfg[rr_key]
    dw_m, dw_lo, dw_hi = cfg["disability_weight"]
    min_age = int(cfg.get("pm25_age_min", cfg.get("no2_age_min", 0)))
    max_age = int(cfg.get("pm25_age_max", cfg.get("no2_age_max", 99)))

    merged["age"] = pd.to_numeric(merged["age"], errors="coerce")
    merged = merged[(merged["age"] >= min_age) & (merged["age"] <= max_age)]

    if merged.empty:
        return pd.DataFrame()

    # Pre-calculate sampling distributions
    beta_morb_central = np.log(RRm)
    RRmort, RRmort_lo, RRmort_hi = RR_MORT.get(pol, (1, 1, 1))
    beta_mort_central = np.log(RRmort)

    esp = {2019: 80, 2030: 83, 2050: 85}.get(annee, 80)
    duration = cfg.get("duration", 1.0)

    # Optimization: Aggregate by age first to minimize loop iterations
    merged["pop_x_delta"] = merged["meandelta"] * merged[pop_col]
    age_agg = merged.groupby("age").agg({
        "pop_x_delta": "sum",
        pop_col: "sum",
        "absolute_incidence_iris": "sum",
        mort_col: "sum"
    })
    age_agg["meandelta"] = (age_agg["pop_x_delta"] / age_agg[pop_col].replace(0, np.nan)).fillna(0)

    rows = []
    for age, data in age_agg.iterrows():
        delta = data["meandelta"]
        b_cases = data["absolute_incidence_iris"]
        b_deaths = data[mort_col]
        rem_life = max(esp - age, 30)
        eff_dur = min(duration, rem_life)

        # Central
        AF_morb_c = 1 - np.exp(-beta_morb_central * delta / 10)
        AF_mort_c = 1 - np.exp(-beta_mort_central * delta / 10)

        # Initialize Monte Carlo accumulators for percentiles
        all_cases_mc = []
        all_deaths_mc = []
        all_yld_mc = []
        all_yll_mc = []

        # Process MC in chunks
        for i in range(0, n_mc, chunk_size):
            current_chunk = min(chunk_size, n_mc - i)

            beta_morb_chunk = rng.normal(np.log(RRm), (np.log(RRm_hi) - np.log(RRm_lo)) / 3.92, current_chunk)
            beta_mort_chunk = rng.normal(np.log(RRmort), (np.log(RRmort_hi) - np.log(RRmort_lo)) / 3.92, current_chunk)
            dw_samples_chunk = np.clip(rng.normal(dw_m, (dw_hi - dw_lo) / 3.92, current_chunk), 0, None)

            AF_morb_chunk = 1 - np.exp(-beta_morb_chunk * delta / 10)
            AF_mort_chunk = 1 - np.exp(-beta_mort_chunk * delta / 10)

            cases_chunk = b_cases * AF_morb_chunk
            deaths_chunk = b_deaths * AF_mort_chunk

            all_cases_mc.append(cases_chunk)
            all_deaths_mc.append(deaths_chunk)
            all_yld_mc.append(cases_chunk * dw_samples_chunk * eff_dur)
            all_yll_mc.append(deaths_chunk * rem_life)

        cases_mc = np.concatenate(all_cases_mc)
        deaths_mc = np.concatenate(all_deaths_mc)
        yld_mc = np.concatenate(all_yld_mc)
        yll_mc = np.concatenate(all_yll_mc)

        rows.append({
            "age": age,
            "avoided_cases_central": b_cases * AF_morb_c,
            "avoided_cases_low": np.percentile(cases_mc, 2.5),
            "avoided_cases_high": np.percentile(cases_mc, 97.5),
            "mortality_avoided_central": b_deaths * AF_mort_c,
            "mortality_avoided_low": np.percentile(deaths_mc, 2.5),
            "mortality_avoided_high": np.percentile(deaths_mc, 97.5),
            "YLD_central": (b_cases * AF_morb_c) * dw_m * eff_dur,
            "YLD_low": np.percentile(yld_mc, 2.5),
            "YLD_high": np.percentile(yld_mc, 97.5),
            "YLL_central": (b_deaths * AF_mort_c) * rem_life,
            "YLL_low": np.percentile(yll_mc, 2.5),
            "YLL_high": np.percentile(yll_mc, 97.5),
        })

    donnees = pd.DataFrame(rows)
    donnees[["disease", "pollutant", "year", "mc_iterations"]] = [disease, pol, annee, n_mc]
    return donnees.sort_values("age").reset_index(drop=True)