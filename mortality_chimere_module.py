import numpy as np
import pandas as pd
import scipy.stats as stats
import logging

#Modified function of health_impact assessment for mortality with Monte Carlo simulations for CIs
## Relative risks
RR_PM25 = 1.15      #1.15
RR_PM25_high = 1.25
RR_PM25_low = 1.05
RR_PM25_Hrapie = 1.10     #exposure range 5-70
RR_PM25_Hrapie_high = 1.13
RR_PM25_Hrapie_low = 1.06
RR_NO2 = 1.023
RR_NO2_high = 1.037
RR_NO2_low = 1.008
RR_NO2_Hrapie = 1.05    #exposure range 10-130
RR_NO2_Hrapie_high = 1.07
RR_NO2_Hrapie_low = 1.03
print(f'loaded defined RR values')

RR_MORTALITY_MAP_Main = {
    "ug_PM25_RH50": (RR_PM25_Hrapie, RR_PM25_Hrapie_low, RR_PM25_Hrapie_high),
    "ug_NO2": (RR_NO2_Hrapie, RR_NO2_Hrapie_low, RR_NO2_Hrapie_high)
}

RR_MORTALITY_MAP_sensitivity = {
    "ug_PM25_RH50": (RR_PM25, RR_PM25_low, RR_PM25_high),
    "ug_NO2": (RR_NO2, RR_NO2_low, RR_NO2_high)
}
# ----------------------------
# Life Table Function
# ----------------------------
def life_expectancy_table(age, m_x, interval=1, a_x=0.5, radix=100000):
    """
    Calculate life expectancy at each age using a standard abridged life table.

    Parameters:
        age (array-like): Ages.
        m_x (array-like): Mortality rates for each age.
        interval (int or float): Age interval length (default 1).
        a_x (float or array-like): Fraction of interval lived by those who die.
            Default is 0.5.
        radix (int or float): Starting population (default 100,000; scaling factor only).

    Returns:
        np.ndarray: Life expectancy at each age (e_x).
    """
    age = np.asarray(age)
    n = len(age)

    # Ensure `a_x` is an array of length n
    if a_x is None:
        a_x = np.full(n, 0.5, dtype=float)
    elif np.isscalar(a_x):
        a_x = np.full(n, float(a_x), dtype=float)
    else:
        a_x = np.asarray(a_x, dtype=float)
        if len(a_x) != n:
            raise ValueError("`a_x` must be a scalar or have the same length as `age`.")

    # Clean mortality rates
    m_x = np.asarray(
        np.nan_to_num(m_x, nan=0.0, posinf=np.finfo(float).max, neginf=0.0),
        dtype=float
    )
    m_x[m_x < 0] = 0.0

    # Standard probability of dying in the interval
    q_x = (interval * m_x) / (1.0 + (interval - a_x) * m_x)
    q_x = np.clip(q_x, 0.0, 1.0)

    # Survivors to exact age x
    l_x = np.zeros(n, dtype=float)
    l_x[0] = float(radix)
    for i in range(1, n):
        l_x[i] = l_x[i - 1] * (1.0 - q_x[i - 1])

    # Person-years lived in each interval
    L_x = np.zeros(n, dtype=float)
    for i in range(n - 1):
        L_x[i] = interval * (l_x[i + 1] + a_x[i] * (l_x[i] - l_x[i + 1]))

    # Open-ended last age interval
    if m_x[-1] > 0:
        L_x[-1] = l_x[-1] / m_x[-1]
    else:
        L_x[-1] = interval * l_x[-1]

    # Cumulative person-years remaining above age x
    T_x = np.zeros(n, dtype=float)
    T_x[-1] = L_x[-1]
    for i in range(n - 2, -1, -1):
        T_x[i] = T_x[i + 1] + L_x[i]

    # Life expectancy at age x
    e_x = T_x / np.where(l_x > 0, l_x, np.nan)

    return e_x

#Function to compute mortality avoided at the commune level
# with robust confidence intervals (CIs) using Monte Carlo simulations.
def mortalite_age_commune_monte_carlo(
        mort_comm_age, donnees_expo, annee, pol,
        num_simulations=1000):
    # --- Step 1: Data Preparation ---
    mort_comm_age = mort_comm_age.copy()
    donnees_expo  = donnees_expo.copy()

    mort_comm_age["iriscod"] = mort_comm_age["iriscod"].astype(str).str.upper().str.zfill(5)
    mort_comm_age["comcod"]  = mort_comm_age["comcod"].astype(str).str.upper().str.zfill(5)
    donnees_expo["iriscod"]  = donnees_expo["iriscod"].astype(str).str.upper().str.zfill(5)
    donnees_expo["comcod"]   = donnees_expo["iriscod"].str[:5]

    mort_comm_age = mort_comm_age[
        (mort_comm_age["age"] >= 30) & (mort_comm_age["age"] <= 99)
    ].copy()

    # Compute per-IRIS population and attach to donnees_expo
    mort_comm_age["pop_age"] = pd.to_numeric(
        mort_comm_age[f"pop{annee}"], errors="coerce"
    ).fillna(0)
    mort_comm_age["morta_age"] = pd.to_numeric(
        mort_comm_age[f"mort{annee}"], errors="coerce"
    ).fillna(0)

    # Attach IRIS population to exposure for pop-weighted aggregation
    pop_iris = (
        mort_comm_age.groupby(['iriscod'])['pop_age'].sum().reset_index()
    )
    donnees_expo = donnees_expo.merge(pop_iris, on='iriscod', how='left')
    donnees_expo['pop_age'] = donnees_expo['pop_age'].fillna(0)

    # --- Step 3: Aggregate exposure to commune level (pop-weighted meandelta)
    expo_commune = (
        donnees_expo
        .groupby("comcod", group_keys=False)
        .apply(
            lambda g: pd.Series({
                "meandelta_commune": (
                    np.average(g["meandelta"], weights=g["pop_age"])
                    if g["pop_age"].sum() > 0 else g["meandelta"].mean()
                )
            }),
            include_groups=False
        )
        .reset_index()
    )

    # --- Step 4: Aggregate mortality to commune × age level
    mort_commune_age = (
        mort_comm_age
        .groupby(["comcod", "age"], as_index=False)
        .agg(
            pop_age   =("pop_age",   "sum"),
            morta_age =("morta_age", "sum"),
        )
    )

    # --- Step 5: Final merge: add commune-level exposure
    merged = mort_commune_age.merge(expo_commune, on="comcod", how="left")
    merged["meandelta_commune"] = merged["meandelta_commune"].fillna(0)

    # --- RR selection ---
    RR_dict_main = {
        "ug_PM25_RH50": (RR_PM25_Hrapie, RR_PM25_Hrapie_low, RR_PM25_Hrapie_high),
        "ug_NO2":       (RR_NO2_Hrapie,  RR_NO2_Hrapie_low,  RR_NO2_Hrapie_high),
    }
    RR, RR_low, RR_high = RR_dict_main[pol]
    lnRR_mean = np.log(RR)
    lnRR_sd   = (np.log(RR_high) - np.log(RR_low)) / (2 * 1.96)

    # --- Step 7: Age-level aggregation (population- and exposure-weighted)
    ages    = np.arange(30, 100)
    results = []
    for age in ages:
        sub = merged[merged["age"] == age]
        pop    = sub["pop_age"].sum()
        deaths = sub["morta_age"].sum()
        meandelta_commune = (
            np.average(sub["meandelta_commune"], weights=sub["pop_age"])
            if pop > 0 else 0.0
        )
        simulated_RRs     = np.random.lognormal(lnRR_mean, lnRR_sd, num_simulations)
        simulated_avoided = deaths * (
            1 - np.exp(-np.log(simulated_RRs) * meandelta_commune / 10)
        )
        mort_LCI = np.percentile(simulated_avoided, 2.5)
        mort_UCI = np.percentile(simulated_avoided, 97.5)
        AF_central   = 1 - np.exp(-np.log(RR) * meandelta_commune / 10)
        mort_central = deaths * AF_central
        m0 = deaths / pop if pop > 0 else 0.0
        mc = (deaths - mort_central) / pop if pop > 0 else 0.0
        results.append({
            "age":               age,
            "pop_age":           pop,
            "morta_age":         deaths,
            "meandelta_commune": meandelta_commune,
            "AF_central":        AF_central,
            "mortpol_age":       mort_central,
            "mortpol_LCI":       mort_LCI,
            "mortpol_UCI":       mort_UCI,
            "taux_initial":      m0,
            "taux_corrige":      mc,
        })

    grouped     = pd.DataFrame(results)
    pop_nonzero = grouped["pop_age"].replace(0, np.nan)

    grouped["taux_corrige_LCI"] = np.clip(
        grouped["taux_initial"] - grouped["mortpol_LCI"] / pop_nonzero, 0, 1
    ).fillna(0)
    grouped["taux_corrige_UCI"] = np.clip(
        grouped["taux_initial"] - grouped["mortpol_UCI"] / pop_nonzero, 0, 1
    ).fillna(0)

    # --- Step 8: Life table and summary stats ---
    le_initial     = life_expectancy_table(grouped["age"].values, grouped["taux_initial"].values,     interval=1, a_x=0.5, radix=1)
    le_corrige     = life_expectancy_table(grouped["age"].values, grouped["taux_corrige"].values,     interval=1, a_x=0.5, radix=1)
    le_corrige_LCI = life_expectancy_table(grouped["age"].values, grouped["taux_corrige_LCI"].values, interval=1, a_x=0.5, radix=1)
    le_corrige_UCI = life_expectancy_table(grouped["age"].values, grouped["taux_corrige_UCI"].values, interval=1, a_x=0.5, radix=1)

    grouped["LifeTable_LE_initial"]     = le_initial
    grouped["LifeTable_LE_corrige"]     = le_corrige
    grouped["LifeTable_LE_corrige_LCI"] = le_corrige_LCI
    grouped["LifeTable_LE_corrige_UCI"] = le_corrige_UCI

    grouped["LifeTable_LEgain_mo"]     = (le_corrige     - le_initial) * 12
    grouped["LifeTable_LEgain_mo_LCI"] = (le_corrige_LCI - le_initial) * 12
    grouped["LifeTable_LEgain_mo_UCI"] = (le_corrige_UCI - le_initial) * 12

    grouped["YLG"]     = grouped["LifeTable_LE_initial"] * grouped["mortpol_age"]
    grouped["YLG_LCI"] = grouped["LifeTable_LE_initial"] * grouped["mortpol_LCI"]
    grouped["YLG_UCI"] = grouped["LifeTable_LE_initial"] * grouped["mortpol_UCI"]

    grouped["YLG_per_avoided_death"] = np.where(
        grouped["mortpol_age"] > 0,
        grouped["YLG"] / grouped["mortpol_age"],
        grouped["LifeTable_LE_initial"]
    )

    total_pop = grouped["pop_age"].sum()
    if total_pop > 0:
        avg_LE_initial     = np.sum(le_initial     * grouped["pop_age"]) / total_pop
        avg_LE_corrige     = np.sum(le_corrige     * grouped["pop_age"]) / total_pop
        avg_LE_corrige_LCI = np.sum(le_corrige_LCI * grouped["pop_age"]) / total_pop
        avg_LE_corrige_UCI = np.sum(le_corrige_UCI * grouped["pop_age"]) / total_pop
    else:
        avg_LE_initial = avg_LE_corrige = avg_LE_corrige_LCI = avg_LE_corrige_UCI = 0.0

    grouped["overall_LifeTable_LEgain_mo"]     = (avg_LE_corrige     - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_LCI"] = (avg_LE_corrige_LCI - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_UCI"] = (avg_LE_corrige_UCI - avg_LE_initial) * 12

    total_avoided_deaths     = grouped["mortpol_age"].sum()
    total_avoided_deaths_LCI = grouped["mortpol_LCI"].sum()
    total_avoided_deaths_UCI = grouped["mortpol_UCI"].sum()
    total_YLG                = grouped["YLG"].sum()
    total_YLG_LCI            = grouped["YLG_LCI"].sum()
    total_YLG_UCI            = grouped["YLG_UCI"].sum()

    grouped["overall_avoided_deaths"]        = total_avoided_deaths
    grouped["overall_avoided_deaths_LCI"]    = total_avoided_deaths_LCI
    grouped["overall_avoided_deaths_UCI"]    = total_avoided_deaths_UCI
    grouped["overall_YLG"]                   = total_YLG
    grouped["overall_YLG_LCI"]               = total_YLG_LCI
    grouped["overall_YLG_UCI"]               = total_YLG_UCI
    grouped["overall_YLG_per_avoided_death"] = (
        total_YLG / total_avoided_deaths if total_avoided_deaths > 0 else 0.0
    )

    return grouped.sort_values("age").reset_index(drop=True)

############################################################################################
# CODE TO APPLY WHO THRESHOLDS FOR SENSTIVITY
############################################################################################
def mortalite_age_commune_monte_carlo_who(
        mort_comm_age, donnees_expo, annee, pol,
        num_simulations=1000):
    # --- Step 1: Data Preparation ---
    mort_comm_age = mort_comm_age.copy()
    donnees_expo  = donnees_expo.copy()

    mort_comm_age["iriscod"] = mort_comm_age["iriscod"].astype(str).str.upper().str.zfill(5)
    mort_comm_age["comcod"]  = mort_comm_age["iriscod"].str[:5]
    mort_comm_age["comcod"] = mort_comm_age["comcod"].astype(str).str.upper().str.zfill(5)
    donnees_expo["iriscod"]  = donnees_expo["iriscod"].astype(str).str.upper().str.zfill(5)
    donnees_expo["comcod"]   = donnees_expo["iriscod"].str[:5]
    donnees_expo["comcod"] = donnees_expo["comcod"].astype(str).str.upper().str.zfill(5)

    mort_comm_age = mort_comm_age[(mort_comm_age["age"] >= 30) & (mort_comm_age["age"] <= 99)]

    mort_comm_age["pop_age"] = pd.to_numeric(
        mort_comm_age[f"pop{annee}"], errors="coerce"
    ).fillna(0)
    mort_comm_age["morta_age"] = pd.to_numeric(
        mort_comm_age[f"mort{annee}"], errors="coerce"
    ).fillna(0)

    # --- Attach IRIS-level population to exposure ---
    pop_iris = mort_comm_age.groupby(['iriscod'])['pop_age'].sum().reset_index()
    donnees_expo = donnees_expo.merge(pop_iris, on='iriscod', how='left')
    donnees_expo['pop_age'] = donnees_expo['pop_age'].fillna(0)

    # --- Step 2: Aggregate pop-weighted exposure delta at the commune level (vs WHO threshold) ---
    SPF_THRESHOLDS = {"ug_PM25_RH50": 5, "ug_NO2": 10}
    threshold = SPF_THRESHOLDS.get(pol, 0.0)

    def weighted_delta_who(g):
        total_pop = g["pop_age"].sum()
        meanconc_weighted = np.average(g["meanconc"], weights=g["pop_age"]) if total_pop > 0 else 0.0
        return pd.Series({
            "meandelta_commune": max(meanconc_weighted - threshold, 0.0),
            "pop_total_commune": total_pop
        })

    expo_commune = (
        donnees_expo
        .groupby("comcod", group_keys=False)
        .apply(weighted_delta_who, include_groups=False)
        .reset_index()
    )

    # --- Step 3: Aggregate mortality to commune x age level ---
    mort_commune_age = (
        mort_comm_age
        .groupby(["comcod", "age"], as_index=False)
        .agg(
            pop_age   =("pop_age", "sum"),
            morta_age =("morta_age", "sum"),
        )
    )

    # --- Step 4: Merge commune exposure ---
    merged = mort_commune_age.merge(expo_commune, on="comcod", how="left")
    merged["meandelta_commune"] = merged["meandelta_commune"].fillna(0)

    # --- RR selection (WHO/Sante Publique values) ---
    RR_dict_who = {
        "ug_PM25_RH50": (RR_PM25, RR_PM25_low, RR_PM25_high),
        "ug_NO2":       (RR_NO2,  RR_NO2_low,  RR_NO2_high),
    }
    RR, RR_low, RR_high = RR_dict_who[pol]
    lnRR_mean = np.log(RR)
    lnRR_sd = (np.log(RR_high) - np.log(RR_low)) / (2 * 1.96)

    # --- Step 5: Age-wise MC simulation and attribution ---
    ages = np.arange(30, 100)
    results = []
    for age in ages:
        sub = merged[merged["age"] == age]
        pop = sub["pop_age"].sum()
        deaths = sub["morta_age"].sum()
        meandelta_commune = np.average(sub["meandelta_commune"], weights=sub["pop_age"]) if pop > 0 else 0.0
        simulated_RRs = np.random.lognormal(lnRR_mean, lnRR_sd, num_simulations)
        simulated_avoided = deaths * (1 - np.exp(-np.log(simulated_RRs) * meandelta_commune / 10))
        mort_LCI = np.percentile(simulated_avoided, 2.5)
        mort_UCI = np.percentile(simulated_avoided, 97.5)
        AF_central = 1 - np.exp(-np.log(RR) * meandelta_commune / 10)
        mort_central = deaths * AF_central
        m0 = deaths / pop if pop > 0 else 0.0
        mc = (deaths - mort_central) / pop if pop > 0 else 0.0
        results.append({
            "age": age,
            "pop_age": pop,
            "morta_age": deaths,
            "meandelta_commune": meandelta_commune,
            "AF_central": AF_central,
            "mortpol_age": mort_central,
            "mortpol_LCI": mort_LCI,
            "mortpol_UCI": mort_UCI,
            "taux_initial": m0,
            "taux_corrige": mc,
        })
    grouped = pd.DataFrame(results)
    pop_nonzero = grouped["pop_age"].replace(0, np.nan)
    grouped["taux_corrige_LCI"] = np.clip(
        grouped["taux_initial"] - grouped["mortpol_LCI"] / pop_nonzero, 0, 1
    ).fillna(0)
    grouped["taux_corrige_UCI"] = np.clip(
        grouped["taux_initial"] - grouped["mortpol_UCI"] / pop_nonzero, 0, 1
    ).fillna(0)

    # --- Step 6: Life table and summary stats ---
    le_initial     = life_expectancy_table(grouped["age"].values, grouped["taux_initial"].values,     interval=1, a_x=0.5, radix=1)
    le_corrige     = life_expectancy_table(grouped["age"].values, grouped["taux_corrige"].values,     interval=1, a_x=0.5, radix=1)
    le_corrige_LCI = life_expectancy_table(grouped["age"].values, grouped["taux_corrige_LCI"].values, interval=1, a_x=0.5, radix=1)
    le_corrige_UCI = life_expectancy_table(grouped["age"].values, grouped["taux_corrige_UCI"].values, interval=1, a_x=0.5, radix=1)

    grouped["LifeTable_LE_initial"]     = le_initial
    grouped["LifeTable_LE_corrige"]     = le_corrige
    grouped["LifeTable_LE_corrige_LCI"] = le_corrige_LCI
    grouped["LifeTable_LE_corrige_UCI"] = le_corrige_UCI

    grouped["LifeTable_LEgain_mo"]     = (le_corrige     - le_initial) * 12
    grouped["LifeTable_LEgain_mo_LCI"] = (le_corrige_LCI - le_initial) * 12
    grouped["LifeTable_LEgain_mo_UCI"] = (le_corrige_UCI - le_initial) * 12

    grouped["YLG"]     = grouped["LifeTable_LE_initial"] * grouped["mortpol_age"]
    grouped["YLG_LCI"] = grouped["LifeTable_LE_initial"] * grouped["mortpol_LCI"]
    grouped["YLG_UCI"] = grouped["LifeTable_LE_initial"] * grouped["mortpol_UCI"]

    grouped["YLG_per_avoided_death"] = np.where(
        grouped["mortpol_age"] > 0,
        grouped["YLG"] / grouped["mortpol_age"],
        grouped["LifeTable_LE_initial"]
    )

    total_pop = grouped["pop_age"].sum()
    if total_pop > 0:
        avg_LE_initial     = np.sum(le_initial     * grouped["pop_age"]) / total_pop
        avg_LE_corrige     = np.sum(le_corrige     * grouped["pop_age"]) / total_pop
        avg_LE_corrige_LCI = np.sum(le_corrige_LCI * grouped["pop_age"]) / total_pop
        avg_LE_corrige_UCI = np.sum(le_corrige_UCI * grouped["pop_age"]) / total_pop
    else:
        avg_LE_initial = avg_LE_corrige = avg_LE_corrige_LCI = avg_LE_corrige_UCI = 0.0

    grouped["overall_LifeTable_LEgain_mo"]     = (avg_LE_corrige     - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_LCI"] = (avg_LE_corrige_LCI - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_UCI"] = (avg_LE_corrige_UCI - avg_LE_initial) * 12

    total_avoided_deaths     = grouped["mortpol_age"].sum()
    total_avoided_deaths_LCI = grouped["mortpol_LCI"].sum()
    total_avoided_deaths_UCI = grouped["mortpol_UCI"].sum()
    total_YLG                = grouped["YLG"].sum()
    total_YLG_LCI            = grouped["YLG_LCI"].sum()
    total_YLG_UCI            = grouped["YLG_UCI"].sum()

    grouped["overall_avoided_deaths"]        = total_avoided_deaths
    grouped["overall_avoided_deaths_LCI"]    = total_avoided_deaths_LCI
    grouped["overall_avoided_deaths_UCI"]    = total_avoided_deaths_UCI
    grouped["overall_YLG"]                   = total_YLG
    grouped["overall_YLG_LCI"]               = total_YLG_LCI
    grouped["overall_YLG_UCI"]               = total_YLG_UCI
    grouped["overall_YLG_per_avoided_death"] = (
        total_YLG / total_avoided_deaths if total_avoided_deaths > 0 else 0.0
    )

    return grouped.sort_values("age").reset_index(drop=True)

