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

RR_MORTALITY_MAP = {
    "ug_PM25_RH50": (RR_PM25, RR_PM25_low, RR_PM25_high),
    "ug_NO2": (RR_NO2, RR_NO2_low, RR_NO2_high)
}

RR_MORTALITY_MAP_Sensitivity = {
    "ug_PM25_RH50": (RR_PM25_Hrapie, RR_PM25_Hrapie_low, RR_PM25_Hrapie_high),
    "ug_NO2": (RR_NO2_Hrapie, RR_NO2_Hrapie_low, RR_NO2_Hrapie_high)
}

# ----------------------------
# Relative risks (assumed predefined)
# ----------------------------
RR_values = {
    "ug_PM25_RH50": RR_PM25,
    "ug_PM25_RH50_low": RR_PM25_low,
    "ug_PM25_RH50_high": RR_PM25_high,
    "ug_NO2": RR_NO2,
    "ug_NO2_low": RR_NO2_low,
    "ug_NO2_high": RR_NO2_high
}
# Sensitivity analysis RR values
RR_values = {
    "ug_PM25_RH50": RR_PM25_Hrapie,
    "ug_PM25_RH50_low": RR_PM25_Hrapie_low,
    "ug_PM25_RH50_high": RR_PM25_Hrapie_high,
    "ug_NO2": RR_NO2_Hrapie,
    "ug_NO2_low": RR_NO2_Hrapie_low,
    "ug_NO2_high": RR_NO2_Hrapie_high
}

# ----------------------------
# Life Table Function
# ----------------------------
def life_expectancy_table(age, m_x, interval=1, a_x=0.5, radix=100000):
    """
    Calculate life expectancy at each age using a simple life table.

    Parameters input:
        age (array-like): Ages.
        m_x (array-like): Mortality rates for each age.
        interval (int): Age interval (default 1).
        a_x : Fraction of interval lived by those who die. Default is 0.5.
        radix (int): Starting population (default 100,000).

    Returns:
        e_x (np.ndarray): Life expectancy at each age.
    """
    n = len(age)

    # Assign default value to `a_x` if it is None
    if a_x is None or isinstance(a_x, (float, int)):
        a_x = np.full(n, a_x)  # Create an array if scalar passed

    # Ensure mortality rate values are valid
    m_x = np.nan_to_num(m_x)  # Replace NaNs with 0
    m_x[m_x < 0] = 0          # No negative mortality rates
    m_x = np.clip(m_x, 0, 1)  # Clip mortality rates to [0, 1]

    # Probability of dying between ages x and x+1
    q_x = (interval * m_x) / (1 + (interval - a_x) * m_x)
    q_x = np.minimum(q_x, 1.0)

    # Compute l_x (probability of surviving to exact age)
    l_x = np.zeros(n)
    l_x[0] = radix
    for i in range(1, n):
        l_x[i] = l_x[i - 1] * (1 - q_x[i - 1])

    # Compute L_x (person-years lived in each interval)
    L_x = np.zeros(n)
    for i in range(n - 1):
        L_x[i] = interval * (l_x[i + 1] + a_x[i] * (l_x[i] - l_x[i + 1]))
    if m_x[-1] > 0:  # Avoid divide by zero for last L_x calculation
        L_x[-1] = l_x[-1] / m_x[-1]
    else:
        L_x[-1] = interval * l_x[-1]

    # Compute T_x (cumulative person-years lived above each age x)
    T_x = np.zeros(n)
    T_x[-1] = L_x[-1]
    for i in range(n - 2, -1, -1):
        T_x[i] = T_x[i + 1] + L_x[i]

    # Compute e_x (life expectancy at age x)
    e_x = T_x / np.where(l_x > 0, l_x, 1e-10)  # Ensure no division by zero

    return e_x

def mortalite_age_commune_monte_carlo(
        mort_comm_age, donnees_expo, annee, pol,
        num_simulations=1000, chunk_size=10000):
    """
    Estimate mortality avoided at the commune level with robust confidence intervals (CIs)
    using Monte Carlo simulations.
    """
    # --- Step 1: Data Preparation ---
    mort_comm_age = mort_comm_age.copy()
    donnees_expo = donnees_expo.copy()

    mort_comm_age["iriscod"] = mort_comm_age["iriscod"].astype(str).str.upper().str.zfill(5)
    mort_comm_age["comcod"] = mort_comm_age["comcod"].astype(str).str.upper().str.zfill(5)
    donnees_expo["iriscod"] = donnees_expo["iriscod"].astype(str).str.upper().str.zfill(5)

    mort_comm_age = mort_comm_age[(mort_comm_age["age"] >= 30) & (mort_comm_age["age"] <= 99)]

    mort_comm_age["pop_age"] = pd.to_numeric(
        mort_comm_age[f"pop{annee}"], errors="coerce"
    ).fillna(0)
    mort_comm_age["morta_age"] = pd.to_numeric(
        mort_comm_age[f"mort{annee}"], errors="coerce"
    ).fillna(0)

    donnees_expo_commune = pd.merge(
        donnees_expo, mort_comm_age, on="iriscod", how="left"
    )

    for col in ["pop_age", "morta_age", "meanconc", "meandelta"]:
        donnees_expo_commune[col] = donnees_expo_commune[col].fillna(0)

    # --- Step 2: Commune-level population-weighted delta ---
    def weighted_mean(row):
        total_pop = row["pop_age"].sum()
        return pd.Series({
            "meandelta_commune": np.average(row["meandelta"], weights=row["pop_age"])
            if total_pop > 0 else 0.0,
            "pop_total_commune": total_pop
        })

    weighted_df = (
        donnees_expo_commune
        .groupby("comcod", group_keys=False)
        .apply(weighted_mean)
        .reset_index()
    )

    donnees_expo_commune = donnees_expo_commune.merge(
        weighted_df[["comcod", "meandelta_commune"]],
        on="comcod",
        how="left"
    )

    # --- RR selection ---
    RR_dict = {
        "ug_PM25_RH50": (RR_PM25, RR_PM25_low, RR_PM25_high),
        "ug_NO2": (RR_NO2, RR_NO2_low, RR_NO2_high),
    } #select the values for main analysis
    RR_dict_sensitivity = {
        "ug_PM25_RH50": (RR_PM25_Hrapie, RR_PM25_Hrapie_low, RR_PM25_Hrapie_high),
        "ug_NO2": (RR_NO2_Hrapie, RR_NO2_Hrapie_low, RR_NO2_Hrapie_high)
    }  #select the values for sensitivity analysis
    RR, RR_low, RR_high = RR_dict_sensitivity[pol] #adjust as per required analysis

    # --- Step 3: Age-specific mortality impacts ---
    ages = np.arange(30, 100)
    results = []

    lnRR_mean = np.log(RR)
    lnRR_sd = (np.log(RR_high) - lnRR_mean) / 1.96

    for age in ages:
        sub = donnees_expo_commune[donnees_expo_commune["age"] == age]

        pop = sub["pop_age"].sum()
        deaths = sub["morta_age"].sum()
        meandelta_commune = sub["meandelta_commune"].mean()

        simulated_RRs = np.random.lognormal(lnRR_mean, lnRR_sd, num_simulations)

        simulated_avoided = deaths * (
            1 - np.exp(-np.log(simulated_RRs) * meandelta_commune / 10)
        )

        mort_LCI = np.percentile(simulated_avoided, 2.5)
        mort_UCI = np.percentile(simulated_avoided, 97.5)

        AF_central = 1 - np.exp(-np.log(RR) * meandelta_commune / 10)
        mort_central = deaths * AF_central

        if pop > 0:
            m0 = deaths / pop
            mc = (deaths - mort_central) / pop
        else:
            m0 = mc = 0.0

        results.append({
            "age": age,
            "pop_age": pop,
            "morta_age": deaths,
            "mortpol_age": mort_central,
            "mortpol_LCI": mort_LCI,
            "mortpol_UCI": mort_UCI,
            "taux_initial": m0,
            "taux_corrige": mc,
        })

    grouped = pd.DataFrame(results)

    grouped["taux_corrige_LCI"] = np.clip(
        grouped["taux_initial"]
        - grouped["mortpol_LCI"] / grouped["pop_age"].replace(0, np.nan),
        0, 1
    )

    grouped["taux_corrige_UCI"] = np.clip(
        grouped["taux_initial"]
        - grouped["mortpol_UCI"] / grouped["pop_age"].replace(0, np.nan),
        0, 1
    )

    # --- Step 4: Life tables ---
    le_initial = life_expectancy_table(
        grouped["age"].values, grouped["taux_initial"].values,
        interval=1, a_x=0.4, radix=100000
    )

    le_corrige = life_expectancy_table(
        grouped["age"].values, grouped["taux_corrige"].values,
        interval=1, a_x=0.4, radix=100000
    )

    le_corrige_LCI = life_expectancy_table(
        grouped["age"].values, grouped["taux_corrige_LCI"].values,
        interval=1, a_x=0.4, radix=100000
    )

    le_corrige_UCI = life_expectancy_table(
        grouped["age"].values, grouped["taux_corrige_UCI"].values,
        interval=1, a_x=0.4, radix=100000
    )

    grouped["LifeTable_LE_initial"] = le_initial
    grouped["LifeTable_LE_corrige"] = le_corrige
    grouped["LifeTable_LE_corrige_LCI"] = le_corrige_LCI
    grouped["LifeTable_LE_corrige_UCI"] = le_corrige_UCI

    # --- Step 5: YLG calculation (life-table based) ---
    total_pop = grouped["pop_age"].sum()

    avg_LE_initial = np.sum(le_initial * grouped["pop_age"]) / total_pop
    avg_LE_corrige = np.sum(le_corrige * grouped["pop_age"]) / total_pop
    avg_LE_corrige_LCI = np.sum(le_corrige_LCI * grouped["pop_age"]) / total_pop
    avg_LE_corrige_UCI = np.sum(le_corrige_UCI * grouped["pop_age"]) / total_pop

    grouped["YLG_total"] = (avg_LE_corrige - avg_LE_initial) * total_pop
    grouped["YLG_total_LCI"] = (avg_LE_corrige_LCI - avg_LE_initial) * total_pop
    grouped["YLG_total_UCI"] = (avg_LE_corrige_UCI - avg_LE_initial) * total_pop

    grouped["overall_LifeTable_LEgain_mo"] = (avg_LE_corrige - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_LCI"] = (avg_LE_corrige_LCI - avg_LE_initial) * 12
    grouped["overall_LifeTable_LEgain_mo_UCI"] = (avg_LE_corrige_UCI - avg_LE_initial) * 12

    return grouped.sort_values("age").reset_index(drop=True)

