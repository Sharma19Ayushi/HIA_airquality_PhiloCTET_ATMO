import numpy as np
import pandas as pd
import scipy.stats as stats
import logging
import geopandas as gpd
import os

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

    if a_x is None:
        a_x = np.full(n, 0.5, dtype=float)
    elif np.isscalar(a_x):
        a_x = np.full(n, float(a_x), dtype=float)
    else:
        a_x = np.asarray(a_x, dtype=float)
        if len(a_x) != n:
            raise ValueError("`a_x` must be a scalar or have the same length as `age`.")

    m_x = np.asarray(
        np.nan_to_num(m_x, nan=0.0, posinf=np.finfo(float).max, neginf=0.0),
        dtype=float
    )
    m_x[m_x < 0] = 0.0

    q_x = (interval * m_x) / (1.0 + (interval - a_x) * m_x)
    q_x = np.clip(q_x, 0.0, 1.0)

    l_x = np.zeros(n, dtype=float)
    l_x[0] = float(radix)
    for i in range(1, n):
        l_x[i] = l_x[i - 1] * (1.0 - q_x[i - 1])

    L_x = np.zeros(n, dtype=float)
    for i in range(n - 1):
        L_x[i] = interval * (l_x[i + 1] + a_x[i] * (l_x[i] - l_x[i + 1]))

    # Open-ended last age interval
    if m_x[-1] > 0:
        L_x[-1] = l_x[-1] / m_x[-1]
    else:
        L_x[-1] = interval * l_x[-1]

    T_x = np.zeros(n, dtype=float)
    T_x[-1] = L_x[-1]
    for i in range(n - 2, -1, -1):
        T_x[i] = T_x[i + 1] + L_x[i]

    e_x = T_x / np.where(l_x > 0, l_x, np.nan)
    return e_x
############################################################################################
# CODE TO APPLY POP-WEIGHTED COMMUNE CONC FOR MORTALITY IMPACTS
############################################################################################
def mortalite_age_commune_monte_carlo(
    mort_comm_age, donnees_expo, donnees_transformed, annee, pol, num_simulations=1000, seed=None
):
    import numpy as np
    import pandas as pd

    if seed is not None:
        np.random.seed(seed)

    mort_comm_age       = mort_comm_age.copy()
    donnees_expo        = donnees_expo.copy()
    donnees_transformed = donnees_transformed.copy()

    def clean_iris(s):
        return s.astype(str).str.strip().str.upper().str.replace(r"\.0$", "", regex=True).str.zfill(9)
    def clean_com(s):
        return s.astype(str).str.strip().str.upper().str.replace(r"\.0$", "", regex=True).str.zfill(5)

    mort_comm_age["iriscod"] = clean_iris(mort_comm_age["iriscod"])
    mort_comm_age["comcod"]  = clean_com(mort_comm_age["comcod"])
    donnees_expo["iriscod"]        = clean_iris(donnees_expo["iriscod"])
    donnees_transformed["iriscod"] = clean_iris(donnees_transformed["iriscod"])
    donnees_transformed["comcod"]  = clean_com(donnees_transformed["comcod"])

    # --- DEDUPLICATION: enforce one row per (iriscod, age) ---
    before = len(mort_comm_age)
    mort_comm_age = mort_comm_age.drop_duplicates(subset=["iriscod", "age"]).reset_index(drop=True)
    after = len(mort_comm_age)
    if before != after:
        print(f"[DEDUP] mort_comm_age: {before:,} → {after:,} rows ({before - after:,} duplicates removed)")

    donnees_expo = donnees_expo.merge(
        donnees_transformed[["iriscod", "comcod"]].drop_duplicates(),
        on="iriscod", how="left"
    )
    n_null = donnees_expo["comcod"].isna().sum()
    if n_null > 0:
        print(f"WARNING: {n_null} IRISes in expo not found in donnees_transformed")
    donnees_expo["comcod"] = donnees_expo["comcod"].fillna(donnees_expo["iriscod"].str[:5])

    mort_comm_age["age"]          = pd.to_numeric(mort_comm_age["age"],          errors="coerce")
    mort_comm_age[f"pop{annee}"]  = pd.to_numeric(mort_comm_age[f"pop{annee}"],  errors="coerce")
    mort_comm_age[f"mort{annee}"] = pd.to_numeric(mort_comm_age[f"mort{annee}"], errors="coerce")
    donnees_expo["meandelta"]     = pd.to_numeric(donnees_expo["meandelta"],      errors="coerce")

    mort_comm_age = mort_comm_age.dropna(subset=["age", f"pop{annee}", f"mort{annee}"])
    mort_comm_age = mort_comm_age[(mort_comm_age["age"] >= 30) & (mort_comm_age["age"] <= 99)].copy()
    mort_comm_age["pop_age"]   = mort_comm_age[f"pop{annee}"].clip(lower=0)
    mort_comm_age["morta_age"] = mort_comm_age[f"mort{annee}"].clip(lower=0)

    # --- Step 1: broadcast delta to every IRIS×age row ---
    mort_comm_age = mort_comm_age.merge(
        donnees_expo[["iriscod", "comcod", "meandelta"]].drop_duplicates("iriscod"),
        on="iriscod", how="left", suffixes=("", "_expo")
    )
    mort_comm_age["comcod"] = mort_comm_age["comcod_expo"].combine_first(mort_comm_age["comcod"])
    mort_comm_age = mort_comm_age.drop(columns=["comcod_expo"], errors="ignore")
    mort_comm_age["meandelta"] = mort_comm_age["meandelta"].fillna(0)

    # --- Step 2: commune×age aggregation ---
    mort_comm_age["pop_x_delta"] = mort_comm_age["pop_age"] * mort_comm_age["meandelta"]

    comm_age = (
        mort_comm_age.groupby(["comcod", "age"], as_index=False)
        .agg(
            pop_age     = ("pop_age",     "sum"),
            morta_age   = ("morta_age",   "sum"),
            pop_x_delta = ("pop_x_delta", "sum"),
        )
    )
    comm_age["meandelta_comm"] = np.where(
        comm_age["pop_age"] > 0,
        comm_age["pop_x_delta"] / comm_age["pop_age"],
        0.0
    )
    comm_age = comm_age.drop(columns=["pop_x_delta"])
    comm_age["comcod"] = clean_com(comm_age["comcod"].astype(str))

    # --- Step 3: Coverage reporting ---
    pop_total    = comm_age["pop_age"].sum()
    mort_total   = comm_age["morta_age"].sum()
    pop_covered  = comm_age.loc[comm_age["meandelta_comm"] > 0,   "pop_age"].sum()
    pop_zero     = comm_age.loc[comm_age["meandelta_comm"] == 0,  "pop_age"].sum()
    pop_excluded = comm_age.loc[comm_age["meandelta_comm"].isna(), "pop_age"].sum()
    mort_covered = comm_age.loc[comm_age["meandelta_comm"] > 0,   "morta_age"].sum()
    mort_excluded= comm_age.loc[comm_age["meandelta_comm"].isna(), "morta_age"].sum()

    n_comm_covered = comm_age.loc[comm_age["meandelta_comm"] > 0,   "comcod"].nunique()
    n_comm_zero    = comm_age.loc[comm_age["meandelta_comm"] == 0,  "comcod"].nunique()
    n_comm_na      = comm_age.loc[comm_age["meandelta_comm"].isna(),"comcod"].nunique()

    print(f"[COVERAGE] Communes with delta > 0   : {n_comm_covered}")
    print(f"[COVERAGE] Communes with delta = 0   : {n_comm_zero}")
    print(f"[COVERAGE] Communes with delta NaN   : {n_comm_na} (excluded)")
    print(f"\n[COVERAGE] Population with delta > 0 : {pop_covered:>12,.0f} ({pop_covered/pop_total*100:.1f}%)")
    print(f"[COVERAGE] Population with delta = 0 : {pop_zero:>12,.0f} ({pop_zero/pop_total*100:.1f}%)")
    print(f"[COVERAGE] Population excluded (NaN) : {pop_excluded:>12,.0f} ({pop_excluded/pop_total*100:.1f}%)")
    print(f"\n[COVERAGE] Deaths with delta > 0     : {mort_covered:>12,.0f} ({mort_covered/mort_total*100:.1f}%)")
    print(f"[COVERAGE] Deaths excluded (NaN)     : {mort_excluded:>12,.0f} ({mort_excluded/mort_total*100:.1f}%)")

    covered_mask = comm_age["meandelta_comm"] > 0
    if covered_mask.any():
        natl_delta = np.average(
            comm_age.loc[covered_mask, "meandelta_comm"],
            weights=comm_age.loc[covered_mask, "pop_age"].clip(lower=1e-9)
        )
        RR_ref = float(RR_PM25_Hrapie if pol == "ug_PM25_RH50" else RR_NO2_Hrapie)
        print(f"\n[SANITY] Pop-weighted mean delta (covered): {natl_delta:.4f} µg/m³")
        print(f"[SANITY] Expected avoided deaths if fully covered: "
              f"{mort_total * (1.0 - np.exp(-np.log(RR_ref) * natl_delta / 10.0)):,.1f}")

    # --- Step 4: RR and central AF ---
    RR_dict_main = {
        "ug_PM25_RH50": (RR_PM25_Hrapie, RR_PM25_Hrapie_low, RR_PM25_Hrapie_high),
        "ug_NO2":       (RR_NO2_Hrapie,  RR_NO2_Hrapie_low,  RR_NO2_Hrapie_high),
    }

    RR, RR_low, RR_high = RR_dict_main[pol]
    RR, RR_low, RR_high = float(RR), float(RR_low), float(RR_high)

    lnRR_mean = np.log(RR)
    lnRR_sd   = max((np.log(RR_high) - np.log(RR_low)) / (2 * 1.96), 1e-12)
    sim_lnRRs = np.random.normal(lnRR_mean, lnRR_sd, size=num_simulations)

    delta_vals     = comm_age["meandelta_comm"].fillna(0).values
    comm_age["AF"] = 1.0 - np.exp(-lnRR_mean * delta_vals / 10.0)
    comm_age["mortpol"] = comm_age["morta_age"] * comm_age["AF"]

    # --- Step 5: aggregate to age level (all communes, for avoided deaths / YLG) ---
    comm_age["mort_x_delta"] = comm_age["morta_age"] * comm_age["meandelta_comm"].fillna(0)

    age_summary = (
        comm_age.groupby("age", as_index=False)
        .agg(
            pop_age      = ("pop_age",      "sum"),
            morta_age    = ("morta_age",    "sum"),
            mortpol_age  = ("mortpol",      "sum"),
            n_communes   = ("comcod",       "nunique"),
        )
    )

    # MC CIs for avoided deaths (all communes)
    age_eff_delta_all = np.nan_to_num(
        np.where(
            age_summary["morta_age"].values > 0,
            comm_age.groupby("age")["mort_x_delta"].sum().reindex(age_summary["age"]).values
            / age_summary["morta_age"].values,
            0.0
        ), nan=0.0
    )
    deaths_age      = age_summary["morta_age"].values[:, np.newaxis]
    delta_age_all   = age_eff_delta_all[:, np.newaxis]
    sim_avoided_age = deaths_age * (1.0 - np.exp(-sim_lnRRs[np.newaxis, :] * delta_age_all / 10.0))

    age_summary["mortpol_LCI"] = np.percentile(sim_avoided_age, 2.5,  axis=1)
    age_summary["mortpol_UCI"] = np.percentile(sim_avoided_age, 97.5, axis=1)

    # --- Step 6: LE rates and life table — EXPOSED POPULATION ONLY (delta > 0) ---
    comm_age_exposed = comm_age[comm_age["meandelta_comm"] > 0].copy()

    age_exposed = (
        comm_age_exposed.groupby("age", as_index=False)
        .agg(
            pop_age_exp   = ("pop_age",       "sum"),
            morta_age_exp = ("morta_age",     "sum"),
            mortpol_exp   = ("mortpol",       "sum"),
            mort_x_delta  = ("mort_x_delta",  "sum"),
        )
    )

    age_summary = age_summary.merge(age_exposed, on="age", how="left")
    age_summary[["pop_age_exp", "morta_age_exp", "mortpol_exp", "mort_x_delta"]] = \
        age_summary[["pop_age_exp", "morta_age_exp", "mortpol_exp", "mort_x_delta"]].fillna(0)

    # Mort-weighted effective delta for exposed pop — use nan_to_num, not .fillna (numpy array)
    age_summary["eff_delta_exp"] = np.nan_to_num(
        np.where(
            age_summary["morta_age_exp"].values > 0,
            age_summary["mort_x_delta"].values / age_summary["morta_age_exp"].values,
            0.0
        ), nan=0.0
    )

    # MC CIs for LE — exposed population only
    deaths_exp      = age_summary["morta_age_exp"].values[:, np.newaxis]
    delta_exp       = age_summary["eff_delta_exp"].values[:, np.newaxis]
    sim_avoided_exp = deaths_exp * (1.0 - np.exp(-sim_lnRRs[np.newaxis, :] * delta_exp / 10.0))

    mortpol_exp_LCI = np.percentile(sim_avoided_exp, 2.5,  axis=1)
    mortpol_exp_UCI = np.percentile(sim_avoided_exp, 97.5, axis=1)

    # Mortality rates for exposed population only
    pop_nz_exp = age_summary["pop_age_exp"].replace(0, np.nan)

    age_summary["taux_initial"]     = (age_summary["morta_age_exp"] / pop_nz_exp).fillna(0)
    age_summary["taux_corrige"]     = ((age_summary["morta_age_exp"] - age_summary["mortpol_exp"])  / pop_nz_exp).clip(0).fillna(0)
    age_summary["taux_corrige_LCI"] = ((age_summary["morta_age_exp"] - mortpol_exp_LCI)             / pop_nz_exp).clip(0, 1).fillna(0)
    age_summary["taux_corrige_UCI"] = ((age_summary["morta_age_exp"] - mortpol_exp_UCI)             / pop_nz_exp).clip(0, 1).fillna(0)

    # Life table on ages 30–99 directly; e_x[0] = remaining LE at age 30 (e_30)
    for rate_col, e30_col, ex_col in [
        ("taux_initial",     "e30_initial",     "LifeTable_LE_initial"),
        ("taux_corrige",     "e30_corrige",     "LifeTable_LE_corrige"),
        ("taux_corrige_LCI", "e30_corrige_LCI", "LifeTable_LE_corrige_LCI"),
        ("taux_corrige_UCI", "e30_corrige_UCI", "LifeTable_LE_corrige_UCI"),
    ]:
        e_x = life_expectancy_table(
            age_summary["age"].values,
            age_summary[rate_col].values,
            interval=1, a_x=0.5, radix=100000
        )
        age_summary[ex_col]  = e_x      # remaining LE at each age 30–99
        age_summary[e30_col] = e_x[0]   # scalar: remaining LE at age 30

    # --- e_30-based LE gain (scalar) ---
    overall_LEgain_mo     = (age_summary["e30_corrige"].iloc[0]     - age_summary["e30_initial"].iloc[0]) * 12
    overall_LEgain_mo_LCI = (age_summary["e30_corrige_LCI"].iloc[0] - age_summary["e30_initial"].iloc[0]) * 12
    overall_LEgain_mo_UCI = (age_summary["e30_corrige_UCI"].iloc[0] - age_summary["e30_initial"].iloc[0]) * 12

    # --- Pop-weighted remaining LE gain (ages 30–99, exposed pop) ---
    total_pop_exp = age_summary["pop_age_exp"].sum()
    if total_pop_exp > 0:
        pw_LE_initial     = np.sum(age_summary["LifeTable_LE_initial"]     * age_summary["pop_age_exp"]) / total_pop_exp
        pw_LE_corrige     = np.sum(age_summary["LifeTable_LE_corrige"]     * age_summary["pop_age_exp"]) / total_pop_exp
        pw_LE_corrige_LCI = np.sum(age_summary["LifeTable_LE_corrige_LCI"] * age_summary["pop_age_exp"]) / total_pop_exp
        pw_LE_corrige_UCI = np.sum(age_summary["LifeTable_LE_corrige_UCI"] * age_summary["pop_age_exp"]) / total_pop_exp
    else:
        pw_LE_initial = pw_LE_corrige = pw_LE_corrige_LCI = pw_LE_corrige_UCI = 0.0

    pw_LEgain_mo     = (pw_LE_corrige     - pw_LE_initial) * 12
    pw_LEgain_mo_LCI = (pw_LE_corrige_LCI - pw_LE_initial) * 12
    pw_LEgain_mo_UCI = (pw_LE_corrige_UCI - pw_LE_initial) * 12

    # --- YLG (uses all-commune mortpol but exposed-pop LE for remaining life) ---
    age_summary["YLG"]     = age_summary["LifeTable_LE_initial"] * age_summary["mortpol_age"]
    age_summary["YLG_LCI"] = age_summary["LifeTable_LE_initial"] * age_summary["mortpol_LCI"]
    age_summary["YLG_UCI"] = age_summary["LifeTable_LE_initial"] * age_summary["mortpol_UCI"]

    total_avoided = age_summary["mortpol_age"].sum()
    total_YLG     = age_summary["YLG"].sum()

    age_summary["overall_avoided_deaths"]        = total_avoided
    age_summary["overall_avoided_deaths_LCI"]    = age_summary["mortpol_LCI"].sum()
    age_summary["overall_avoided_deaths_UCI"]    = age_summary["mortpol_UCI"].sum()
    age_summary["overall_YLG"]                   = total_YLG
    age_summary["overall_YLG_LCI"]               = age_summary["YLG_LCI"].sum()
    age_summary["overall_YLG_UCI"]               = age_summary["YLG_UCI"].sum()
    age_summary["overall_YLG_per_avoided_death"] = total_YLG / total_avoided if total_avoided > 0 else 0.0

    age_summary["overall_LifeTable_LEgain_mo"]     = overall_LEgain_mo
    age_summary["overall_LifeTable_LEgain_mo_LCI"] = overall_LEgain_mo_LCI
    age_summary["overall_LifeTable_LEgain_mo_UCI"] = overall_LEgain_mo_UCI

    age_summary["overall_pw_LEgain_mo"]     = pw_LEgain_mo
    age_summary["overall_pw_LEgain_mo_LCI"] = pw_LEgain_mo_LCI
    age_summary["overall_pw_LEgain_mo_UCI"] = pw_LEgain_mo_UCI
    # --- Full-coverage scaling: adjust LE gains for unexposed population ---
    # Observed gains are restricted to exposed pop (83%). To estimate national
    # full-coverage equivalent, scale by coverage fraction.
    coverage_fraction = pop_covered / pop_total if pop_total > 0 else 1.0
    # Scale pop-weighted gain by coverage fraction (unexposed contribute 0 gain)
    pw_LEgain_mo_fullcov     = pw_LEgain_mo     * coverage_fraction
    pw_LEgain_mo_fullcov_LCI = pw_LEgain_mo_LCI * coverage_fraction
    pw_LEgain_mo_fullcov_UCI = pw_LEgain_mo_UCI * coverage_fraction
    # e30 gain scaled similarly — assumes unexposed have same age structure
    e30_LEgain_mo_fullcov     = overall_LEgain_mo     * coverage_fraction
    e30_LEgain_mo_fullcov_LCI = overall_LEgain_mo_LCI * coverage_fraction
    e30_LEgain_mo_fullcov_UCI = overall_LEgain_mo_UCI * coverage_fraction

    age_summary["coverage_fraction"]                  = coverage_fraction
    age_summary["overall_pw_LEgain_mo_fullcov"]       = pw_LEgain_mo_fullcov
    age_summary["overall_pw_LEgain_mo_fullcov_LCI"]   = pw_LEgain_mo_fullcov_LCI
    age_summary["overall_pw_LEgain_mo_fullcov_UCI"]   = pw_LEgain_mo_fullcov_UCI
    age_summary["overall_e30_LEgain_mo_fullcov"]      = e30_LEgain_mo_fullcov
    age_summary["overall_e30_LEgain_mo_fullcov_LCI"]  = e30_LEgain_mo_fullcov_LCI
    age_summary["overall_e30_LEgain_mo_fullcov_UCI"]  = e30_LEgain_mo_fullcov_UCI

    print(f"\n[COVERAGE ADJUSTMENT]")
    print(f"[COVERAGE ADJUSTMENT] Exposed population fraction : {coverage_fraction*100:.1f}%")
    print(f"[COVERAGE ADJUSTMENT] LE gain e30 - exposed only  : {overall_LEgain_mo:>8.2f}  [{overall_LEgain_mo_LCI:.2f}, {overall_LEgain_mo_UCI:.2f}] mo")
    print(f"[COVERAGE ADJUSTMENT] LE gain e30 - full nat. pop : {e30_LEgain_mo_fullcov:>8.2f}  [{e30_LEgain_mo_fullcov_LCI:.2f}, {e30_LEgain_mo_fullcov_UCI:.2f}] mo")
    print(f"[COVERAGE ADJUSTMENT] LE gain pw  - exposed only  : {pw_LEgain_mo:>8.2f}  [{pw_LEgain_mo_LCI:.2f}, {pw_LEgain_mo_UCI:.2f}] mo")
    print(f"[COVERAGE ADJUSTMENT] LE gain pw  - full nat. pop : {pw_LEgain_mo_fullcov:>8.2f}  [{pw_LEgain_mo_fullcov_LCI:.2f}, {pw_LEgain_mo_fullcov_UCI:.2f}] mo")

    print(f"\n[RESULT] Total avoided deaths       : {total_avoided:>12,.1f}")
    print(f"[RESULT] Total YLG                  : {total_YLG:>12,.1f}")
    print(f"[RESULT] LE gain e30 (mo)           : {overall_LEgain_mo:>12,.2f}  [{overall_LEgain_mo_LCI:.2f}, {overall_LEgain_mo_UCI:.2f}]")
    print(f"[RESULT] LE gain pop-weighted (mo)  : {pw_LEgain_mo:>12,.2f}  [{pw_LEgain_mo_LCI:.2f}, {pw_LEgain_mo_UCI:.2f}]")

    return age_summary.sort_values("age").reset_index(drop=True)