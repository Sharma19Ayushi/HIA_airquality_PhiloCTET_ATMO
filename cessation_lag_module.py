import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

## Relative risks
RR_PM25 = 1.15      #1.15
RR_PM25_high = 1.25
RR_PM25_low = 1.05
RR_NO2 = 1.023
RR_NO2_high = 1.037
RR_NO2_low = 1.008
RR_PM25_Hrapie = 1.10     #exposure range 5-70
RR_PM25_Hrapie_high = 1.13
RR_PM25_Hrapie_low = 1.06
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

def _resolve_rr_key(pollutant_key, rr_map):
    k = str(pollutant_key).strip()
    if k in rr_map:
        return k

    kl = k.lower()
    for kk in rr_map.keys():
        if str(kk).strip().lower() == kl:
            return kk

    raise KeyError(f"{pollutant_key} not found in RR_MORTALITY_MAP. Available: {list(rr_map.keys())}")
#Functions for cessation lag weights
def comeap_delay_weights(
        lag_years=20,
        first_year_share=0.30,
        years_2_5_share=0.50,
        years_6_20_share=0.20,
):
    """
    COMEAP Cessation Lag Weights with flexible full-effect realization over `lag_years`:
    - Year 1: `first_year_share`
    - Years 2-5: `years_2_5_share` evenly split (as long as those years exist)
    - Years 6..lag_years: remaining share evenly split (so full effect is realized by `lag_years`)
    """
    lag_years = int(lag_years)
    if lag_years <= 0:
        raise ValueError("lag_years must be > 0.")

    first_year_share = float(first_year_share)
    years_2_5_share = float(years_2_5_share)
    years_6_20_share = float(years_6_20_share)

    total = first_year_share + years_2_5_share + years_6_20_share
    if total <= 0:
        raise ValueError("Sum of shares must be > 0.")
    first_year_share /= total
    years_2_5_share /= total
    years_6_20_share /= total

    w = np.zeros(lag_years, dtype=float)

    # Year 1
    w[0] = first_year_share

    # Years 2-5 (if present)
    if lag_years >= 2:
        end_2_5 = min(lag_years, 5)
        n_2_5 = max(0, end_2_5 - 1)
        if n_2_5 > 0:
            w[1:end_2_5] = years_2_5_share / float(n_2_5)

    # Years 6..lag_years (if present)
    if lag_years >= 6:
        n_6_L = lag_years - 5
        if n_6_L > 0:
            w[5:lag_years] = years_6_20_share / float(n_6_L)

    s = float(w.sum())
    return w / s if s > 0 else w

def linear_delay_weights(lag_years=20):
    """
    Linear lag weights (uniform realization over lag_years).
    """
    lag_years = int(lag_years)
    if lag_years <= 0:
        raise ValueError("lag_years must be > 0.")
    return np.ones(lag_years, dtype=float) / float(lag_years)

def realized_effect_fraction(eval_years, base_year, delay_weights):
    """
    Compute cumulative realization fractions using delay weights relative to base_year.
    """
    eval_years = np.asarray(eval_years, dtype=int)
    L = len(delay_weights)
    out = np.zeros_like(eval_years, dtype=float)
    for i, y in enumerate(eval_years):
        difference = y - int(base_year)
        if difference < 0:
            out[i] = 0.0
        else:
            out[i] = float(delay_weights[: min(L, difference + 1)].sum())
    return out

def generate_cessation_lag_fractions(base_year, eval_years, lag_years=20):
    """
    Generate cessation lag fractions using:
    - COMEAP methodology (flexible full realization over `lag_years`)
    - Linear (uniform) lag
    """
    delay_weights_comeap = comeap_delay_weights(lag_years=lag_years)
    delay_weights_linear = linear_delay_weights(lag_years=lag_years)

    frac_comeap = realized_effect_fraction(eval_years, base_year, delay_weights_comeap)
    frac_linear = realized_effect_fraction(eval_years, base_year, delay_weights_linear)

    return delay_weights_comeap, frac_comeap, delay_weights_linear, frac_linear

def adjusted_rr_with_cessation_lag(rr_full_effect, effect_fraction):
    """
    Adjust RR values using cessation lag fractions:
         RR_adjusted = 1.0 + (RR_full_effect - 1.0) * effect_fraction
    """
    rr_full_effect = np.asarray(rr_full_effect, dtype=float)
    effect_fraction = np.asarray(effect_fraction, dtype=float)
    return 1.0 + (rr_full_effect - 1.0) * effect_fraction

def _pick_exposure_series(donnees_expo, pollutant):
    if pollutant in donnees_expo.columns:
        s = pd.to_numeric(donnees_expo[pollutant], errors="coerce")
        if s.notna().any():
            return s

    numeric_cols = donnees_expo.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) == 0:
        s = pd.to_numeric(donnees_expo.select_dtypes(exclude=[object]).iloc[:, 0], errors="coerce")
        return s

    s = pd.to_numeric(donnees_expo[numeric_cols[0]], errors="coerce")
    return s

def exposure_mean_from_cache(expo_cache, pollutant):
    means = {}
    for y, df in expo_cache.items():
        s = _pick_exposure_series(df, pollutant)
        means[int(y)] = float(np.nanmean(s.to_numpy(dtype=float)))
    return means

ANCHOR_YEARS = [2019, 2030, 2050]

def _validate_anchor_years(anchor_years):
    anchor_years = sorted(set(int(y) for y in anchor_years))
    if len(anchor_years) < 2:
        raise ValueError("Need at least two anchor years to interpolate.")
    return anchor_years
ANCHOR_YEARS = _validate_anchor_years(ANCHOR_YEARS)

def interpolate_exposure_means_linear(eval_years, anchor_means):
    anchor_years = _validate_anchor_years(anchor_means.keys())
    x = np.array(anchor_years, dtype=float)
    y = np.array([float(anchor_means[yr]) for yr in anchor_years], dtype=float)
    f = interp1d(x, y, kind="linear", fill_value="extrapolate")
    return np.asarray(f(np.asarray(eval_years, dtype=float)), dtype=float)

def rr_from_exposure_change(exposure_year, exposure_base, rr_per_unit, unit_change=10.0):
    exposure_year = np.asarray(exposure_year, dtype=float)
    rr_per_unit = float(rr_per_unit)
    unit_change = float(unit_change)
    if unit_change <= 0:
        raise ValueError("unit_change must be > 0.")
    # If exposure_base is None: compute absolute RR relative to 0 exposure (RR reflects full exposure level)
    if exposure_base is None:
        return np.power(rr_per_unit, exposure_year / unit_change)
    exposure_base = float(exposure_base)
    delta = exposure_year - exposure_base
    return np.power(rr_per_unit, delta / unit_change)

print("Successfully loaded Cessation Lag functions")

#helper function to compute the overall excel sheet with data from all files:
import re
import os

# -------------------------
# Helper Functions
# -------------------------
import os
import re

os.environ["NUMEXPR_MAX_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

def _normalize_outcome_text(outcome_text: str) -> str:
    if pd.isna(outcome_text):
        return ""
    cleaned = str(outcome_text).strip().lower()
    cleaned = re.sub(r"\([^)]*\)", "", cleaned)
    cleaned = cleaned.replace("–", "-").replace("—", "-")
    cleaned = re.sub(r"[^0-9a-z\s]", " ", cleaned)
    return re.sub(r"\s+", " ", cleaned).strip()


_MORT_STEMS = [
    "multiplier_lag_v_no_lag_AVOIDED__s4_ug_PM25_RH50_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s4_ug_NO2_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s3_ug_PM25_RH50_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s3_ug_NO2_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s2_ug_PM25_RH50_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s2_ug_NO2_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s1_ug_PM25_RH50_2010_2070",
    "multiplier_lag_v_no_lag_AVOIDED__s1_ug_NO2_2010_2070",
]

_MORB_STEMS_BY_POL = {
    "ug_PM25_RH50": "morbidity_multiplier_lag_v_no_lag_AVOIDED__ug_PM25_RH50_2019_2070",
    "ug_NO2": "morbidity_multiplier_lag_v_no_lag_AVOIDED__ug_NO2_2019_2070",
}

_MULTIPLIER_COL = "cumulative_multiplier_lag_over_no_lag"
CESSATION_LAGS_DIR = r"data\2-output-data\Cessation_lags"

_MORB_MULTIPLIER_CACHE = {}
_MORB_CONFIG_CACHE = {}


def _to_valid_year(year, default):
    year_numeric = pd.to_numeric(year, errors="coerce")
    if not np.isfinite(year_numeric):
        return None, default
    return int(year_numeric), None

def _read_cached_csv(csv_path, cache=None):
    if not os.path.exists(csv_path):
        return None
    if cache is not None and csv_path in cache:
        return cache[csv_path].copy()
    df = pd.read_csv(csv_path)
    if cache is not None:
        cache[csv_path] = df.copy()
    return df.copy()

def _pick_multiplier_from_csv(df: pd.DataFrame, year_int: int, default: float = 1.0) -> float:
    if df is None or df.empty:
        return default

    df = df.copy()
    df.columns = df.columns.str.strip()

    if _MULTIPLIER_COL not in df.columns:
        return default

    year_col = next((col for col in ("year", "Year", "annee", "Annee") if col in df.columns), None)

    if year_col is not None:
        df["_year"] = pd.to_numeric(df[year_col], errors="coerce")
        match = df[df["_year"] == year_int]
        if match.empty:
            match = df[df["_year"] <= year_int].sort_values("_year").tail(1)
    else:
        period_start_col = next(
            (col for col in ("period_start", "start_year", "start", "from_year", "from") if col in df.columns),
            None,
        )
        period_end_col = next(
            (col for col in ("period_end", "end_year", "end", "to_year", "to") if col in df.columns),
            None,
        )

        if period_start_col is None or period_end_col is None:
            return default

        df["_period_start"] = pd.to_numeric(df[period_start_col], errors="coerce")
        df["_period_end"] = pd.to_numeric(df[period_end_col], errors="coerce")
        match = df[(df["_period_start"] <= year_int) & (df["_period_end"] >= year_int)]

    if match.empty:
        return default

    value = pd.to_numeric(match[_MULTIPLIER_COL].iloc[0], errors="coerce")
    return float(value) if np.isfinite(value) else default

def _pick_morbidity_multiplier(
        health_outcome,
        pol,
        year,
        disease=None,
        cfg=None,
        base_dir=CESSATION_LAGS_DIR,
        default=1.0,
):
    year_int, fallback = _to_valid_year(year, default)
    if fallback is not None:
        return fallback

    pol_key = str(pol).strip()
    stem = _MORB_STEMS_BY_POL.get(pol_key)
    if not stem:
        return default

    csv_path = os.path.join(base_dir, f"{stem}.csv")
    cache_key = f"{pol_key}|{csv_path}"
    if cache_key in _MORB_MULTIPLIER_CACHE:
        df = _MORB_MULTIPLIER_CACHE[cache_key].copy()
    else:
        df = _read_cached_csv(csv_path)
        if df is None:
            return default
        _MORB_MULTIPLIER_CACHE[cache_key] = df.copy()

    if df.empty:
        return default

    df.columns = df.columns.str.strip()
    outcome_col = find_col(df, "endpoint_short", "endpoint", "disease")
    if not outcome_col:
        return _pick_multiplier_from_csv(df, year_int, default=default)

    candidates = []

    if health_outcome is not None:
        candidates.append(_normalize_outcome_text(health_outcome))
    if disease is not None:
        candidates.append(_normalize_outcome_text(disease))

    if isinstance(cfg, dict):
        for key in ("endpoint_short", "endpoint", "disease"):
            value = cfg.get(key)
            if value is not None:
                candidates.append(_normalize_outcome_text(value))

    candidates = list(dict.fromkeys(candidate for candidate in candidates if candidate))
    if not candidates:
        return _pick_multiplier_from_csv(df, year_int, default=default)

    df["_outcome_norm"] = df[outcome_col].map(_normalize_outcome_text)

    df_match = df[df["_outcome_norm"].isin(candidates)]
    if df_match.empty:
        for candidate in candidates:
            df_match = df[df["_outcome_norm"].str.contains(re.escape(candidate), na=False)]
            if not df_match.empty:
                break

    if not df_match.empty:
        return _pick_multiplier_from_csv(df_match, year_int, default=default)

    return _pick_multiplier_from_csv(df, year_int, default=default)

def cessation_lag(health_outcome, pol, year, scenario=None, base_dir=CESSATION_LAGS_DIR, default=1.0):
    year_int, fallback = _to_valid_year(year, default)
    if fallback is not None:
        return fallback

    outcome_normalized = _normalize_outcome_text(health_outcome)
    pollutant_normalized = str(pol).strip()
    is_mortality = "mortality" in outcome_normalized

    if is_mortality:
        if scenario is None:
            return default

        sc = str(scenario).strip().lower()
        stem = next(
            (item for item in _MORT_STEMS if f"__{sc}_" in item and f"_{pollutant_normalized}_" in item),
            None,
        )
        if stem is None:
            return default

        csv_path = os.path.join(base_dir, f"{stem}.csv")
        df = _read_cached_csv(csv_path)
        if df is None:
            print(f"[Warning] Mortality lag file not found: {csv_path}")
            return default

        return _pick_multiplier_from_csv(df, year_int, default=default)

    stem = _MORB_STEMS_BY_POL.get(pollutant_normalized)
    if not stem:
        return default

    csv_path = os.path.join(base_dir, f"{stem}.csv")
    df = _read_cached_csv(csv_path)
    if df is None:
        return default

    return _pick_multiplier_from_csv(df, year_int, default=default)

def find_col(df, *candidates):
    cols = list(df.columns)
    lower_map = {col.lower(): col for col in cols}

    for candidate in candidates:
        if candidate:
            match = lower_map.get(candidate.lower())
            if match:
                return match

    for candidate in candidates:
        if candidate:
            candidate_lower = candidate.lower()
            for col in cols:
                if candidate_lower in col.lower():
                    return col

    return None

def safe_sum(df, col):
    if df is None or df.empty or col is None or col not in df.columns:
        return 0.0

    series = pd.to_numeric(df[col], errors="coerce").fillna(0)

    if str(col).lower() in {"ylg", "ylg_lci", "ylg_uci"}:
        age_col = next((c for c in df.columns if "age" in str(c).lower()), None)
        if age_col is not None:
            age_series = pd.to_numeric(df[age_col], errors="coerce")
            mask = age_series.between(30, 98, inclusive="both")
            series = series[mask]

    return float(series.sum())

def safe_series_sum(df, col):
    return safe_sum(df, col)

def safe_first(df, col):
    if df is None or df.empty or col is None or col not in df.columns:
        return 0.0
    series = pd.to_numeric(df[col], errors="coerce").dropna()
    if series.empty:
        return 0.0
    return float(series.iloc[0])



