from __future__ import annotations

import numpy as np
import pandas as pd


def compute_effect_scale(
    exposure: pd.DataFrame,
    primary_location_ids,
    primary_mean: float,
    primary_lower: float,
    primary_upper: float,
    bootstrap_median: float,
    bootstrap_lower: float,
    bootstrap_upper: float,
) -> pd.DataFrame:
    """Translate the registered cumulative log-rate contrast to reader-facing scales.

    Parameters are supplied in memory. ``exposure`` is the GBD High alcohol use SEV
    table; ``primary_location_ids`` is the frozen 204-unit set. The returned ratios
    are model-implied AC-to-OC rate contrasts, not causal or individual risk ratios.
    """
    frame = exposure.copy()
    integer_columns = ("location_id", "year", "rei_id", "measure_id", "metric_id", "age_id", "sex_id")
    for column in integer_columns:
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(int)

    ids = {int(x) for x in primary_location_ids}
    if len(ids) != 204:
        raise ValueError("Expected 204 unique primary location identifiers")
    active = frame.loc[
        frame["location_id"].isin(ids)
        & frame["year"].between(2000, 2023)
        & frame["rei_id"].eq(102)
        & frame["measure_id"].eq(29)
        & frame["metric_id"].eq(3)
        & frame["age_id"].eq(22)
        & frame["sex_id"].eq(3)
    ]
    values = pd.to_numeric(active["val"], errors="raise")
    if len(values) != 204 * 24 or not np.isfinite(values).all():
        raise ValueError("Unexpected registered exposure coverage")

    source_sd = float(values.std(ddof=1))
    q25, q75 = (float(values.quantile(q)) for q in (0.25, 0.75))
    source_iqr = q75 - q25
    if source_sd <= 0 or source_iqr <= 0:
        raise ValueError("Exposure scale must be positive")

    rows = []
    for scenario, raw_shift, standardized_shift in (
        ("one_sd_equal_weight_lag_contrast", source_sd, 1.0),
        ("raw_exposure_iqr_equal_weight_lag_contrast", source_iqr, source_iqr / source_sd),
    ):
        rows.append(
            {
                "scenario": scenario,
                "raw_exposure_shift": raw_shift,
                "standardized_exposure_shift": standardized_shift,
                "model_implied_ratio_mean": float(np.exp(primary_mean * standardized_shift)),
                "posterior_ratio_lower_95": float(np.exp(primary_lower * standardized_shift)),
                "posterior_ratio_upper_95": float(np.exp(primary_upper * standardized_shift)),
                "bootstrap_ratio_median": float(np.exp(bootstrap_median * standardized_shift)),
                "bootstrap_ratio_lower_95": float(np.exp(bootstrap_lower * standardized_shift)),
                "bootstrap_ratio_upper_95": float(np.exp(bootstrap_upper * standardized_shift)),
            }
        )
    return pd.DataFrame(rows)
