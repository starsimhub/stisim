"""
Shared utilities for STIsim interventions.

Coverage parsing, stratified targeting, and related helpers used across
intervention classes (ART, VMMC, etc.).
"""

import warnings
import starsim as ss
import numpy as np
import pandas as pd
import sciris as sc


# %% Coverage parsing

def parse_coverage(data, valid_names=None, yearvec=None, smoothness=0, **interp_kw):
    """
    Parse coverage data into a per-timestep array.

    Returns (coverage_array, coverage_format, age_bins, sex_keys) where:
        - coverage_array: 1D array (len=yearvec) for aggregate, or dict of arrays for stratified
        - coverage_format: 'n' (absolute numbers) or 'p' (proportion)
        - age_bins: list of age bin strings if stratified, else None
        - sex_keys: list of sex values if stratified, else None

    Examples::

        # Constant proportion
        parse_coverage(0.9, yearvec=yearvec)

        # Time-varying dict (format auto-detected from values)
        parse_coverage({'year': [2000, 2020], 'value': [0, 0.9]}, yearvec=yearvec)

        # Time-varying dict with explicit format
        parse_coverage({'year': [2000, 2020], 'value': [1000, 50000], 'format': 'n'}, yearvec=yearvec)

        # Single-column DataFrame
        parse_coverage(pd.DataFrame(index=years, data={'n_art': vals}), valid_names=['n_art'], yearvec=yearvec)

        # Age/sex stratified DataFrame (Gender column is optional)
        parse_coverage(strat_df, valid_names=['p_art'], yearvec=yearvec)

        # Smooth interpolation between data points
        parse_coverage(data, yearvec=yearvec, smoothness=5)

    Args:
        data:        coverage input in any of the formats above, or None
        valid_names: list of valid column names, e.g. ['n_art', 'p_art']
        yearvec:     simulation year vector for interpolation
        smoothness:  interpolation smoothness (0=linear, higher=smoother S-curves);
                     passed to sc.smoothinterp
        **interp_kw: additional keyword arguments passed to sc.smoothinterp
                     (e.g. method='nearest', growth=0.1)
    """
    if valid_names is None:
        valid_names = []

    # None → no coverage
    if data is None:
        return None, None, None, None

    # Scalar → constant proportion
    if sc.isnumber(data):
        arr = np.full(len(yearvec), float(data))
        return arr, 'p', None, None

    # Dict with year + value keys → interpolate
    if isinstance(data, dict) and 'year' in data and 'value' in data:
        years  = np.array(data['year'], dtype=float)
        values = np.array(data['value'], dtype=float)
        fmt = data.get('format', None)
        if fmt is None:
            fmt = 'p' if np.all(values <= 1.0) else 'n'
        arr = sc.smoothinterp(yearvec, years, values, smoothness=smoothness, **interp_kw)
        return arr, fmt, None, None

    # DataFrame — check for stratified vs single-column
    if isinstance(data, pd.DataFrame):
        return _parse_coverage_df(data, valid_names, yearvec, smoothness=smoothness, **interp_kw)

    errormsg = f'Coverage data format not recognized: {type(data)}. Expected None, number, dict, or DataFrame.'
    raise ValueError(errormsg)


def _parse_coverage_df(data, valid_names, yearvec, smoothness=0, **interp_kw):
    """
    Parse a DataFrame of coverage data.

    Handles two cases:
        1. Single-column: index=years, column in valid_names → 1D interpolated array
        2. Stratified: columns include Year, AgeBin, and optionally Gender/Sex → dict of arrays by stratum
    """

    # Check for stratified format — normalize column names for flexible matching
    col_lower = {c.lower(): c for c in data.columns}
    has_year   = 'year' in col_lower
    has_agebin = 'agebin' in col_lower or 'age_bin' in col_lower or 'age' in col_lower

    if has_year and has_agebin:  # Gender/Sex column is optional
        return _parse_stratified_df(data, yearvec, smoothness=smoothness, **interp_kw)

    # Single-column format: index=years, column in valid_names (e.g. n_art or p_art)
    if len(data.columns) == 1 and data.columns[0] in valid_names:
        colname = data.columns[0]
        fmt = 'n' if colname.startswith('n_') else 'p'
        arr = sc.smoothinterp(yearvec, data.index.values, data[colname].values, smoothness=smoothness, **interp_kw)
        return arr, fmt, None, None

    # Try to find a valid column
    for col in data.columns:
        if col in valid_names:
            fmt = 'n' if col.startswith('n_') else 'p'
            if data.index.name in ['year', 'Year'] or np.all(data.index > 1900):
                arr = sc.smoothinterp(yearvec, data.index.values, data[col].values, smoothness=smoothness, **interp_kw)
            else:
                arr = sc.smoothinterp(yearvec, np.arange(len(data)), data[col].values, smoothness=smoothness, **interp_kw)
            return arr, fmt, None, None

    errormsg = f'DataFrame columns {list(data.columns)} do not match any valid names {valid_names}.'
    raise ValueError(errormsg)


def _normalize_sex(val):
    """
    Normalize a gender/sex value to integer form.
    Convention: 0 = female, 1 = male.

    Accepts: 0, 1, 'f', 'm', 'female', 'male' (case-insensitive).
    """
    if isinstance(val, str):
        v = val.strip().lower()
        if v in ('f', 'female'): return 0
        if v in ('m', 'male'):   return 1
        raise ValueError(f'Cannot parse sex value: {val!r}. Expected 0/1, f/m, or female/male.')
    return int(val)


def _normalize_stratified_cols(data):
    """
    Normalize column names and gender values in a stratified DataFrame.
    Supports: Year/year, Gender/Sex/gender/sex, AgeBin/age_bin/agebin/age,
    Count/count, lb/LB, ub/UB.

    Gender values are normalized to integers: 0 = female, 1 = male.
    Accepts: 0/1, 'f'/'m', 'female'/'male' (case-insensitive).

    Returns a copy of the DataFrame with canonical column names and normalized
    gender values.
    """
    rename = {}
    col_lower = {c.lower(): c for c in data.columns}

    # Year
    for alias in ['year']:
        if alias in col_lower:
            rename[col_lower[alias]] = 'Year'
            break
    else:
        raise ValueError(f'Stratified coverage DataFrame must have a "Year" column. Found: {list(data.columns)}')

    # Sex/Gender (optional — age-only stratification is valid, e.g. for VMMC)
    has_sex = False
    for alias in ['gender', 'sex']:
        if alias in col_lower:
            rename[col_lower[alias]] = 'Gender'
            has_sex = True
            break

    # Age bin
    for alias in ['agebin', 'age_bin', 'age']:
        if alias in col_lower:
            rename[col_lower[alias]] = 'AgeBin'
            break
    else:
        raise ValueError(f'Stratified coverage DataFrame must have an "AgeBin" or "Age" column. Found: {list(data.columns)}')

    # Optional columns
    for alias, canonical in [('count', 'Count'), ('lb', 'lb'), ('ub', 'ub')]:
        if alias in col_lower:
            rename[col_lower[alias]] = canonical

    df = data.rename(columns=rename)

    # Normalize gender values to integers: 0=female, 1=male
    if has_sex:
        df['Gender'] = df['Gender'].map(_normalize_sex)

    return df


def _parse_stratified_df(data, yearvec, smoothness=0, **interp_kw):
    """
    Parse a stratified coverage DataFrame.

    Accepts flexible column names (Year/year, Gender/Sex, AgeBin/age_bin/age).
    Gender/Sex is optional — if absent, data is stratified by age only and the
    coverage dict is keyed by age_bin string; otherwise by (age_bin, sex) tuple.

    The value column is detected as the first numeric column that's not a
    metadata column (Year, Gender, AgeBin, Count, lb, ub).

    Returns (coverage_dict, format, age_bins, sex_keys) where sex_keys is None
    if no Gender column is present.
    """
    data = _normalize_stratified_cols(data)
    has_sex = 'Gender' in data.columns

    skip_cols = {'Year', 'Gender', 'AgeBin', 'Count', 'lb', 'ub'}
    val_col = None
    for col in data.columns:
        if col not in skip_cols and pd.api.types.is_numeric_dtype(data[col]):
            val_col = col
            break
    if val_col is None:
        raise ValueError(f'Could not find a numeric value column in stratified coverage DataFrame. Columns: {list(data.columns)}')

    fmt = 'p' if data[val_col].dropna().max() <= 1.0 else 'n'

    # Extract unique age bins (sort numerically, not lexicographically)
    age_bins = sorted(data['AgeBin'].unique(), key=lambda ab: ss.parse_age_range(ab)[0])
    sex_keys = sorted(data['Gender'].unique()) if has_sex else None

    # Build interpolated arrays for each stratum
    coverage = {}
    for ab in age_bins:
        for sex in (sex_keys or [None]):
            # Filter data for this stratum
            ab_mask = data['AgeBin'] == ab
            if sex is not None:
                subset = data[ab_mask & (data['Gender'] == sex)].sort_values('Year')
                key = (ab, sex)
            else:
                subset = data[ab_mask].sort_values('Year')
                key = ab

            if len(subset) == 0:
                coverage[key] = np.zeros(len(yearvec))
            else:
                years = subset['Year'].values.astype(float)
                vals  = subset[val_col].values.astype(float)
                # Fill NaN gaps before interpolating to yearvec
                mask = ~np.isnan(vals)
                if mask.any():
                    vals = sc.smoothinterp(years, years[mask], vals[mask], smoothness=0)
                coverage[key] = sc.smoothinterp(yearvec, years, vals, smoothness=smoothness, **interp_kw)

    return coverage, fmt, age_bins, sex_keys


def _handle_deprecated_coverage(coverage, coverage_data, kwargs):
    """
    Handle deprecated coverage_data and future_coverage kwargs.
    Returns the normalized coverage input and any remaining future_coverage.
    """
    future_coverage = kwargs.pop('future_coverage', None)

    if coverage_data is not None and coverage is not None:
        raise ValueError('Cannot specify both coverage and coverage_data. Use coverage only.')

    if coverage_data is not None:
        warnings.warn('coverage_data is deprecated; use coverage instead', FutureWarning, stacklevel=3)
        coverage = coverage_data

    if future_coverage is not None:
        warnings.warn('future_coverage is deprecated; extend coverage data to cover the full simulation period instead', FutureWarning, stacklevel=3)

    return coverage, future_coverage


# %% Coverage targeting helpers

def age_sex_mask(age_bin, sex, people):
    """
    Return a BoolArr mask for agents matching an age/sex stratum.

    Args:
        age_bin: age range string (e.g. '[15,25)')
        sex: 0 (female), 1 (male), or None (no sex filter)
        people: sim.people

    Returns:
        BoolArr: True for agents in the stratum
    """
    mask = ss.apply_age_range(age_bin, people.age)
    if sex is not None:
        sex_mask = people.female if sex == 0 else people.male
        mask = mask & sex_mask
    return mask


def coverage_to_number(cov_val, coverage_format, pop_scale=None, n_eligible=None):
    """
    Convert a coverage value to a target count.

    For 'n' (absolute) format, scales down by pop_scale to match the
    model population. For 'p' (proportion) format, multiplies by the
    number of eligible agents.

    Args:
        cov_val: coverage value (absolute number or proportion)
        coverage_format: 'n' (absolute) or 'p' (proportion)
        pop_scale: population scale factor (required for 'n' format)
        n_eligible: number of eligible agents (required for 'p' format)

    Returns:
        int: target count for the model population
    """
    if coverage_format == 'n':
        return int(cov_val / pop_scale)
    else:
        return int(cov_val * n_eligible)


def compute_coverage_target(coverage, coverage_format, age_bins, sex_keys,
                            ti, eligible_uids, sim, future_coverage=None):
    """
    Compute the target number of people to cover this timestep.

    Handles aggregate, stratified, and legacy future_coverage modes.
    Returns None if no coverage constraint, int otherwise.

    Shared by ART, VMMC, and potentially other coverage-targeted interventions.
    """
    # Legacy future_coverage mode (deprecated, will be removed)
    if future_coverage is not None and sim.t.now('year') >= future_coverage['year']:
        return int(future_coverage['prop'] * len(eligible_uids))

    # No coverage data — no constraint
    if coverage is None:
        return None

    # Stratified coverage: loop over age/sex strata and sum targets
    if isinstance(coverage, dict):
        total = 0
        for ab in age_bins:
            for sex in (sex_keys or [None]):
                # Look up the coverage value for this stratum
                key = (ab, sex) if sex is not None else ab
                cov = coverage.get(key, np.zeros(1))
                cov_val = cov[ti] if len(cov) > ti else cov[-1]

                # Count eligible agents and convert to a target number
                n = (age_sex_mask(ab, sex, sim.people) & eligible_uids).count()
                total += coverage_to_number(cov_val, coverage_format,
                                            pop_scale=sim.pars.pop_scale, n_eligible=n)
        return total

    # Aggregate coverage: single value for the whole population
    cov_val = coverage[ti]
    n_eligible = len(eligible_uids)
    return coverage_to_number(cov_val, coverage_format,
                              pop_scale=sim.pars.pop_scale, n_eligible=n_eligible)
