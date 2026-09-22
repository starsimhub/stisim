"""Data-preparation utilities for stisim country studies.

Currently: HIV-deleted background mortality construction, adapted from
Adam Akullian's ``mortality_construction`` in the ``hivsim_eswatini``
repository (github.com/starsimhub/hivsim_eswatini). The construction,
the argument for why it is needed, and the interpolation methods are
Adam's work. This module drops the Eswatini-specific file plumbing and
exposes the reusable core so any country study can adopt it.

HIV-deleted background mortality
================================

Problem
-------
All-cause mortality data (UN WPP and similar) includes AIDS deaths. When
stisim feeds these rates to ``ss.Deaths``, agents die from AIDS-attributable
causes via the demographics module and separately via the HIV module
(``p_hiv_death`` and the ``ti_zero`` AIDS pathway) — double-counting HIV
mortality and inflating apparent impact.

``dedup_deaths`` constructs a non-AIDS counterfactual by log-linear
interpolation between a pre-epidemic anchor year and a post-epidemic
anchor year, per age and sex. The interpolated line carries the genuine
secular trend (falling child mortality and so on); only the observed
excess above it is treated as AIDS and removed.

Stated assumptions
------------------
- The base year is AIDS-free (very nearly true for anchor years like 1985
  in most high-HIV-burden African countries; HIV prevalence was still low).
- The end year is AIDS-free (not quite — residual AIDS mortality persists,
  so this slightly under-deletes).
- Non-AIDS mortality moved smoothly between the endpoints. Where it did
  not, the residual gets misattributed to AIDS.
- Only years strictly between the endpoints are modified. Rates at and
  beyond the end year (typically projections with no AIDS hump) are left
  untouched.

Typical usage
-------------
::

    import pandas as pd
    from stisim.data.data_utils import dedup_deaths, deleted_fraction

    all_cause = pd.read_csv('data/<country>_deaths.csv')
    hiv_deleted = dedup_deaths(all_cause, base_year=1985, end_year=2025)
    diagnostic = deleted_fraction(all_cause, hiv_deleted)

    # Write the HIV-deleted version under the name stisim reads; preserve
    # the original all-cause file under a name stisim will not pick up.
    all_cause.to_csv('data/<country>_deaths_all_cause.csv', index=False)
    hiv_deleted.to_csv('data/<country>_deaths.csv', index=False)

The ``diagnostic`` DataFrame carries per-cell AIDS shares — a structural
sanity check that does not require running the model. In practice the
AIDS share of all-cause mortality should peak in the mid-30s (roughly
70–85% depending on country and sex) and fall to essentially zero above
age 80. Sanity-checked against Zimbabwe (2005 peak year: ~82% female age
31, ~71% male age 36) and Eswatini (2005 peak year: ~84% mid-30s).

Input schema
------------
The input DataFrame must have columns ``Time``, ``Sex``, ``AgeStart``,
``Value`` — the shape used by stisim's demographic loaders. Rows at the
base year and end year define the anchors; rows strictly between them
are candidates for deletion. Cells where the observed rate is already at
or below the counterfactual are left untouched.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

__all__ = [
    'dedup_deaths',
    'deleted_fraction',
    'make_datafolder',
    'percentiles_to_pars',
    'logn_percentiles_to_pars',
]

DEFAULT_BASE_YEAR = 1985
DEFAULT_END_YEAR = 2025


def _warp(u: np.ndarray, method: str, par: float) -> np.ndarray:
    """Map elapsed fraction ``u`` in [0,1] to interpolation weight in [0,1].

    The weight decides how much of the decline has happened by a given year,
    and therefore how much of the observed rate is called AIDS. It is the
    single most consequential assumption in the construction, and it is
    *not* identifiable from this data — there are only two AIDS-free
    anchors, with every year between them contaminated.

    Methods
    -------
    - ``loglinear`` / ``linear``: ``w = u``, constant proportional or
      absolute change.
    - ``power``: ``w = u ** par``. ``par > 1`` delays the decline (the
      counterfactual stays high mid-period, so LESS is attributed to AIDS);
      ``par < 1`` front-loads it.
    - ``sigmoid``: symmetric S-curve, slow-fast-slow. A no-op at the
      midpoint year — a symmetric sigmoid passes through the midpoint at
      ``u = 0.5``, exactly where linear does. It changes the shoulders,
      not the peak.
    """
    if method in ('loglinear', 'linear'):
        return u
    if method == 'power':
        return u ** par
    if method == 'sigmoid':
        f = lambda x: 1.0 / (1.0 + np.exp(-par * (x - 0.5)))
        return (f(u) - f(0.0)) / (f(1.0) - f(0.0))
    raise ValueError(f'unknown method: {method}')


def dedup_deaths(deaths_df: pd.DataFrame,
                 base_year: int = DEFAULT_BASE_YEAR,
                 end_year: int = DEFAULT_END_YEAR,
                 method: str = 'loglinear',
                 par: float = 1.0) -> pd.DataFrame:
    """Return a copy of ``deaths_df`` with the AIDS hump removed.

    Rates strictly between ``base_year`` and ``end_year`` are replaced by
    an interpolation of the two anchors, per ``(Sex, AgeStart)``. Rates are
    only ever lowered — if the observed rate is already at or below the
    counterfactual it is kept, so bins AIDS never touched are unchanged.

    ``method='loglinear'`` (the default) interpolates geometrically:
    constant proportional change (exponential decay). ``method='linear'``
    interpolates arithmetically, which sits above the geometric curve
    mid-period and so attributes less to AIDS.
    """
    df = deaths_df.copy()
    anchors = df[df.Time.isin([base_year, end_year])]
    lo = anchors[anchors.Time == base_year].set_index(['Sex', 'AgeStart'])['Value']
    hi = anchors[anchors.Time == end_year].set_index(['Sex', 'AgeStart'])['Value']

    mid = df.Time.between(base_year, end_year, inclusive='neither')
    idx = pd.MultiIndex.from_frame(df.loc[mid, ['Sex', 'AgeStart']])
    u = ((df.loc[mid, 'Time'] - base_year) / (end_year - base_year)).values
    w = _warp(u, method, par)

    a = lo.reindex(idx).values
    b = hi.reindex(idx).values
    with np.errstate(divide='ignore', invalid='ignore'):
        if method == 'linear':
            cf = (1 - w) * a + w * b
        else:
            cf = np.exp((1 - w) * np.log(a) + w * np.log(b))
    cf = np.where(np.isfinite(cf), cf, df.loc[mid, 'Value'].values)

    # Never raise a rate: only the excess above the counterfactual is AIDS.
    df.loc[mid, 'Value'] = np.minimum(df.loc[mid, 'Value'].values, cf)
    return df


def deleted_fraction(deaths_df: pd.DataFrame,
                     hiv_deleted: pd.DataFrame) -> pd.DataFrame:
    """Per row, how much mortality was removed — the audit trail.

    Returns a DataFrame with the original and HIV-deleted rates side by
    side, plus a ``deleted_rate`` column (what was removed) and an
    ``aids_share`` column (deleted rate as a fraction of all-cause). Bins
    with a large deleted fraction are the ones the construction is doing
    the most work on and where its assumptions are most consequential.
    """
    m = deaths_df.merge(hiv_deleted, on=['Time', 'Sex', 'AgeStart'],
                        suffixes=('_all_cause', '_non_aids'))
    m['deleted_rate'] = m['Value_all_cause'] - m['Value_non_aids']
    m['aids_share'] = np.where(m['Value_all_cause'] > 0,
                               m['deleted_rate'] / m['Value_all_cause'], 0.0)
    return m


def make_datafolder(src: Path, dest: Path, hiv_deleted: pd.DataFrame,
                    deaths_filename: str) -> Path:
    """Create an alternate datafolder with HIV-deleted mortality.

    Copies every CSV from ``src`` so stisim finds whatever demographic
    input it asks for, then overwrites the deaths file with the HIV-deleted
    version. Useful when the user wants to run with and without HIV-
    deletion side by side without editing the primary data directory.
    """
    dest.mkdir(parents=True, exist_ok=True)
    for f in src.glob('*.csv'):
        shutil.copy2(f, dest / f.name)
    hiv_deleted.to_csv(dest / deaths_filename, index=False)
    return dest


# ---------------------------------------------------------------------------
# Distribution fitting from empirical percentiles
# ---------------------------------------------------------------------------
#
# Sexual-behaviour data from surveys (DHS, PHIA, IBBS) is often reported as
# cumulative fractions by exact age — e.g. "fraction of women who have had
# first sexual intercourse by age 15, 18, 20, 22, 25". Two such (value,
# fraction) points identify a two-parameter distribution, so the pair can
# be inverted to distribution parameters that stisim can consume.


def percentiles_to_pars(x1, p1, x2, p2):
    """Find (location, scale) of a normal distribution given two quantiles.

    Solves ``P(X < x1) = p1`` and ``P(X < x2) = p2`` for a normal
    distribution. Returns ``(location, scale)`` — the mean and standard
    deviation of the fitted normal. Feed those to ``scipy.stats.norm``
    or, if the distribution is not skewed, directly to a normal-shaped
    stisim distribution.
    """
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    location = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = (x2 - x1) / (p2ppf - p1ppf)
    return location, scale


def logn_percentiles_to_pars(x1, p1, x2, p2):
    """Find (s, scale) of a lognormal distribution given two quantiles.

    Solves ``P(X < x1) = p1`` and ``P(X < x2) = p2`` for a lognormal
    distribution. Returns ``(s, scale)`` — the shape (standard deviation
    of the underlying normal) and scale (``exp`` of its mean) parameters
    used by ``scipy.stats.lognorm``.

    To convert to the ``(mean, std)`` form used by ``ss.lognorm_ex``::

        s, scale = logn_percentiles_to_pars(x1, p1, x2, p2)
        mean = scale * np.exp(s ** 2 / 2)
        std = mean * np.sqrt(np.exp(s ** 2) - 1)
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    s = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = np.exp(mean)
    return s, scale
