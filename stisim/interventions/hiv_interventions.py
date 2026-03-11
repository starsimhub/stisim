"""
Define HIV interventions for STIsim
By default, these all have units of a year and timesteps of 1/12
"""

import warnings
import starsim as ss
import numpy as np
import pandas as pd
import sciris as sc
from stisim.interventions.base_interventions import STITest


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


def parse_coverage(data, valid_names=None, yearvec=None):
    """
    Parse coverage data into a per-timestep array.

    Accepts multiple input formats and normalizes them into
    (coverage_array, coverage_format, age_bins, sex_keys) where:
        - coverage_array: 1D array (len=yearvec) for aggregate, or dict of arrays for stratified
        - coverage_format: 'n' (absolute numbers) or 'p' (proportion)
        - age_bins: list of (lo, hi) tuples if stratified, else None
        - sex_keys: list of sex values if stratified, else None

    Supported input formats:
        None                      → (None, None, None, None)
        0.9                       → constant proportion
        {'year': [...], 'value': [...]}                          → interpolated, infer format
        {'year': [...], 'value': [...], 'format': 'n'}          → interpolated, explicit format
        pd.DataFrame(index=years, columns=['n_art'])             → legacy single-column
        pd.DataFrame(columns=['Year','Gender','AgeBin',...])     → age/sex stratified

    Args:
        data:        coverage input in any supported format
        valid_names: list of valid column names, e.g. ['n_art', 'p_art']
        yearvec:     simulation year vector for interpolation
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
        arr = sc.smoothinterp(yearvec, years, values, smoothness=0)
        return arr, fmt, None, None

    # DataFrame — check for stratified vs legacy
    if isinstance(data, pd.DataFrame):
        return _parse_coverage_df(data, valid_names, yearvec)

    errormsg = f'Coverage data format not recognized: {type(data)}. Expected None, number, dict, or DataFrame.'
    raise ValueError(errormsg)


def _parse_coverage_df(data, valid_names, yearvec):
    """
    Parse a DataFrame of coverage data.

    Handles two cases:
        1. Legacy: index=years, single column in valid_names → 1D interpolated array
        2. Stratified: columns include Year, Gender/Sex, AgeBin, and a value column → dict of arrays by (age, sex)
    """

    # Check for stratified format — normalize column names for flexible matching
    col_lower = {c.lower(): c for c in data.columns}
    has_year   = 'year' in col_lower
    has_sex    = 'gender' in col_lower or 'sex' in col_lower
    has_agebin = 'agebin' in col_lower or 'age_bin' in col_lower or 'age' in col_lower

    if has_year and has_sex and has_agebin:
        return _parse_stratified_df(data, yearvec)

    # Legacy single-column format: index=years, column in valid_names
    if len(data.columns) == 1 and data.columns[0] in valid_names:
        colname = data.columns[0]
        fmt = 'n' if colname.startswith('n_') else 'p'
        arr = sc.smoothinterp(yearvec, data.index.values, data[colname].values, smoothness=0)
        return arr, fmt, None, None

    # Try to find a valid column
    for col in data.columns:
        if col in valid_names:
            fmt = 'n' if col.startswith('n_') else 'p'
            if data.index.name in ['year', 'Year'] or np.all(data.index > 1900):
                arr = sc.smoothinterp(yearvec, data.index.values, data[col].values, smoothness=0)
            else:
                arr = sc.smoothinterp(yearvec, np.arange(len(data)), data[col].values, smoothness=0)
            return arr, fmt, None, None

    errormsg = f'DataFrame columns {list(data.columns)} do not match any valid names {valid_names}.'
    raise ValueError(errormsg)


def _normalize_stratified_cols(data):
    """
    Normalize column names in a stratified DataFrame to canonical form.
    Supports: Year/year, Gender/Sex/gender/sex, AgeBin/age_bin/agebin/age,
    Count/count, lb/LB, ub/UB.

    Returns a copy of the DataFrame with canonical column names, and raises
    ValueError if required columns (year, sex/gender, age bin) are missing.
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

    # Sex/Gender
    for alias in ['gender', 'sex']:
        if alias in col_lower:
            rename[col_lower[alias]] = 'Gender'
            break
    else:
        raise ValueError(f'Stratified coverage DataFrame must have a "Gender" or "Sex" column. Found: {list(data.columns)}')

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
    return df


def _parse_stratified_df(data, yearvec):
    """
    Parse a stratified coverage DataFrame.

    Accepts flexible column names (Year/year, Gender/Sex, AgeBin/age_bin/age).
    The value column is detected as the first numeric column that's not a
    metadata column (Year, Gender, AgeBin, Count, lb, ub).

    Returns a dict keyed by (age_bin, sex) with interpolated arrays, plus
    the list of age_bins and sex_keys for iteration.
    """
    data = _normalize_stratified_cols(data)

    skip_cols = {'Year', 'Gender', 'AgeBin', 'Count', 'lb', 'ub'}
    val_col = None
    for col in data.columns:
        if col not in skip_cols and pd.api.types.is_numeric_dtype(data[col]):
            val_col = col
            break
    if val_col is None:
        raise ValueError(f'Could not find a numeric value column in stratified coverage DataFrame. Columns: {list(data.columns)}')

    fmt = 'p' if data[val_col].dropna().max() <= 1.0 else 'n'

    # Extract unique age bins and sex keys
    age_bins = sorted(data['AgeBin'].unique())
    sex_keys = sorted(data['Gender'].unique())

    # Build interpolated arrays for each (age_bin, sex) combination
    coverage = {}
    for ab in age_bins:
        for sex in sex_keys:
            subset = data[(data['AgeBin'] == ab) & (data['Gender'] == sex)].sort_values('Year')
            if len(subset) == 0:
                coverage[(ab, sex)] = np.zeros(len(yearvec))
            else:
                years = subset['Year'].values.astype(float)
                vals  = subset[val_col].values.astype(float)
                mask = ~np.isnan(vals)
                if mask.any():
                    vals = np.interp(years, years[mask], vals[mask])
                coverage[(ab, sex)] = sc.smoothinterp(yearvec, years, vals, smoothness=0)

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


# %% HIV classes
__all__ = ["HIVDx", "HIVTest", "ART", "VMMC", "Prep"]


class HIVDx(ss.Product):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.result_list = ['positive', 'negative']

    def administer(self, sim, uids):
        outcomes = {r: ss.uids() for r in self.result_list}
        outcomes['positive'] = sim.diseases.hiv.infected.uids.intersect(uids)
        outcomes['negative'] = sim.diseases.hiv.susceptible.uids.intersect(uids)
        return outcomes


class HIVTest(STITest):
    """
    Base class for HIV testing
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None, eligibility=None, name=None, label=None, **kwargs):
        if product is None: product = HIVDx(name=f'HIVDx_{name}')
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years, start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = self.ti
        return outcomes


class ART(ss.Intervention):
    """
    ART-treatment intervention.

    Coverage can be specified as:
        - A scalar (constant proportion of infected)
        - A dict with 'year' and 'value' keys (interpolated over time)
        - A DataFrame with index=years, single column 'n_art' or 'p_art' (legacy)
        - A DataFrame with Year/Gender/AgeBin columns (age/sex stratified)

    Args:
        coverage:         coverage data in any supported format
        art_initiation:   probability that a newly diagnosed person initiates ART (ss.bernoulli)
        coverage_data:    deprecated, use coverage instead
    """

    def __init__(self, pars=None, coverage=None, coverage_data=None, **kwargs):
        super().__init__()

        # Handle deprecated kwargs
        coverage, future_coverage = _handle_deprecated_coverage(coverage, coverage_data, kwargs)

        # Handle deprecated init_prob → art_initiation
        init_prob = kwargs.pop('init_prob', None)
        if init_prob is not None:
            warnings.warn('init_prob is deprecated; use art_initiation instead', FutureWarning, stacklevel=2)

        self.define_pars(
            art_initiation=init_prob if init_prob is not None else ss.bernoulli(p=0.9),
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._future_coverage = future_coverage  # Legacy compat, removed once deprecated
        self.coverage         = None  # Set in init_pre
        self.coverage_format  = None  # 'n' or 'p'
        self.age_bins         = None  # For stratified coverage
        self.sex_keys         = None  # For stratified coverage
        return

    def init_pre(self, sim):
        super().init_pre(sim)

        # Warn if no HIV testing intervention found
        has_hiv_test = any(isinstance(m, HIVTest) for m in sim.interventions())
        if not has_hiv_test:
            ss.warn('ART intervention added without an HIV testing intervention; diagnosed agents may never be identified')

        # Parse coverage data
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parse_coverage(
            self._raw_coverage, valid_names=['n_art', 'p_art'], yearvec=self.t.yearvec,
        )
        self.initialized = True
        return

    def _get_n_to_treat(self, eligible_uids):
        """
        Get the target number of people on ART this timestep.

        Returns None if no coverage data is provided (meaning: no capacity
        constraint, just treat everyone who initiates). Returns an int otherwise.
        """
        # Legacy future_coverage mode
        if self._future_coverage is not None and self.t.now('year') >= self._future_coverage['year']:
            p_cov = self._future_coverage['prop']
            return int(p_cov * len(eligible_uids))

        if self.coverage is None:
            return None

        # Stratified coverage
        if isinstance(self.coverage, dict):
            return self._get_n_to_treat_stratified(eligible_uids)

        # Aggregate coverage
        if self.coverage_format == 'n':
            return int(self.coverage[self.ti] / self.sim.pars.pop_scale)
        else:
            return int(self.coverage[self.ti] * len(eligible_uids))

    def _get_n_to_treat_stratified(self, eligible_uids):
        """
        Compute target ART numbers by age/sex stratum.

        Returns the total target across all strata.
        """
        sim = self.sim
        total = 0

        for ab in self.age_bins:
            lo, hi = _parse_age_bin(ab)
            for sex in self.sex_keys:
                cov = self.coverage.get((ab, sex), np.zeros(1))
                cov_val = cov[self.ti] if len(cov) > self.ti else cov[-1]

                # Find eligible agents in this stratum
                age_mask = (sim.people.age >= lo) & (sim.people.age < hi)
                sex_mask = sim.people.female if sex == 1 else sim.people.male
                stratum_uids = eligible_uids[age_mask[eligible_uids] & sex_mask[eligible_uids]]

                if self.coverage_format == 'n':
                    total += int(cov_val / sim.pars.pop_scale)
                else:
                    total += int(cov_val * len(stratum_uids))

        return total

    def step(self):
        """
        Apply ART at each timestep: stop ART for those scheduled, initiate for newly diagnosed,
        and correct overall coverage to match targets.
        """
        sim = self.sim
        hiv = sim.diseases.hiv
        inf_uids = hiv.infected.uids

        # Determine treatment target (None = no capacity constraint)
        n_to_treat = self._get_n_to_treat(inf_uids)

        # Check who is stopping ART
        if hiv.on_art.any():
            stopping = hiv.on_art & (hiv.ti_stop_art <= self.ti)
            if stopping.any():
                try:
                    hiv.stop_art(stopping.uids)
                except:
                    errormsg = f'Error stopping ART for {stopping.uids}'
                    raise ValueError(errormsg)

        # Initiate ART for newly diagnosed
        diagnosed = hiv.ti_diagnosed == self.ti
        if len(diagnosed.uids):
            dx_to_treat = self.pars.art_initiation.filter(diagnosed.uids)

            if n_to_treat is None:
                # No coverage target — treat all who initiate
                hiv.start_art(dx_to_treat)
            else:
                # Coverage target — only treat if spots available
                on_art = hiv.on_art
                n_available_spots = n_to_treat - len(on_art.uids)
                if n_available_spots > 0:
                    self.prioritize_art(sim, n=n_available_spots, awaiting_art_uids=dx_to_treat)

        # Correct coverage to match target (only if target is set)
        if n_to_treat is not None:
            self.art_coverage_correction(sim, target_coverage=n_to_treat)

        # Adjust rel_sus for protected unborn agents (only if pregnancy is modeled)
        if hasattr(sim.people, 'pregnancy') and hasattr(sim.networks, 'maternalnet'):
            if hiv.on_art[sim.people.pregnancy.pregnant].any():
                mother_uids = (hiv.on_art & sim.people.pregnancy.pregnant).uids
                infants = sim.networks.maternalnet.find_contacts(mother_uids)
                hiv.rel_sus[ss.uids(infants)] = 0

        return

    def prioritize_art(self, sim, n=None, awaiting_art_uids=None):
        """ Prioritize ART to n agents among those awaiting treatment """
        hiv = sim.diseases.hiv
        if awaiting_art_uids is None:
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids

        # Enough spots for everyone
        if n > len(awaiting_art_uids):
            start_uids = awaiting_art_uids

        # Not enough spots — prioritize by CD4 and care seeking
        else:
            cd4_counts   = hiv.cd4[awaiting_art_uids]
            care_seeking = hiv.care_seeking[awaiting_art_uids]
            weights = cd4_counts * (1 / care_seeking)
            choices = np.argsort(weights)[:n]
            start_uids = awaiting_art_uids[choices]

        hiv.start_art(start_uids)

        return

    def art_coverage_correction(self, sim, target_coverage=None):
        """ Adjust ART coverage to match data """
        hiv = sim.diseases.hiv
        on_art = hiv.on_art

        # Too many on treatment → remove
        if len(on_art.uids) > target_coverage:
            n_to_stop    = int(len(on_art.uids) - target_coverage)
            on_art_uids  = on_art.uids
            cd4_counts   = hiv.cd4[on_art_uids]
            care_seeking = hiv.care_seeking[on_art_uids]
            weights  = cd4_counts / care_seeking
            choices  = np.argsort(-weights)[:n_to_stop]
            stop_uids = on_art_uids[choices]
            hiv.ti_stop_art[stop_uids] = self.ti
            hiv.stop_art(stop_uids)

        # Not enough on treatment → add
        elif len(on_art.uids) < target_coverage:
            n_to_add = target_coverage - len(on_art.uids)
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids
            self.prioritize_art(sim, n=n_to_add, awaiting_art_uids=awaiting_art_uids)

        return


class VMMC(ss.Intervention):
    """
    Voluntary medical male circumcision intervention.

    Coverage can be specified as:
        - A scalar (constant proportion of males)
        - A dict with 'year' and 'value' keys (interpolated over time)
        - A DataFrame with index=years, single column 'n_vmmc' or 'p_vmmc' (legacy)
        - A DataFrame with Year/Gender/AgeBin columns (age stratified)

    Args:
        coverage:    coverage data in any supported format
        eff_circ:    efficacy of circumcision (reduction in HIV acquisition risk)
        eligibility: optional eligibility function or filter
        coverage_data: deprecated, use coverage instead
    """

    def __init__(self, pars=None, coverage=None, coverage_data=None, eligibility=None, **kwargs):
        super().__init__()

        # Handle deprecated kwargs
        coverage, future_coverage = _handle_deprecated_coverage(coverage, coverage_data, kwargs)

        self.define_pars(
            eff_circ=0.6,
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._future_coverage = future_coverage
        self.coverage         = None
        self.coverage_format  = None
        self.age_bins         = None
        self.sex_keys         = None
        self.eligibility      = eligibility

        # States
        self.willingness     = ss.FloatArr('willingness', default=ss.random())
        self.circumcised     = ss.BoolArr('circumcised', default=False)
        self.ti_circumcised  = ss.FloatArr('ti_circumcised')

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parse_coverage(
            self._raw_coverage, valid_names=['n_vmmc', 'p_vmmc'], yearvec=self.t.yearvec,
        )
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_circumcisions', dtype=int, label='New circumcisions', auto_plot=False),
            ss.Result('n_circumcised',     dtype=int, label='Number circumcised', auto_plot=False),
        )
        return

    def _get_n_to_circ(self, eligible_uids):
        """ Get the target number of circumcisions this timestep """
        # Legacy future_coverage mode
        if self._future_coverage is not None and self.t.now('year') >= self._future_coverage['year']:
            return int(self._future_coverage['prop'] * len(eligible_uids))

        if self.coverage is None:
            return None

        # Stratified coverage
        if isinstance(self.coverage, dict):
            return self._get_n_stratified(eligible_uids)

        # Aggregate coverage
        if self.coverage_format == 'n':
            return int(self.coverage[self.ti] / self.sim.pars.pop_scale)
        else:
            return int(self.coverage[self.ti] * len(eligible_uids))

    def _get_n_stratified(self, eligible_uids):
        """ Compute target circumcisions by age stratum """
        sim = self.sim
        total = 0
        for ab in self.age_bins:
            lo, hi = _parse_age_bin(ab)
            for sex in (self.sex_keys or [0]):  # VMMC is males only, but handle data gracefully
                cov = self.coverage.get((ab, sex), np.zeros(1))
                cov_val = cov[self.ti] if len(cov) > self.ti else cov[-1]
                age_mask = (sim.people.age >= lo) & (sim.people.age < hi)
                stratum_uids = eligible_uids[age_mask[eligible_uids]]
                if self.coverage_format == 'n':
                    total += int(cov_val / sim.pars.pop_scale)
                else:
                    total += int(cov_val * len(stratum_uids))
        return total

    def step(self):
        sim = self.sim
        m_uids = sim.people.male.uids

        n_to_circ = self._get_n_to_circ(m_uids)

        if n_to_circ is not None and n_to_circ > 0:
            eligible_uids = (sim.people.male & ~self.circumcised).uids
            weights = self.willingness[eligible_uids]
            choices = np.argsort(-weights)[:n_to_circ]
            new_circs = eligible_uids[choices]

            self.circumcised[new_circs] = True
            self.ti_circumcised[new_circs] = self.ti

        self.results['new_circumcisions'][self.ti] = n_to_circ or 0
        self.results['n_circumcised'][self.ti] = count(self.circumcised)

        # Reduce rel_sus
        sim.diseases.hiv.rel_sus[self.circumcised] *= 1 - self.pars.eff_circ

        return


class Prep(ss.Intervention):
    """ PrEP for FSW """
    def __init__(self, pars=None, eligibility=None, **kwargs):
        super().__init__()
        self.define_pars(
            coverage_dist=ss.bernoulli(p=0),
            coverage=[0, 0.01, 0.5, 0.8],
            years=[2004, 2005, 2015, 2025],
            eff_prep=0.8,
        )
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        self.define_states(
            ss.BoolArr('on_prep', label='On PrEP'),
        )
        return

    def step(self):
        sim = self.sim
        self.coverage = np.interp(self.t.yearvec, self.pars.years, self.pars.coverage)
        if self.coverage[self.ti] > 0:
            self.pars.coverage_dist.set(p=self.coverage[self.ti])
            el_fsw = self.sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected & ~self.on_prep
            fsw_on_prep = self.pars.coverage_dist.filter(el_fsw)
            self.sim.diseases.hiv.rel_sus[fsw_on_prep] *= 1 - self.pars.eff_prep

        return


# %% Utility functions
def _parse_age_bin(ab_str):
    """
    Parse age bin string like "[15,25)" into (lo, hi).
    Also handles numeric tuples and plain strings like "15-25".
    """
    if isinstance(ab_str, (tuple, list)):
        return float(ab_str[0]), float(ab_str[1])
    s = str(ab_str).strip('[]() ')
    if ',' in s:
        parts = s.split(',')
    elif '-' in s:
        parts = s.split('-')
    else:
        raise ValueError(f'Cannot parse age bin: {ab_str}')
    return float(parts[0].strip()), float(parts[1].strip().rstrip(')'))
