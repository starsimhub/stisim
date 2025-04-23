""" STIsim utilities """

import sciris as sc
import pandas as pd
import numpy as np
import starsim as ss

__all__ = ['TimeSeries', 'make_init_prev_fn', "finalize_results", "shrink_calib", "standardize_data"]


class TimeSeries:
    """
    Class to store time-series data

    Internally values are stored as lists rather than numpy arrays because
    insert/remove operations on lists tend to be faster (and working with sparse
    data is a key role of TimeSeries objects). Note that methods like :meth:`interpolate()`
    return numpy arrays, so the output types from such functions should generally match up
    with what is required by the calling function.

    :param t: Optionally specify a scalar, list, or array of time values
    :param vals: Optionally specify a scalar, list, or array of values (must be same size as ``t``)
    :param units: Optionally specify units (as a string)
    :param assumption: Optionally specify a scalar assumption
    :param sigma: Optionally specify a scalar uncertainty

    """

    # Use slots here to guarantee that __deepcopy__() and __eq__() only have to check these
    # specific fields - otherwise, would need to do a more complex recursive dict comparison
    __slots__ = ["t", "vals", "units", "assumption", "sigma", "_sampled"]

    def __init__(self, t=None, vals=None, units: str = None, assumption: float = None, sigma: float = None):

        self.t = []  #: Sorted array of time points. Normally interacted with via methods like :meth:`insert()`
        self.vals = []  #: Time-specific values - indices correspond to ``self.t``
        self.units = units  #: The units of the quantity
        self.assumption = assumption  #: The time-independent scalar assumption
        self.sigma = sigma  #: Uncertainty value, assumed to be a standard deviation
        self._sampled = False  #: Flag to indicate whether sampling has been performed. Once sampling has been performed, cannot sample again

        # Using insert() means that array/list inputs containing None or duplicate entries will
        # be sanitized via insert()
        self.insert(t, vals)

    def __repr__(self):
        output = sc.prepr(self)
        return output

    def __eq__(self, other):
        """
        Check TimeSeries equality

        Two TimeSeries instances are equal if all of their attributes are equal. This is easy to
        implement because `==` is directly defined for all of the attribute types (lists and scalars)
        and due to `__slots__` there are guaranteed not to be any other attributes

        :param other:
        :return:
        """

        return all(getattr(self, x) == getattr(other, x) for x in self.__slots__)

    def __deepcopy__(self, memodict={}):
        new = TimeSeries.__new__(TimeSeries)
        new.t = self.t.copy()
        new.vals = self.vals.copy()
        new.units = self.units
        new.assumption = self.assumption
        new.sigma = self.sigma
        new._sampled = self._sampled
        return new

    def __getstate__(self):
        return dict([(k, getattr(self, k, None)) for k in self.__slots__])

    def __setstate__(self, data):

        if "format" in data:
            # 'format' was changed to 'units' but the attribute was not dropped, however now this is a
            # hard error because of the switch to __slots__ so we need to make sure it gets removed.
            # This can't be done as a Migration because a Migration expects an instance of the
            data = sc.dcp(data)
            if "units" not in data:
                data["units"] = data["format"]
            del data["format"]

        for k, v in data.items():
            setattr(self, k, v)

    def copy(self):
        """
        Return a copy of the ``TimeSeries``

        :return: An independent copy of the ``TimeSeries``
        """

        return self.__deepcopy__(self)

    @property
    def has_data(self) -> bool:
        """
        Check if any data has been provided

        :return: ``True`` if any data has been entered (assumption or time-specific)

        """
        return self.assumption is not None or self.has_time_data

    @property
    def has_time_data(self) -> bool:
        """
        Check if time-specific data has been provided

        Unlike ``has_data``, this will return ``False`` if only an assumption has been entered

        :return: ``True`` if any time-specific data has been entered

        """
        # Returns true if any time-specific data has been entered (not just an assumption)
        return len(self.t) > 0

    def insert(self, t, v) -> None:
        """
        Insert a value or list of at a particular time

        If the value already exists in the ``TimeSeries``, it will be overwritten/updated.
        The arrays are internally sorted by time value, and this order will be maintained.

        :param t: Time value to insert or update. If ``None``, the value will be assigned to the assumption
        :param v: Value to insert. If ``None``, this function will return immediately without doing anything

        """

        # Check if inputs are iterable
        iterable_input = True
        try:
            assert len(t) == len(v), "Cannot insert non-matching lengths or types of time and values %s and %s" % (t, v)
        except TypeError:
            iterable_input = False

        # If inputs are iterable, call insert() for each zipped item
        if iterable_input:
            for ti, vi in zip(t, v):
                self.insert(ti, vi)
            return

        if v is None:  # Can't cast a None to a float, so just skip it
            return

        v = float(v)  # Convert input to float

        if t is None:  # Store the value in the assumption
            self.assumption = v
            return

        idx = bisect_left(self.t, t)
        if idx < len(self.t) and self.t[idx] == t:
            # Overwrite an existing entry
            self.vals[idx] = v
        else:
            self.t.insert(idx, t)
            self.vals.insert(idx, v)

    def get(self, t) -> float:
        """
        Retrieve value at a particular time

        This function will automatically retrieve the value of the assumption if
        no time specific values have been provided, or if any time specific values
        are provided, will return the value entered at that time. If time specific
        values have been entered and the requested time is not explicitly present,
        an error will be raised.

        This function may be deprecated in future because generally it is more useful
        to either call ``TimeSeries.interpolate()`` if interested in getting values at
        arbitrary times, or ``TimeSeries.get_arrays()`` if interested in retrieving
        values that have been entered.

        :param t: A time value. If ``None``, will return assumption regardless of whether
                  time data has been entered or not
        :return: The value at the corresponding time. Returns None if the value no value present
        """

        if t is None or len(self.t) == 0:
            return self.assumption
        elif t in self.t:
            return self.vals[self.t.index(t)]
        else:
            return None

    def get_arrays(self):
        """
        Return arrays with the contents of this TimeSeries

        The TimeSeries instance may have time values, or may simply have
        an assumption. If obtaining raw arrays is desired, this function will
        return arrays with values extracted from the appropriate attribute of the
        TimeSeries. However, in general, it is usually `.interpolate()` that is
        desired, rather than `.get_arrays()`

        :return: Tuple with two arrays - the first item is times (with a single NaN if
                 the TimeSeries only has an assumption) and the second item is values

        """
        if len(self.t) == 0:
            t = np.array([np.nan])
            v = np.array([self.assumption])
        else:
            t = np.array(self.t)
            v = np.array(self.vals)
        return t, v

    def remove(self, t) -> None:
        """
        Remove single time point

        :param t: Time value to remove. Set to ``None`` to remove the assumption

        """
        # To remove the assumption, set t=None
        if t is None:
            self.assumption = None
        elif t in self.t:
            idx = self.t.index(t)
            del self.t[idx]
            del self.vals[idx]
        else:
            raise Exception("Item not found")

    def remove_before(self, t_remove) -> None:
        """
        Remove times from start

        :param tval: Remove times up to but not including this time

        """

        for tval in sc.dcp(self.t):
            if tval < t_remove:
                self.remove(tval)

    def remove_after(self, t_remove) -> None:
        """
        Remove times from start

        :param tval: Remove times up to but not including this time

        """

        for tval in sc.dcp(self.t):
            if tval > t_remove:
                self.remove(tval)

    def remove_between(self, t_remove) -> None:
        """
        Remove a range of times

        Note that the endpoints are not included

        :param t_remove: two element iterable e.g. array, with [min,max] times

        """

        for tval in sc.dcp(self.t):
            if t_remove[0] < tval < t_remove[1]:
                self.remove(tval)

    def interpolate(self, t2: np.array, method="linear", **kwargs) -> np.array:
        """
        Return interpolated values

        This method returns interpolated values from the time series at time points `t2`
        according to a given interpolation method. There are 4 possibilities for the method

        - 'linear' - normal linear interpolation (with constant, zero-gradient extrapolation)
        - 'pchip' - legacy interpolation with some curvature between points (with constant, zero-gradient extrapolation)
        - 'previous' - stepped interpolation, maintain value until the next timepoint is reached (with constant, zero-gradient extrapolation)
        - Interpolation class or generator function

        That final option allows the use of arbitrary interpolation methods. The underlying call will be::

            c = method(t1, v1, **kwargs)
            return c(t2)

        so for example, if you wanted to use the base Scipy pchip method with no extrapolation, then could pass in::

            TimeSeries.interpolate(...,method=scipy.interpolate.PchipInterpolator)

        Note that the following special behaviours apply:

        - If there is no data at all, this function will return ``np.nan`` for all requested time points
        - If only an assumption exists, this assumption will be returned for all requested time points
        - Otherwise, arrays will be formed with all finite time values
        
            - If no finite time values remain, an error will be raised (in general, a TimeSeries should not store such values anyway)
            - If only one finite time value remains, then that value will be returned for all requested time points
            - Otherwise, the specified interpolation method will be used

        :param t2: float, list, or array, with times
        :param method: A string 'linear', 'pchip' or 'previous' OR a callable item that returns an Interpolator
        :return: array the same length as t2, with interpolated values

        """

        t2 = sc.promotetoarray(t2)  # Deal with case where user prompts for single time point

        # Deal with not having time-specific data first
        if not self.has_data:
            return np.full(t2.shape, np.nan)
        elif not self.has_time_data:
            return np.full(t2.shape, self.assumption)

        # Then, deal with having only 0 or 1 valid time points
        t1, v1 = self.get_arrays()
        idx = ~np.isnan(t1) & ~np.isnan(v1)
        t1, v1 = t1[idx], v1[idx]
        if t1.size == 0:
            raise Exception("No time points remained after removing NaNs from the TimeSeries")
        elif t1.size == 1:
            return np.full(t2.shape, v1[0])

        # # Finally, perform interpolation
        if sc.isstring(method):
            if method == "linear":
                # Default linear interpolation
                return np.interp(t2, t1, v1, left=v1[0], right=v1[-1])
            elif method == "pchip":
                # Legacy pchip interpolation
                f = scipy.interpolate.PchipInterpolator(t1, v1, axis=0, extrapolate=False)
                y2 = np.zeros(t2.shape)
                y2[(t2 >= t1[0]) & (t2 <= t1[-1])] = f(t2[(t2 >= t1[0]) & (t2 <= t1[-1])])
                y2[t2 < t1[0]] = v1[0]
                y2[t2 > t1[-1]] = v1[-1]
                return y2
            elif method == "previous":
                return scipy.interpolate.interp1d(t1, v1, kind="previous", copy=False, assume_sorted=True, bounds_error=False, fill_value=(v1[0], v1[-1]))(t2)
            else:
                raise Exception('Unknown interpolation type - must be one of "linear", "pchip", or "previous"')

        # Otherwise, `method` is a callable (class instance e.g. `scipy.interpolate.PchipInterpolator` or generating function) that
        # produces a callable function representation of the interpolation. This function is then called with the new time points
        interpolator = method(t1, v1, **kwargs)
        return interpolator(t2)

    def sample(self, constant=True):
        """
        Return a sampled copy of the TimeSeries

        This method returns a copy of the TimeSeries in which the values have been
        perturbed based on the uncertainty value.

        :param constant: If True, time series will be perturbed by a single constant offset. If False,
                         an different perturbation will be applied to each time specific value independently.
        :return: A copied ``TimeSeries`` with perturbed values

        """

        if self._sampled:
            raise Exception("Sampling has already been performed - can only sample once")

        new = self.copy()
        if self.sigma is not None:
            delta = self.sigma * np.random.randn(1)[0]
            if self.assumption is not None:
                new.assumption += delta

            if constant:
                # Use the same delta for all data points
                new.vals = [v + delta for v in new.vals]
            else:
                # Sample again for each data point
                for i, (v, delta) in enumerate(zip(new.vals, self.sigma * np.random.randn(len(new.vals)))):
                    new.vals[i] = v + delta

        # Sampling flag only needs to be set if the TimeSeries had data to change
        if new.has_data:
            new._sampled = True

        return new


def make_init_prev_fn(module, sim, uids, data=None, active=False):
    """ Initialize prevalence by sex and risk group """

    if data is None: data = module.init_prev_data

    if sc.isnumber(data):
        init_prev = data

    elif isinstance(data, pd.DataFrame):

        init_prev = pd.Series(index=uids)
        df = data

        nw = sim.networks.structuredsexual
        n_risk_groups = nw.pars.n_risk_groups
        for rg in range(n_risk_groups):
            for sex in ['female', 'male']:
                for sw in [0, 1]:
                    thisdf = df.loc[(df.risk_group==rg) & (df.sex==sex) & (df.sw==sw)]
                    conditions = sim.people[sex] & (nw.risk_group==rg)
                    if active:
                        conditions = conditions & nw.active(sim.people)
                    if sw:
                        if sex == 'female': conditions = conditions & sim.networks.structuredsexual.fsw
                        if sex == 'male':   conditions = conditions & sim.networks.structuredsexual.client
                    init_prev[conditions[uids]] = thisdf.init_prev.values[0]

    else:
        errormsg = 'Format of init_prev_data must be float or dataframe.'
        raise ValueError(errormsg)

    # Scale and validate
    init_prev = init_prev * module.pars.rel_init_prev
    init_prev = np.clip(init_prev, a_min=0, a_max=1)

    return init_prev


def finalize_results(sim, modules_to_drop=None):
    modules_to_drop = sc.promotetolist(modules_to_drop)
    modules_to_drop += ['yearvec']

    # Drop modules
    for res in modules_to_drop:
        if res in sim.results: del sim.results[res]

    # Make new results
    newres = ss.Results(module=None)
    newres.yearvec = np.arange(sim.pars.start, sim.pars.stop)  # np.unique(sim.yearvec.astype(int))[:-1]

    periods = int(1/sim.pars.dt)
    for rname, res in sim.results.items():
        if isinstance(res, ss.Results):
            newres[rname] = ss.Results(module=res)
            for mrname, modres in res.items():
                if mrname.startswith('new_'):
                    newres[rname][mrname] = np.sum(modres[:-1].reshape(-1, periods), axis=1)
                elif 'prev' in mrname or 'inci' in mrname or 'rel_' in mrname or mrname.startswith('n_'):
                    newres[rname][mrname] = np.mean(modres[:-1].reshape(-1, periods), axis=1)
        else:
            if rname.startswith('new_'):
                newres[rname] = np.sum(res[:-1].reshape(-1, periods), axis=1)
            elif 'prev' in rname or 'inci' in rname or 'rel_' in rname or rname.startswith('n_'):
                newres[rname] = np.mean(res[:-1].reshape(-1, periods), axis=1)

    def flatten_results(d, prefix=''):
        flat = {}
        for key, val in d.items():
            if isinstance(val, dict):
                flat.update(flatten_results(val, prefix=prefix+key+'.'))
            else:
                flat[prefix+key] = val
        return flat

    resdict = flatten_results(newres)
    df = sc.dataframe.from_dict(resdict).set_index('yearvec')

    return df


def shrink_calib(calib, n_results=100):
    cal = sc.objdict()
    plot_indices = calib.df.iloc[0:n_results, 0].values
    cal.sim_results = [calib.sim_results[i] for i in plot_indices]
    cal.analyzer_results = [calib.analyzer_results[i] for i in plot_indices]
    cal.target_data = calib.target_data
    cal.df = calib.df.iloc[0:n_results, ]
    return cal


def standardize_data(data=None, metadata=None, min_year=1800, out_of_range=0, default_age=0, default_year=2024):
    """
    Standardize formats of input data

    Input data can arrive in many different forms. This function accepts a variety of data
    structures, and converts them into a Pandas Series containing one variable, based on
    specified metadata, or an ``ss.Dist`` if the data is already an ``ss.Dist`` object.

    The metadata is a dictionary that defines columns of the dataframe or keys
    of the dictionary to use as indices in the output Series. It should contain:

    - ``metadata['data_cols']['value']`` specifying the name of the column/key to draw values from
    - ``metadata['data_cols']['year']`` optionally specifying the column containing year values; otherwise the default year will be used
    - ``metadata['data_cols']['age']`` optionally specifying the column containing age values; otherwise the default age will be used
    - ``metadata['data_cols'][<arbitrary>]`` optionally specifying any other columns to use as indices. These will form part of the multi-index for the standardized Series output.

    If a ``sex`` column is part of the index, the metadata can also optionally specify a string mapping to convert
    the sex labels in the input data into the 'm'/'f' labels used by Starsim. In that case, the metadata can contain
    an additional key like ``metadata['sex_keys'] = {'Female':'f','Male':'m'}`` which in this case would map the strings
    'Female' and 'Male' in the original data into 'm'/'f' for Starsim.

    Args:
        data (pandas.DataFrame, pandas.Series, dict, int, float): An associative array  or a number, with the input data to be standardized.
        metadata (dict): Dictionary specifiying index columns, the value column, and optionally mapping for sex labels
        min_year (float): Optionally specify a minimum year allowed in the data. Default is 1800.
        out_of_range (float): Value to use for negative ages - typically 0 is a reasonable choice but other values (e.g., np.inf or np.nan) may be useful depending on the calculation. This will automatically be added to the dataframe with an age of ``-np.inf``

    Returns:

        - A `pd.Series` for all supported formats of `data` *except* an ``ss.Dist``. This series will contain index columns for 'year'
          and 'age' (in that order) and then subsequent index columns for any other variables specified in the metadata, in the order
          they appeared in the metadata (except for year and age appearing first).
        - An ``ss.Dist`` instance - if the ``data`` input is an ``ss.Dist``, that same object will be returned by this function
    """
    # It's a format that can be used directly: return immediately
    if sc.isnumber(data) or isinstance(data, (ss.Dist, ss.TimePar)):
        return data

    # Convert series and dataframe inputs into dicts
    if isinstance(data, pd.Series) or isinstance(data, pd.DataFrame):
        data = data.reset_index().to_dict(orient='list')

    # Check that the input is now a dict (scalar types have already been handled above
    assert isinstance(data, dict), 'Supported inputs are ss.Dict, scalar numbers, DataFrames, Series, or dictionaries'

    # Extract the values and index columns
    assert 'value' in metadata['data_cols'], 'The metadata is missing a column name for "value", which must be provided if the input data is in the form of a DataFrame, Series, or dict'
    values = sc.promotetoarray(data[metadata['data_cols']['value']])
    index = sc.objdict()
    for k, col in metadata['data_cols'].items():
        if k != 'value':
            index[k] = data[col]

    # Add defaults for year and age
    if 'year' not in index:
        index['year'] = np.full(values.shape, fill_value=default_year, dtype=ss.dtypes.float)
    if 'age' not in index:
        index['age'] = np.full(values.shape, fill_value=default_age, dtype=ss.dtypes.float)

    # Reorder the index so that it starts with age first
    index.insert(0, 'age', index.pop('age'))
    index.insert(0, 'year', index.pop('year'))

    # Map sex values
    if 'sex' in index:
        if 'sex_keys' in metadata:
            index['sex'] = [metadata['sex_keys'][x] for x in index['sex']]
        assert set(index['sex']) == {'f','m'}, 'If standardized data contains a "sex" column, it should use "m" and "f" to specify the sex.'

    # Create the series
    output = pd.Series(data=values, index=pd.MultiIndex.from_arrays(index.values(), names=index.keys()))

    # Add an entry for negative ages
    new_entries = output.index.set_codes([0] * len(output), level=1).set_levels([-np.inf], level=1).unique()
    new = pd.Series(data=out_of_range, index=new_entries)
    output = pd.concat([output, new])

    # Truncate years
    output = output.loc[output.index.get_level_values('year')>min_year]

    # Perform a final sort to prevent "indexing past lexsort depth" warnings
    output = output.sort_index()

    return output

