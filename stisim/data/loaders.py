'''
Load data
'''

#%% Housekeeping
import pandas as pd
import numpy as np
import sciris as sc
import starsim as ss


__all__ = ['get_age_distribution', 'get_rates', 'load', 'load_local', 'DataLoader']


thisdir = sc.thispath(__file__)
filesdir = thisdir / 'files'
files = sc.objdict()
files.age_dist_sex = 'populations_by_sex.obj'
files.birth = 'births.csv'
files.asfr = 'asfr.csv'
files.death = 'deaths.csv'
files.migration = 'migration.csv'


def get_age_distribution(location=None, year=None, datafolder=None):
    """ Load age distribution for a given location & year"""
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename= f'{location}_age_{year}.csv', folder=datafolder)
    try:
        raw_df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for age distribution.'
        raise ValueError(errormsg) from E

    # Remap column names
    raw_df = raw_df.rename(columns={'Time': 'year', 'AgeStart': 'age', 'Value': 'value'})

    # Add year column if not present (some files encode year in the filename only)
    if 'year' not in raw_df.columns:
        raw_df['year'] = year
    raw_df = raw_df[['year', 'age', 'value']]

    return raw_df


def get_rates(location, which, datafolder=None):
    """ Load birth/death/fertility/migration rates for a given location """
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename=f'{location}_{files[which]}', folder=datafolder)
    try:
        df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for {which}.'
        raise ValueError(errormsg) from E
    return df


def load(location):
    """
    Resolve the data path for a location.

    Checks for a local cache of data files. In the future, will auto-download
    from the stisim_data GitHub repository if not found locally.

    Args:
        location (str): Location name (e.g. 'zimbabwe')

    Returns:
        Path: Path to the location's data directory
    """
    location = location.lower().strip()
    local_path = filesdir / location
    if local_path.exists():
        return local_path
    # TODO: auto-download from https://github.com/starsimhub/stisim_data
    raise NotImplementedError(
        f'Automatic data download for "{location}" is not yet implemented. '
        f'Please provide data manually via data_path="path/to/folder".'
    )


def load_local(path):
    """
    Validate and return a local data folder path.

    Args:
        path (str/Path): Path to folder containing CSV data files

    Returns:
        Path: Validated path to the data folder
    """
    path = sc.path(path)
    if not path.exists():
        raise FileNotFoundError(f'Data folder not found: {path}')
    return path


class DataLoader:
    """
    Load location-specific data and create configured STI simulation modules.

    Loads CSV data files (initial prevalence, condom use, ART/VMMC coverage, etc.)
    from a data directory and creates disease, network, and intervention module
    instances from the loaded data. Parameter assignment (e.g. ``sti_pars``,
    ``nw_pars``) is handled by the Sim, not by the DataLoader.

    Args:
        data_path (str/Path): Path to directory containing CSV data files
        location (str): Location name (used for file naming conventions, e.g.
            to find ``{location}_hiv_data.csv``)
        diseases (str/list): Disease name(s) to load data for

    Examples:
        >>> dl = sti.DataLoader(data_path='path/to/data', diseases='hiv')
        >>> dl.load()
        >>> modules = dl.make_modules()

        >>> # Or let Sim handle it automatically:
        >>> sim = sti.Sim(diseases='hiv', location='zimbabwe',
        ...               data_path='path/to/data', sti_pars=sti_pars)
    """

    def __init__(self, data_path=None, location=None, diseases=None):
        self.location = location
        self.data_path = sc.path(data_path) if data_path is not None else None
        self.diseases = sc.tolist(diseases) if diseases is not None else []
        self.data = sc.objdict()

    def load(self):
        """
        Load CSV data files from the data directory.

        Populates ``self.data`` with DataFrames for each available CSV:
          - ``data.diseases.{name}.init_prev`` from ``init_prev_{name}.csv``
          - ``data.condom_use`` from ``condom_use.csv``
          - ``data.art_coverage`` from ``art_coverage.csv``
          - ``data.vmmc_coverage`` from ``vmmc_coverage.csv``
          - ``data.hiv_data`` from ``{location}_hiv_data.csv``

        Returns:
            self: For method chaining (e.g. ``dl.load().make_modules()``)
        """
        if self.data_path is None:
            raise ValueError('data_path is required for DataLoader.load()')
        if not self.data_path.exists():
            raise FileNotFoundError(f'Data path not found: {self.data_path}')

        # Load disease-specific data
        self.data.location = self.location
        self.data.diseases = sc.objdict()
        for disease in self.diseases:
            disease_data = sc.objdict()
            init_prev = self._load_csv(f'init_prev_{disease.lower()}.csv')
            if init_prev is not None:
                disease_data.init_prev = init_prev
            self.data.diseases[disease.lower()] = disease_data

        # Load behavioral and intervention data
        for key in ['condom_use', 'art_coverage', 'vmmc_coverage']:
            df = self._load_csv(f'{key}.csv')
            if df is not None:
                self.data[key] = df

        # Load calibration/comparison data for plotting overlay
        if self.location is not None:
            hiv_data = self._load_csv(f'{self.location}_hiv_data.csv')
            if hiv_data is not None:
                self.data.hiv_data = hiv_data

        return self

    def make_modules(self):
        """
        Create module instances from loaded data.

        Creates disease instances with loaded init_prev data, network instances
        with condom data, and intervention instances (ART, VMMC) from coverage
        CSVs. Does NOT apply user parameters — that is handled by the Sim's
        ``process_stis()`` and ``process_networks()`` methods.

        Returns:
            sc.objdict: With keys ``diseases``, ``networks``,
                ``interventions``, and ``data``.
        """
        import stisim as sti

        result = sc.objdict()

        # Create diseases with loaded data
        diseases = []
        for disease in self.diseases:
            d_lower = disease.lower()
            d_kwargs = {}
            if d_lower in self.data.diseases and 'init_prev' in self.data.diseases[d_lower]:
                d_kwargs['init_prev_data'] = self.data.diseases[d_lower].init_prev
            if d_lower == 'hiv':
                diseases.append(sti.HIV(**d_kwargs))
            else:
                diseases.append(sti.make_sti(d_lower, pars=d_kwargs))
        result.diseases = diseases

        # Create networks with loaded condom data
        nw_kwargs = {}
        if 'condom_use' in self.data:
            nw_kwargs['condom_data'] = self.data.condom_use
        result.networks = [sti.StructuredSexual(**nw_kwargs), ss.MaternalNet()]

        # Create interventions from coverage data
        interventions = []
        if 'art_coverage' in self.data:
            art_df = self.data.art_coverage.set_index('year')
            interventions.append(sti.ART(coverage_data=art_df))
        if 'vmmc_coverage' in self.data:
            vmmc_df = self.data.vmmc_coverage.set_index('year')
            interventions.append(sti.VMMC(coverage_data=vmmc_df))
        result.interventions = interventions

        # Calibration/comparison data for plotting
        result.data = self.data.get('hiv_data', None)

        return result

    def _load_csv(self, filename):
        """Load a CSV file from the data path, returning None if not found."""
        filepath = self.data_path / filename
        if filepath.exists():
            return pd.read_csv(filepath)
        return None
