"""
Data loading utilities for location-specific STIsim examples.
"""
import pandas as pd
import sciris as sc

__all__ = ['load_location_data', 'list_locations', 'LOCATIONS']

# Registry of available locations
LOCATIONS = {
    'zimbabwe': 'Zimbabwe',
    'kenya': 'Kenya',
    'demo': 'Demo/Generic location',
}


def list_locations():
    """
    List available locations.

    Returns:
        dict: Dictionary mapping location names to descriptions
    """
    return LOCATIONS.copy()


def load_location_data(location, diseases=None):
    """
    Load location-specific data for specified diseases.

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        diseases (str/list): Disease(s) to load data for

    Returns:
        sc.objdict: Nested dict with disease parameters and behavioral data
            - location (str): Location name
            - diseases (dict): Disease-specific data keyed by disease name
            - condom_use (DataFrame): Condom use data (if available)
            - art_coverage (DataFrame): ART coverage data (if available)
            - vmmc_coverage (DataFrame): VMMC coverage data (if available)

    Example:
        >>> data = load_location_data('zimbabwe', diseases='hiv')
        >>> print(data.diseases.hiv.init_prev.head())
        >>> print(data.condom_use.head())
    """
    # Validate location
    if location not in LOCATIONS:
        available = ', '.join(LOCATIONS.keys())
        raise ValueError(f"Location '{location}' not found. Available locations: {available}")

    # Normalize diseases to list
    if isinstance(diseases, str):
        diseases = [diseases]
    elif diseases is None:
        diseases = []

    # Get location path
    loc_path = sc.thispath(__file__) / location
    if not loc_path.exists():
        raise FileNotFoundError(f"Location directory not found: {loc_path}")

    # Initialize data structure
    data = sc.objdict()
    data.location = location
    data.diseases = sc.objdict()

    # Load disease-specific data
    for disease in diseases:
        disease_data = sc.objdict()
        disease_lower = disease.lower()

        # Load initial prevalence
        init_prev_file = loc_path / f'init_prev_{disease_lower}.csv'
        if init_prev_file.exists():
            disease_data.init_prev = pd.read_csv(init_prev_file)

        # Store disease data (even if empty, for consistency)
        data.diseases[disease_lower] = disease_data

    # Load behavioral data (shared across diseases)
    condom_file = loc_path / 'condom_use.csv'
    if condom_file.exists():
        data.condom_use = pd.read_csv(condom_file)

    # Load intervention data
    art_file = loc_path / 'art_coverage.csv'
    if art_file.exists():
        data.art_coverage = pd.read_csv(art_file)

    vmmc_file = loc_path / 'vmmc_coverage.csv'
    if vmmc_file.exists():
        data.vmmc_coverage = pd.read_csv(vmmc_file)

    # Load calibration/comparison data (for plotting overlay)
    hiv_data_file = loc_path / f'{location}_hiv_data.csv'
    if hiv_data_file.exists():
        data.hiv_data = pd.read_csv(hiv_data_file)

    return data


def load_csv_if_exists(filepath, required=False):
    """
    Load a CSV file if it exists.

    Args:
        filepath (Path/str): Path to CSV file
        required (bool): If True, raise error if file doesn't exist

    Returns:
        DataFrame or None: Loaded data or None if file doesn't exist
    """
    filepath = sc.path(filepath)
    if filepath.exists():
        return pd.read_csv(filepath)
    elif required:
        raise FileNotFoundError(f"Required data file not found: {filepath}")
    else:
        return None
