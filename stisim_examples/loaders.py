"""
Data loading utilities for location-specific STIsim examples.
"""
import sciris as sc
import stisim as sti

__all__ = ['load_location_data', 'list_locations', 'LOCATIONS']

# Registry of available locations
LOCATIONS = {
    'zimbabwe': 'Zimbabwe',
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

    Delegates to ``sti.DataLoader`` for the actual CSV loading.

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        diseases (str/list): Disease(s) to load data for

    Returns:
        sc.objdict: Nested dict with disease parameters and behavioral data

    Example:
        >>> data = load_location_data('zimbabwe', diseases='hiv')
        >>> print(data.diseases.hiv.init_prev.head())
    """
    if location not in LOCATIONS:
        available = ', '.join(LOCATIONS.keys())
        raise ValueError(f"Location '{location}' not found. Available locations: {available}")

    data_path = sc.thispath(__file__) / location
    loader = sti.DataLoader(data_path=data_path, location=location, diseases=diseases)
    loader.load()
    return loader.data


def load_csv_if_exists(filepath, required=False):
    """
    Load a CSV file if it exists.

    Args:
        filepath (Path/str): Path to CSV file
        required (bool): If True, raise error if file doesn't exist

    Returns:
        DataFrame or None: Loaded data or None if file doesn't exist
    """
    import pandas as pd
    filepath = sc.path(filepath)
    if filepath.exists():
        return pd.read_csv(filepath)
    elif required:
        raise FileNotFoundError(f"Required data file not found: {filepath}")
    else:
        return None
