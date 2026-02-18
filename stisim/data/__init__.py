from .loaders import *
from .downloaders import *

import sciris as sc

_cache_dir = sc.thispath(__file__) / 'files'


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
    local_path = _cache_dir / location
    if local_path.exists():
        return local_path
    # TODO: auto-download from https://github.com/starsimhub/stisim_data
    raise NotImplementedError(
        f'Automatic data download for "{location}" is not yet implemented. '
        f'Please provide data manually via data="path/to/folder".'
    )


def download(location, folder=None):
    """
    Download location data from the stisim_data repository.

    Not yet implemented.

    Args:
        location (str): Location name (e.g. 'kenya')
        folder (str/Path): Local folder to save data to. If None, uses default cache.
    """
    raise NotImplementedError(
        'Automatic data download is not yet implemented.'
    )


def load_local(path):
    """
    Load data from a local folder.

    Args:
        path (str/Path): Path to folder containing CSV data files

    Returns:
        Path: Validated path to the data folder
    """
    path = sc.path(path)
    if not path.exists():
        raise FileNotFoundError(f'Data folder not found: {path}')
    return path
