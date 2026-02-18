"""
HIVsim Examples - Pre-configured location-specific HIV simulations.

Usage:
    >>> import hivsim_examples as hx
    >>> sim = hx.Sim(location='zimbabwe')
    >>> sim.run()
"""
from stisim_examples.loaders import list_locations, LOCATIONS

__all__ = ['Sim', 'list_locations', 'LOCATIONS']


def Sim(location=None, **kwargs):
    """
    Create a pre-configured HIV simulation for a given location.

    Imports the location's make_sim() from stisim_examples and calls it.

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        **kwargs: Additional parameters passed to make_sim()
    """
    import importlib
    if location is None:
        raise ValueError("location parameter is required. Use one of: " + ", ".join(LOCATIONS.keys()))
    loc_module = importlib.import_module(f'stisim_examples.{location}.sim')
    return loc_module.make_sim(**kwargs)
