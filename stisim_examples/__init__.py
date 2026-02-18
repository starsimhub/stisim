"""
STIsim Examples - Pre-configured location-specific simulations.

This package provides factory functions to create STIsim simulations with
location-specific data (demographics, disease parameters, behavioral data).

Usage:
    >>> import stisim_examples as stx
    >>> sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
    >>> sim.run()

    >>> # HIV-specific convenience
    >>> sim = stx.HIVSim(location='zimbabwe')
    >>> sim.run()
"""
import stisim as sti
import sciris as sc
from .loaders import load_location_data, list_locations, LOCATIONS

__all__ = ['Sim', 'HIVSim', 'MultiSim', 'list_locations', 'LOCATIONS']


def Sim(demographics=None, diseases=None, **kwargs):
    """
    Factory function to create pre-configured stisim simulations.

    This function loads location-specific data (demographics, disease parameters,
    behavioral data) and creates a fully configured stisim.Sim instance.

    Args:
        demographics (str): Location name for demographics (zimbabwe, kenya, demo)
        diseases (str/list): Disease(s) to include (e.g., 'hiv', ['hiv', 'syphilis'])
        **kwargs: Additional parameters passed to stisim.Sim

    Returns:
        stisim.Sim: Configured simulation instance

    Examples:
        >>> # Zimbabwe HIV model
        >>> sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
        >>> sim.run()

        >>> # Zimbabwe HIV-Syphilis model
        >>> sim = stx.Sim(demographics='zimbabwe', diseases=['hiv', 'syphilis'])
        >>> sim.run()

        >>> # Kenya HIV model with custom parameters
        >>> sim = stx.Sim(demographics='kenya', diseases='hiv', n_agents=5000, dur=30)
        >>> sim.run()
    """
    if demographics is None:
        raise ValueError("demographics parameter is required. Use one of: " + ", ".join(LOCATIONS.keys()))

    if diseases is None:
        raise ValueError("diseases parameter is required (e.g., 'hiv' or ['hiv', 'syphilis'])")

    # Load location-specific data
    loc_data = load_location_data(demographics, diseases)

    # Import location-specific helper to configure parameters
    configured_pars = {}
    try:
        location_module = __import__(f'stisim_examples.{demographics}.sim', fromlist=['configure_sim_pars'])
        configure_func = getattr(location_module, 'configure_sim_pars', None)
        if configure_func is not None:
            configured_pars = configure_func(loc_data, diseases, kwargs)
    except (ImportError, AttributeError):
        # If no location-specific configuration exists, proceed with basic setup
        # This allows locations to be added incrementally
        pass

    # Use configured disease instances if provided, otherwise fall back to string names
    if 'diseases' in configured_pars:
        kwargs['diseases'] = configured_pars.pop('diseases')
    else:
        kwargs['diseases'] = diseases

    # Use configured networks; user-provided networks replace defaults
    if 'networks' in configured_pars:
        if 'networks' not in kwargs:
            kwargs['networks'] = configured_pars.pop('networks')
        else:
            configured_pars.pop('networks')

    # Configured interventions always included; user-provided interventions are additive
    if 'interventions' in configured_pars:
        base_intvs = configured_pars.pop('interventions')
        user_intvs = sc.tolist(kwargs.pop('interventions', []))
        kwargs['interventions'] = base_intvs + user_intvs

    # Extract data for plotting overlay (not a valid sti.Sim kwarg)
    data = configured_pars.pop('data', kwargs.pop('data', None))

    # Merge remaining configured params (e.g., disease_pars, sim_pars)
    kwargs = sc.mergedicts(kwargs, configured_pars)

    # Always set demographics for demographic data loading
    kwargs['demographics'] = demographics

    # Create the sim and attach data for plot overlay
    sim = sti.Sim(**kwargs)
    if data is not None:
        sim.data = data

    return sim


def HIVSim(location=None, **kwargs):
    """
    Convenience function for HIV-only simulations.

    This is a shorthand for Sim(demographics=location, diseases='hiv', **kwargs).

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        **kwargs: Additional parameters passed to Sim()

    Returns:
        stisim.Sim: Configured HIV simulation

    Examples:
        >>> # Zimbabwe HIV model
        >>> sim = stx.HIVSim(location='zimbabwe')
        >>> sim.run()

        >>> # Kenya HIV model with custom duration
        >>> sim = stx.HIVSim(location='kenya', dur=50)
        >>> sim.run()
    """
    if location is None:
        raise ValueError("location parameter is required. Use one of: " + ", ".join(LOCATIONS.keys()))

    return Sim(demographics=location, diseases='hiv', **kwargs)


