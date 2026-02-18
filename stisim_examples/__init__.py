"""
STIsim Examples - Pre-configured location-specific simulations.

Each location subdirectory (e.g. zimbabwe/, demo/) contains a sim.py with:
  - Parameter dicts (sim_pars, sti_pars, nw_pars)
  - A make_sim() function that creates a configured sti.Sim

Usage:
    >>> from stisim_examples.zimbabwe.sim import make_sim
    >>> sim = make_sim()
    >>> sim.run()
"""
from .loaders import load_location_data, list_locations, LOCATIONS

__all__ = ['list_locations', 'LOCATIONS', 'load_location_data']
