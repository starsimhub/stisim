"""
HIVsim Examples - Pre-configured location-specific HIV simulations.

This is a convenience wrapper around stisim_examples.HIVSim for HIV-specific workflows.

Usage:
    >>> import hivsim_examples as hx
    >>> sim = hx.Sim(location='zimbabwe')
    >>> sim.run()
"""
import stisim_examples as stx

# Alias HIVSim as Sim for cleaner syntax
Sim = stx.HIVSim
list_locations = stx.list_locations
LOCATIONS = stx.LOCATIONS

__all__ = ['Sim', 'list_locations', 'LOCATIONS']
