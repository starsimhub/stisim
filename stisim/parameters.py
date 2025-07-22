"""
Set parameters
"""
import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
import hpvsim as hpv


__all__ = ['SimPars', 'make_sim_pars', 'NetworkPars', 'make_network_pars']


class SimPars(ss.SimPars):
    """
    Subclass of Starsim's SimPars with defaults for STI simulations. Refer to
    Starsim's SimPars for more information on the parameters.
    """
    def __init__(self, **kwargs):

        # Initialize the parent class
        super().__init__()

        # General parameters
        self.label   = ''  # The label of the simulation
        self.verbose = ss.options.verbose  # Whether or not to display information during the run -- options are 0 (silent), 0.1 (some; default), 1 (default), 2 (everything)

        # Population parameters
        self.n_agents  = 1e3   # Number of agents
        self.total_pop = None  # If defined, used for calculating the scale factor
        self.pop_scale = None  # How much to scale the population

        # Simulation parameters
        self.unit      = 'year' # The time unit to use; options are 'year' (default), 'day', 'week', 'month', or 'none'
        self.start     = 2000   # Start of the simulation
        self.stop      = 2030   # End of the simulation
        self.dur       = None   # Duration of time to run, if stop isn't specified
        self.dt        = 1/12   # Timestep (in units of self.unit)
        self.rand_seed = 1      # Random seed; if None, don't reset
        self.verbose = 1/12

        # Demographic parameters
        self.birth_rate = None
        self.death_rate = None
        self.use_aging  = True

        # Update with any supplied parameter values and generate things that need to be generated
        self.update(kwargs)
        return


def make_sim_pars(**kwargs):
    """ Shortcut for making a new instance of SimPars """
    return SimPars(**kwargs)


def make_network_pars(**kwargs):
    """ Shortcut for making a new instance of NetworkPars """
    return NetworkPars(**kwargs)

