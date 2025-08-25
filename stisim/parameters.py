"""
Set parameters
"""
import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti


__all__ = ['SimPars', 'sti_aliases', 'sti_register', 'merged_sti_pars', 'make_sti']


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
        self.start     = ss.years(2000)   # Start of the simulation
        self.stop      = ss.years(2030)   # End of the simulation
        self.dur       = None   # Duration of time to run, if stop isn't specified
        self.dt        = ss.years(1/12)   # Timestep
        self.rand_seed = 1      # Random seed; if None, don't reset
        self.verbose = 1/12

        # Demographic parameters
        self.birth_rate = None
        self.death_rate = None
        self.use_aging  = False
        self.use_pregnancy = True
        self.use_migration = False

        # Update with any supplied parameter values and generate things that need to be generated
        self.update(kwargs)
        return


def sti_aliases():
    """
    Define aliases for STIs
    """
    # List of choices available
    choices = {
        'ng':    ['gonorrhea', 'gonorrhoea', 'ng', 'gon'],
        'ct':    ['chlamydia', 'ct', 'chlam'],
        'tv':    ['trichomoniasis', 'trichomonas', 'trich', 'tv'],
        'bv':    ['bacterial vaginosis', 'bv'],
        'syph':  ['syphilis', 'syph'],
        'gud':   ['genital ulcerative disease', 'gud'],
        'hiv':   ['hiv', 'human immunodeficiency virus']
    }
    mapping = {name: key for key,synonyms in choices.items() for name in synonyms} # Flip from key:value to value:key
    return choices, mapping


def sti_register(key=None):
    """
    Registry of STI names linked to disease classes
    """
    sti_dict = sc.objdict(
        bv=sti.BV,
        ct=sti.Chlamydia,
        ctbl=sti.ChlamydiaBL,
        gud=sti.GUDPlaceholder,  # Placeholder for GUD
        ng=sti.Gonorrhea,
        hiv=sti.HIV,
        syph=sti.Syphilis,
        tv=sti.Trichomoniasis,
    )
    return sti_dict if key is None else sti_dict[key]


def merged_sti_pars():
    """ Merge all the parameters from the STI disease modules """
    # Initialize an empty dictionary to hold the merged parameters
    merged_pars = {}
    # Get the parameters from each disease module
    for dname, disease in sti_register().items():
        disease_pars = disease().pars
        # Merge the parameters into the merged_pars dictionary
        merged_pars.update(disease_pars)
    return merged_pars


def make_sti(name, pars=None):
    """
    Create an STI disease module based on the name and parameters provided.

    Args:
        name (str): Name of the STI disease module to create.
        pars (dict, optional): Parameters to initialize the disease module with.

    Returns:
        ss.Disease: An instance of the specified STI disease module.
    """
    sti_dict = sti_register()
    if name not in sti_dict:
        raise ValueError(f'Unknown STI name: {name}. Available options are: {list(sti_dict.keys())}')

    # Create an instance of the disease class
    disease_class = sti_dict[name]
    return disease_class(pars=pars)
