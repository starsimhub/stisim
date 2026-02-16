import sciris as sc
import starsim as ss
import stisim as sti

from stisim.diseases.hiv import HIV, HIVPars
from stisim.interventions.hiv_interventions import HIVTest, ART, VMMC, Prep

__all__ = ['Sim', 'HIV', 'HIVPars', 'HIVTest', 'ART', 'VMMC', 'Prep', 'ss', 'sti']

class Sim(sti.Sim):
    """
    A subclass of stisim.Sim that is specifically designed for HIV simulations.

    This parses input parameters among the sim and the HIV module,
    and adds default demographics (pregnancy and deaths), networks (sexual and maternal),
    and interventions (testing, ART, VMMC, and PrEP).

    Args:
        pars (dict): Parameters for the simulation
        sim_pars (dict): Parameters specific to the simulation
        hiv_pars (dict): Parameters specific to HIV
        location (str): Location name to load demographic data (e.g., 'zimbabwe', 'kenya')
        datafolder (str/Path): Custom path to data folder (defaults to stisim/data/files)
        **kwargs: Additional keyword arguments passed to stisim.Sim

    Example:
        >>> # Basic HIV sim
        >>> sim = hs.Sim()

        >>> # HIV sim with location-based demographics
        >>> sim = hs.Sim(location='zimbabwe')

        >>> # HIV sim with custom data folder
        >>> sim = hs.Sim(location='zimbabwe', datafolder='/path/to/data')
    """
    def __init__(self, pars=None, sim_pars=None, hiv_pars=None, location=None, datafolder=None, **kwargs):

        # Handle location parameter - pass through as demographics to parent
        if location is not None:
            kwargs['demographics'] = location

        # Handle datafolder parameter - pass through to parent
        if datafolder is not None:
            kwargs['datafolder'] = datafolder

        # Initialize parameters
        pars = sc.mergedicts(pars, kwargs)
        sim_pars = sc.mergedicts(sim_pars)
        hiv_pars = sc.mergedicts(hiv_pars)
        default_sim_keys = ss.SimPars().keys()
        default_hiv_keys = HIVPars().keys()

        # Pull modules out for special processing
        modules = sc.objdict()
        for mod_type in ['diseases', 'networks', 'demographics', 'interventions']:
            # Special handling for demographics: if it's a string (location name),
            # don't pop it - let it pass through to parent stisim.Sim
            if mod_type == 'demographics' and isinstance(pars.get('demographics'), str):
                modules[mod_type] = []  # Empty list, will use default
                # Don't pop demographics - let it go to parent
            else:
                modules[mod_type] = sc.mergelists(pars.pop(mod_type, None)) # Remove from kwargs and turn into a list

        # Handle diseases -- first, figure out what parameters belong in HIV
        for key in list(pars.keys()):
            if key in default_hiv_keys:
                if key in default_sim_keys: # If the key is in both, copy
                    val = pars[key]
                else: # Else, pop
                    val = pars.pop(key)
                hiv_pars[key] = val

        hiv = sti.HIV(pars=hiv_pars)
        modules.diseases.insert(0, hiv)

        # Handle demographics
        # Only add defaults if no location was specified and no demographics were provided
        if not modules.demographics and location is None:
            modules.demographics = [ss.Pregnancy(), ss.Deaths()]

        # Handle networks
        if not modules.networks:
            modules.networks = [sti.StructuredSexual(), ss.MaternalNet()]

        # Handle interventions
        if not modules.interventions:
            modules.interventions = [HIVTest(), ART(), VMMC(), Prep()]

        # Handle interventions
        # If location was specified, don't pass demographics modules - let parent process the location string
        if location is None:
            super().__init__(
                pars          = sim_pars,
                demographics  = modules.demographics,
                networks      = modules.networks,
                diseases      = modules.diseases,
                interventions = modules.interventions,
                **pars
            )
        else:
            super().__init__(
                pars          = sim_pars,
                networks      = modules.networks,
                diseases      = modules.diseases,
                interventions = modules.interventions,
                **pars
            )
        return