import sciris as sc
import starsim as ss
import stisim as sti

from stisim.diseases.hiv import HIV, HIVPars
from stisim.interventions.hiv_interventions import HIVTest, ART, VMMC, Prep

__all__ = ['Sim', 'HIV', 'HIVPars', 'HIVTest', 'ART', 'VMMC', 'Prep', 'ss', 'sti']

class Sim(sti.Sim):
    """
    A subclass of stisim.Sim that is specifically designed for HIV simulations.

    Currently, this simply parses input parameters among the sim and the HIV module,
    and adds default demographics (pregnancy and deaths), networks (sexual and maternal),
    and interventions (testing, ART, VMMC, and PrEP).

    In future this will support location data and other features.
    """
    def __init__(self, pars=None, sim_pars=None, hiv_pars=None, location=None, **kwargs):

        if location is not None:
            raise NotImplementedError('Location-based sim creation is not implemented yet')

        # Initialize parameters
        pars = sc.mergedicts(pars, kwargs)
        sim_pars = sc.mergedicts(sim_pars)
        hiv_pars = sc.mergedicts(hiv_pars)
        default_sim_keys = ss.SimPars().keys()
        default_hiv_keys = HIVPars().keys()

        # Pull modules out for special processing
        modules = sc.objdict()
        for mod_type in ['diseases', 'networks', 'demographics', 'interventions']:
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
        if not modules.demographics:
            modules.demographics = [ss.Pregnancy(), ss.Deaths()]

        # Handle networks
        if not modules.networks:
            modules.networks = [sti.StructuredSexual(), ss.MaternalNet()]

        # Handle interventions
        if not modules.interventions:
            modules.interventions = [HIVTest(), ART(), VMMC(), Prep()]

        # Handle interventions
        super().__init__(
            pars          = sim_pars,
            demographics  = modules.demographics,
            networks      = modules.networks,
            diseases      = modules.diseases,
            interventions = modules.interventions,
            **pars
        )
        return