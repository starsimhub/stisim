import sciris as sc
import starsim as ss
import stisim as sti

from stisim.diseases.hiv import HIV, HIVPars
from stisim.interventions.hiv_interventions import HIVTest, ART, VMMC, Prep

__all__ = ['Sim', 'HIV', 'HIVPars', 'HIVTest', 'ART', 'VMMC', 'Prep']

class Sim(sti.Sim):
    """
    A subclass of stisim.Sim that is specifically designed for HIV simulations.
    """
    def __init__(self, sim_pars=None, hiv_pars=None, location=None, **kwargs):

        if location is not None:
            raise NotImplementedError('Location-based sim creation is not implemented yet')

        # Initialize parameters
        sim_pars = sc.mergedicts(sim_pars)
        hiv_pars = sc.mergedicts(hiv_pars)
        default_sim_keys = ss.SimPars().keys()
        default_hiv_keys = sti.diseases.hiv().pars.keys()

        # Pull things out for special processing
        diseases = sc.mergelists(kwargs.pop('diseases'))
        networks = sc.mergelists(kwargs.pop('networks'))
        demographics = sc.mergelists(kwargs.pop('demographics'))
        interventions = sc.mergelists(kwargs.pop('interventions'))

        # Handle diseases -- first, figure out what parameters belong in HIV
        for key in kwargs.keys():
            if key in default_hiv_keys:
                if key in default_sim_keys: # If the key is in both, copy
                    val = kwargs[key]
                else: # Else, pop
                    val = kwargs.pop(key)
                hiv_pars[key] = val
        
        hiv = sti.HIV(pars=hiv_pars)
        diseases.insert(0, hiv)

        # Handle demographics
        if not demographics:
            demographics = [sti.Pregnancy(), ss.Deaths()]

        # Handle networks
        if not networks:
            networks = [sti.StructuredSexual(), ss.MaternalNet()]

        # Handle interventions 
        if not interventions:
            interventions = [HIVTest(), ART(), VMMC(), Prep()]

        # Handle interventions
        super().__init__(
            pars=sim_pars, 
            demographics=demographics, 
            networks=networks, 
            diseases=diseases, 
            interventions=interventions,
            **kwargs
        )
        return