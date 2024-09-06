"""
Define syphilis interventions for STIsim
"""

import starsim as ss
import numpy as np
from stisim.interventions.base_interventions import STIDx, STITreatment


__all__ = ["SyphDx", "SyphTx", "NewbornTreatment"]


class SyphDx(STIDx):
    def __init__(self, df, *args, **kwargs):
        super().__init__(df, 'syphilis', *args, **kwargs)
        return


class SyphTx(STITreatment):
    """
    Treat a fixed number of people each timestep.
    """

    def __init__(self, pars=None, max_capacity=None, years=None, eligibility=None, **kwargs):
        super().__init__(disease='syphilis', eligibility=eligibility, years=years, max_capacity=max_capacity)
        self.default_pars(
            rel_treat_prob=1,
            treat_prob=ss.bernoulli(p=0.9),
            treat_eff=ss.bernoulli(p=0.95),
            fetus_age_cutoff_treat_eff=-0.25,  # Reduced treatment efficacy for fetuses in the last trimester
            treat_eff_reduced=ss.bernoulli(p=0.2)  # Reduced efficacy for fetuses older than cut off
        )
        self.update_pars(pars, **kwargs)
        return

    def administer(self, sim, uids, return_format='dict'):
        """ Administer treatment, keeping track of unnecessarily treated individuals """

        inf = sim.diseases.syphilis.infected
        sus = sim.diseases.syphilis.susceptible
        inf_uids = uids[inf[uids]]
        sus_uids = uids[sus[uids]]

        successful = self.pars.treat_eff.filter(inf_uids)
        unsuccessful = np.setdiff1d(inf_uids, successful)
        unnecessary = sus_uids

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
        elif return_format == 'array':
            output = successful

        return output

    def change_states(self, sim, treat_succ):
        """ Change the states of people who are treated """
        sim.diseases.syphilis.primary[treat_succ] = False
        sim.diseases.syphilis.secondary[treat_succ] = False
        sim.diseases.syphilis.latent[treat_succ] = False
        sim.diseases.syphilis.tertiary[treat_succ] = False
        sim.diseases.syphilis.ti_primary[treat_succ] = np.nan
        sim.diseases.syphilis.ti_secondary[treat_succ] = np.nan
        sim.diseases.syphilis.ti_latent[treat_succ] = np.nan
        sim.diseases.syphilis.ti_tertiary[treat_succ] = np.nan
        sim.diseases.syphilis.susceptible[treat_succ] = True
        sim.diseases.syphilis.infected[treat_succ] = False

    def treat_fetus(self, sim, mother_uids):
        """
        Treat fetuses of successfully treated mothers.
        Birth outcomes of successfully treated fetuses will be updated. If fetus is not yet infected, the treatment doesn't do anything.
        Treatment success depends on when in the pregnancy the mother gets treated
        with reduced treatment efficacy above the cut off age
        """
        # Get Fetuses
        children_uids = ss.uids(sim.networks.maternalnet.find_contacts(mother_uids))  # Children and Fetuses of treated mothers
        fetus_gt_cutoff = children_uids[sim.people.age[children_uids] > self.pars.fetus_age_cutoff_treat_eff]  # Fetuses above cutoff age
        fetus_st_cutoff = children_uids[sim.people.age[children_uids] <= self.pars.fetus_age_cutoff_treat_eff]  # Fetuses below cutoff age

        # Treat fetuses above cutoff age (e.g. in the last trimester) with a reduced treatment efficacy
        treat_gt_cutoff = self.pars.treat_eff_reduced.filter(fetus_gt_cutoff)

        # Treat fetuses below cutoff age (e.g. in the first two trimesters)
        treat_st_cutoff = self.pars.treat_eff.filter(fetus_st_cutoff)

        # Combine
        treat_uids = treat_gt_cutoff.concat(treat_st_cutoff)

        # Store results - Failure
        sim.diseases.syphilis.results['new_fetus_treated_failure'][sim.ti] = len(fetus_gt_cutoff) + len(fetus_st_cutoff) - len(treat_uids)
        # Store results - Unnecessary. Subtract successful uids below
        sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][sim.ti] += len(treat_uids)

        # Change birth outcomes for successfully treated unborn babies
        # If fetus is not yet infected, this will do nothing
        if len(treat_uids):
            # Birth outcomes must be modified to add probability of susceptible birth
            birth_outcomes = sim.diseases.syphilis.pars.birth_outcome_keys

            # Change states
            for oi, outcome in enumerate(birth_outcomes):
                ti_outcome = f'ti_{outcome}'
                vals = getattr(sim.diseases.syphilis, ti_outcome)
                successful_uids = treat_uids & vals.notnan.uids
                vals[successful_uids] = np.nan
                setattr(sim.diseases.syphilis, ti_outcome, vals)

                # Store results - Success
                sim.diseases.syphilis.results['new_fetus_treated_success'][sim.ti] += len(successful_uids)
                # Store results - Unnecessary. Subtract the uids that were successful
                sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][sim.ti] -= len(successful_uids)

        return

    def apply(self, sim):
        """
        Apply treatment
        """
        treat_uids = super().apply(sim)
        # Treat unborn babies of successfully treated mothers
        treat_pregnant_uids = sim.people.pregnancy.pregnant.uids & self.outcomes['successful']
        if len(treat_pregnant_uids):
            self.treat_fetus(sim, mother_uids=treat_pregnant_uids)
        return treat_uids


class NewbornTreatment(STITreatment):

    def change_states(self, sim, treat_succ):
        """ Change states of congenital cases """
        sim.diseases.syphilis.congenital[treat_succ] = False
        sim.diseases.syphilis.ti_congenital[treat_succ] = np.nan
        # sim.diseases.syphilis.susceptible[treat_succ] = True  # Leave this out for now

    def administer(self, sim, uids, return_format='dict'):
        """ Administer treatment to newborns """
        sus = sim.diseases.syphilis.susceptible
        con = sim.diseases.syphilis.congenital
        sus_uids = uids[sus[uids]]
        con_uids = uids[con[uids]]

        successful = self.pars.treat_eff.filter(con_uids)
        unsuccessful = np.setdiff1d(con_uids, successful)
        unnecessary = sus_uids

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
        elif return_format == 'array':
            output = successful

        return output



