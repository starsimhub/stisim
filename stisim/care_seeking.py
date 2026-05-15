"""
Define care-seeking propensity module for STIsim.

Provides a cross-disease individual-level care-seeking propensity that
can be used by testing interventions across HIV, syphilis, and other diseases.
"""

import numpy as np
import starsim as ss
import sciris as sc

__all__ = ['CareSeeking']


class CareSeeking(ss.Module):
    """
    Cross-disease care-seeking propensity for each individual.

    Each agent gets a baseline propensity drawn from a distribution at birth.
    Men have lower propensity by default (controlled by rel_care_m2f).
    The propensity may be temporarily modified by factors like pregnancy.
    """

    def __init__(self, pars=None, **kwargs):
        super().__init__(name='care_seeking')
        self.define_pars(
            propensity_dist=ss.lognorm_ex(1, 0.3),  # Mean=1, moderate variation
            rel_care_m2f=0.7,  # Men are less likely to seek care than women
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.FloatArr('baseline_propensity', default=1.0),  # Fixed at birth, doesn't change
            ss.FloatArr('propensity', default=1.0),           # Working value, may be modified by pregnancy etc.
        )

    def init_post(self):
        super().init_post()
        self.set_care_seeking_states()

    def set_care_seeking_states(self, upper_age=None):
        """
        Set baseline care-seeking propensity for each agent.
        Called at init (upper_age=None for all agents) and each step
        (upper_age=dt to initialize newborns).
        """
        ppl = self.sim.people
        if upper_age is None:
            uids = ppl.auids
        else:
            uids = (ppl.age <= upper_age).uids

        if len(uids) == 0:
            return

        # Draw base propensity from distribution
        base = self.pars.propensity_dist.rvs(uids)

        # Apply sex differential: men have lower care-seeking
        male_uids = ppl.male[uids].uids
        base[male_uids] *= self.pars.rel_care_m2f

        self.baseline_propensity[uids] = base
        self.propensity[uids] = base

    def step(self):
        """
        Initialize propensity for newborns, then reset all propensities
        to baseline before applying temporary modifiers.
        """
        # Set propensity for newly born agents
        self.set_care_seeking_states(upper_age=self.t.dt_year)

        # Reset to baseline (undoes any temporary modifiers from previous step)
        self.propensity[:] = self.baseline_propensity[:]

        # Future: apply temporary modifiers, e.g.:
        # if 'pregnancy' in self.sim.demographics:
        #     pregnant = self.sim.demographics.pregnancy.pregnant
        #     self.propensity[pregnant] *= pregnancy_scale
