"""
Gonorrhea disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
from stisim.diseases.sti import SEIS, by_sex, lognorm_by_sex

__all__ = ['Gonorrhea']


class Gonorrhea(SEIS):
    """
    Neisseria gonorrhoeae disease module.

    Models gonorrhea as an SEIS infection (see `sti.SEIS`) with no latent period, sex-specific
    symptoms and care seeking, and natural clearance. Drug resistance is tracked by the
    `sti.GonorrheaTreatment` intervention.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'ng'``)
        kwargs (dict): parameter values (e.g. ``beta_m2f=0.08``), and the arguments of `sti.SEIS`
    """
    def __init__(self, pars=None, name='ng', **kwargs):
        defaults = dict(
            dur_exp = ss.constant(0), # Infectious at once

            # Symptoms
            p_symp = by_sex(0.35, 0.65), # https://doi.org/10.1371/journal.pone.0143304
            dur_presymp = lognorm_by_sex([ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(1)]), # Men: symptoms should appear within days

            # Care seeking
            p_symp_care = by_sex(0.66, 0.83), # See Table 2: https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
            dur_symp2care = lognorm_by_sex([ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]),

            # Clearance
            dur_asymp2clear = lognorm_by_sex([ss.months(8), ss.months(2)], [ss.months(6), ss.months(3)]),
            dur_symp2clear = lognorm_by_sex([ss.months(9), ss.months(2)], [ss.months(6), ss.months(3)]), # Assumption

            # PID
            p_pid = 0.0, # TODO
            dur_prepid = ss.lognorm_ex(ss.months(1.5), ss.months(3)),

            # Transmission and initial conditions
            beta_m2f = 0.06,
            init_prev = 0.01,
            eff_condom = 0.9,
        )
        super().__init__(name=name, **sc.mergedicts(defaults, pars, kwargs))
        return

    def init_results(self):
        super().init_results()
        self.define_results(ss.Result('rel_treat', dtype=float, label='Drug resistance', auto_plot=False)) # Written by sti.GonorrheaTreatment
        return

    def step(self):
        """ Transmission; then the people infected take on the response to treatment of whoever infected them (see `sti.GonorrheaTreatment`) """
        new, sources, networks = super().step()
        if 'ng_tx' in self.sim.interventions:
            known = sources >= 0
            rel_treat = self.sim.interventions.ng_tx.rel_treat
            rel_treat[new[known]] = rel_treat[ss.uids(sources[known])]
        return new, sources, networks
