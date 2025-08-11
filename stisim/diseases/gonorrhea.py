"""
Gonorrhea disease module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import SEIS, STIPars

__all__ = ['Gonorrhea', 'NGPars']


class NGPars(STIPars):
    def __init__(self, **kwargs):
        super().__init__()
        self.dur_exp = ss.constant(0)  # How long after exposure before you can infect others

        # Symptoms
        self.p_symp = [0.35, 0.65]  # https://doi.org/10.1371/journal.pone.0143304
        self.dur_presymp = [  # For those who develop symptoms, how long before symptoms appear
            [ss.dur(1, 'week'), ss.dur(12, 'week')],  # Women:
            [ss.dur(0.25, 'week'), ss.dur(1, 'week')],  # Men: symptoms should appear within days
        ]

        # Care seeking
        self.p_symp_care = [0.66, 0.83]  # See Table 2: https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
        self.dur_symp2care = [  # For those who test, how long before they seek care
            [ss.dur(2, 'month'), ss.dur(1, 'month')],  # Women
            [ss.dur(1, 'week'), ss.dur(2, 'week')],  # Men
        ]

        # Clearance: lognormal distribution
        self.dur_asymp2clear = [
            [ss.dur(8, 'month'), ss.dur(2, 'month')],  # Women
            [ss.dur(6, 'month'), ss.dur(3, 'month')],  # Men
        ]
        self.dur_symp2clear = [
            [ss.dur(9, 'month'), ss.dur(2, 'month')],  # Assumption
            [ss.dur(6, 'month'), ss.dur(3, 'month')],  # Assumption
        ]

        self.p_pid = ss.bernoulli(p=0.0)  # TODO
        self.dur_prepid = ss.lognorm_ex(ss.dur(1.5, 'month'), ss.dur(3, 'month'))

        # Initial conditions
        self.init_prev = ss.bernoulli(p=0.01)
        self.eff_condom = 0.9
        self.update(kwargs)
        return


class Gonorrhea(SEIS):

    def __init__(self, name='ng', pars=None, init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)
        default_pars = NGPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('rel_treat', dtype=float, label='Drug resistance')
        )
        return

    def set_prognoses(self, uids, source_uids=None):
        super().set_prognoses(uids, source_uids)
        # Also pass on the relative treatability
        if 'ng_tx' in self.sim.interventions and source_uids is not None:
            self.sim.people.ng_tx.rel_treat[uids] = self.sim.people.ng_tx.rel_treat[source_uids]


