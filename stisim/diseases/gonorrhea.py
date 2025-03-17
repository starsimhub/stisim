"""
Gonorrhea disease module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import SEIS

__all__ = ['Gonorrhea']


class Gonorrhea(SEIS):

    def __init__(self, name='ng', init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        self.define_pars(
            dur_exp=ss.constant(0),  # Initial latent period: how long after exposure before you can infect others

            # Symptoms
            p_symp=[0.35, 0.65],  # https://doi.org/10.1371/journal.pone.0143304
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                [ss.weeks(1), ss.weeks(12)],  # Women:
                [ss.weeks(0.25), ss.weeks(1)],  # Men: symptoms should appear within days
            ],

            # Care seeking
            p_symp_care=[0.66, 0.83],  # See Table 2: https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
            dur_symp2care=[  # For those who test, how long before they seek care
                [ss.months(1), ss.months(1)],  # Women
                [ss.weeks(1), ss.weeks(2)],  # Men
            ],

            # Clearance: lognormal distribution
            dur_asymp2clear=[
                [ss.months(7), ss.months(2)],  # Women
                [ss.months(5), ss.months(3)],  # Men
            ],
            dur_symp2clear=[
                [ss.months(7), ss.months(2)],  # Assumption
                [ss.months(5), ss.months(3)],  # Assumption
            ],

            p_pid=ss.bernoulli(p=0.2),  # TODO
            dur_prepid=ss.lognorm_ex(ss.months(1.5), ss.months(3)),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01),
            eff_condom=0.0,
        )
        self.update_pars(**kwargs)

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


