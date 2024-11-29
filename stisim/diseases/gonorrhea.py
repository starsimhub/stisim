"""
Gonorrhea disease module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import SEIS

__all__ = ['Gonorrhea']


class Gonorrhea(SEIS):

    def __init__(self, name='ng', pars=None, init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        self.define_pars(
            dur_exp=ss.constant(0),  # Initial latent period: how long after exposure before you can infect others
            p_symp=[
                ss.bernoulli(p=0.35),  # Women: https://doi.org/10.1371/journal.pone.0143304
                ss.bernoulli(p=0.65),  # Men:   https://doi.org/10.1371/journal.pone.0143304
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(12, 'week')),  # Women:
                ss.lognorm_ex(ss.dur(0.25, 'week'), ss.dur(1, 'week')),  # Men: symptoms should appear within days
            ],
            p_symp_clear=[  # NOT USED
                ss.bernoulli(p=0.0),
                ss.bernoulli(p=0.0),
            ],
            p_symp_care=[  # See Table 2: https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
                ss.bernoulli(p=0.66),
                ss.bernoulli(p=0.83),
            ],
            dur_symp=[  # NOT USED
                ss.constant(ss.dur(1)),  # Duration of symptoms: not used, as symptoms don't self-resolve
                ss.constant(ss.dur(1)),  # Duration of symptoms
            ],
            dur_asymp2clear=[
                ss.lognorm_ex(ss.dur(7, 'month'), ss.dur(2, 'month')),  # Women
                ss.lognorm_ex(ss.dur(5, 'month'), ss.dur(3, 'month')),  # Men
            ],
            dur_symp2clear=[
                ss.lognorm_ex(ss.dur(7, 'month'), ss.dur(2, 'month')),  # Assumption
                ss.lognorm_ex(ss.dur(5, 'month'), ss.dur(3, 'month')),  # Assumption
            ],
            dur_postsymp2clear=[  # NOT USED
                ss.constant(ss.dur(1)),  # Duration of infection after symptom clearance: not used, as symptoms don't self-resolve
                ss.constant(ss.dur(1)),  # As above
            ],
            dur_symp2care=[  # For those who test, how long before they seek care
                ss.lognorm_ex(ss.dur(1, 'month'), ss.dur(1, 'month')),  # Women
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(2, 'week')),  # Men
            ],
            p_pid=ss.bernoulli(p=0.2),  # TODO
            dur_prepid=ss.lognorm_ex(ss.dur(1.5, 'month'), ss.dur(3, 'month')),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01),
            eff_condom=0.0,
        )
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


