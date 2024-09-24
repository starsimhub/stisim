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

        self.default_pars(
            dur_exp=ss.constant(0),  # Initial latent period: how long after exposure before you can infect others
            p_symp=[
                ss.bernoulli(p=0.3),  # Women: https://doi.org/10.1371/journal.pone.0143304
                ss.bernoulli(p=0.6),  # Men:   https://doi.org/10.1371/journal.pone.0143304
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(1/52, 12/52),  # Women:
                ss.lognorm_ex(0.25/52, 1/52),  # Men: symptoms should appear within days
            ],
            p_symp_clear=[
                ss.bernoulli(p=0.0),
                ss.bernoulli(p=0.0),
            ],
            p_symp_care=[
                ss.bernoulli(p=0.5),
                ss.bernoulli(p=0.8),
            ],
            dur_symp=[
                ss.lognorm_ex(1/12, 1/12),  # Women
                ss.lognorm_ex(1/12, 1/12),  # Men
            ],
            dur_asymp2clear=[
                ss.lognorm_ex(7/12, 2/12),  # Women
                ss.lognorm_ex(5/12, 3/12),  # Men
            ],
            dur_symp2clear=[
                ss.lognorm_ex(7/12, 2/12),  # Assumption
                ss.lognorm_ex(5/12, 3/12),  # Assumption
            ],
            dur_postsymp2clear=[
                ss.lognorm_ex(7/12, 1.5/12),  # Women
                ss.lognorm_ex(4/12, 1.5/12),  # Men
            ],
            dur_symp2care=[  # For those who test, how long before they seek care
                ss.lognorm_ex(1/12, 1/12),  # Women
                ss.lognorm_ex(1/52, 1/52),  # Men
            ],
            p_pid=ss.bernoulli(p=0.2),  # TODO
            dur_prepid=ss.lognorm_ex(1.5/12, 3/12),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01),
            eff_condom=0.0,
        )
        self.update_pars(pars, **kwargs)

        return

    def init_results(self):
        super().init_results()
        self.results += [
            ss.Result(self.name, 'rel_treat', self.sim.npts, dtype=float, scale=False),
        ]
        return

    def set_prognoses(self, uids, source_uids=None):
        super().set_prognoses(uids, source_uids)
        # Also pass on the relative treatability
        if 'ng_tx' in self.sim.interventions and source_uids is not None:
            self.sim.people.ng_tx.rel_treat[uids] = self.sim.people.ng_tx.rel_treat[source_uids]


