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
            dur_exp2inf=ss.constant(0),
            dur_inf2clear=[
                ss.lognorm_ex(6/12, 1.2/12),  # Women
                ss.lognorm_ex(4/12, 1.2/12),  # Men
            ],
            p_symp=[
                ss.bernoulli(p=0.35),  # Women: https://doi.org/10.1371/journal.pone.0143304
                ss.bernoulli(p=0.65),  # Men:   https://doi.org/10.1371/journal.pone.0143304
            ],
            p_pid=ss.bernoulli(p=0.2),  # TODO
            dur_inf2pid=ss.lognorm_ex(1.5/12, 1/12),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01),
            eff_condom=0.4,
        )
        self.update_pars(pars, **kwargs)

        return

    def set_prognoses(self, uids, source_uids=None):
        super().set_prognoses(uids, source_uids)
        # Also pass on the relative treatability
        if 'gonorrheatreatment' in self.sim.interventions:
            self.sim.interventions.gonorrheatreatment.rel_treat[uids] = self.sim.interventions.gonorrheatreatment.rel_treat[source_uids]


