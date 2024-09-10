"""
Chlamydia trachomatis disease module
"""

import numpy as np
import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['DischargingSTI']


class DischargingSTI(SEIS):

    def __init__(self, pars=None, name='bv', init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        self.default_pars(
            p_symp=[
                ss.bernoulli(p=0.8),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            p_symp_care=[
                ss.bernoulli(p=0.3),
                ss.bernoulli(p=0.0),
            ],
            p_pid=ss.bernoulli(p=0),
            init_prev=ss.bernoulli(p=0.025),
            eff_condom=0.6,
        )
        self.update_pars(pars, **kwargs)

        return
 
