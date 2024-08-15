"""
Gonorrhea disease module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.seis import SEIS

__all__ = ['Gonorrhea']


class Gonorrhea(SEIS):

    def __init__(self, pars=None, **kwargs):
        super().__init__()

        self.default_pars(
            dur_exp=ss.lognorm_ex(1/52, 1/52),
            dur_inf=[
                ss.lognorm_ex(26/52, 5/52),  # Women
                ss.lognorm_ex(13/52, 5/52),  # Men
            ],
            p_symp=[
                ss.bernoulli(p=0.35),  # Women: https://doi.org/10.1371/journal.pone.0143304
                ss.bernoulli(p=0.65),  # Men:   https://doi.org/10.1371/journal.pone.0143304
            ],
            p_pid=ss.bernoulli(p=0.2),  # TODO
            dur_prepid=ss.lognorm_ex(3/52, 6/52),  # TODO

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        return


