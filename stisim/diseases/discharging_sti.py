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
                # Amsel score or Nugent score, clinical diagnosis for BV, cal be asymptomatic
                # Asymptomatic still has the same risk for STIs/HIV
                #
                ss.bernoulli(p=0.1),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            p_symp_care=[
                ss.bernoulli(p=0.4),
                ss.bernoulli(p=0.0),
            ],
            dur_asymp2clear=[  # Duration of untreated asymptomatic infection (excl initial latent)
                ss.uniform(1/52, 20/52),  # Women
                ss.constant(100),  # Men
            ],
            dur_symp2clear=[  # Duration of untreated symptomatic infection (excl initial latent)
                ss.uniform(1/52, 20/52),  # Women
                ss.constant(100),  # Men
            ],
            # Care-seeking based on partner dynamics - if their partner notices changes
            dur_symp2care=[  # For those who test, how long before they seek care
                ss.uniform(1/52, 5/52),  # Women
                ss.constant(100),  # Men
            ],
            p_pid=ss.bernoulli(p=0),
            init_prev=ss.bernoulli(p=0.025),
            eff_condom=0.0,
        )
        self.update_pars(pars, **kwargs)

        return
 
