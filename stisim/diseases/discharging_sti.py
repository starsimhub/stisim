"""
Chlamydia trachomatis disease module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import SEIS

__all__ = ['DischargingSTI']


class DischargingSTI(SEIS):

    def __init__(self, pars=None, name='bv', init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        self.define_pars(
            unit='month',
            p_symp=[
                # Amsel score or Nugent score, clinical diagnosis for BV, can be asymptomatic
                # Asymptomatic still has the same risk for STIs/HIV
                ss.bernoulli(p=0.1),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            p_symp_care=[
                ss.bernoulli(p=0.4),
                ss.bernoulli(p=0.0),
            ],
            dur_asymp2clear=[  # Duration of untreated asymptomatic infection (excl initial latent)
                ss.uniform(ss.dur(1, 'week'), ss.dur(18, 'week')),  # Women
                ss.constant(ss.dur(100, 'year')),  # Men
            ],
            dur_symp2clear=[  # Duration of untreated symptomatic infection (excl initial latent)
                ss.uniform(ss.dur(1, 'week'), ss.dur(18, 'week')),  # Women
                ss.constant(ss.dur(100, 'year')),  # Men
            ],
            # Care-seeking based on partner dynamics - if their partner notices changes
            dur_symp2care=[  # For those who test, how long before they seek care
                ss.uniform(ss.dur(1, 'week'), ss.dur(18, 'week')),  # Women
                ss.constant(ss.dur(100, 'year')),  # Men
            ],
            p_pid=ss.bernoulli(p=0),

            # Transmission parameters
            init_prev=ss.bernoulli(p=0.025),
            eff_condom=0.0,
            rel_beta_f2m=0,

            # Spontaneous occurrence parameters. These will be used within a logistic regression
            # model to calculate the probability of spontaneous occurrence. The model is flexible
            # but should always include an intercept term.
            p_bv=ss.bernoulli(p=0.01),  # Probability of BV in the general population. Overwritten by the model below
            p_spontaneous=sc.objdict(
                intercept=-np.log(99),   # Set so that 1/(1+np.exp(-intercept)) = 0.01
                douching=3,             # OR of BV for douching
                prior_bv=2,             # OR of BV for prior episodes of BV
                n_partners_12m=2,       # OR of BV for each additional partner in the previous 12 months
                poor_menstrual_hygiene=2,    # OR of BV for poor menstrual hygiene
            )
        )
        self.update_pars(pars, **kwargs)

        self.define_states(
            ss.BoolArr('douching'),
            ss.BoolArr('prior_bv'),
            ss.FloatArr('n_partners_12m', 0),
            ss.BoolArr('poor_menstrual_hygiene'),
        )

        return

    def spontaneous(self, uids):
        """
        Create new cases via spontaneous occurrence
        """
        # Set intercept
        p = sc.dcp(self.pars.p_spontaneous)
        intercept = p.pop('intercept')
        rhs = np.full_like(uids, fill_value=intercept, dtype=float)

        # Add all covariates
        for term, val in p.items():
            rhs += val * getattr(self, term)[uids]

        # Calculate the probability of spontaneous occurrence of BV
        p_bv = 1 / (1+np.exp(-rhs))

        return p_bv

    def step(self):
        """
        Create new cases (via both sexual transmission and spontaneous occurence), and set prognoses
        """
        # Create new cases via sexual transmission
        new_cases, sources, networks = self.infect()

        # Create new cases via spontaneous occurrence
        f_uids = (self.sim.people.female & (self.sim.people.age > 15)).uids
        p_bv = self.spontaneous(f_uids)
        self.pars.p_bv.set(p_bv)
        bv_cases = self.pars.p_bv.filter(f_uids)
        new_cases = new_cases | bv_cases

        # Set prognoses
        if len(new_cases):
            self.set_prognoses(new_cases)

        return new_cases, sources, networks
