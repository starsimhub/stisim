"""
BV module
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.diseases.sti import SEIS

__all__ = ['BV']


class BV(SEIS):

    def __init__(self, pars=None, name='bv', init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)

        self.define_pars(
            unit='month',
            dur_exp=ss.constant(0),  # No exposure period
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
            p_douching=ss.bernoulli(p=0.1),  # Share of population douching
            p_poor_menstrual_hygiene=ss.bernoulli(p=0.1),  # Share of population with poor menstrual hygiene
            p_base=0.1,                 # Used to calculate the baseline (intercept) probability of spontaneous occurrence
            p_spontaneous=sc.objdict(
                douching=3,             # OR of BV for douching
                n_partners_12m=2,       # OR of BV for each additional partner in the past 12 months - not working yet
                poor_menstrual_hygiene=2,    # OR of BV for poor menstrual hygiene
            )
        )
        self.update_pars(pars, **kwargs)

        # States that elevate risk of BV
        self.define_states(
            ss.BoolArr('douching'),
            ss.FloatArr('n_partners_12m', 0),
            ss.BoolArr('poor_menstrual_hygiene'),
        )

        return

    def _get_uids(self, upper_age=None):
        """ Get uids of females younger than upper_age """
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        within_age = people.age < upper_age
        return (within_age & people.female).uids

    def bv_sus(self):
        return self.sim.people.female & (self.sim.people.age > 15) & self.susceptible

    def set_hygiene_states(self, upper_age=None):
        """ Set vaginal hygiene states """
        f_uids = self._get_uids(upper_age=upper_age)
        self.douching[f_uids] = self.pars.p_douching.rvs(f_uids)
        self.poor_menstrual_hygiene[f_uids] = self.pars.p_poor_menstrual_hygiene.rvs(f_uids)
        return

    def init_post(self):
        """ Initialize with sim properties """
        for state in self.states:
            if not state.initialized:
                state.init_vals()
        self.initialized = True

        # Set hygiene states, but don't create any initial infections
        self.set_hygiene_states()

        return

    def spontaneous(self, uids):
        """ Create new cases via spontaneous occurrence """
        # Set intercept
        p = sc.dcp(self.pars.p_spontaneous)
        intercept = -np.log(1/self.pars.p_base-1)  # Use a transformation consistent with the logistic regression
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
        # Update VMB-relevant states
        self.set_hygiene_states()

        # Create new cases via spontaneous occurrence
        uids = self.bv_sus().uids
        p_bv = self.spontaneous(uids)
        self.pars.p_bv.set(p_bv)
        bv_cases = self.pars.p_bv.filter(uids)

        # Set prognoses
        if len(bv_cases):
            self.set_prognoses(bv_cases)

        return
