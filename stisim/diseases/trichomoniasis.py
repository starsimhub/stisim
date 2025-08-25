"""
Trichomoniasis disease module
"""

import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['Trichomoniasis']


class Trichomoniasis(SEIS):

    def __init__(self, pars=None, name='tv', init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)
        self.requires = 'structuredsexual'

        self.define_pars(
            dt='month',
            dur_exp=ss.constant(0),

            # Symptoms
            p_symp=[0.4, 0.5],  # https://sti.bmj.com/content/76/4/248
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                [ss.weeks(1), ss.weeks(12)],  # Women:
                [ss.weeks(0.25), ss.weeks(1)],  # Men: symptoms should appear within days
            ],

            # Care seeking
            p_symp_care=[0.39, 0.27],
            dur_symp2care=[  # For those who test, how long before they seek care
                [ss.months(2), ss.months(1)],  # Women
                [ss.weeks(1), ss.weeks(2)],  # Men
            ],

            # Clearance
            dur_asymp2clear=[
                # Average duration of infection in women is at least 3–5 years and approximately 4 months for men
                # Source: https://sti.bmj.com/content/76/4/248
                [ss.months(48), ss.months(6)],  # Women
                [ss.weeks(26), ss.weeks(4)],  # Men
            ],
            dur_symp2clear=[
                [ss.weeks(20), ss.weeks(4)],  # Women - assumptions
                [ss.weeks(18), ss.weeks(4)],  # Men - assumptions
            ],
            p_clear=ss.bernoulli(p=0.1),  # Most women do not spontaneously clear, men do (https://sti.bmj.com/content/76/4/248)
            dur_persist=ss.years(100),

            p_pid=ss.bernoulli(p=0.00),
            dur_prepid=ss.lognorm_ex(ss.weeks(6), ss.weeks(4)),
            eff_condom=0.0,

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        return

    def set_duration(self, p, symp, asymp, pid):
        """ Overwrite duration setting with persistence for female asymptomatic infection """
        # Firstly, set clearance time as per base class
        super().set_duration(p, symp, asymp, pid)

        # Next, overwrite time of clearance for a subset of asymptomatic women
        potential_persist = asymp[self.sim.people.female[asymp]]
        _, f_persist = p.p_clear.split(potential_persist)
        self.ti_clearance[f_persist] = self.ti_infected[f_persist] + self.pars.dur_persist
        return
