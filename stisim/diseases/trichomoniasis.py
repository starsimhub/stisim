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
            unit='month',
            dur_exp=ss.constant(0),

            # Symptoms
            p_symp=[0.4, 0.5],  # https://sti.bmj.com/content/76/4/248
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                [ss.dur(1, 'week'), ss.dur(12, 'week')],  # Women:
                [ss.dur(0.25, 'week'), ss.dur(1, 'week')],  # Men: symptoms should appear within days
            ],

            # Care seeking
            p_symp_care=[0.39, 0.27],
            dur_symp2care=[  # For those who test, how long before they seek care
                [ss.dur(2, 'month'), ss.dur(1, 'month')],  # Women
                [ss.dur(1, 'week'), ss.dur(2, 'week')],  # Men
            ],

            # Clearance
            dur_asymp2clear=[
                # Average duration of infection in women is at least 3–5 years and approximately 4 months for men
                # Source: https://sti.bmj.com/content/76/4/248
                [ss.dur(48, 'month'), ss.dur(6, 'month')],  # Women
                [ss.dur(26, 'week'), ss.dur(4, 'week')],  # Men
            ],
            dur_symp2clear=[
                [ss.dur(20, 'week'), ss.dur(4, 'week')],  # Women - assumptions
                [ss.dur(18, 'week'), ss.dur(4, 'week')],  # Men - assumptions
            ],
            p_clear=ss.bernoulli(p=0.1),  # Most women do not spontaneously clear, men do (https://sti.bmj.com/content/76/4/248)
            dur_persist=ss.years(100),

            p_pid=ss.bernoulli(p=0.025),
            dur_prepid=ss.lognorm_ex(ss.dur(6, 'week'), ss.dur(4, 'week')),
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
