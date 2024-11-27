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
            p_symp=[
                ss.bernoulli(p=0.4),  # Women: https://sti.bmj.com/content/76/4/248
                ss.bernoulli(p=0.5),  # Men: https://sti.bmj.com/content/76/4/248
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(12, 'week')),  # Women:
                ss.lognorm_ex(ss.dur(0.25, 'week'), ss.dur(1, 'week')),  # Men: symptoms should appear within days
            ],
            p_symp_clear=[
                ss.bernoulli(p=0.0),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            p_symp_care=[
                ss.bernoulli(p=0.39),  # Women
                ss.bernoulli(p=0.27),  # Men
            ],
            dur_symp=[
                ss.lognorm_ex(ss.dur(20, 'week'), ss.dur(5, 'week')),  # Women
                ss.lognorm_ex(ss.dur(10, 'week'), ss.dur(4, 'week')),  # Men
            ],
            dur_asymp2clear=[
                # Average duration of infection in women is at least 3–5 years and approximately 4 months for men
                # Source: https://sti.bmj.com/content/76/4/248
                ss.normal(ss.dur(15, 'week'), ss.dur(5, 'week')),  # Women
                ss.normal(ss.dur(26, 'week'), ss.dur(4, 'week')),  # Men
            ],
            dur_symp2clear=[
                # Assumptions...
                ss.lognorm_ex(ss.dur(20, 'week'), ss.dur(4, 'week')),  # Women
                ss.lognorm_ex(ss.dur(18, 'week'), ss.dur(4, 'week')),  # Men
            ],
            p_clear=[
                ss.bernoulli(p=0.1),  # Most women do not spontaneously clear
                ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
            ],
            dur_persist=ss.constant(ss.dur(100, 'year')),
            p_pid=ss.bernoulli(p=0.025),
            dur_prepid=ss.lognorm_ex(ss.dur(6, 'week'), ss.dur(4, 'week')),
            eff_condom=0.0,

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        return

    def set_duration(self, p, f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist, f_asymp, m_asymp, pid):
        """ Overwrite duration setting with persistence for female asymptomatic infection """
        # Firstly, set clearance time as per base class
        super().set_duration(p, f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist, f_asymp, m_asymp, pid)

        # Next, overwrite time of clearance for a subset of asymptomatic and postsymptomatic women
        potential_persist = f_symp_clear | f_asymp
        _, f_persist = p.p_clear[0].split(potential_persist)
        self.ti_clearance[f_persist] = self.ti_infected[f_persist] + self.pars.dur_persist.rvs(f_persist)
        return
