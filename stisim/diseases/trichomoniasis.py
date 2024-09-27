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

        self.default_pars(
            dur_exp=ss.constant(0),
            p_symp=[
                ss.bernoulli(p=0.4),  # Women: https://sti.bmj.com/content/76/4/248
                ss.bernoulli(p=0.5),  # Men: https://sti.bmj.com/content/76/4/248
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(1/52, 12/52),  # Women:
                ss.lognorm_ex(0.25/52, 1/52),  # Men: symptoms should appear within days
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
                ss.lognorm_ex(20/52, 5/52),  # Women
                ss.lognorm_ex(10/52, 4/52),  # Men
            ],
            dur_asymp2clear=[
                # Average duration of infection in women is at least 3–5 years and approximately 4 months for men
                # Source: https://sti.bmj.com/content/76/4/248
                ss.normal(150/52, 20/52),  # Women
                ss.normal(26/52, 4/52),  # Men
            ],
            dur_symp2clear=[
                # Assumptions...
                ss.lognorm_ex(20/52, 4/52),  # Women
                ss.lognorm_ex(15/52, 5/52),  # Men
            ],
            p_clear=[
                ss.bernoulli(p=0.1),  # Most women do not spontaneously clear
                ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
            ],
            p_pid=ss.bernoulli(p=0.025),
            dur_prepid=ss.lognorm_ex(1.5/12, 1/12),
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
        dt = self.sim.dt
        _, f_persist = p.p_clear[0].split(potential_persist)
        self.ti_clearance[f_persist] = self.ti_infected[f_persist] + 1e2/dt
        return
