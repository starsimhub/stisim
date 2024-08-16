"""
Trichomoniasis disease module
"""

import starsim as ss
from stisim.diseases.seis import SEIS

__all__ = ['Trichomoniasis']


class Trichomoniasis(SEIS):

    def __init__(self, pars=None, name='tv', **kwargs):
        super().__init__(name=name)
        self.requires = 'structuredsexual'

        self.default_pars(
            dur_exp2inf=ss.lognorm_ex(1/52, 1/52),
            dur_inf2clear=[
                # Average duration of infection in women is at least 3–5 years and approximately 4 months for men
                # Source: https://sti.bmj.com/content/76/4/248
                ss.lognorm_ex(20/52, 4/52),  # Men
                ss.lognorm_ex(104/52, 15/52),  # Women
            ],
            p_symp=[
                ss.bernoulli(p=0.4),  # Women: https://sti.bmj.com/content/76/4/248
                ss.bernoulli(p=0.5),  # Men: https://sti.bmj.com/content/76/4/248
            ],
            p_clear=[
                ss.bernoulli(p=0.1),  # Most women do not spontaneously clear
                ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
            ],
            p_pid=ss.bernoulli(p=0.025),
            dur_inf2pid=ss.lognorm_ex(1.5/12, 1/12),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        return

    def set_duration(self, p, f_uids, m_uids):
        """ Set duration of infection"""
        f_clear, f_persist = p.p_clear[0].split(f_uids)
        m_clear = p.p_clear[1].filter(m_uids)
        dur_inf_f = p.dur_inf2clear[0].rvs(f_clear)
        dur_inf_m = p.dur_inf2clear[1].rvs(m_clear)
        self.ti_clearance[f_clear] = self.ti_infected[f_clear] + dur_inf_f/self.sim.dt
        self.ti_clearance[m_clear] = self.ti_infected[m_clear] + dur_inf_m/self.sim.dt
        self.dur_inf[f_persist] = 100
        self.dur_inf[f_clear] = dur_inf_f
        self.dur_inf[m_clear] = dur_inf_m

    # def set_prognoses(self, uids, source_uids=None):
    #     """
    #     Set initial prognoses for adults newly infected
    #     """
    #     ppl = self.sim.people
    #     p = self.pars
    #     m_uids = ppl.male.uids.intersect(uids)
    #     f_uids = ppl.female.uids.intersect(uids)
    #
    #     self.set_exposure(uids)
    #     f_symp, m_symp = self.set_symptoms(p, f_uids, m_uids)
    #     self.set_duration(p, f_uids, m_uids)
    #     self.set_care_seeking(p, f_symp, m_symp)
    #     self.set_pid(p, f_uids)
    #     self.set_pid_care_seeking(p, f_uids)
    #
    #     # Determine when people recover
    #     self.ti_clearance[uids] = self.ti_infected[uids] + self.dur_inf[uids]/self.sim.dt
    #
    #     return
    #

