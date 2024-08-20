"""
Chlamydia trachomatis disease module
"""

import numpy as np
import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['Chlamydia', 'ChlamydiaBL']


class Chlamydia(SEIS):

    def __init__(self, pars=None, name='ct', **kwargs):
        super().__init__(name=name)

        self.default_pars(
            dur_exp2inf=ss.constant(1/52),
            dur_inf2clear=[
                ss.lognorm_ex(60/52, 5/52),  # Women: 433 days (https://doi.org/10.1016/j.epidem.2010.04.002)
                ss.lognorm_ex(60/52, 5/52),  # Men: as above
            ],
            p_symp=[
                ss.bernoulli(p=0.375),  # Women: 62.5% asymptomatic (https://doi.org/10.1016/j.epidem.2010.04.002)
                ss.bernoulli(p=0.375),  # Men: as above
            ],
            p_pid=ss.bernoulli(p=0.2),  # Assumption used in https://doi.org/10.1086/598983, based on https://doi.org/10.1016/s0029-7844(02)02118-x
            dur_inf2pid=ss.lognorm_ex(1.5/12, 1/12),
            init_prev=ss.bernoulli(p=0.01),
            eff_condom=0.6,  # doi:10.1001/archpedi.159.6.536
        )
        self.update_pars(pars, **kwargs)

        return


class ChlamydiaBL(Chlamydia):

    def __init__(self, pars=None, **kwargs):
        super().__init__()

        self.default_pars(
            # Bacterial load dynamics
            init_load=1,
            peak_load=10e7,
            time_to_peak=8/52,
            half_life=ss.lognorm_ex(2.5/52, 0.5/52),
            ct_beta=0.5,  # Growth rate in logistic function mapping CT load to rel_trans
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            # Bacterial load
            ss.FloatArr('ct_load'),
            ss.FloatArr('ct_peak_time'),
            ss.FloatArr('ct_growth_rate'),
            ss.FloatArr('ct_decay_rate'),
            ss.FloatArr('ct_half_life'),
        )

        return

    def update_ct_load(self, uids):
        ct_incr_uids = (self.ct_peak_time > self.sim.ti).uids.intersect(uids)
        ct_decr_uids = (self.ct_peak_time <= self.sim.ti).uids.intersect(uids)
        if len(ct_incr_uids):
            self.ct_load[ct_incr_uids] *= np.exp(self.ct_growth_rate[ct_incr_uids])
        if len(ct_decr_uids):
            self.ct_load[ct_decr_uids] *= np.exp(-self.ct_decay_rate[ct_decr_uids])
        return

    def update_pre(self):
        """ Updates prior to interventions """
        super().update_pre()
        # Update CT load and scale rel_trans
        self.update_ct_load(self.infected.uids)
        self.rel_trans[self.infected] = 2 / (1 + np.exp(-self.pars.ct_beta*np.log10(self.ct_load[self.infected]))) - 1

        return

    def set_ct_load(self, uids):
        """ Bacterial load dynamics """
        ti = self.sim.ti
        dt = self.sim.dt
        p = self.pars
        timesteps_to_peak = p.time_to_peak/dt
        if timesteps_to_peak < 1:
            self.ct_load[uids] = p.peak_load
        else:
            self.ct_load[uids] = p.init_load

        self.ct_peak_time[uids] = ti + timesteps_to_peak
        self.ct_half_life[uids] = p.half_life.rvs(uids)
        self.ct_growth_rate[uids] = np.log(p.peak_load/p.init_load)/timesteps_to_peak
        self.ct_decay_rate[uids] = np.log(2) / (self.ct_half_life[uids]/dt)
        dur_inf = (-np.log(p.init_load/p.peak_load)/self.ct_decay_rate[uids])*dt

        return dur_inf

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected
        """
        ppl = self.sim.people
        p = self.pars
        f_uids = ppl.female.uids.intersect(uids)
        m_uids = ppl.male.uids.intersect(uids)

        self.set_exposure(uids)
        self.set_symptoms(p, f_uids, m_uids)
        dur_inf = self.set_ct_load(uids)
        self.dur_inf[uids] = dur_inf
        self.set_pid(p, f_uids)

        # Determine when people recover
        self.ti_clearance[uids] = self.ti_infected[uids] + self.dur_inf[uids]/self.sim.dt

        return


