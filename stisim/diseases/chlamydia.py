"""
Chlamydia trachomatis disease module
"""

import numpy as np
import starsim as ss

__all__ = ['Chlamydia']


class Chlamydia(ss.Infection):

    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.requires = 'structuredsexual'

        self.default_pars(
            # Bacterial load dynamics
            init_load=1,
            peak_load=10e7,
            time_to_peak=8/52,
            half_life=ss.lognorm_ex(2.5/52, 0.5/52),
            ct_beta=0.5,  # Growth rate in logistic function mapping CT load to rel_trans
            model_ct_load=True,  # Whether to model CT load
            dur_inf=ss.lognorm_ex(60/52, 5/52),  # Duration of infection - only used if model_ct_load is False

            # Transmission
            beta=1.0,  # Placeholder
            beta_m2f=None,
            beta_f2m=None,
            beta_m2c=None,

            # Symptoms
            p_symp=[
                ss.bernoulli(p=0.375),
                ss.bernoulli(p=0.375),
            ],
            p_pid_symp=ss.bernoulli(p=0.25),
            p_pid_asymp=ss.bernoulli(p=0.25),
            dur_presymp=ss.lognorm_ex(2/52, 1/52),
            dur_prepid=ss.lognorm_ex(3/52, 6/52),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            # Bacterial load
            ss.FloatArr('ct_load'),
            ss.FloatArr('ct_peak_time'),
            ss.FloatArr('ct_growth_rate'),
            ss.FloatArr('ct_decay_rate'),
            ss.FloatArr('ct_half_life'),

            # Natural history
            ss.BoolArr('asymptomatic'),
            ss.BoolArr('symptomatic'),
            ss.BoolArr('pid'),
            ss.FloatArr('dur_inf'),
            ss.FloatArr('ti_symptomatic'),
            ss.FloatArr('ti_pid'),
            ss.FloatArr('ti_clearance'),

            # Immunity
            ss.FloatArr('immunity', default=0.0),
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

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.beta_m2f is not None:
            self.pars.beta['structuredsexual'][0] *= self.pars.beta_m2f
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        return

    def init_post(self):
        """ Make initial cases """
        super().init_post()
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        return

    def clear_infection(self, uids):
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.pid[uids] = False
        self.susceptible[uids] = True
        self.ti_clearance[uids] = self.sim.ti

    def update_pre(self):
        """ Updates prior to interventions """
        ti = self.sim.ti
        dt = self.sim.dt

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Progress symptoms
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        self.asymptomatic[new_symptomatic] = False
        self.symptomatic[new_symptomatic] = True
        self.ti_symptomatic[new_symptomatic] = ti

        # Progress PID
        new_pid = (self.infected & (self.ti_pid <= ti)).uids
        self.pid[new_pid] = False

        # Update CT load and scale res_sus
        if self.pars.model_ct_load:
            self.update_ct_load(self.infected.uids)
            self.rel_trans[self.infected] = 2 / (1 + np.exp(-self.pars.ct_beta*np.log10(self.ct_load[self.infected]))) - 1

        return

    def update_results(self):
        super().update_results()
        return

    def finalize_results(self):
        super().finalize_results()
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
        Set initial prognoses for adults newly infected with syphilis
        """
        super().set_prognoses(uids, source_uids)

        ti = self.sim.ti
        dt = self.sim.dt
        ppl = self.sim.people
        p = self.pars

        self.susceptible[uids] = False
        self.infected[uids] = True
        self.asymptomatic[uids] = True
        self.ti_infected[uids] = ti

        # Calculate duration of infection
        if self.pars.model_ct_load:
            dur_inf = self.set_ct_load(uids)
        else:
            dur_inf = self.pars.dur_inf.rvs(uids)
        self.dur_inf[uids] = dur_inf

        # Symptoms
        m_uids = ppl.male.uids.intersect(uids)
        f_uids = ppl.female.uids.intersect(uids)
        m_symp = p.p_symp[0].filter(m_uids)
        f_symp, f_asymp = p.p_symp[1].split(f_uids)
        dur_presymp_m = np.minimum(p.dur_presymp.rvs(m_symp), self.dur_inf[m_symp])
        dur_presymp_f = np.minimum(p.dur_presymp.rvs(f_symp), self.dur_inf[f_symp])
        self.ti_symptomatic[m_symp] = ti + dur_presymp_m/dt
        self.ti_symptomatic[f_symp] = ti + dur_presymp_f/dt

        pid_symp = p.p_pid_symp.filter(f_symp)
        pid_asymp = p.p_pid_asymp.filter(f_asymp)
        dur_prepid_symp = np.minimum(p.dur_prepid.rvs(pid_symp), self.dur_inf[pid_symp])
        dur_prepid_asymp = np.minimum(p.dur_prepid.rvs(pid_asymp), self.dur_inf[pid_asymp])
        self.ti_pid[pid_symp] = ti + dur_prepid_symp/dt
        self.ti_pid[pid_asymp] = ti + dur_prepid_asymp/dt

        # Determine when people recover
        self.ti_clearance[uids] = ti + dur_inf/dt

        return


