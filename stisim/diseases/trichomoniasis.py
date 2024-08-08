"""
Trichomoniasis disease module
"""

import numpy as np
import starsim as ss
 
__all__ = ['Trichomoniasis']


class Trichomoniasis(ss.Infection):

    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.requires = 'structuredsexual'

        self.default_pars(
            # Transmission
            beta=1.0,  # Placeholder
            beta_m2f=None,
            beta_f2m=None,

            # Natural history
            dur_inf=ss.lognorm_ex(60/52, 5/52),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
        )

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.beta_m2f is not None:
            self.pars.beta['structuredsexual'][0] *= self.pars.beta_m2f
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        return

    def update_pre(self):
        """ Updates prior to interventions """
        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Clear infections
        newly_cleared = (self.infected & (self.ti_clearance <= self.sim.ti)).uids
        self.infected[newly_cleared] = False
        self.susceptible[newly_cleared] = True

        return

    def update_results(self):
        super().update_results()
        return

    def finalize_results(self):
        super().finalize_results()
        return

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected with syphilis
        """
        super().set_prognoses(uids, source_uids)

        ti = self.sim.ti
        dt = self.sim.dt

        self.susceptible[uids] = False
        self.infected[uids] = True
        self.ti_infected[uids] = ti

        # Calculate duration of infection
        dur_inf = self.pars.dur_inf.rvs(uids)

        # Determine when people recover
        self.ti_clearance[uids] = ti + dur_inf/dt

        return


