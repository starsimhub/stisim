"""
Define gonorrhea interventions for STIsim
"""

import starsim as ss
import numpy as np
from stisim.interventions.base_interventions import STITreatment


# %% Gonorrhea classes
__all__ = ["GonorrheaTreatment", "UpdateDrugs"]


class GonorrheaTreatment(STITreatment):
    """Treatment for gonorrhea infection with antimicrobial resistance tracking.

    Successful treatment clears infection immediately. Unsuccessful treatment
    reduces infection duration and relative transmissibility but lowers the
    agent's future treatment response (``rel_treat``). Unnecessary treatment
    also reduces ``rel_treat``, modeling resistance accumulation.

    Args:
        pars (dict): Parameter overrides (e.g. ``base_treat_eff``,
            ``rel_treat_unsucc``, ``rel_treat_unneed``).
        eligibility (func): Function returning eligible UIDs.
        name (str): Intervention name.
        **kwargs: Additional parameter overrides.
    """
    def __init__(self, pars=None, eligibility=None, name=None, *args, **kwargs):
        super().__init__(diseases='ng', name=name, eligibility=eligibility, *args)
        self.requires = ['ng', 'structuredsexual']
        self.define_pars(
            base_treat_eff=0.96,
            rel_treat_unsucc=0.01,
            rel_treat_unneed=0.005,
        )
        self.update_pars(pars, **kwargs)

        # Treatment outcome: probability is set each step from base_treat_eff and the agent's rel_treat
        self.treat_eff = ss.bernoulli(p=0)

        # States
        self.define_states(
            ss.FloatArr('rel_treat', default=1),  # How well a person will respond to treatment
        )

    def get_treat_eff(self, uids):
        """ Set the probability of treatment success from base_treat_eff and each agent's rel_treat """
        rt = self.rel_treat[uids]
        rt = np.where(np.isnan(rt), 1.0, rt)  # default for agents whose slot was inactive at init
        new_treat_eff = rt * self.pars.base_treat_eff
        self.treat_eff.set(new_treat_eff)
        return self.treat_eff

    def step(self):
        """
        Apply treatment to eligible people who accept it, and update their treatment response.
        """
        treat_uids = super().step()

        # Change treatment resistance for those unsuccessfully treated
        treat_unsucc = self.outcomes['ng']['unsuccessful']
        if len(treat_unsucc):
            self.rel_treat[treat_unsucc] *= (1 - self.pars.rel_treat_unsucc)
        treat_unneed = self.outcomes['ng']['unnecessary']
        if len(treat_unneed):
            self.rel_treat[treat_unneed] *= (1 - self.pars.rel_treat_unneed)

        return treat_uids

    def record_events(self):
        super().record_events()
        self.sim.diseases.ng.results['rel_treat'][self.ti] = np.mean(self.rel_treat)
        return


class UpdateDrugs(ss.Intervention):
    """Intervention that switches the gonorrhea drug regimen when resistance is high.

    Monitors the population-level mean ``rel_treat`` and, when it drops below
    ``threshold_amr``, resets all agents' treatment response to 1.0 (simulating
    a switch to a new drug). The previous ``rel_treat`` values are stored for
    reference.

    Args:
        pars (dict): Parameter overrides; key parameter is ``threshold_amr``
            (default 0.05).
        eligibility (func): Function returning UIDs to monitor.
        **kwargs: Additional parameter overrides.
    """
    def __init__(self, pars=None, eligibility=None, *args, **kwargs):
        super().__init__(*args, who=eligibility)
        self.requires = ['ng', 'gonorrheatreatment']
        self.define_pars(
            threshold_amr=0.05
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.FloatArr('rel_treat_prev'),  # Store a copy of AMR to the previous regimen
        )
        self.change_time = None

    def init_post(self):
        super().init_post()
        self.results += [
        ]

    def step(self):
        tx = self.sim.interventions.gonorrheatreatment
        target_uids = self.check_eligibility()
        pop_rel_treat = np.mean(tx.rel_treat[target_uids])
        if pop_rel_treat < self.pars.threshold_amr:
            self.rel_treat_prev[:] = tx.rel_treat[:]
            tx.rel_treat[:] = 1  # Reset
            self.change_time = self.t.now('year')
        return

