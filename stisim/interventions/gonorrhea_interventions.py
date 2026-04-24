"""
Define gonorrhea interventions for STIsim
"""

import starsim as ss
import numpy as np
import sciris as sc
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
        max_capacity (int): Maximum treatments per timestep.
        years (list): Calendar years during which the intervention is active.
        name (str): Intervention name.
        **kwargs: Additional parameter overrides.
    """
    def __init__(self, pars=None, eligibility=None, max_capacity=None, years=None, name=None, *args, **kwargs):
        super().__init__(diseases='ng', name=name, eligibility=eligibility, max_capacity=max_capacity, years=years, *args)
        self.requires = ['ng', 'structuredsexual']
        self.define_pars(
            base_treat_eff=0.96,
            treat_eff=ss.bernoulli(p=0),  # Reset each time step depending on base_treat_eff and population AMR
            rel_treat_unsucc=0.01,
            rel_treat_unneed=0.005,
        )
        self.update_pars(pars, **kwargs)

        # States
        self.define_states(
            ss.FloatArr('rel_treat', default=1),  # How well a person will respond to treatment
        )

    def set_treat_eff(self, uids):
        new_treat_eff = self.rel_treat[uids] * self.pars.base_treat_eff
        self.pars.treat_eff.set(new_treat_eff)
        return

    def step(self):
        """
        Apply treatment. On each timestep, this method will add eligible people who are willing to accept treatment to a
        queue, and then will treat as many people in the queue as there is capacity for.
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

    def update_results(self):
        super().update_results()
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
        years (list): Calendar years during which the intervention is active.
        **kwargs: Additional parameter overrides.
    """
    def __init__(self, pars=None, eligibility=None, years=None, *args, **kwargs):
        super().__init__(*args)
        self.requires = ['ng', 'gonorrheatreatment']
        self.define_pars(
            threshold_amr=0.05
        )
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        self.years = years
        self.define_states(
            ss.FloatArr('rel_treat_prev'),  # Store a copy of AMR to the previous regimen
        )
        self.change_time = None

    def init_post(self):
        super().init_post()
        self.results += [
        ]

    def step(self):
        sim = self.sim
        target_uids = self.check_eligibility(sim)
        pop_rel_treat = np.mean(self.sim.people.rel_treat[target_uids])
        if pop_rel_treat < self.pars.threshold_amr:
            self.sim.people.rel_treat_prev = sc.dcp(self.sim.people.rel_treat)
            self.sim.people.rel_treat[:] = 1  # Reset
            self.change_time = self.sim.year
        return

