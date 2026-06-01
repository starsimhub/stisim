"""
Pregnancy behaviour modifier.

Pregnant agents in many sexual-health models temporarily reduce sexual
risk-taking (drop out of sex work, exit high-risk concurrency patterns,
have fewer partners). ``PregnancyRiskReduction`` switches those network
attributes during pregnancy and restores them postpartum (defined as
breastfeeding & no longer pregnant).

Ported from the syph_dx_zim study (``interventions.pregnancy_risk_reduction``)
into the core stisim package so other localizations can reuse it.
"""

import numpy as np
import starsim as ss

__all__ = ['PregnancyRiskReduction']


class PregnancyRiskReduction(ss.Intervention):
    """Pregnancy-driven sexual-behaviour modifier.

    During pregnancy, optionally:

    - Drop FSW status (``fsw_redux``)
    - Drop ``high_risk_group`` membership down to ``default_risk_group``
      (``high_risk_redux``)
    - Drop concurrency to 0 (``concurrency_redux``)

    Postpartum (``breastfeeding & ~pregnant``), restore each agent's
    previous FSW status, risk group, and concurrency (using the per-agent
    high-water marks recorded during pregnancy).

    Designed to be layered on top of ``sti.StructuredSexual`` (or any
    network that exposes ``fsw``, ``risk_group``, ``concurrency``) and an
    ``ss.Pregnancy`` demographics module. Requires both to be present in
    the sim.

    Args:
        fsw_redux (ss.bernoulli): probability per pregnant FSW that she
            exits FSW during pregnancy (default 1.0 = all of them).
        high_risk_redux (ss.bernoulli): probability per pregnant agent in
            ``high_risk_group`` that her risk group is reduced to
            ``default_risk_group`` (default 1.0).
        concurrency_redux (ss.bernoulli): probability per pregnant agent
            that her concurrency drops to 0 (default 1.0).
        network_name (str): name of the network on the sim whose
            ``fsw`` / ``risk_group`` / ``concurrency`` attributes are
            modified (default ``'structuredsexual'``).
        high_risk_group (int): integer risk-group label that is "reduced"
            during pregnancy (default 2 — matches the 3-tier 0/1/2 scheme
            in ``StructuredSexual``).
        default_risk_group (int): risk group to which reduced agents are
            moved (default 0).

    Example:

        sim = sti.Sim(
            diseases='hiv',
            networks=[sti.StructuredSexual()],
            demographics=[ss.Pregnancy(), ss.Deaths()],
            interventions=[sti.PregnancyRiskReduction()],
        )

    Notes:
        - Recovery of pre-pregnancy behaviour happens at ``breastfeeding``,
          NOT at birth. If you don't include a ``BreastfeedingNet`` /
          breastfeeding state, restoration never fires; ensure ``Pregnancy``
          publishes a ``breastfeeding`` boolean.
        - The intervention takes high-water marks (``ever_fsw``,
          ``ever_high_risk``, ``default_concurrency``) so that even agents
          who joined a high-risk group AFTER the intervention started are
          correctly restored postpartum.
    """

    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.define_pars(
            fsw_redux=ss.bernoulli(p=1.0),
            high_risk_redux=ss.bernoulli(p=1.0),
            concurrency_redux=ss.bernoulli(p=1.0),
            network_name='structuredsexual',
            high_risk_group=2,
            default_risk_group=0,
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('ever_fsw'),
            ss.BoolState('ever_high_risk'),
            ss.FloatArr('default_concurrency', default=0.0),
        )
        return

    def _network(self):
        return self.sim.networks[self.pars.network_name]

    def init_post(self):
        super().init_post()
        nw = self._network()
        self.ever_fsw[:] = nw.fsw
        self.ever_high_risk[:] = (nw.risk_group == self.pars.high_risk_group)
        self.default_concurrency[:] = nw.concurrency
        return

    def step(self):
        nw = self._network()
        preg = self.sim.demographics.pregnancy

        # High-water-mark updates: anyone who has been FSW / high-risk at any
        # point gets restored to that status postpartum.
        self.ever_fsw[:] = self.ever_fsw[:] | nw.fsw
        self.ever_high_risk[:] = self.ever_high_risk[:] | (nw.risk_group == self.pars.high_risk_group)
        self.default_concurrency[:] = np.maximum(nw.concurrency, self.default_concurrency)

        # Postpartum restore (breastfeeding & no longer pregnant).
        is_postpartum = preg.breastfeeding & ~preg.pregnant
        nw.fsw[self.ever_fsw & is_postpartum] = True
        nw.risk_group[self.ever_high_risk & is_postpartum] = self.pars.high_risk_group
        nw.concurrency[is_postpartum] = self.default_concurrency[is_postpartum]

        # Pregnancy-time reductions, sampled per-agent.
        is_preg = preg.pregnant
        preg_fsw = (is_preg & nw.fsw).uids
        preg_hr  = (is_preg & (nw.risk_group == self.pars.high_risk_group)).uids
        fsw_quitters  = self.pars.fsw_redux.filter(preg_fsw)
        hr_quitters   = self.pars.high_risk_redux.filter(preg_hr)
        conc_reducers = self.pars.concurrency_redux.filter(is_preg)

        nw.fsw[fsw_quitters] = False
        nw.risk_group[hr_quitters] = self.pars.default_risk_group
        nw.concurrency[conc_reducers] = 0
        return
