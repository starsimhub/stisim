"""
Pregnancy behaviour modifier.

Pregnant agents in many sexual-health models temporarily reduce sexual
risk-taking (drop out of sex work, exit high-risk concurrency patterns,
have fewer partners). ``PregnancyRiskReduction`` switches those network
attributes during pregnancy and restores them as soon as the pregnancy
ends.

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

    As soon as the pregnancy ends (``~pregnant`` for an agent that has been
    pregnant at any prior step), restore each agent's previous FSW status,
    risk group, and concurrency from per-agent high-water marks. No
    dependency on a breastfeeding state.

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
        - The intervention takes per-agent high-water marks
          (``was_fsw``, ``was_high_risk``, ``default_concurrency``,
          ``ever_pregnant``) so that agents who first entered a high-risk
          group or first became pregnant AFTER the intervention started
          are still correctly tracked and restored.
        - FSW status is toggled by mutating SWNetwork's reversible
          ``paused`` BoolArr (read by the ``fsw`` property). The
          lifetime ``ever_fsw`` flag and the SW window bounds
          (``age_sw_start``, ``dur_sw``) are NEVER modified by this
          intervention.
        - Restoration is gated on ``ever_pregnant`` so the intervention is
          a no-op for agents who have never been pregnant. It does not
          rewrite the network state of unrelated agents.
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
        # Renamed off the network's own ``ever_fsw`` to avoid collisions.
        self.define_states(
            ss.BoolState('was_fsw'),
            ss.BoolState('was_high_risk'),
            ss.BoolState('ever_pregnant'),
            ss.FloatArr('default_concurrency', default=0.0),
        )
        return

    def _network(self):
        return self.sim.networks[self.pars.network_name]

    def init_post(self):
        super().init_post()
        nw = self._network()
        # Tracks who is "ever" FSW in the network's lifetime sense. ``ever_fsw``
        # is the writable BoolArr on SWNetwork; for an MFNetwork without sex
        # work (no ever_fsw / paused / fsw at all), fall back to ``fsw`` if
        # the network exposes it, else mark everyone as never-FSW so the
        # intervention is a no-op on that axis.
        if hasattr(nw, 'ever_fsw'):
            self.was_fsw[:] = nw.ever_fsw
        elif hasattr(nw, 'fsw'):
            self.was_fsw[:] = nw.fsw
        self.was_high_risk[:] = (nw.risk_group == self.pars.high_risk_group)
        self.default_concurrency[:] = nw.concurrency
        return

    def step(self):
        nw = self._network()
        preg = self.sim.demographics.pregnancy

        # FSW status is reversibly toggled via SWNetwork's ``paused`` flag,
        # which the ``fsw`` property reads. ``ever_fsw`` and the SW window
        # bounds (age_sw_start, dur_sw) are never mutated by this
        # intervention. risk_group and concurrency are real arrays on
        # MFNetwork and can be written directly.
        has_pause = hasattr(nw, 'paused')

        # High-water-mark updates: track anyone who has ever been FSW,
        # high-risk, or pregnant, so post-pregnancy restoration only fires
        # for agents this intervention has actually modified.
        if hasattr(nw, 'ever_fsw'):
            self.was_fsw[:] = self.was_fsw[:] | nw.ever_fsw
        self.was_high_risk[:] = self.was_high_risk[:] | (nw.risk_group == self.pars.high_risk_group)
        self.default_concurrency[:] = np.maximum(nw.concurrency, self.default_concurrency)
        self.ever_pregnant[:] = self.ever_pregnant[:] | preg.pregnant

        # End-of-pregnancy restore: anyone who has been pregnant at some
        # prior step and is no longer pregnant gets ``paused`` cleared and
        # risk_group / concurrency restored. Clearing ``paused`` puts the
        # agent back into the derived-``fsw`` set provided she is still
        # inside her SW window; if she has aged out during pregnancy,
        # ``nw.fsw`` correctly stays False.
        is_post = ~preg.pregnant & self.ever_pregnant
        if has_pause:
            nw.paused[is_post] = False
        nw.risk_group[self.was_high_risk & is_post] = self.pars.high_risk_group
        nw.concurrency[is_post] = self.default_concurrency[is_post]

        # Pregnancy-time reductions, sampled per-agent.
        is_preg = preg.pregnant
        preg_fsw = (is_preg & nw.fsw).uids
        preg_hr  = (is_preg & (nw.risk_group == self.pars.high_risk_group)).uids
        fsw_quitters  = self.pars.fsw_redux.filter(preg_fsw)
        hr_quitters   = self.pars.high_risk_redux.filter(preg_hr)
        conc_reducers = self.pars.concurrency_redux.filter(is_preg)

        if has_pause:
            nw.paused[fsw_quitters] = True
        nw.risk_group[hr_quitters] = self.pars.default_risk_group
        nw.concurrency[conc_reducers] = 0
        return
