"""
Tests for the PregnancyRiskReduction intervention (ported from
syph_dx_zim, issue #482).

The intervention should:
1. Set fsw=False for pregnant FSWs (per fsw_redux probability)
2. Set risk_group=default_risk_group for pregnant high-risk agents
   (per high_risk_redux probability)
3. Set concurrency=0 for pregnant agents (per concurrency_redux probability)
4. Restore those attributes postpartum (breastfeeding & ~pregnant), using
   the per-agent high-water mark
"""

import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti


def _make_sim(intervention_pars=None, n_agents=3_000, dur=15, rand_seed=1):
    return sti.Sim(
        diseases=[sti.HIV(init_prev=0.05, beta_m2f=0.05)],
        networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
        demographics=[ss.Pregnancy(), ss.Deaths()],
        interventions=[sti.PregnancyRiskReduction(intervention_pars or {})],
        n_agents=n_agents, dur=dur, start=2005, verbose=-1, rand_seed=rand_seed,
    )


@sc.timer()
def test_pregnant_agents_are_not_fsw():
    """At most a handful of pregnant FSWs should remain (one-step lag for
    agents who newly become pregnant + FSW within the final timestep is
    expected); the bulk should be cleared.
    """
    sim = _make_sim()
    sim.run()
    preg = sim.demographics.pregnancy.pregnant
    fsw  = sim.networks.structuredsexual.fsw
    n_preg_fsw = (preg & fsw).count()
    n_preg = preg.count()
    # Allow a small one-step-lag residual but require the intervention
    # caught the vast majority — <5% of pregnant agents.
    assert n_preg == 0 or n_preg_fsw / n_preg < 0.05, \
        f'{n_preg_fsw}/{n_preg} pregnant agents still FSW — fsw_redux=1.0 should drop them'


@sc.timer()
def test_pregnant_agents_have_zero_concurrency():
    """Pregnant agents should have concurrency 0 when concurrency_redux=1.0."""
    sim = _make_sim()
    sim.run()
    preg = sim.demographics.pregnancy.pregnant
    conc = sim.networks.structuredsexual.concurrency
    if preg.count() == 0:
        return  # no pregnant agents in this run; nothing to assert
    assert np.all(conc[preg.uids] == 0), \
        f'Found pregnant agents with non-zero concurrency: max={conc[preg.uids].max()}'


@sc.timer()
def test_pregnant_agents_not_in_high_risk_group():
    """Pregnant agents should not be in high_risk_group when high_risk_redux=1.0."""
    sim = _make_sim()
    sim.run()
    preg = sim.demographics.pregnancy.pregnant
    rg   = sim.networks.structuredsexual.risk_group
    high = sim.interventions.pregnancyriskreduction.pars.high_risk_group
    assert (preg & (rg == high)).count() == 0, \
        f'Found pregnant agents still in risk group {high}'


@sc.timer()
def test_redux_zero_means_no_change():
    """With all redux probabilities = 0, the intervention should be a no-op."""
    sim_off = _make_sim(intervention_pars=dict(
        fsw_redux=ss.bernoulli(p=0.0),
        high_risk_redux=ss.bernoulli(p=0.0),
        concurrency_redux=ss.bernoulli(p=0.0),
    ))
    sim_off.run()
    preg = sim_off.demographics.pregnancy.pregnant
    nw = sim_off.networks.structuredsexual
    intv = sim_off.interventions.pregnancyriskreduction
    if preg.count() == 0:
        return
    # ever_fsw / ever_high_risk are high-water marks that include "current"
    # status, so a pregnant agent who started high-risk should still be
    # flagged ever_high_risk regardless of redux. The redux-only check is
    # that the LIVE network attributes weren't zeroed by the intervention:
    # if any pregnant agent currently has concurrency==0 AND default_concurrency>0,
    # the intervention erroneously fired.
    live_conc = nw.concurrency[preg.uids]
    cached_conc = intv.default_concurrency[preg.uids]
    erroneously_zeroed = (live_conc == 0) & (cached_conc > 0)
    assert erroneously_zeroed.sum() == 0, \
        f'{erroneously_zeroed.sum()} pregnant agents had concurrency zeroed despite redux=0'


@sc.timer()
def test_postpartum_restoration():
    """Postpartum agents (breastfeeding & not pregnant) should have their
    ever_fsw and ever_high_risk statuses restored by the intervention.
    """
    sim = _make_sim()
    sim.run()
    nw = sim.networks.structuredsexual
    preg = sim.demographics.pregnancy
    intv = sim.interventions.pregnancyriskreduction
    postpartum = preg.breastfeeding & ~preg.pregnant
    if postpartum.count() == 0:
        return  # too few agents/years to observe postpartum window
    # Anyone flagged ever_fsw who is currently postpartum should be FSW now.
    pp_ever_fsw = (postpartum & intv.ever_fsw).uids
    if len(pp_ever_fsw):
        assert nw.fsw[pp_ever_fsw].all(), \
            f'Postpartum agents flagged ever_fsw are not restored to FSW: ' \
            f'{(~nw.fsw[pp_ever_fsw]).sum()}/{len(pp_ever_fsw)}'


# %% Run as a script
if __name__ == '__main__':
    test_pregnant_agents_are_not_fsw()
    test_pregnant_agents_have_zero_concurrency()
    test_pregnant_agents_not_in_high_risk_group()
    test_redux_zero_means_no_change()
    test_postpartum_restoration()
