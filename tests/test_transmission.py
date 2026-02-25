"""
Test HIV transmission dynamics (rel_trans by stage)
"""
import numpy as np
import starsim as ss
import stisim as sti


class track_rel_trans(ss.Analyzer):
    """Track rel_trans values by HIV stage on each timestep."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.acute_vals = []   # rel_trans for prior-timestep acute agents
        self.latent_vals = []  # rel_trans for latent agents
        self.aids_vals = []    # rel_trans for AIDS agents (cd4 < 200)

    def step(self):
        hiv = self.sim.diseases.hiv
        ti = self.sim.ti
        if ti == 0:
            return

        # Acute agents infected on the previous timestep
        acute_uids = hiv.acute.uids
        if len(acute_uids):
            prev = acute_uids[np.asarray(hiv.ti_infected[acute_uids]) == ti - 1]
            if len(prev):
                self.acute_vals.extend(np.asarray(hiv.rel_trans[prev]).tolist())

        # Latent agents
        latent_uids = hiv.latent.uids
        if len(latent_uids):
            self.latent_vals.extend(np.asarray(hiv.rel_trans[latent_uids]).tolist())

        # AIDS agents (cd4 < 200)
        aids_uids = hiv.aids.uids
        if len(aids_uids):
            self.aids_vals.extend(np.asarray(hiv.rel_trans[aids_uids]).tolist())


def test_rel_trans():
    """Check that rel_trans is correct for acute, latent, and AIDS stages."""
    tracker = track_rel_trans()
    sim = sti.Sim(
        diseases=sti.HIV(init_prev=0.15),
        n_agents=2000,
        start=1990,
        dur=10,
        analyzers=[tracker],
    )
    sim.run()

    tracker = sim.analyzers[0]

    # Acute: should be ~6x for agents infected on the prior timestep
    assert len(tracker.acute_vals) > 0, 'Expected some prior-timestep acute agents'
    acute_mean = np.mean(tracker.acute_vals)
    acute_min = np.min(tracker.acute_vals)
    print(f'Acute (prior-ts) rel_trans: mean={acute_mean:.2f}, min={acute_min:.2f} (n={len(tracker.acute_vals)})')
    assert acute_mean > 5, f'Acute mean rel_trans {acute_mean:.2f} should be >5 (expected ~6)'
    assert acute_min > 1, f'Acute min rel_trans {acute_min:.2f} should be >1'

    # Latent: should be ~1
    assert len(tracker.latent_vals) > 0, 'Expected some latent agents'
    latent_mean = np.mean(tracker.latent_vals)
    print(f'Latent rel_trans:           mean={latent_mean:.2f} (n={len(tracker.latent_vals)})')
    assert np.isclose(latent_mean, 1.0, atol=0.01), f'Latent mean rel_trans {latent_mean:.2f} should be ~1'

    # AIDS: should be ~8x
    assert len(tracker.aids_vals) > 0, 'Expected some AIDS agents'
    aids_mean = np.mean(tracker.aids_vals)
    aids_min = np.min(tracker.aids_vals)
    print(f'AIDS rel_trans:             mean={aids_mean:.2f}, min={aids_min:.2f} (n={len(tracker.aids_vals)})')
    assert aids_mean > 6, f'AIDS mean rel_trans {aids_mean:.2f} should be >6 (expected ~8)'
    assert aids_min > 1, f'AIDS min rel_trans {aids_min:.2f} should be >1'

    return tracker


if __name__ == '__main__':
    tracker = test_rel_trans()
    print('\nDone.')
