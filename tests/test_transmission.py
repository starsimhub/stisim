"""
Test HIV transmission dynamics
"""
import numpy as np
import starsim as ss
import stisim as sti


class track_rel_trans(ss.Analyzer):
    """Track rel_trans statistics by HIV stage on each timestep."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stages = ['acute', 'latent', 'falling']
        self.stats = ['min', 'max', 'mean', 'median']
        self.data = {stage: {stat: [] for stat in self.stats} for stage in self.stages}
        self.counts = {stage: [] for stage in self.stages}
        self.tvec = []

    def step(self):
        hiv = self.sim.diseases.hiv
        self.tvec.append(self.sim.ti)
        for stage in self.stages:
            uids = getattr(hiv, stage).uids
            self.counts[stage].append(len(uids))
            if len(uids) > 0:
                vals = np.asarray(hiv.rel_trans[uids])
                self.data[stage]['min'].append(float(np.min(vals)))
                self.data[stage]['max'].append(float(np.max(vals)))
                self.data[stage]['mean'].append(float(np.mean(vals)))
                self.data[stage]['median'].append(float(np.median(vals)))
            else:
                for stat in self.stats:
                    self.data[stage][stat].append(np.nan)

    def report(self):
        """Print a summary table of rel_trans by stage."""
        for stage in self.stages:
            n = np.nanmax(self.counts[stage]) if self.counts[stage] else 0
            if n > 0:
                means = [v for v in self.data[stage]['mean'] if not np.isnan(v)]
                mins = [v for v in self.data[stage]['min'] if not np.isnan(v)]
                maxs = [v for v in self.data[stage]['max'] if not np.isnan(v)]
                print(f"\n{stage.upper()} (max n={int(n)}):")
                print(f"  mean of means: {np.mean(means):.3f}")
                print(f"  overall min:   {np.min(mins):.3f}")
                print(f"  overall max:   {np.max(maxs):.3f}")
            else:
                print(f"\n{stage.upper()}: no agents in this stage")


def test_rel_trans():
    """Check rel_trans values across HIV stages at monthly dt."""
    tracker = track_rel_trans()
    sim = sti.Sim(
        diseases=sti.HIV(init_prev=0.15),
        networks=sti.StructuredSexual(),
        n_agents=2000,
        start=1990,
        dur=10,
        analyzers=[tracker],
    )
    sim.run()

    tracker = sim.analyzers[0]
    tracker.report()

    acute_means = [v for v in tracker.data['acute']['mean'] if not np.isnan(v)]
    latent_means = [v for v in tracker.data['latent']['mean'] if not np.isnan(v)]
    falling_means = [v for v in tracker.data['falling']['mean'] if not np.isnan(v)]

    assert len(acute_means) > 0, 'Expected some agents in acute stage'
    assert len(latent_means) > 0, 'Expected some agents in latent stage'

    if acute_means:
        print(f"\nAcute mean rel_trans:   {np.mean(acute_means):.2f} (expected ~6)")
    if latent_means:
        print(f"Latent mean rel_trans:  {np.mean(latent_means):.2f} (expected ~1)")
    if falling_means:
        print(f"Falling mean rel_trans: {np.mean(falling_means):.2f} (expected ~8)")

    return tracker


if __name__ == '__main__':
    tracker = test_rel_trans()
    print('\nDone.')
