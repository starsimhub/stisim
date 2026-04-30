"""
Quick comparison of HIV transmission patterns under different acute-stage
rel_trans assumptions (HIVsim default ~6× vs Hollingsworth et al. ~26×).

Adds a minimal analyzer that partitions new transmissions each timestep
by the source's HIV stage: acute, chronic (latent/falling with CD4>=200),
or AIDS (CD4<200). At the end, compares the two scenarios on the share
of cumulative transmissions originating from the acute stage.
"""
import numpy as np
import starsim as ss
import hivsim
import matplotlib.pyplot as plt


class StageTransmission(ss.Analyzer):
    """Partition new HIV transmissions by source stage each timestep."""

    def init_pre(self, 
    sim):
        super().init_pre(sim)
        self.define_results(
            ss.Result('n_from_acute',   dtype=int, scale=True),
            ss.Result('n_from_chronic', dtype=int, scale=True),
            ss.Result('n_from_aids',    dtype=int, scale=True),
        )

    def step(self):
        hiv = self.sim.diseases.hiv
        ti = self.ti
        n = hiv.new_transmissions  # per-source count of new infections this ti

        self.results['n_from_acute'][ti]   = n[hiv.acute].sum()
        self.results['n_from_aids'][ti]    = n[hiv.aids].sum()
        chronic = hiv.infected & ~hiv.acute & ~hiv.aids
        self.results['n_from_chronic'][ti] = n[chronic].sum()


def make_sim(rel_trans_acute, label):
    sim = hivsim.demo(
        'zimbabwe', run=False,
        n_agents=5_000, dur=30, rand_seed=1, verbose=0,
        rel_trans_acute=rel_trans_acute,
        dur_acute=ss.year,
        analyzers=[StageTransmission()],
    )
    sim.label = label
    return sim


s_lo = make_sim(6,  'HIVsim default (6×)')
s_hi = make_sim(26, 'Hollingsworth (26×)')
msim = ss.parallel(s_lo, s_hi)

# Aggregate results
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Panel 1: stacked fraction of cumulative transmissions by source stage
ax = axes[0]
labels, acute, chronic, aids = [], [], [], []
for sim in msim.sims:
    r = sim.analyzers.stagetransmission.results
    tot = r.n_from_acute.sum() + r.n_from_chronic.sum() + r.n_from_aids.sum()
    tot = max(tot, 1)
    labels.append(sim.label)
    acute.append(r.n_from_acute.sum() / tot)
    chronic.append(r.n_from_chronic.sum() / tot)
    aids.append(r.n_from_aids.sum() / tot)
acute, chronic, aids = map(np.array, (acute, chronic, aids))
ax.bar(labels, acute,   label='Acute',   color='#d62728')
ax.bar(labels, chronic, label='Chronic', color='#7f7f7f', bottom=acute)
ax.bar(labels, aids,    label='AIDS',    color="#cb8062", bottom=acute + chronic)
ax.set_ylabel('Share of cumulative transmissions')
ax.set_ylim(0, 1)
ax.set_title('Transmissions by source stage')
ax.legend(loc='upper right')

# Panel 2: fraction from acute over time (rolling-summed for readability)
ax = axes[1]
for sim, color in zip(msim.sims, ['#888', '#d62728']):
    r = sim.analyzers.stagetransmission.results
    total = r.n_from_acute + r.n_from_chronic + r.n_from_aids
    # Rolling window to smooth noise
    w = 6
    total_s = np.convolve(total, np.ones(w), mode='same')
    acute_s = np.convolve(r.n_from_acute, np.ones(w), mode='same')
    frac = np.where(total_s > 0, acute_s / np.maximum(total_s, 1), np.nan)
    ax.plot(sim.timevec, frac, label=sim.label, color=color, lw=2)
ax.set_xlabel('Year')
ax.set_ylabel('Fraction of transmissions from acute stage')
ax.set_ylim(0, 1)
ax.grid(alpha=0.3)
ax.legend()

fig.tight_layout()
plt.show()
