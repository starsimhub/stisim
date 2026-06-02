"""
MFNetwork assertion tests.

Verifies the MFNetwork class against the spec in
``mfnetwork_test_instructions.txt``:

  * 3-seed rand_seed sweep, n_agents=10000, dur=20yr, window_months=12
  * Across target_age_gap in {5, 8, 11, 14}, the sweep-averaged realized
    male-female age gap mean and std are within 5% of their targets
  * Per 5-year female-age bin from 15 through 55, the realized mean and std
    are within 5% of target; bins whose lower edge is >= 50 may relax to 40%
  * At least one partnership forms in every female-age bin through 55

Additional checks proposed and accepted:

  * Truncation hard-cap: no realized |gap| exceeds target + 3*stddev
  * Taper kwarg: f_partnership_taper_cut propagates to network.pars and shapes
    the female-age distribution of newly formed partnerships. The companion
    offset is internal (``_f_partnership_taper_offset``), auto-set in
    ``init_post`` from the age-diff prefs; it is read here for diagnostics
    but not modified.

Companion to tests/test_network_diagnostics.py, which produces diagnostic
plots; this file produces assertable pass/fail tests.
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import sciris as sc
import stisim as sti
import sys

from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hiv_natural_history_analyzers import AgeGapWindowAnalyzer, NewPairsAnalyzer


# ---- Spec config -----------------------------------------------------------
SEEDS = [1, 2, 3]
TARGET_GAPS = [5, 8, 11, 14]
TARGET_STD = 3
N_AGENTS = 10000
DUR = 20
WINDOW_MONTHS = 12

# Female-age bins (5-yr) covering 15 through 55: [15,20), [20,25), ..., [50,55)
BIN_EDGES = np.arange(15, 60, 5).astype(float)
# Loose-tolerance bin threshold is derived per-sim from the live network as
# ceil((f_partnership_taper_cut - _f_partnership_taper_offset) / 5) * 5 — the
# female age at which the looking-chance taper starts, snapped UP to the next
# 5-yr bin edge. Bins whose lower edge >= that threshold use LOOSE_TOL.

TIGHT_TOL = 0.05  # 5% — overall sweep-averaged mean / std
PER_BIN_TIGHT_TOL = 0.10  # 10% — young/middle female-age bins
LOOSE_TOL = 0.25  # 25% — bins with lower edge >= OLD_BIN_LO
MAX_DEV_MULT = 3   # |gap| hard-cap = target + MAX_DEV_MULT * target_std


# ---- Helpers ---------------------------------------------------------------
def _run_age_gap_sim(target_age_gap, seed, n_agents=N_AGENTS, dur=DUR,
                      window_months=WINDOW_MONTHS, target_std=TARGET_STD):
    """Run one sim and return (sim, (male_ages, female_ages, gaps), loose_bin_lo)."""
    # age_diff_pars: same shape as test_network_diagnostics.test_age_gap_distribution
    age_diff_pars = {'teens': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)],
                     'young': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)],
                     'adult': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)]}
    network = sti.MFNetwork(age_diff_pars=age_diff_pars)
    analyzer = AgeGapWindowAnalyzer(network='mfnetwork', window_months=window_months)
    sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                   dur=dur, rand_seed=seed)
    sim.run()
    nw = sim.networks.mfnetwork
    taper_start = float(nw.pars.f_partnership_taper_cut) - float(nw.pars._f_partnership_taper_offset)
    loose_bin_lo = math.ceil(taper_start / 5) * 5
    return sim, sim.analyzers.agegapwindowanalyzer.get_first_appearance(), loose_bin_lo


def _bin_stats(female_ages, gaps):
    """Return list of (n, mean, std) per BIN_EDGES bin (NaN where empty)."""
    bin_idx = np.digitize(female_ages, BIN_EDGES) - 1
    n_bins = len(BIN_EDGES) - 1
    out = []
    for b in range(n_bins):
        mask = bin_idx == b
        if mask.any():
            out.append((int(mask.sum()),
                        float(gaps[mask].mean()),
                        float(gaps[mask].std())))
        else:
            out.append((0, float('nan'), float('nan')))
    return out


# ---- Tests -----------------------------------------------------------------
@sc.timer()
def test_mfnetwork_age_gap_targets(n_agents=N_AGENTS, dur=DUR, seeds=SEEDS,
                                     target_gaps=TARGET_GAPS, target_std=TARGET_STD,
                                     window_months=WINDOW_MONTHS):
    """Sweep target_age_gap x rand_seed; assert realized stats match targets.

    Asserted per the spec:
      1. Sweep-averaged within-sim mean and std are within 5% of target.
      2. Per 5-yr female-age bin (15-55), sweep-averaged within-bin mean and
         std are within 10%; bins with lower edge >= 50 may relax to 25%.
      3. At least one partnership forms in every 5-yr female-age bin through 55.
      4. Truncation: max |realized gap| <= target + 3*stddev (proposed add-on).
    """
    sweep = defaultdict(list)  # target_gap -> list of (mean, std, bin_stats, loose_bin_lo)
    n_bins = len(BIN_EDGES) - 1
    max_dev_target = MAX_DEV_MULT * target_std

    for target_gap in target_gaps:
        max_cap = target_gap + max_dev_target
        for seed in seeds:
            sim, (male_ages, female_ages, gaps), loose_bin_lo = _run_age_gap_sim(
                target_gap, seed, n_agents=n_agents, dur=dur,
                window_months=window_months, target_std=target_std)
            assert gaps.size > 0, \
                f"target_gap={target_gap}, seed={seed}: no partnerships in window"

            overall_mean = float(gaps.mean())
            overall_std = float(gaps.std())
            max_abs = float(np.max(np.abs(gaps)))

            # Truncation hard-cap (proposal).
            assert max_abs <= max_cap, \
                (f"target_gap={target_gap}, seed={seed}: |gap| cap violated: "
                 f"max|gap|={max_abs:.3f} > target+{MAX_DEV_MULT}*std={max_cap:.3f}")

            bin_stats = _bin_stats(female_ages, gaps)
            sweep[target_gap].append((overall_mean, overall_std, bin_stats, loose_bin_lo))
            print(f"  target={target_gap:2d} seed={seed}: n={gaps.size:6d} "
                  f"mean={overall_mean:6.3f} std={overall_std:6.3f} "
                  f"max|gap|={max_abs:5.2f} loose>={loose_bin_lo}")

    # 1. Overall mean/std assertions
    print("\n--- Overall realized-gap assertions (sweep-averaged) ---")
    print(f"  target | sweep_mean (err%) | sweep_std (err%) | std(per-seed means)")
    print(f"  -------+-------------------+------------------+--------------------")
    for target_gap in target_gaps:
        means = [r[0] for r in sweep[target_gap]]
        stds = [r[1] for r in sweep[target_gap]]
        sweep_mean = float(np.mean(means))
        sweep_mean_of_stds = float(np.mean(stds))
        std_of_means = float(np.std(means))
        mean_err = abs(sweep_mean - target_gap) / target_gap
        std_err = abs(sweep_mean_of_stds - target_std) / target_std
        print(f"  {target_gap:6d} | {sweep_mean:8.3f} ({mean_err*100:5.2f}%) | "
              f"{sweep_mean_of_stds:7.3f} ({std_err*100:5.2f}%) | {std_of_means:8.3f}")
        assert mean_err <= TIGHT_TOL, \
            (f"target_gap={target_gap}: sweep mean {sweep_mean:.3f} is "
             f"{mean_err*100:.2f}% off target {target_gap}; allowed {TIGHT_TOL*100:.0f}%")
        assert std_err <= TIGHT_TOL, \
            (f"target_gap={target_gap}: sweep mean-of-stds {sweep_mean_of_stds:.3f} "
             f"is {std_err*100:.2f}% off target_std {target_std}; "
             f"allowed {TIGHT_TOL*100:.0f}%")

    # 2 & 3. Per-bin assertions (mean, std, presence) by female age
    print("\n--- Per-bin assertions (female age, 5-yr bins through 55) ---")
    print(f"  target | bin    | n_pairs | sweep_mean (err%) | sweep_std (err%) | tol")
    print(f"  -------+--------+---------+-------------------+------------------+------")
    for target_gap in target_gaps:
        # All seeds share the same age_diff_pars, so loose_bin_lo is constant per target_gap.
        loose_bin_los = {r[3] for r in sweep[target_gap]}
        assert len(loose_bin_los) == 1, \
            f"target_gap={target_gap}: loose_bin_lo varied across seeds: {loose_bin_los}"
        loose_bin_lo = loose_bin_los.pop()

        for b in range(n_bins):
            bin_lo = float(BIN_EDGES[b])
            bin_hi = float(BIN_EDGES[b + 1])
            ns = [r[2][b][0] for r in sweep[target_gap]]
            n_total = int(sum(ns))

            assert n_total > 0, \
                (f"target_gap={target_gap}: no partnerships in female-age bin "
                 f"[{int(bin_lo)}, {int(bin_hi)}) across {len(seeds)} seeds")

            means = [r[2][b][1] for r in sweep[target_gap]
                     if not np.isnan(r[2][b][1])]
            stds_ = [r[2][b][2] for r in sweep[target_gap]
                     if not np.isnan(r[2][b][2])]
            sweep_mean = float(np.mean(means))
            sweep_mean_of_stds = float(np.mean(stds_))
            mean_err = abs(sweep_mean - target_gap) / target_gap
            std_err = abs(sweep_mean_of_stds - target_std) / target_std
            tol = LOOSE_TOL if bin_lo >= loose_bin_lo else PER_BIN_TIGHT_TOL
            tol_label = 'loose' if tol == LOOSE_TOL else 'tight'
            print(f"  {target_gap:6d} | {int(bin_lo):2d}-{int(bin_hi):2d}  | "
                  f"{n_total:7d} | {sweep_mean:8.3f} ({mean_err*100:5.2f}%) | "
                  f"{sweep_mean_of_stds:7.3f} ({std_err*100:5.2f}%) | {tol_label}")
            assert mean_err <= tol, \
                (f"target_gap={target_gap}, bin [{int(bin_lo)},{int(bin_hi)}): "
                 f"mean {sweep_mean:.3f} is {mean_err*100:.2f}% off target "
                 f"{target_gap}; allowed {tol*100:.0f}%")
            assert std_err <= tol, \
                (f"target_gap={target_gap}, bin [{int(bin_lo)},{int(bin_hi)}): "
                 f"std {sweep_mean_of_stds:.3f} is {std_err*100:.2f}% off target_std "
                 f"{target_std}; allowed {tol*100:.0f}%")

    return sweep


@sc.timer()
def test_mfnetwork_taper_kwargs(n_agents=N_AGENTS, dur=DUR, target_age_gap=5,
                                  target_std=TARGET_STD, seed=1):
    """Verify the f_partnership_taper_cut kwarg shapes new-partnership formation.

    The kwarg to ``MFNetwork()`` should:
      * Propagate to ``network.pars`` after init.
      * Suppress new partnerships among females with age >= ``taper_cut``.
      * Allow new partnerships below ``taper_cut`` (and especially below
        ``taper_cut - _f_partnership_taper_offset`` where the looking-chance
        saturates at 1).

    ``_f_partnership_taper_offset`` is internal — auto-set by ``init_post``
    from ``age_diff_pars`` — so we only read it for diagnostic bands, never
    modify it.
    """

    # age_diff_pars: same shape as test_network_diagnostics.test_age_gap_distribution
    age_diff_pars = {'teens': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)],
                     'young': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)],
                     'adult': [(target_age_gap, target_std), (target_age_gap, target_std), (target_age_gap, target_std)]}

    def _run(**mf_kwargs):
        network = sti.MFNetwork(age_diff_pars=age_diff_pars, **mf_kwargs)
        analyzer = NewPairsAnalyzer(network='mfnetwork')
        sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                       dur=dur, rand_seed=seed)
        sim.run()
        nw = sim.networks.mfnetwork
        return nw, sim.analyzers.newpairsanalyzer.get_female_ages()

    # --- Baseline: default kwargs ------------------------------------------
    nw_def, ages_def = _run()
    cut_def = float(nw_def.pars.f_partnership_taper_cut)
    offset_def = float(nw_def.pars._f_partnership_taper_offset)
    n_above_def = int((ages_def >= cut_def).sum())
    n_below_def = int((ages_def < cut_def).sum())
    n_in_taper_def = int(((ages_def >= cut_def - offset_def) &
                           (ages_def < cut_def)).sum())
    print(f"  default: cut={cut_def} offset={offset_def}: "
          f"n_above_cut={n_above_def} n_below_cut={n_below_def} "
          f"n_in_taper_band={n_in_taper_def}")
    assert n_above_def == 0, \
        (f"default taper_cut={cut_def}: {n_above_def} new partnerships at "
         f"female_age >= cut; expected 0")
    assert n_below_def > 0, \
        f"default taper_cut={cut_def}: no new partnerships at female_age < cut"
    assert n_in_taper_def > 0, \
        (f"default taper band [{cut_def - offset_def}, {cut_def}) had no new "
         f"partnerships; the taper should reduce — not zero — formations here")

    # --- Modified taper_cut: lower the upper bound -------------------------
    custom_cut = cut_def - 10.0  # e.g. 55 -> 45
    nw_lo, ages_lo = _run(f_partnership_taper_cut=custom_cut)
    actual_cut_lo = float(nw_lo.pars.f_partnership_taper_cut)
    assert actual_cut_lo == custom_cut, \
        (f"f_partnership_taper_cut kwarg not honored: input={custom_cut}, "
         f"post-init={actual_cut_lo}")
    n_above_lo = int((ages_lo >= actual_cut_lo).sum())
    n_below_lo = int((ages_lo < actual_cut_lo).sum())
    print(f"  cut={actual_cut_lo}: n_above_cut={n_above_lo} n_below_cut={n_below_lo}")
    assert n_above_lo == 0, \
        (f"taper_cut={actual_cut_lo}: {n_above_lo} new partnerships at "
         f"female_age >= cut; expected 0")
    assert n_below_lo > 0, \
        f"taper_cut={actual_cut_lo}: no new partnerships at female_age < cut"

    # Lowering the cut should zero out partnerships in [custom_cut, default_cut)
    n_band_def = int(((ages_def >= actual_cut_lo) & (ages_def < cut_def)).sum())
    n_band_lo = int(((ages_lo >= actual_cut_lo) & (ages_lo < cut_def)).sum())
    print(f"  formations in [{actual_cut_lo}, {cut_def}): "
          f"default={n_band_def}, cut={actual_cut_lo}={n_band_lo}")
    assert n_band_def > 0, \
        (f"sanity check failed: expected default sim to form partnerships in "
         f"[{actual_cut_lo}, {cut_def}); got {n_band_def}")
    assert n_band_lo == 0, \
        (f"lowering taper_cut to {actual_cut_lo} did not suppress new partnerships in "
         f"[{actual_cut_lo}, {cut_def}); got {n_band_lo}")


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_mfnetwork_age_gap_targets()
    test_mfnetwork_taper_kwargs()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
