"""
Analyzers for PFA comparison: partners-in-last-year and pair-age heatmap.

These live in devtests, not in stisim.analyzers, because they're for this
benchmark only. They piggy-back on the edge history kept by MFNetwork.
"""
import starsim as ss


def _steps_per_year(sim):
    """ Steps per year for the sim's timeline (fallback 12 for monthly). """
    try:
        return int(round(1.0 / sim.t.dt_year))
    except Exception:
        return 12


class PartnersLastYearAnalyzer(ss.Analyzer):
    """At sim end, count unique partners per agent over the last year, by sex and rel_type.

    Output stored on ``self.records``: list of dicts with keys
    {'uid', 'age', 'sex', 'rel_type', 'n_partners'} where rel_type is one of
    'stable', 'casual', 'onetime', 'sw'.
    """

    def __init__(self):
        super().__init__()
        self.records = []
        return

    def step(self):
        pass

    def finalize_results(self):
        super().finalize_results()
        sim = self.sim
        ppl = sim.people
        nets = [n for n in sim.networks() if hasattr(n, 'edge_types')]
        if not nets:
            return
        end_ti = sim.t.ti
        one_year = _steps_per_year(sim)
        cutoff_ti = end_ti - one_year
        per_agent = {}  # uid -> {rel_type: set(partner_uid)}
        for net in nets:
            rel_durs = getattr(net, 'relationship_durs', {})
            inv_types = {v: k for k, v in net.edge_types.items()}
            for (a, b), events in rel_durs.items():
                for ev in events:
                    start = ev['start']
                    end = start + ev['dur']
                    if end < cutoff_ti:
                        continue
                    rt = inv_types.get(int(ev['edge_type']), 'unknown')
                    per_agent.setdefault(int(a), {}).setdefault(rt, set()).add(int(b))
                    per_agent.setdefault(int(b), {}).setdefault(rt, set()).add(int(a))
        ages = ppl.age.values
        is_female = ppl.female.values
        n_alive = len(ages)
        for uid, by_type in per_agent.items():
            if uid >= n_alive:  # agent died and was removed; skip
                continue
            for rt, partners in by_type.items():
                self.records.append(dict(
                    uid=int(uid),
                    age=float(ages[uid]),
                    sex='F' if is_female[uid] else 'M',
                    rel_type=rt,
                    n_partners=len(partners),
                ))
        return


class PairAgeHeatmapAnalyzer(ss.Analyzer):
    """At sim end, collect (age_p1, age_p2, rel_type) tuples from edges formed in the last year.

    Output stored on ``self.records``: list of dicts.
    """

    def __init__(self):
        super().__init__()
        self.records = []
        return

    def step(self):
        pass

    def finalize_results(self):
        super().finalize_results()
        sim = self.sim
        end_ti = sim.t.ti
        one_year = _steps_per_year(sim)
        cutoff_ti = end_ti - one_year
        nets = [n for n in sim.networks() if hasattr(n, 'edge_types')]
        ages = sim.people.age.values
        n_alive = len(ages)
        for net in nets:
            inv_types = {v: k for k, v in net.edge_types.items()}
            rel_durs = getattr(net, 'relationship_durs', {})
            for (a, b), events in rel_durs.items():
                if a >= n_alive or b >= n_alive:  # one party died; skip
                    continue
                for ev in events:
                    if ev['start'] < cutoff_ti:
                        continue
                    rt = inv_types.get(int(ev['edge_type']), 'unknown')
                    self.records.append(dict(
                        age_p1=float(ages[int(a)]),
                        age_p2=float(ages[int(b)]),
                        rel_type=rt,
                    ))
        return
