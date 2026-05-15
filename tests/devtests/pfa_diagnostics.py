"""
Analyzers for PFA comparison.

These live in devtests, not in stisim.analyzers, because they're for this benchmark
only. They piggy-back on the edge history kept by MFNetwork.

- PartnersLastYearAnalyzer: count of unique partners per agent over the last year,
  by sex × rel_type. Used for the concurrency-rate figure.
- PairFormationAgesAnalyzer: (age_male, age_female) at every pair formation over the
  full sim. Captures the algorithmic age-mixing imposed by the PFA before any
  differential survival/dissolution. Reads ``edges.age_p1`` / ``edges.age_p2`` which
  the network stores at formation (StructuredSexual convention: p1=male, p2=female).
- PairPrevalenceAnalyzer: snapshot of currently-active edges at sim end with the
  *current* ages of partnered agents. Captures the age-mixing of surviving
  long-lived partnerships.
"""
import starsim as ss


def _steps_per_year(sim):
    """Steps per year for the sim's timeline (fallback 12 for monthly)."""
    try:
        return int(round(1.0 / sim.t.dt_year))
    except Exception:
        return 12


class PartnersLastYearAnalyzer(ss.Analyzer):
    """At sim end, count unique partners per agent over the last year, by sex and rel_type.

    Records: list of dicts with keys ``{uid, age, sex, rel_type, n_partners}``.
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
        per_agent = {}
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
            if uid >= n_alive:
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


class PairFormationAgesAnalyzer(ss.Analyzer):
    """Records (age_male, age_female, rel_type) at every pair formation throughout the sim.

    Hooks ``step()`` and reads each MFNetwork-derived network's current edges. New
    edges are identified by fingerprint (p1, p2, edge_type, age_p1, age_p2) so the
    same pair re-forming after dissolution produces a new record. Uses the edge's
    stored ``age_p1`` / ``age_p2`` (set at formation in ``add_pairs``), which by
    convention are male / female ages respectively.

    Records: list of dicts with keys ``{age_male, age_female, rel_type}``.
    """

    def __init__(self):
        super().__init__()
        self.records = []
        self._seen = set()
        return

    def step(self):
        for net in self.sim.networks():
            if not hasattr(net, 'edge_types'):
                continue
            inv_types = {v: k for k, v in net.edge_types.items()}
            edges = net.edges
            # If edges aren't initialised or empty, skip
            if not hasattr(edges, 'p1') or len(edges.p1) == 0:
                continue
            for i in range(len(edges.p1)):
                key = (
                    int(edges.p1[i]),
                    int(edges.p2[i]),
                    int(edges.edge_type[i]),
                    float(edges.age_p1[i]),
                    float(edges.age_p2[i]),
                )
                if key in self._seen:
                    continue
                self._seen.add(key)
                self.records.append(dict(
                    age_male=float(edges.age_p1[i]),
                    age_female=float(edges.age_p2[i]),
                    rel_type=inv_types.get(int(edges.edge_type[i]), 'unknown'),
                ))
        return


class PairPrevalenceAnalyzer(ss.Analyzer):
    """At sim end, snapshot of currently-active edges with *current* partner ages.

    Records: list of dicts with keys ``{age_male, age_female, rel_type}``. Uses
    ``ppl.age[p1]`` / ``ppl.age[p2]`` for current ages; the StructuredSexual
    convention p1=male, p2=female still applies.
    """

    def __init__(self):
        super().__init__()
        self.records = []
        return

    def step(self):
        pass

    def finalize_results(self):
        super().finalize_results()
        ppl = self.sim.people
        ages = ppl.age.values
        n_alive = len(ages)
        for net in self.sim.networks():
            if not hasattr(net, 'edge_types'):
                continue
            inv_types = {v: k for k, v in net.edge_types.items()}
            edges = net.edges
            if not hasattr(edges, 'p1') or len(edges.p1) == 0:
                continue
            for i in range(len(edges.p1)):
                p1 = int(edges.p1[i])
                p2 = int(edges.p2[i])
                if p1 >= n_alive or p2 >= n_alive:
                    continue
                self.records.append(dict(
                    age_male=float(ages[p1]),
                    age_female=float(ages[p2]),
                    rel_type=inv_types.get(int(edges.edge_type[i]), 'unknown'),
                ))
        return
