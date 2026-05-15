"""
Pair-formation algorithm variants for MFNetwork.

Each class is a thin subclass of stisim.MFNetwork that overrides
``match_pairs`` only. Used for benchmarking and comparison; see
``tests/devtests/devtest_pfa_comparison.py``.
"""
import numpy as np
import scipy.optimize as spo
import scipy.spatial as spsp
import starsim as ss
from .networks import MFNetwork, NoPartnersFound

__all__ = [
    'MFNetwork_LSA',
    'MFNetwork_SortPair',
    'MFNetwork_DesiredAgeBucket',
    'MFNetwork_SortBisect',
    'MFNetwork_GreedyOldEnough',
    'MFNetwork_KDTreeNN',
    'MFNetwork_BandMatch',
]


class MFNetwork_LSA(MFNetwork):
    """Linear sum assignment on the full age-distance matrix. O(n^3)."""

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        if len(m_ages) == 0 or len(desired_ages) == 0:
            raise NoPartnersFound()
        dist_mat = spsp.distance_matrix(m_ages[:, np.newaxis], desired_ages[:, np.newaxis])
        ind_m, ind_f = spo.linear_sum_assignment(dist_mat)
        p1 = m_eligible.uids[ind_m]
        p2 = f_looking[ind_f]
        return p1, p2


class MFNetwork_SortPair(MFNetwork):
    """argsort both groups, zip, truncate to min(len). No tail trim."""

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        ind_m = np.argsort(m_ages, stable=True)
        ind_f = np.argsort(desired_ages, stable=True)
        maxlen = min(len(ind_m), len(ind_f))
        if maxlen == 0:
            raise NoPartnersFound()
        ind_m = ind_m[:maxlen]
        ind_f = ind_f[:maxlen]
        p1 = m_eligible.uids[ind_m]
        p2 = f_looking[ind_f]
        return p1, p2


class MFNetwork_DesiredAgeBucket(MFNetwork):
    """Bucket women by integer desired-age; sample men in that bucket (with
    replacement if M < W). Post-filter by male concurrency cap (processing women
    in random order) and by risk-group-conditional relationship acceptance.

    This duplicates logic that currently sits in ``add_pairs``; the duplication
    is accepted for cleanliness of the variant.
    """

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        m_uids = m_eligible.uids
        f_uids = f_looking

        # Integer-bucket desired ages and male ages.
        desired_int = np.floor(desired_ages).astype(int)
        m_int = np.floor(m_ages).astype(int)

        # Group male indices by age bucket.
        buckets = {}
        for i, age in enumerate(m_int):
            buckets.setdefault(int(age), []).append(i)

        # Match each woman to a man in her bucket (with replacement if needed).
        p1_list, p2_list = [], []
        rng = np.random.default_rng(self.sim.pars.rand_seed + self.ti)
        order = rng.permutation(len(f_uids))
        for idx in order:
            bucket = buckets.get(int(desired_int[idx]))
            if not bucket:
                continue
            choice = bucket[rng.integers(0, len(bucket))]
            p1_list.append(int(m_uids[choice]))
            p2_list.append(int(f_uids[idx]))

        if not p1_list:
            raise NoPartnersFound()

        p1 = ss.uids(np.array(p1_list, dtype=np.int64))
        p2 = ss.uids(np.array(p2_list, dtype=np.int64))

        # Post-filter 1: enforce male concurrency cap.
        male_count = {}
        keep = np.ones(len(p1), dtype=bool)
        for i, uid in enumerate(p1):
            male_count[int(uid)] = male_count.get(int(uid), 0) + 1
            if male_count[int(uid)] > int(self.concurrency[uid]):
                keep[i] = False
        p1 = p1[keep]
        p2 = p2[keep]

        if len(p1) == 0:
            raise NoPartnersFound()

        # Post-filter 2: relationship-acceptance Bernoulli.
        # Mirrors logic from MFNetwork.add_pairs.
        matched_risk = (self.risk_group[p1] == self.risk_group[p2])
        mismatched_risk = ~matched_risk
        p_match = np.zeros(len(p1), dtype=float)
        for rg in range(self.pars.n_risk_groups):
            p_match[matched_risk & (self.risk_group[p1] == rg)] = self.pars.p_matched_stable[rg]
            p_match[mismatched_risk & (self.risk_group[p2] == rg)] = self.pars.p_mismatched_casual[rg]
        accept = rng.random(len(p1)) < p_match
        p1 = p1[accept]
        p2 = p2[accept]

        if len(p1) == 0:
            raise NoPartnersFound()

        return p1, p2


class MFNetwork_SortBisect(MFNetwork):
    """Current production: argsort + bisect-trim + subsample. Inherits match_pairs unchanged."""
    pass


class MFNetwork_GreedyOldEnough(MFNetwork):
    """Sort women by desired age ascending; for each, take youngest available
    male with age >= desired_age. No replacement.
    """

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        m_uids = m_eligible.uids
        f_uids = f_looking

        order_f = np.argsort(desired_ages)
        order_m = np.argsort(m_ages)
        sorted_m_ages = m_ages[order_m]
        sorted_m_uids = m_uids[order_m]
        available = np.ones(len(sorted_m_uids), dtype=bool)

        p1_list, p2_list = [], []
        cursor = 0
        for fi in order_f:
            target = desired_ages[fi]
            while cursor < len(sorted_m_ages) and (not available[cursor] or sorted_m_ages[cursor] < target):
                cursor += 1
            if cursor >= len(sorted_m_ages):
                break
            p1_list.append(int(sorted_m_uids[cursor]))
            p2_list.append(int(f_uids[fi]))
            available[cursor] = False
            cursor += 1

        if not p1_list:
            raise NoPartnersFound()

        p1 = ss.uids(np.array(p1_list, dtype=np.int64))
        p2 = ss.uids(np.array(p2_list, dtype=np.int64))
        return p1, p2


class MFNetwork_KDTreeNN(MFNetwork):
    """Build KDTree on male ages; each woman queries k=1; resolve collisions
    by giving each contested man to the closest-by-age woman.
    """

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        m_uids = m_eligible.uids
        f_uids = f_looking

        tree = spsp.KDTree(m_ages[:, np.newaxis])
        dists, idxs = tree.query(desired_ages[:, np.newaxis], k=1)
        dists = dists.ravel()
        idxs = idxs.ravel()

        order = np.argsort(dists)
        taken_m = set()
        keep = np.zeros(len(f_uids), dtype=bool)
        for ranked in order:
            mi = int(idxs[ranked])
            if mi in taken_m:
                continue
            taken_m.add(mi)
            keep[ranked] = True

        if not keep.any():
            raise NoPartnersFound()

        p1 = m_uids[idxs[keep]]
        p2 = f_uids[keep]
        return p1, p2


class MFNetwork_BandMatch(MFNetwork):
    """Bucket both groups into 5-year age bands; shuffle and zip within band."""

    band_width = 5

    def match_pairs(self):
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
        m_ages = ppl.age[m_eligible]
        m_uids = m_eligible.uids
        f_uids = f_looking

        f_band = (desired_ages // self.band_width).astype(int)
        m_band = (m_ages // self.band_width).astype(int)

        rng = np.random.default_rng(self.sim.pars.rand_seed + self.ti)
        p1_list, p2_list = [], []
        bands = np.unique(np.concatenate([f_band, m_band]))
        for b in bands:
            mi = np.where(m_band == b)[0]
            fi = np.where(f_band == b)[0]
            if len(mi) == 0 or len(fi) == 0:
                continue
            rng.shuffle(mi)
            rng.shuffle(fi)
            k = min(len(mi), len(fi))
            p1_list.extend(m_uids[mi[:k]].tolist())
            p2_list.extend(f_uids[fi[:k]].tolist())

        if not p1_list:
            raise NoPartnersFound()

        p1 = ss.uids(np.array(p1_list, dtype=np.int64))
        p2 = ss.uids(np.array(p2_list, dtype=np.int64))
        return p1, p2
