"""
Pair-formation algorithm variants for MFNetwork.

Each class is a thin subclass of stisim.MFNetwork that overrides
``match_pairs`` only. Used for benchmarking and comparison; see
``tests/devtests/devtest_pfa_comparison.py``.
"""
import numpy as np
import scipy.optimize as spo
import scipy.spatial as spsp
from .networks import MFNetwork, NoPartnersFound


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
