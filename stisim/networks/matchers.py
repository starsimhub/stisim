"""Pair-formation algorithms (matchers) used by MFNetwork.match_pairs.

Each matcher is a module-level function with signature::

    def matcher(net) -> (p1, p2)

where ``net`` is an MFNetwork instance (so the matcher can call
``net._get_eligible()`` and ``net._sample_desired_ages()``) and ``p1``,
``p2`` are ``ss.uids`` arrays of equal length giving male/female partners.

Raise :class:`NoPartnersFound` to signal an empty match.

The registry ``MATCHERS`` maps string names (used as the ``match_method``
parameter on ``MFNetwork``) to functions. Callers may also pass a callable
directly as ``match_method``.
"""
import numpy as np
import scipy.optimize as spo
import scipy.spatial as spsp
import starsim as ss
from bisect import bisect_left

from .base import NoPartnersFound

__all__ = ['MATCHERS', 'BAND_MATCH_WIDTH',
           'sort_bisect', 'sort_pair', 'kdtree_nn',
           'desired_age_bucket', 'greedy_old_enough', 'band_match', 'lsa']


BAND_MATCH_WIDTH = 5  # Year-band width for the ``band_match`` matcher.


def sort_bisect(net):
    """Current production: sort by age + bisect-trim of support tails + subsample.

    Does NOT honour ``age_diff_pars`` — the bisect trim only corrects the
    support boundaries, not the actual age gap. Kept for backward compat.
    """
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    m_ages = ppl.age[m_eligible]
    ind_m = np.argsort(m_ages, stable=True)
    ind_f = np.argsort(desired_ages, stable=True)

    if len(ind_m) == 0 or len(ind_f) == 0:
        raise NoPartnersFound()

    youngest_preferred_male_age = desired_ages[ind_f[0]]
    youngest_male_age = m_ages[ind_m[0]]
    if youngest_male_age < youngest_preferred_male_age:
        cutoff_index = bisect_left(m_ages[ind_m], youngest_preferred_male_age)
        ind_m = ind_m[cutoff_index:]
    elif youngest_preferred_male_age < youngest_male_age:
        cutoff_index = bisect_left(desired_ages[ind_f], youngest_male_age)
        ind_f = ind_f[cutoff_index:]
    if len(ind_m) == 0 or len(ind_f) == 0:
        raise NoPartnersFound()

    oldest_preferred_male_age = desired_ages[ind_f[-1]]
    oldest_male_age = m_ages[ind_m[-1]]
    if oldest_male_age > oldest_preferred_male_age:
        cutoff_index = bisect_left(m_ages[ind_m], oldest_preferred_male_age)
        ind_m = ind_m[:cutoff_index]
    elif oldest_preferred_male_age > oldest_male_age:
        cutoff_index = bisect_left(desired_ages[ind_f], oldest_male_age)
        ind_f = ind_f[:cutoff_index]

    if len(ind_m) < len(ind_f):
        ind_f_subset = np.random.choice(len(ind_f), size=len(ind_m), replace=False)
        ind_f_subset.sort()
        ind_f = ind_f[ind_f_subset]
    elif len(ind_f) < len(ind_m):
        ind_m_subset = np.random.choice(len(ind_m), size=len(ind_f), replace=False)
        ind_m_subset.sort()
        ind_m = ind_m[ind_m_subset]

    if len(ind_m) == 0 or len(ind_f) == 0:
        raise NoPartnersFound()

    p1 = m_eligible.uids[ind_m]
    p2 = f_looking[ind_f]
    return p1, p2


def sort_pair(net):
    """argsort both groups, zip, truncate to min(len). No tail trim."""
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    m_ages = ppl.age[m_eligible]
    ind_m = np.argsort(m_ages, stable=True)
    ind_f = np.argsort(desired_ages, stable=True)
    maxlen = min(len(ind_m), len(ind_f))
    if maxlen == 0:
        raise NoPartnersFound()
    ind_m = ind_m[:maxlen]
    ind_f = ind_f[:maxlen]
    return m_eligible.uids[ind_m], f_looking[ind_f]


def kdtree_nn(net):
    """Build KDTree on male ages; each woman queries 1-NN. Resolve collisions
    by giving each contested man to the closest-by-age woman.
    """
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
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


def desired_age_bucket(net):
    """Bucket women by integer desired-age; sample men in that bucket (with
    replacement if M < W). Post-filter by male concurrency cap and by the
    risk-group-conditional relationship-acceptance Bernoulli.

    Duplicates the acceptance logic that MFNetwork.add_pairs runs after
    matching — accepted technical debt; see spec.
    """
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    m_ages = ppl.age[m_eligible]
    m_uids = m_eligible.uids
    f_uids = f_looking

    desired_int = np.floor(desired_ages).astype(int)
    m_int = np.floor(m_ages).astype(int)

    buckets = {}
    for i, age in enumerate(m_int):
        buckets.setdefault(int(age), []).append(i)

    p1_list, p2_list = [], []
    rng = np.random.default_rng(net.sim.pars.rand_seed + net.ti)
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

    male_count = {}
    keep = np.ones(len(p1), dtype=bool)
    for i, uid in enumerate(p1):
        male_count[int(uid)] = male_count.get(int(uid), 0) + 1
        if male_count[int(uid)] > int(net.concurrency[uid]):
            keep[i] = False
    p1 = p1[keep]
    p2 = p2[keep]

    if len(p1) == 0:
        raise NoPartnersFound()

    matched_risk = (net.risk_group[p1] == net.risk_group[p2])
    mismatched_risk = ~matched_risk
    p_match = np.zeros(len(p1), dtype=float)
    for rg in range(net.pars.n_risk_groups):
        p_match[matched_risk & (net.risk_group[p1] == rg)] = net.pars.p_matched_stable[rg]
        p_match[mismatched_risk & (net.risk_group[p2] == rg)] = net.pars.p_mismatched_casual[rg]
    accept = rng.random(len(p1)) < p_match
    p1 = p1[accept]
    p2 = p2[accept]

    if len(p1) == 0:
        raise NoPartnersFound()

    return p1, p2


def greedy_old_enough(net):
    """Sort women by desired age ascending; for each, take youngest available
    male with age >= desired_age. No replacement.
    """
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
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


def band_match(net):
    """Bucket both groups into BAND_MATCH_WIDTH-year age bands;
    shuffle and zip within band.
    """
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    m_ages = ppl.age[m_eligible]
    m_uids = m_eligible.uids
    f_uids = f_looking

    f_band = (desired_ages // BAND_MATCH_WIDTH).astype(int)
    m_band = (m_ages // BAND_MATCH_WIDTH).astype(int)

    rng = np.random.default_rng(net.sim.pars.rand_seed + net.ti)
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


def lsa(net):
    """Linear sum assignment on the full age-distance matrix. O(n^3); reference only."""
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    m_ages = ppl.age[m_eligible]
    if len(m_ages) == 0 or len(desired_ages) == 0:
        raise NoPartnersFound()
    dist_mat = spsp.distance_matrix(m_ages[:, np.newaxis], desired_ages[:, np.newaxis])
    ind_m, ind_f = spo.linear_sum_assignment(dist_mat)
    return m_eligible.uids[ind_m], f_looking[ind_f]


MATCHERS = {
    'kdtree_nn':          kdtree_nn,
    'sort_bisect':        sort_bisect,
    'sort_pair':          sort_pair,
    'desired_age_bucket': desired_age_bucket,
    'greedy_old_enough':  greedy_old_enough,
    'band_match':         band_match,
    'lsa':                lsa,
}
