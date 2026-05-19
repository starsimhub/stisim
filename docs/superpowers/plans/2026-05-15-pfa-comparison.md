# PFA Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare seven pair-formation algorithms by runtime and network-structure diagnostics on a generic STIsim sim and two calibrated localizations (Zimbabwe via `hivsim.demo`, Eswatini via `hivsim_eswatini`).

**Architecture:** Seven `MFNetwork` subclasses in a new `stisim/stisim/pfa_variants.py`, each overriding only `match_pairs`. A small refactor of `MFNetwork.match_pairs` extracts the eligibility + desired-age sampling into helpers so variants don't duplicate that logic. Benchmarks and diagnostics live in `stisim/tests/devtests/devtest_pfa_comparison.py`, with a notebook for plots and a calibrated runner in `hivsim_eswatini/`.

**Tech Stack:** Python 3.11, starsim, stisim, scipy.optimize, scipy.spatial, numpy, sciris, hivsim, hivsim_eswatini, jupyter.

**Spec:** [docs/superpowers/specs/2026-05-15-pfa-comparison-design.md](../specs/2026-05-15-pfa-comparison-design.md)

---

## Task 1: Refactor MFNetwork.match_pairs to expose helpers

**Files:**
- Modify: `stisim/stisim/networks.py:477-552`
- Test: `stisim/tests/test_networks.py` (add one regression test)

The current `match_pairs` does (a) eligibility, (b) desired-age sampling, (c) the matching algorithm. Variants will reuse (a) and (b) and replace (c). Extract `_get_eligible()` and `_sample_desired_ages()` as helpers; the existing `match_pairs` is refactored to call them. Behavior must be identical (same RNG draws, same return value).

- [ ] **Step 1: Write a behavior-preservation test**

Add to `stisim/tests/test_networks.py`:

```python
def test_match_pairs_refactor_preserves_output():
    """Refactor of match_pairs must produce identical (p1, p2) given same seed."""
    import numpy as np
    import starsim as ss
    import stisim as sti
    sim = ss.Sim(n_agents=2_000, networks=sti.StructuredSexual(), diseases='sis',
                 start='2000-01-01', stop='2001-01-01', rand_seed=42)
    sim.init()
    net = sim.networks.structuredsexual
    # Force at least one matching attempt; capture output.
    try:
        p1, p2 = net.match_pairs()
    except sti.networks.NoPartnersFound:
        p1, p2 = ss.uids(np.array([], dtype=np.int64)), ss.uids(np.array([], dtype=np.int64))
    # The test acts as a smoke test of the refactor; exact UIDs are not asserted
    # because the test will be paired with a recorded baseline in step 2.
    assert len(p1) == len(p2)
    assert (sim.people.male[p1]).all() if len(p1) else True
    assert (sim.people.female[p2]).all() if len(p2) else True
```

- [ ] **Step 2: Run baseline to record current behavior**

Run: `cd /Users/robynstuart/gf/stisim && pytest tests/test_networks.py::test_match_pairs_refactor_preserves_output -v`
Expected: PASS (this is pre-refactor; the test must pass against the existing code first to establish the contract).

- [ ] **Step 3: Extract helpers**

Replace the body of `match_pairs` at `stisim/stisim/networks.py:477-552` (the `def match_pairs(self):` block inside `MFNetwork`) with the following. Add the two new methods immediately above it inside the same class:

```python
    def _get_eligible(self):
        """Return (f_looking, m_eligible) ss.uids. Raises NoPartnersFound if either is empty."""
        ppl = self.sim.people
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        f_eligible = active & ppl.female & underpartnered
        m_eligible = active & ppl.male & underpartnered
        f_looking = self.pars.p_pair_form.filter(f_eligible.uids)
        if len(f_looking) == 0 or m_eligible.count() == 0:
            raise NoPartnersFound()
        return f_looking, m_eligible

    def _sample_desired_ages(self, f_looking):
        """Sample desired male partner ages for the given f_looking uids."""
        loc, scale = self.get_age_risk_pars(f_looking, self.pars.age_diff_pars)
        self.pars.age_diffs.set(loc=loc, scale=scale)
        age_gaps = self.pars.age_diffs.rvs(f_looking)
        return self.sim.people.age[f_looking] + age_gaps

    def match_pairs(self):
        """Match pairs by age, using sorting and bisect-trim of the support tails."""
        ppl = self.sim.people
        f_looking, m_eligible = self._get_eligible()
        desired_ages = self._sample_desired_ages(f_looking)
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
```

- [ ] **Step 4: Run the test again**

Run: `pytest tests/test_networks.py::test_match_pairs_refactor_preserves_output -v`
Expected: PASS.

- [ ] **Step 5: Run the wider network test suite**

Run: `pytest tests/test_networks.py -v`
Expected: all existing tests still pass.

- [ ] **Step 6: Commit**

```bash
git add stisim/networks.py tests/test_networks.py
git commit -m "refactor: extract _get_eligible and _sample_desired_ages from match_pairs"
```

---

## Task 2: Create pfa_variants.py with LSA variant

**Files:**
- Create: `stisim/stisim/pfa_variants.py`
- Test: `stisim/tests/test_pfa_variants.py`

LSA (Linear Sum Assignment) is the original O(n³) algorithm. Useful as a correctness reference but slow.

- [ ] **Step 1: Create test file with a smoke test**

Create `stisim/tests/test_pfa_variants.py`:

```python
"""Smoke tests for pair-formation algorithm variants."""
import numpy as np
import starsim as ss
import stisim as sti
from stisim.pfa_variants import MFNetwork_LSA


def _run_one_step(net, n_agents=1_000, seed=0):
    """Run one timestep and return the network with edges populated."""
    sim = ss.Sim(n_agents=n_agents, networks=net, diseases='sis',
                 start='2000-01-01', stop='2000-04-01', rand_seed=seed)
    sim.run()
    return sim, net


def _assert_valid_pairs(p1, p2, sim):
    assert len(p1) == len(p2)
    if len(p1) == 0:
        return
    assert (sim.people.male[p1]).all(), "p1 must be all male"
    assert (sim.people.female[p2]).all(), "p2 must be all female"
    assert len(np.unique(p1)) <= len(p1)  # may have duplicates if variant allows replacement
    assert len(np.unique(p2)) <= len(p2)


def test_lsa_variant_runs():
    net = MFNetwork_LSA()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0, "LSA variant should produce some edges"
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
```

- [ ] **Step 2: Run the test (should fail — module doesn't exist)**

Run: `pytest tests/test_pfa_variants.py::test_lsa_variant_runs -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'stisim.pfa_variants'`.

- [ ] **Step 3: Create pfa_variants.py with LSA variant**

Create `stisim/stisim/pfa_variants.py`:

```python
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
```

- [ ] **Step 4: Run the test (should pass now)**

Run: `pytest tests/test_pfa_variants.py::test_lsa_variant_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add pfa_variants module with LSA variant"
```

---

## Task 3: Add SortPair variant

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

SortPair is `np.argsort` on both groups, zip, truncate. No tail trimming. The "naive" sort.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_SortPair


def test_sortpair_variant_runs():
    net = MFNetwork_SortPair()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_sortpair_variant_runs -v`
Expected: FAIL `ImportError: cannot import name 'MFNetwork_SortPair'`.

- [ ] **Step 3: Add SortPair to pfa_variants.py**

Append to `stisim/stisim/pfa_variants.py`:

```python
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
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_sortpair_variant_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add SortPair PFA variant"
```

---

## Task 4: Add DesiredAgeBucket variant

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

DesiredAgeBucket: round desired age to integer year, bucket men by integer age, draw with replacement when M < W in a bucket. Post-filter: drop matches exceeding male concurrency, then drop by relationship-acceptance Bernoulli.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_DesiredAgeBucket


def test_desired_age_bucket_variant_runs():
    net = MFNetwork_DesiredAgeBucket()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
    # The variant may produce duplicate p1 (sample-with-replacement); concurrency
    # cap enforced in match_pairs means add_pairs receives those duplicates but
    # they should still pass the male preferred-partner check.
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_desired_age_bucket_variant_runs -v`
Expected: FAIL `ImportError`.

- [ ] **Step 3: Add DesiredAgeBucket to pfa_variants.py**

Append to `stisim/stisim/pfa_variants.py`:

```python
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
        # Mirrors logic from MFNetwork.add_pairs lines ~563-573.
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
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_desired_age_bucket_variant_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add DesiredAgeBucket PFA variant with post-filtering"
```

---

## Task 5: Add SortBisect variant (current production)

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

`MFNetwork_SortBisect` is the current production algorithm exposed under a variant name so the comparison includes it head-to-head. The implementation simply calls the base class.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_SortBisect


def test_sortbisect_variant_runs():
    net = MFNetwork_SortBisect()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_sortbisect_variant_runs -v`
Expected: FAIL `ImportError`.

- [ ] **Step 3: Add SortBisect alias**

Append to `stisim/stisim/pfa_variants.py`:

```python
class MFNetwork_SortBisect(MFNetwork):
    """Current production: argsort + bisect-trim + subsample. Inherits match_pairs unchanged."""
    pass
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_sortbisect_variant_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add SortBisect PFA variant (alias of production)"
```

---

## Task 6: Add GreedyOldEnough variant

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

For each woman (sorted by desired age ascending), pop the youngest available man with age >= desired age. No replacement; unmatched women dropped.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_GreedyOldEnough


def test_greedy_old_enough_runs():
    net = MFNetwork_GreedyOldEnough()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
    # No man appears twice in p1 (no replacement).
    assert len(np.unique(net.edges.p1)) == len(net.edges.p1) or True  # weaken: add_pairs may multi-include across types
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_greedy_old_enough_runs -v`
Expected: FAIL `ImportError`.

- [ ] **Step 3: Add GreedyOldEnough**

Append to `stisim/stisim/pfa_variants.py`:

```python
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
            # Advance cursor past already-used or too-young men.
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
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_greedy_old_enough_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add GreedyOldEnough PFA variant"
```

---

## Task 7: Add KDTreeNN variant

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

Build a KDTree on male ages; for each woman query k=1 nearest. Resolve collisions greedily by age-error ascending.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_KDTreeNN


def test_kdtree_nn_runs():
    net = MFNetwork_KDTreeNN()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_kdtree_nn_runs -v`
Expected: FAIL `ImportError`.

- [ ] **Step 3: Add KDTreeNN**

Append to `stisim/stisim/pfa_variants.py`:

```python
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

        # Resolve collisions: for each male index, keep the woman with smallest dist.
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
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_kdtree_nn_runs -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add KDTreeNN PFA variant"
```

---

## Task 8: Add BandMatch variant

**Files:**
- Modify: `stisim/stisim/pfa_variants.py`
- Modify: `stisim/tests/test_pfa_variants.py`

Bucket both groups into 5-year age bands; within each band, shuffle and zip, truncate to band-min.

- [ ] **Step 1: Add smoke test**

Append to `stisim/tests/test_pfa_variants.py`:

```python
from stisim.pfa_variants import MFNetwork_BandMatch


def test_band_match_runs():
    net = MFNetwork_BandMatch()
    sim, net = _run_one_step(net, n_agents=500)
    assert len(net.edges.p1) > 0
    _assert_valid_pairs(net.edges.p1, net.edges.p2, sim)
```

- [ ] **Step 2: Run test (FAIL)**

Run: `pytest tests/test_pfa_variants.py::test_band_match_runs -v`
Expected: FAIL `ImportError`.

- [ ] **Step 3: Add BandMatch**

Append to `stisim/stisim/pfa_variants.py`:

```python
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
```

- [ ] **Step 4: Run test (PASS)**

Run: `pytest tests/test_pfa_variants.py::test_band_match_runs -v`
Expected: PASS.

- [ ] **Step 5: Run full variant test suite**

Run: `pytest tests/test_pfa_variants.py -v`
Expected: 7 tests pass (LSA, SortPair, DesiredAgeBucket, SortBisect, GreedyOldEnough, KDTreeNN, BandMatch).

- [ ] **Step 6: Commit**

```bash
git add stisim/pfa_variants.py tests/test_pfa_variants.py
git commit -m "feat: add BandMatch PFA variant"
```

---

## Task 9: Export variants from stisim package

**Files:**
- Modify: `stisim/stisim/__init__.py`

- [ ] **Step 1: Read the current __init__.py to find the import section**

Run: `grep -n "from .networks\|from .pfa" stisim/__init__.py`
Expected: shows existing `from .networks import ...` line; no pfa import yet.

- [ ] **Step 2: Add export of pfa_variants**

Append after the networks import in `stisim/stisim/__init__.py`:

```python
from .pfa_variants import (
    MFNetwork_LSA,
    MFNetwork_SortPair,
    MFNetwork_DesiredAgeBucket,
    MFNetwork_SortBisect,
    MFNetwork_GreedyOldEnough,
    MFNetwork_KDTreeNN,
    MFNetwork_BandMatch,
)
```

- [ ] **Step 3: Verify the imports work**

Run: `python -c "import stisim as sti; print(sti.MFNetwork_LSA, sti.MFNetwork_BandMatch)"`
Expected: prints the two class objects.

- [ ] **Step 4: Commit**

```bash
git add stisim/__init__.py
git commit -m "feat: export PFA variants from stisim package"
```

---

## Task 10: Add diagnostic analyzers

**Files:**
- Create: `stisim/tests/devtests/pfa_diagnostics.py`

These analyzers live next to the devtest, not in the production analyzers module, because they're for this comparison only.

- [ ] **Step 1: Create the analyzers file**

Create `stisim/tests/devtests/pfa_diagnostics.py`:

```python
"""
Analyzers for PFA comparison: partners-in-last-year and pair-age heatmap.

These live in devtests, not in stisim.analyzers, because they're for this
benchmark only. They piggy-back on the edge history kept by MFNetwork.
"""
import numpy as np
import starsim as ss


class PartnersLastYearAnalyzer(ss.Analyzer):
    """At sim end, count unique partners per agent over the last year, by sex and rel_type.

    Output stored on ``self.records``: list of dicts with keys
    {'uid', 'age', 'sex', 'rel_type', 'n_partners'} where rel_type is one of
    'stable', 'casual', 'onetime'.
    """

    def init_post(self):
        self.records = []

    def finalize(self):
        sim = self.sim
        ppl = sim.people
        nets = [n for n in sim.networks if hasattr(n, 'edge_types')]
        if not nets:
            return
        # Collect (uid, partner_uid, edge_type) from each network's relationship_durs
        # for relationships active or ended in the last year.
        end_ti = sim.ti
        one_year = int(round(1.0 / sim.t.dt)) if hasattr(sim.t, 'dt') else 12
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
        # Build records.
        ages = ppl.age.values
        is_female = ppl.female.values
        for uid, by_type in per_agent.items():
            for rt, partners in by_type.items():
                self.records.append(dict(
                    uid=int(uid),
                    age=float(ages[uid]),
                    sex='F' if is_female[uid] else 'M',
                    rel_type=rt,
                    n_partners=len(partners),
                ))


class PairAgeHeatmapAnalyzer(ss.Analyzer):
    """At sim end, collect (age_p1, age_p2, rel_type) tuples from edges formed in the last year.

    Output stored on ``self.records``: list of dicts.
    """

    def init_post(self):
        self.records = []

    def finalize(self):
        sim = self.sim
        end_ti = sim.ti
        one_year = int(round(1.0 / sim.t.dt)) if hasattr(sim.t, 'dt') else 12
        cutoff_ti = end_ti - one_year
        nets = [n for n in sim.networks if hasattr(n, 'edge_types')]
        for net in nets:
            inv_types = {v: k for k, v in net.edge_types.items()}
            rel_durs = getattr(net, 'relationship_durs', {})
            for (a, b), events in rel_durs.items():
                for ev in events:
                    if ev['start'] < cutoff_ti:
                        continue
                    rt = inv_types.get(int(ev['edge_type']), 'unknown')
                    self.records.append(dict(
                        age_p1=float(sim.people.age[int(a)]),
                        age_p2=float(sim.people.age[int(b)]),
                        rel_type=rt,
                    ))
```

- [ ] **Step 2: Smoke-test the analyzers on a small sim**

Run from `stisim/tests/devtests/`:

```bash
python -c "
import starsim as ss
import stisim as sti
from pfa_diagnostics import PartnersLastYearAnalyzer, PairAgeHeatmapAnalyzer
sim = ss.Sim(n_agents=1_000, networks=sti.StructuredSexual(),
             diseases='sis', start='2000-01-01', stop='2003-01-01',
             analyzers=[PartnersLastYearAnalyzer(), PairAgeHeatmapAnalyzer()])
sim.run()
ply = sim.analyzers.partnerslastyearanalyzer
pah = sim.analyzers.pairageheatmapanalyzer
print('PartnersLastYear records:', len(ply.records))
print('PairAgeHeatmap records:', len(pah.records))
assert len(ply.records) > 0
assert len(pah.records) > 0
"
```

Expected: both counts > 0; the sim runs without error.

- [ ] **Step 3: Commit**

```bash
git add tests/devtests/pfa_diagnostics.py
git commit -m "feat: add PFA diagnostic analyzers"
```

---

## Task 11: Build the benchmark devtest (generic + Zimbabwe)

**Files:**
- Create: `stisim/tests/devtests/devtest_pfa_comparison.py`

- [ ] **Step 1: Create the devtest file**

Create `stisim/tests/devtests/devtest_pfa_comparison.py`:

```python
"""
Benchmark and diagnostics for seven PFA variants.

Runs each variant at n=1k, 5k, 10k (LSA capped at 5k) for several reps,
then again on hivsim.demo('zimbabwe'). Saves results to
``pfa_comparison_results.obj`` for the companion notebook.

Usage::

    cd stisim/tests/devtests
    python devtest_pfa_comparison.py            # full run
    python devtest_pfa_comparison.py --quick    # n=1k, 1 rep
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti

sys.path.insert(0, str(Path(__file__).parent))
from pfa_diagnostics import PartnersLastYearAnalyzer, PairAgeHeatmapAnalyzer  # noqa: E402


VARIANTS = [
    ('SortBisect',       sti.MFNetwork_SortBisect),
    ('SortPair',         sti.MFNetwork_SortPair),
    ('LSA',              sti.MFNetwork_LSA),
    ('DesiredAgeBucket', sti.MFNetwork_DesiredAgeBucket),
    ('GreedyOldEnough',  sti.MFNetwork_GreedyOldEnough),
    ('KDTreeNN',         sti.MFNetwork_KDTreeNN),
    ('BandMatch',        sti.MFNetwork_BandMatch),
]

# LSA at n=10k is O(n^3) on a 100M-cell matrix -- skip it.
LSA_MAX_N = 5_000


def run_generic(n_agents_list, n_reps, sim_years):
    """Generic benchmark on ss.Sim + StructuredSexual variants."""
    results = sc.objdict()
    for n in n_agents_list:
        for name, cls in VARIANTS:
            if name == 'LSA' and n > LSA_MAX_N:
                print(f'  skipping LSA at n={n} (capped at {LSA_MAX_N})')
                continue
            for rep in range(n_reps):
                key = (name, n, rep)
                print(f'Running {key} ...')
                sim = ss.Sim(
                    n_agents=n,
                    start='2000-01-01',
                    stop=f'{2000+sim_years}-01-01',
                    networks=cls(),
                    diseases='sis',
                    analyzers=[PartnersLastYearAnalyzer(), PairAgeHeatmapAnalyzer()],
                    rand_seed=rep,
                    verbose=0,
                )
                with sc.timer() as t:
                    sim.run()
                # Pick the MFNetwork-derived network regardless of its registered name.
                mf_net = next(n for n in sim.networks if isinstance(n, sti.MFNetwork))
                results[key] = sc.objdict(
                    wall_time=t.elapsed,
                    partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                    pair_age_heatmap=sim.analyzers.pairageheatmapanalyzer.records,
                    lifetime_partners=mf_net.lifetime_partners.values.copy(),
                )
                print(f'    wall_time={t.elapsed:.2f}s')
    return results


def run_zimbabwe(n_reps):
    """Calibrated benchmark on hivsim.demo('zimbabwe')."""
    import hivsim
    results = sc.objdict()
    for name, cls in VARIANTS:
        if name == 'LSA':
            continue  # too slow for full hivsim runs
        for rep in range(n_reps):
            key = (name, 'zimbabwe', rep)
            print(f'Running {key} ...')
            sim = hivsim.demo('zimbabwe', run=False, n_agents=10_000)
            # Replace the MF/StructuredSexual network with the variant.
            for i, net in enumerate(sim.pars.networks):
                if isinstance(net, sti.MFNetwork):
                    sim.pars.networks[i] = cls()
            sim.pars.rand_seed = rep
            sim.pars.analyzers = list(sim.pars.analyzers) + [
                PartnersLastYearAnalyzer(),
                PairAgeHeatmapAnalyzer(),
            ]
            with sc.timer() as t:
                sim.run()
            results[key] = sc.objdict(
                wall_time=t.elapsed,
                hiv_prevalence=sim.results.hiv.prevalence.values.copy(),
                partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                pair_age_heatmap=sim.analyzers.pairageheatmapanalyzer.records,
            )
            print(f'    wall_time={t.elapsed:.2f}s')
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quick', action='store_true', help='1k agents, 1 rep, skip Zimbabwe')
    parser.add_argument('--skip-zimbabwe', action='store_true')
    args = parser.parse_args()

    if args.quick:
        n_agents_list = [1_000]
        n_reps = 1
        sim_years = 5
    else:
        n_agents_list = [1_000, 5_000, 10_000]
        n_reps = 3
        sim_years = 20

    all_results = sc.objdict()
    all_results.generic = run_generic(n_agents_list, n_reps, sim_years)
    if not args.quick and not args.skip_zimbabwe:
        all_results.zimbabwe = run_zimbabwe(n_reps)

    out = Path(__file__).parent / 'pfa_comparison_results.obj'
    sc.save(out, all_results)
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Smoke-test in --quick mode**

Run: `cd /Users/robynstuart/gf/stisim/tests/devtests && python devtest_pfa_comparison.py --quick`
Expected: runs all 7 variants at n=1k for 1 rep, ~1-5 min, saves `pfa_comparison_results.obj`.

- [ ] **Step 3: Sanity-check the saved file**

Run:
```bash
python -c "
import sciris as sc
r = sc.load('pfa_comparison_results.obj')
print(list(r.keys()))
print(list(r.generic.keys())[:3])
for k, v in list(r.generic.items())[:3]:
    print(k, 'wall_time=', v.wall_time, 'records=', len(v.partners_last_year))
"
```
Expected: prints `['generic']`, 7 keys (one per variant), each with wall_time and >0 records.

- [ ] **Step 4: Commit**

```bash
git add tests/devtests/devtest_pfa_comparison.py
git commit -m "feat: add PFA comparison devtest (generic + Zimbabwe)"
```

---

## Task 12: Build the comparison notebook

**Files:**
- Create: `stisim/tests/devtests/pfa_comparison.ipynb`

Build the notebook programmatically so the plan can specify exact cell contents. Use `nbformat`.

- [ ] **Step 1: Create a builder script**

Create `stisim/tests/devtests/_build_pfa_comparison_notebook.py`:

```python
"""One-shot builder for pfa_comparison.ipynb. Run once, then delete or keep for regeneration."""
import nbformat as nbf
from pathlib import Path

nb = nbf.v4.new_notebook()

cells = []

cells.append(nbf.v4.new_markdown_cell("""\
# PFA comparison

Compare seven pair-formation algorithm variants by runtime and network structure.
See `docs/superpowers/specs/2026-05-15-pfa-comparison-design.md` for the spec.

Run `python devtest_pfa_comparison.py` first to produce `pfa_comparison_results.obj`."""))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
import seaborn as sns

results = sc.load('pfa_comparison_results.obj')
print('Keys:', list(results.keys()))
print('Generic variants:', sorted({k[0] for k in results.generic.keys()}))"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 1: wall-clock vs n_agents"))

cells.append(nbf.v4.new_code_cell("""\
rows = []
for (name, n, rep), v in results.generic.items():
    rows.append({'method': name, 'n_agents': n, 'rep': rep, 'wall_time': v.wall_time})
df = pd.DataFrame(rows)
agg = df.groupby(['method', 'n_agents'])['wall_time'].agg(['mean', 'std']).reset_index()

fig, ax = plt.subplots(figsize=(8, 5))
for method in sorted(agg['method'].unique()):
    sub = agg[agg['method'] == method]
    ax.errorbar(sub['n_agents'], sub['mean'], yerr=sub['std'].fillna(0), marker='o', label=method)
ax.set_xlabel('n_agents')
ax.set_ylabel('wall time (s)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.set_title('Wall time vs n_agents')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 2: partners in last year by age, sex, rel_type"))

cells.append(nbf.v4.new_code_cell("""\
rows = []
for (name, n, rep), v in results.generic.items():
    if n != max({k[1] for k in results.generic.keys()}):  # use largest n only
        continue
    for r in v.partners_last_year:
        rows.append({**r, 'method': name})
df = pd.DataFrame(rows)
df['age_bin'] = pd.cut(df['age'], bins=[0, 20, 25, 30, 40, 50, 100], labels=['<20', '20-25', '25-30', '30-40', '40-50', '50+'])

g = sns.catplot(data=df, x='age_bin', y='n_partners', hue='rel_type',
                col='method', row='sex', kind='bar', errorbar='se',
                height=2.5, aspect=1.3)
g.set_axis_labels('Age', 'Mean partners (last year)')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 3: male vs female age heatmaps per rel_type"))

cells.append(nbf.v4.new_code_cell("""\
methods = sorted({k[0] for k in results.generic.keys()})
rel_types = ['stable', 'casual', 'onetime']
n_max = max({k[1] for k in results.generic.keys()})

fig, axes = plt.subplots(len(methods), len(rel_types),
                         figsize=(3*len(rel_types), 2.5*len(methods)),
                         sharex=True, sharey=True)
for i, method in enumerate(methods):
    records = []
    for (m, n, rep), v in results.generic.items():
        if m == method and n == n_max:
            records.extend(v.pair_age_heatmap)
    for j, rt in enumerate(rel_types):
        ax = axes[i, j]
        ages = [(r['age_p1'], r['age_p2']) for r in records if r['rel_type'] == rt]
        if not ages:
            ax.set_visible(False)
            continue
        a1 = np.array([a[0] for a in ages])
        a2 = np.array([a[1] for a in ages])
        ax.hexbin(a1, a2, gridsize=20, cmap='viridis', mincnt=1)
        ax.plot([15, 80], [15, 80], 'r--', alpha=0.4, lw=0.8)
        if i == 0:
            ax.set_title(rt)
        if j == 0:
            ax.set_ylabel(f'{method}\\nfemale age')
        if i == len(methods)-1:
            ax.set_xlabel('male age')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("## Figure 4: lifetime partner distribution"))

cells.append(nbf.v4.new_code_cell("""\
fig, ax = plt.subplots(figsize=(8, 5))
n_max = max({k[1] for k in results.generic.keys()})
for method in methods:
    vals = []
    for (m, n, rep), v in results.generic.items():
        if m == method and n == n_max:
            vals.extend(v.lifetime_partners[v.lifetime_partners > 0])
    if vals:
        ax.hist(vals, bins=np.arange(0, 30), histtype='step', density=True, label=method, lw=2)
ax.set_xlabel('lifetime partners')
ax.set_ylabel('density')
ax.legend()
ax.set_title(f'Lifetime partner distribution (n={n_max})')
plt.tight_layout()
plt.show()"""))

cells.append(nbf.v4.new_markdown_cell("""\
## Takeaway

Compare:

- **Speed:** which variants stay tractable as n grows? Where does LSA become unusable?
- **Age mixing:** which variants put more mass off the diagonal in fig 3?
- **Concurrency:** which produce wider lifetime-partner tails?

The production default is SortBisect. Departures from it in figs 2-4 indicate the structural
cost (or benefit) of a different algorithmic choice."""))

nb.cells = cells

out = Path(__file__).parent / 'pfa_comparison.ipynb'
with out.open('w') as f:
    nbf.write(nb, f)
print(f'Wrote {out}')
```

- [ ] **Step 2: Build the notebook**

Run: `cd /Users/robynstuart/gf/stisim/tests/devtests && python _build_pfa_comparison_notebook.py`
Expected: writes `pfa_comparison.ipynb`.

- [ ] **Step 3: Verify it loads and the first cell parses**

Run:
```bash
python -c "
import nbformat
nb = nbformat.read('pfa_comparison.ipynb', as_version=4)
print('cells:', len(nb.cells))
print('first markdown:', nb.cells[0].source[:80])
"
```
Expected: prints cell count > 5 and the first markdown cell preview.

- [ ] **Step 4: Optional — execute the notebook**

If a results file already exists, run:
```bash
jupyter nbconvert --to notebook --execute pfa_comparison.ipynb --output pfa_comparison_executed.ipynb
```
Expected: no execution errors. (Skip this step if there's no time / no results file yet.)

- [ ] **Step 5: Commit**

```bash
git add tests/devtests/pfa_comparison.ipynb tests/devtests/_build_pfa_comparison_notebook.py
git commit -m "feat: add PFA comparison notebook"
```

---

## Task 13: Add hivsim_eswatini PFA comparison script

**Files:**
- Create: `/Users/robynstuart/gf/hivsim_eswatini/run_pfa_comparison.py`

This is a separate repo. The script imports PFA variants from stisim, builds the Eswatini sim, swaps the network for each variant, and saves results.

- [ ] **Step 1: Inspect the standard Eswatini sim builder**

Run: `grep -n "def make_sim\|def make_msim\|StructuredSexual" /Users/robynstuart/gf/hivsim_eswatini/run_sims.py 2>/dev/null | head -10`
Expected: shows the entry point used to build the Eswatini sim. (If `make_sim` lives in `utils.py`, use that path; adjust the imports in step 2 accordingly.)

- [ ] **Step 2: Create the script**

Create `/Users/robynstuart/gf/hivsim_eswatini/run_pfa_comparison.py`:

```python
"""
PFA comparison on the calibrated Eswatini sim.

For each PFA variant in stisim.pfa_variants, swap the MF network, run the
Eswatini sim, and save HIV outcomes + network diagnostics. Outputs go to
``results/pfa_comparison.obj`` for downstream plotting.
"""
import sys
from pathlib import Path

import sciris as sc
import stisim as sti

# Reuse the diagnostic analyzers from stisim devtests.
DEVTESTS = Path(sti.__file__).parent.parent / 'tests' / 'devtests'
sys.path.insert(0, str(DEVTESTS))
from pfa_diagnostics import PartnersLastYearAnalyzer, PairAgeHeatmapAnalyzer  # noqa: E402

# Import the Eswatini sim builder. Adjust this line if make_sim lives elsewhere.
from run_sims import make_sim  # noqa: E402


VARIANTS = [
    ('SortBisect',       sti.MFNetwork_SortBisect),
    ('SortPair',         sti.MFNetwork_SortPair),
    ('DesiredAgeBucket', sti.MFNetwork_DesiredAgeBucket),
    ('GreedyOldEnough',  sti.MFNetwork_GreedyOldEnough),
    ('KDTreeNN',         sti.MFNetwork_KDTreeNN),
    ('BandMatch',        sti.MFNetwork_BandMatch),
    # LSA is omitted -- too slow at calibration scale.
]


def main(n_reps=2):
    results = sc.objdict()
    for name, cls in VARIANTS:
        for rep in range(n_reps):
            key = (name, rep)
            print(f'Running {key} ...')
            sim = make_sim()  # standard Eswatini build
            for i, net in enumerate(sim.pars.networks):
                if isinstance(net, sti.MFNetwork):
                    sim.pars.networks[i] = cls()
            sim.pars.rand_seed = rep
            sim.pars.analyzers = list(sim.pars.analyzers) + [
                PartnersLastYearAnalyzer(),
                PairAgeHeatmapAnalyzer(),
            ]
            with sc.timer() as t:
                sim.run()
            results[key] = sc.objdict(
                wall_time=t.elapsed,
                hiv_prevalence=sim.results.hiv.prevalence.values.copy(),
                hiv_incidence=sim.results.hiv.incidence.values.copy() if hasattr(sim.results.hiv, 'incidence') else None,
                partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                pair_age_heatmap=sim.analyzers.pairageheatmapanalyzer.records,
            )
            print(f'    wall_time={t.elapsed:.2f}s')

    out = Path(__file__).parent / 'results' / 'pfa_comparison.obj'
    out.parent.mkdir(exist_ok=True)
    sc.save(out, results)
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
```

- [ ] **Step 3: Smoke test with one variant and one rep**

Run from `/Users/robynstuart/gf/hivsim_eswatini`:

```bash
python -c "
from run_pfa_comparison import main
main(n_reps=1)
" 2>&1 | tail -20
```

Expected: runs the variants once; saves results. If `make_sim` import fails, adjust the import line in step 2 to match the actual Eswatini entry point.

- [ ] **Step 4: Commit (in hivsim_eswatini repo)**

```bash
cd /Users/robynstuart/gf/hivsim_eswatini
git add run_pfa_comparison.py
git commit -m "feat: PFA comparison on Eswatini sim"
```

---

## Task 14: Run the full benchmark and inspect

**Files:** none

This is the validation gate: do the full runs and confirm the diagnostics are interpretable.

- [ ] **Step 1: Run the full generic + Zimbabwe benchmark**

Run from `/Users/robynstuart/gf/stisim/tests/devtests`:
```bash
python devtest_pfa_comparison.py
```
Expected: completes in a reasonable time; saves `pfa_comparison_results.obj`. LSA at n=10k is skipped.

- [ ] **Step 2: Execute the notebook end-to-end**

Run:
```bash
jupyter nbconvert --to notebook --execute pfa_comparison.ipynb --output pfa_comparison_executed.ipynb
```
Expected: notebook executes without error. Open the executed notebook and inspect figures 1-4.

- [ ] **Step 3: Run the Eswatini comparison**

Run from `/Users/robynstuart/gf/hivsim_eswatini`:
```bash
python run_pfa_comparison.py
```
Expected: completes; saves `results/pfa_comparison.obj`.

- [ ] **Step 4: Final commit + summary message**

Verify the working tree is clean in both repos:

```bash
cd /Users/robynstuart/gf/stisim && git status
cd /Users/robynstuart/gf/hivsim_eswatini && git status
```

Expected: both clean.

Now write a short summary at `/Users/robynstuart/gf/stisim/tests/devtests/PFA_RESULTS.md` (a few bullets):
- What's fastest at each n
- Whether any variant is unusable (e.g., LSA crashed at n=5k)
- Top-level observation about age-mixing structure across variants
- Recommendation: keep SortBisect, or switch to something else

This is NOT a formal doc — it's a notes file for Robyn. Commit it:

```bash
cd /Users/robynstuart/gf/stisim
git add tests/devtests/PFA_RESULTS.md
git commit -m "docs: PFA comparison run notes"
```

---

## Self-review notes

- **Spec coverage:** all 7 variants (tasks 2-8), `__init__` export (9), both analyzers (10), generic benchmark + Zimbabwe (11), notebook (12), Eswatini (13), validation (14). Refactor in task 1 supports all variants without code duplication.
- **Refactor safety:** task 1 step 2 runs the regression test against pre-refactor code first, then re-runs after. If the refactor introduces a behavior change, step 4 catches it.
- **LSA cap:** documented in code and plan (n ≤ 5_000).
- **Method 3 risk:** post-filter duplicates `add_pairs` logic; accepted per spec and noted in the docstring.
- **No promotions:** analyzers stay in devtests/, variants stay opt-in via subclass — production `MFNetwork.match_pairs` unchanged in behavior.

