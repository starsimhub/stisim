# `match_method` API + `networks/` package layout — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganise `stisim/networks.py` into a `networks/` sub-package mirroring `diseases/` and `interventions/`, then replace the seven `MFNetwork_*` subclasses with a string `match_method` parameter on `MFNetwork` that dispatches through a registry of matcher functions. Default matcher changes from current production (`sort_bisect`) to `kdtree_nn`.

**Architecture:** Two sequential commits on `feat/pfa-comparison`:
1. **Refactor** — pure file move; `networks.py` → `networks/` folder split across `base.py`, `mf.py`, `fsw.py`, `msm.py`, `layered_networks.py`. Re-exports in `networks/__init__.py` keep the flat public API intact (`sti.MFNetwork`, `sti.StructuredSexual`, etc.).
2. **Feature** — new `networks/matchers.py` with module-level matcher functions and a `MATCHERS` dict registry. `MFPars` gains `match_method = 'kdtree_nn'`. `MFNetwork.match_pairs` becomes a dispatcher accepting string OR callable. `stisim/pfa_variants.py` and its seven `MFNetwork_*` subclasses are deleted; tests + devtest move to string-based selection.

**Tech Stack:** Python 3.11, starsim, stisim, numpy, scipy.spatial, scipy.optimize, sciris, pytest.

**Spec:** [docs/superpowers/specs/2026-05-18-match-method-api-design.md](../specs/2026-05-18-match-method-api-design.md)

---

## Commit A: refactor `networks.py` → `networks/` package

### Task A1: Audit downstream `from stisim.networks import …` paths

**Files:**
- Read-only audit across all working trees under `/Users/robynstuart/gf/`

The refactor preserves the flat API (`from stisim.networks import MFNetwork` keeps working). But code that reaches into the implementation file directly (`from stisim.networks.networks import …`) would break. Verify no such code exists.

- [ ] **Step 1: Run the grep**

```bash
grep -rn "from stisim\.networks\.networks\|stisim\.networks\.networks\b" /Users/robynstuart/gf/ 2>/dev/null \
  | grep -v "\.pyc\|__pycache__\|/results/\|/\.git/"
```

Expected: no matches. If any matches appear, note them and patch alongside the refactor in Task A6.

- [ ] **Step 2: Also confirm no external `from stisim.pfa_variants` references**

```bash
grep -rn "from stisim\.pfa_variants\|import stisim\.pfa_variants\b" /Users/robynstuart/gf/ 2>/dev/null \
  | grep -v "\.pyc\|__pycache__\|/\.git/"
```

Expected: only `stisim/tests/test_pfa_variants.py:5` (which is rewritten in Commit B). If anything else appears, add it to the Commit B cleanup list.

### Task A2: Create `stisim/networks/` skeleton

**Files:**
- Create: `stisim/networks/__init__.py`
- Create: `stisim/networks/README.md`
- Move: `stisim/networks.py` → `stisim/networks/_legacy_full.py` (temporary; deleted at end of commit)

The strategy: keep the old `networks.py` content available as `_legacy_full.py` during the split so we can reference it from the new module files, then delete the temp file once everything is migrated. Git renames will detect both halves and produce a clean `diff -M`.

- [ ] **Step 1: Move the existing file under the new package**

```bash
mkdir -p stisim/networks
git mv stisim/networks.py stisim/networks/_legacy_full.py
```

- [ ] **Step 2: Create an empty `__init__.py` so the package imports**

```python
# stisim/networks/__init__.py
"""STIsim networks package — split across base, mf, fsw, msm, layered_networks, matchers."""

# Temporary during refactor; all symbols still come from _legacy_full
from ._legacy_full import *  # noqa: F401, F403
```

- [ ] **Step 3: Verify the package imports cleanly**

```bash
cd /Users/robynstuart/gf/stisim && python -c "import stisim; print(stisim.MFNetwork)"
```

Expected: `<class 'stisim.networks._legacy_full.MFNetwork'>` (or similar). No ImportError.

- [ ] **Step 4: Run existing tests once to record a green baseline**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_pfa_variants.py -x -q
```

Expected: all pass (pre-refactor; pfa_variants still imports the legacy subclasses).

- [ ] **Step 5: Commit the move**

```bash
git add stisim/networks/__init__.py stisim/networks/_legacy_full.py
git commit -m "refactor: move networks.py under stisim/networks/ package (1/N)"
```

### Task A3: Extract `base.py`

**Files:**
- Create: `stisim/networks/base.py`
- Modify: `stisim/networks/__init__.py`
- Modify: `stisim/networks/_legacy_full.py` (delete extracted classes)

Extract `NoPartnersFound`, `BasePars`, `NetworkPars`, `BaseNetwork`, plus the module-level `ss_float`/`ss_int` aliases. `MFPars` and `SWPars` stay in `_legacy_full` for now and move with their networks in later tasks. `NetworkPars` references `BasePars`, `MFPars`, `SWPars` — keep its dependencies straight by importing `MFPars`/`SWPars` from `_legacy_full` for the moment.

- [ ] **Step 1: Create `stisim/networks/base.py` with the extracted code**

```python
"""Base classes and shared parameters for stisim networks."""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

ss_float = ss.dtypes.float
ss_int = ss.dtypes.int

__all__ = ['NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork']


class NoPartnersFound(Exception):
    """Raise if the matching algorithm wasn't able to match any partners."""
    pass


class BasePars(ss.Pars):
    """Defaults shared by every sexual network (debut, acts, condom data, recall_prior)."""
    def __init__(self, **kwargs):
        super().__init__()
        self.recall_prior = False
        self.debut_f = ss.lognorm_ex(20, 3)
        self.debut_m = ss.lognorm_ex(21, 3)
        self.acts = ss.lognorm_ex(ss.freqperyear(80), ss.freqperyear(30))
        self.condom_data = None
        self.condom_smoothness = None
        self.update(kwargs)
        return


class NetworkPars(ss.Pars):
    """Combined defaults — union of base + MF + SW pars.

    Used by :class:`sti.Sim` to filter network kwargs into the right bucket.
    """
    def __init__(self, **kwargs):
        # Local import to avoid circular dependency at module load time.
        from .mf import MFPars
        from .fsw import SWPars
        super().__init__()
        for src in (BasePars(), MFPars(), SWPars()):
            self.update(dict(src), create=True)
        self.update(kwargs)
        return


class BaseNetwork(ss.SexualNetwork):
    """Shared infrastructure for the heterosexual and sex-work networks.

    Provides the common edge meta layout, condom-data processing, debut
    setting, end-of-step edge accounting, and the orchestrating ``step``
    method. Subclasses (``MFNetwork``, ``SWNetwork``) define their own
    parameters, states, edge types, ``set_network_states``, ``add_pairs``,
    and ``set_condom_use``.
    """

    def __init__(self, name=None, **kwargs):
        super().__init__(name=name)
        self.meta.condoms = ss_float
        self.meta.age_p1 = ss_float
        self.meta.age_p2 = ss_float
        self.meta.edge_type = ss_float
        self.define_pars(**BasePars())
        self.define_states(
            ss.BoolArr('participant', default=True),
            ss.FloatArr('debut', default=0),
            reset=True,
        )

    # === Paste the full body of BaseNetwork from _legacy_full here, ===
    # === unchanged. Methods to include: process_condom_data,        ===
    # === init_pre, init_post, set_network_states, over_debut,        ===
    # === _get_uids, set_debut, net_beta, step, add_pairs,            ===
    # === set_condom_use, end_pairs, _on_edge_dissolution.             ===
```

Open `stisim/networks/_legacy_full.py` and copy the body of `BaseNetwork` (lines 220-376 of the original file) into the placeholder above. Strip the duplicate `__init__` (already shown above).

- [ ] **Step 2: Delete the extracted code from `_legacy_full.py`**

In `_legacy_full.py`, delete:
- The `NoPartnersFound` class definition (around line 31).
- The `BasePars` class (around lines 77-88).
- The `NetworkPars` class (around lines 197-210).
- The `BaseNetwork` class (around lines 220-375).
- The `ss_float`/`ss_int` aliases at top (lines 22-23) — they live in `base.py` now, but `_legacy_full` still needs them for the remaining classes. Replace with `from .base import ss_float, ss_int`.

Add at the top of `_legacy_full.py` (just after the existing imports):

```python
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork, ss_float, ss_int
```

- [ ] **Step 3: Update `networks/__init__.py` to re-export from both**

```python
"""STIsim networks package — split across base, mf, fsw, msm, layered_networks, matchers."""

from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from ._legacy_full import *  # noqa: F401, F403 (remaining classes during refactor)
```

- [ ] **Step 4: Run the existing tests**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_pfa_variants.py tests/test_networks.py -x -q
```

Expected: all pass. If `tests/test_networks.py` doesn't exist or is empty, omit it from this command; the smoke test is `tests/test_pfa_variants.py`.

- [ ] **Step 5: Confirm `sti.NoPartnersFound`, `sti.BaseNetwork` etc. still resolve**

```bash
cd /Users/robynstuart/gf/stisim && python -c "
import stisim as sti
assert sti.NoPartnersFound
assert sti.BaseNetwork
assert sti.BasePars
assert sti.NetworkPars
print('OK')
"
```

Expected: prints `OK`.

- [ ] **Step 6: Commit**

```bash
git add stisim/networks/
git commit -m "refactor: extract base.py from networks/_legacy_full (2/N)"
```

### Task A4: Extract `mf.py`

**Files:**
- Create: `stisim/networks/mf.py`
- Modify: `stisim/networks/__init__.py`
- Modify: `stisim/networks/_legacy_full.py`

Extract `MFPars`, `_mf_states`, and `MFNetwork` into `mf.py`. Helpers (`_mf_states`) sit at module scope; `MFNetwork` references them by relative import.

- [ ] **Step 1: Create `stisim/networks/mf.py`**

```python
"""Heterosexual (MF) sexual network and its parameters."""
import numpy as np
import starsim as ss
from collections import defaultdict

from .base import BaseNetwork, NoPartnersFound, ss_float

__all__ = ['MFPars', 'MFNetwork']


def _mf_states():
    """States specific to the heterosexual (MF) network (excludes shared)."""
    return [
        ss.FloatArr('risk_group'),
        ss.FloatArr('concurrency'),
        ss.FloatArr('partners', default=0),
        ss.FloatArr('partners_12', default=0),
        ss.FloatArr('lifetime_partners', default=0),
        ss.FloatArr('casual_partners', default=0),
        ss.FloatArr('stable_partners', default=0),
        ss.FloatArr('onetime_partners', default=0),
        ss.FloatArr('lifetime_casual_partners', default=0),
        ss.FloatArr('lifetime_stable_partners', default=0),
        ss.FloatArr('lifetime_onetime_partners', default=0),
    ]


class MFPars(ss.Pars):
    # === Paste the full body of MFPars from _legacy_full here (lines 91-161). ===


class MFNetwork(BaseNetwork):
    # === Paste the full body of MFNetwork from _legacy_full here (lines 380-681). ===
```

Paste verbatim from `_legacy_full.py` for the bodies of `MFPars` and `MFNetwork`. `_mf_states` moves into this file; references from inside `MFNetwork.__init__` (`self.define_states(*_mf_states())`) work unchanged because `_mf_states` is now in scope at module level.

- [ ] **Step 2: Delete extracted code from `_legacy_full.py`**

Remove:
- `_mf_states` function (lines 38-52).
- `MFPars` class (lines 91-161).
- `MFNetwork` class (lines 380-681).

Add to the top imports section of `_legacy_full.py`:

```python
from .mf import MFPars, MFNetwork, _mf_states
```

`StructuredSexual` (which inherits `MFNetwork`) and `AgeMatchedMSM`/`AgeApproxMSM` still live in `_legacy_full.py` for now; they continue to find `MFNetwork` via this import.

- [ ] **Step 3: Update `networks/__init__.py`**

```python
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from .mf import MFPars, MFNetwork
from ._legacy_full import *  # noqa: F401, F403 (remaining classes during refactor)
```

- [ ] **Step 4: Run tests**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_pfa_variants.py -x -q
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add stisim/networks/
git commit -m "refactor: extract mf.py from networks/_legacy_full (3/N)"
```

### Task A5: Extract `fsw.py`

**Files:**
- Create: `stisim/networks/fsw.py`
- Modify: `stisim/networks/__init__.py`
- Modify: `stisim/networks/_legacy_full.py`

Extract `SWPars`, `_sw_states`, `SWNetwork` into `fsw.py`. The class name stays `SWNetwork`; only the file moves.

- [ ] **Step 1: Create `stisim/networks/fsw.py`**

```python
"""Sex-work (FSW–client) sexual network and its parameters.

The class is named ``SWNetwork`` (unchanged); the file is named ``fsw.py``
for discoverability.
"""
import numpy as np
import starsim as ss

from .base import BaseNetwork, NoPartnersFound, ss_float

__all__ = ['SWPars', 'SWNetwork']


def _sw_states():
    """States specific to the sex-work (SW) network (excludes shared)."""
    return [
        ss.BoolArr('ever_fsw'),
        ss.BoolArr('ever_client'),
        ss.FloatArr('age_sw_start'),
        ss.FloatArr('dur_sw'),
        ss.FloatArr('age_client_start'),
        ss.FloatArr('dur_client'),
        ss.FloatArr('sw_intensity'),
        ss.FloatArr('sw_partners', default=0),
        ss.FloatArr('lifetime_sw_partners', default=0),
    ]


class SWPars(ss.Pars):
    # === Paste the full body of SWPars from _legacy_full here (lines 164-194). ===


class SWNetwork(BaseNetwork):
    # === Paste the full body of SWNetwork from _legacy_full here (lines 683-843). ===
```

- [ ] **Step 2: Delete extracted code from `_legacy_full.py`**

Remove:
- `_sw_states` function (lines 55-67).
- `SWPars` class (lines 164-194).
- `SWNetwork` class (lines 683-843).

Add to the top imports section of `_legacy_full.py`:

```python
from .fsw import SWPars, SWNetwork, _sw_states
```

- [ ] **Step 3: Update `networks/__init__.py`**

```python
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from .mf import MFPars, MFNetwork
from .fsw import SWPars, SWNetwork
from ._legacy_full import *  # noqa: F401, F403
```

- [ ] **Step 4: Run tests**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_pfa_variants.py -x -q
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add stisim/networks/
git commit -m "refactor: extract fsw.py from networks/_legacy_full (4/N)"
```

### Task A6: Extract `msm.py` and `layered_networks.py`; remove `_legacy_full.py`

**Files:**
- Create: `stisim/networks/msm.py`
- Create: `stisim/networks/layered_networks.py`
- Delete: `stisim/networks/_legacy_full.py`
- Modify: `stisim/networks/__init__.py`

Final cleanup: extract the remaining classes and dispose of the temporary `_legacy_full.py`.

- [ ] **Step 1: Create `stisim/networks/msm.py`**

```python
"""Men-who-have-sex-with-men (MSM) sexual networks."""
import numpy as np
import starsim as ss

from .base import NoPartnersFound
from .mf import MFNetwork

__all__ = ['AgeMatchedMSM', 'AgeApproxMSM']


class AgeMatchedMSM(MFNetwork):
    # === Paste the full body of AgeMatchedMSM from _legacy_full (lines 936-996). ===


class AgeApproxMSM(MFNetwork):
    # === Paste the full body of AgeApproxMSM from _legacy_full (lines 999-1039). ===
```

- [ ] **Step 2: Create `stisim/networks/layered_networks.py`**

```python
"""Layered networks — composite networks that combine multiple network modules.

A "layered network" is a single class that bundles two or more network
behaviours so a sim can model their combined dynamics without juggling
multiple network objects. The canonical example is :class:`StructuredSexual`,
which bundles heterosexual partnerships (:class:`MFNetwork`) and sex-work
partnerships (:class:`SWNetwork`) on a single edge list.

:class:`PriorPartners` is a supplementary network you add *alongside* an MF
network when ``recall_prior=True`` — it stores recently-dissolved edges so
partner-notification interventions can trace them.
"""
import numpy as np
import starsim as ss

from .base import BaseNetwork, NoPartnersFound, ss_float
from .mf import MFNetwork, _mf_states
from .fsw import SWNetwork, SWPars, _sw_states

__all__ = ['StructuredSexual', 'PriorPartners']


class StructuredSexual(MFNetwork):
    # === Paste body from _legacy_full (lines 846-896). ===


class PriorPartners(ss.DynamicNetwork):
    # === Paste body from _legacy_full (lines 901-933). ===
```

- [ ] **Step 3: Delete `_legacy_full.py`**

At this point all classes have been extracted. `_legacy_full.py` should contain only imports and the original module docstring/`__all__`. Remove the file entirely:

```bash
git rm stisim/networks/_legacy_full.py
```

- [ ] **Step 4: Replace `networks/__init__.py` with the final form**

```python
"""STIsim networks package.

See ``README.md`` for the layout. Public API exports are flat: e.g.
``stisim.MFNetwork`` works, as does ``stisim.networks.MFNetwork``.

For direct access to matcher functions used by ``MFNetwork.match_pairs``,
use ``stisim.networks.matchers`` (see Commit B).
"""
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from .mf import MFPars, MFNetwork
from .fsw import SWPars, SWNetwork
from .msm import AgeMatchedMSM, AgeApproxMSM
from .layered_networks import StructuredSexual, PriorPartners

__all__ = [
    'NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork',
    'MFPars', 'MFNetwork',
    'SWPars', 'SWNetwork',
    'AgeMatchedMSM', 'AgeApproxMSM',
    'StructuredSexual', 'PriorPartners',
]
```

- [ ] **Step 5: Run the full test suite**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/ -x -q
```

Expected: all pass. If `tests/test_pfa_variants.py` still imports from `stisim.pfa_variants`, that's fine — that module hasn't been touched yet (it's deleted in Commit B).

- [ ] **Step 6: Smoke-check all public symbols**

```bash
cd /Users/robynstuart/gf/stisim && python -c "
import stisim as sti
for name in ['NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork',
             'MFPars', 'MFNetwork', 'SWPars', 'SWNetwork',
             'AgeMatchedMSM', 'AgeApproxMSM',
             'StructuredSexual', 'PriorPartners']:
    assert hasattr(sti, name), f'missing {name}'
print('OK')
"
```

Expected: prints `OK`.

- [ ] **Step 7: Commit**

```bash
git add stisim/networks/
git commit -m "refactor: complete networks/ package split (5/N)"
```

### Task A7: Write `stisim/networks/README.md`

**Files:**
- Create: `stisim/networks/README.md`

Mirror the style of `stisim/diseases/README.md` (look at it first to match tone). Make the layered-networks concept clear.

- [ ] **Step 1: Read `stisim/diseases/README.md` for style reference**

```bash
cd /Users/robynstuart/gf/stisim && cat stisim/diseases/README.md | head -60
```

- [ ] **Step 2: Write `stisim/networks/README.md`**

```markdown
# stisim/networks

Sexual contact networks for STIsim models. Each network is a module within
this package; the table below shows where each lives.

## Layout

| File                    | Provides                                                      |
|-------------------------|---------------------------------------------------------------|
| `base.py`               | `BaseNetwork`, `BasePars`, `NetworkPars`, `NoPartnersFound`   |
| `mf.py`                 | `MFNetwork`, `MFPars` — heterosexual partnerships             |
| `fsw.py`                | `SWNetwork`, `SWPars` — sex-work (FSW–client) partnerships    |
| `msm.py`                | `AgeMatchedMSM`, `AgeApproxMSM` — men who have sex with men   |
| `layered_networks.py`   | `StructuredSexual`, `PriorPartners` — composite networks      |
| `matchers.py`           | Pair-formation algorithms used by `MFNetwork.match_pairs`     |

## Layered networks

A *layered network* bundles two or more network behaviours into a single
class so a sim can model their combined dynamics without juggling multiple
network objects.

`StructuredSexual` is the canonical example: it layers heterosexual
partnerships (from `MFNetwork`) with sex-work partnerships (from
`SWNetwork`) on a single edge list. Most STIsim models in practice use
`StructuredSexual` because real populations have both regular and
transactional partnerships, and the diseases transmit across both.

`PriorPartners` is a supplementary network you add *alongside* an MF
network when `recall_prior=True`. It stores recently-dissolved edges so
partner-notification interventions can trace them.

## Pair-formation algorithms

`MFNetwork.match_pairs` dispatches through a registry in `matchers.py`,
selected by the `match_method` parameter (string or callable):

```python
import stisim as sti

# Default (kdtree_nn — nearest-neighbour by age)
net = sti.MFNetwork()

# Explicit choice
net = sti.MFNetwork(match_method='desired_age_bucket')

# Custom matcher
def my_matcher(net):
    ...
    return p1, p2  # both ss.uids
net = sti.MFNetwork(match_method=my_matcher)
```

See `matchers.py` for the seven built-in methods and their tradeoffs.
```

- [ ] **Step 3: Commit**

```bash
git add stisim/networks/README.md
git commit -m "docs: add networks/ package README (6/N)"
```

### Task A8: Squash commits A2–A7 into a single refactor commit (optional)

**Files:** none directly; `git rebase -i` is interactive and **NOT supported in this plan**. Skip this task — keep the six small commits as a clean refactor sequence. Reviewers can read them in order.

---

## Commit B: `match_method` parameter

### Task B1: Create `stisim/networks/matchers.py` with the seven matchers + registry

**Files:**
- Create: `stisim/networks/matchers.py`
- Modify: `stisim/networks/__init__.py`

Extract the matcher logic from each of the seven `MFNetwork_*` subclasses currently in `stisim/pfa_variants.py` into module-level functions in `matchers.py`. Each function takes `net` (an `MFNetwork` instance) and returns `(p1, p2)` as `ss.uids`. Behavior must match the existing subclasses exactly — same RNG calls, same return shape.

- [ ] **Step 1: Read the source matchers**

```bash
cd /Users/robynstuart/gf/stisim && cat stisim/pfa_variants.py
```

- [ ] **Step 2: Create `stisim/networks/matchers.py`**

```python
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
    """argsort both groups, zip, truncate to ``min(len)``. No tail trim."""
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

    Duplicates the acceptance logic that ``MFNetwork.add_pairs`` runs after
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
    male with ``age >= desired_age``. No replacement.
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
    """Bucket both groups into ``BAND_MATCH_WIDTH``-year age bands;
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
```

- [ ] **Step 3: Update `networks/__init__.py` to expose `matchers`**

```python
"""STIsim networks package.

See ``README.md`` for the layout.
"""
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from .mf import MFPars, MFNetwork
from .fsw import SWPars, SWNetwork
from .msm import AgeMatchedMSM, AgeApproxMSM
from .layered_networks import StructuredSexual, PriorPartners
from . import matchers
from .matchers import MATCHERS

__all__ = [
    'NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork',
    'MFPars', 'MFNetwork',
    'SWPars', 'SWNetwork',
    'AgeMatchedMSM', 'AgeApproxMSM',
    'StructuredSexual', 'PriorPartners',
    'MATCHERS',
]
```

- [ ] **Step 4: Confirm matchers module is reachable**

```bash
cd /Users/robynstuart/gf/stisim && python -c "
from stisim.networks import matchers
print(sorted(matchers.MATCHERS.keys()))
"
```

Expected:

```
['band_match', 'desired_age_bucket', 'greedy_old_enough', 'kdtree_nn', 'lsa', 'sort_bisect', 'sort_pair']
```

- [ ] **Step 5: Commit**

```bash
git add stisim/networks/matchers.py stisim/networks/__init__.py
git commit -m "feat: add networks/matchers.py with seven PFAs and MATCHERS registry"
```

### Task B2: Pin pre-existing behaviour with a regression test

**Files:**
- Create: `stisim/tests/test_match_methods.py`

Lock down the current production behaviour (`sort_bisect`) by capturing the matched pairs at a fixed seed. The test must pass against the current code (still in `MFNetwork.match_pairs` from the legacy path), then continue to pass once we switch the dispatcher to call `matchers.sort_bisect` (Task B3).

- [ ] **Step 1: Create `stisim/tests/test_match_methods.py`**

```python
"""Tests for the match_method API on MFNetwork."""
import hashlib
import numpy as np
import pytest
import starsim as ss
import stisim as sti
from stisim.networks import matchers


def _make_sim(match_method=None, n_agents=2_000, seed=0):
    kw = {} if match_method is None else {'match_method': match_method}
    net = sti.MFNetwork(**kw)
    sim = ss.Sim(n_agents=n_agents, networks=net, diseases=sti.HIV(),
                 start='2000-01-01', stop='2001-01-01', rand_seed=seed)
    sim.init()
    return sim


def _fingerprint(p1, p2):
    """Stable hash over the matched pairs for cross-version comparison."""
    a = np.asarray(p1, dtype=np.int64)
    b = np.asarray(p2, dtype=np.int64)
    raw = a.tobytes() + b'|' + b.tobytes()
    return hashlib.sha256(raw).hexdigest()


def test_sort_bisect_dispatch_matches_legacy():
    """The 'sort_bisect' matcher must produce the same pairs as the legacy
    MFNetwork.match_pairs at the same seed.
    """
    sim_legacy = _make_sim(match_method=None, seed=42)
    net_legacy = sim_legacy.networks[0]
    np.random.seed(0)
    p1_legacy, p2_legacy = net_legacy.match_pairs()

    sim_new = _make_sim(match_method='sort_bisect', seed=42)
    net_new = sim_new.networks[0]
    np.random.seed(0)
    p1_new, p2_new = net_new.match_pairs()

    assert _fingerprint(p1_legacy, p2_legacy) == _fingerprint(p1_new, p2_new)


@pytest.mark.parametrize('method', sorted(matchers.MATCHERS.keys()))
def test_each_matcher_produces_valid_pairs(method):
    """Smoke test: each registered matcher runs and produces valid (p1, p2)."""
    sim = _make_sim(match_method=method, n_agents=500, seed=method.__hash__() % 1000)
    net = sim.networks[0]
    try:
        p1, p2 = net.match_pairs()
    except sti.networks.NoPartnersFound:
        # 'desired_age_bucket' can post-filter to empty at small n; that's fine.
        return
    assert len(p1) == len(p2)
    if len(p1):
        assert sim.people.male[p1].all(), 'p1 must be all male'
        assert sim.people.female[p2].all(), 'p2 must be all female'


def test_callable_match_method_is_invoked():
    """Passing a callable as match_method bypasses the registry."""
    sentinel = {'called': False}

    def custom(net):
        sentinel['called'] = True
        empty = ss.uids(np.array([], dtype=np.int64))
        return empty, empty

    sim = _make_sim(match_method=custom, n_agents=200, seed=0)
    net = sim.networks[0]
    try:
        net.match_pairs()
    except sti.networks.NoPartnersFound:
        pass
    assert sentinel['called'], 'custom matcher must be called'


def test_unknown_string_raises_keyerror():
    """Unknown method strings raise KeyError with a helpful keyset."""
    sim = _make_sim(match_method='no_such_method', n_agents=200, seed=0)
    net = sim.networks[0]
    with pytest.raises(KeyError):
        net.match_pairs()
```

- [ ] **Step 2: Run the regression test against the *current* code (pre-dispatch refactor)**

The `match_method` parameter doesn't exist yet, so most tests will fail. The point of this step is to confirm `test_sort_bisect_dispatch_matches_legacy` correctly produces the legacy fingerprint. Run just that test before adding `MFPars.match_method`:

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_match_methods.py::test_sort_bisect_dispatch_matches_legacy -v 2>&1 | head -30
```

Expected: FAIL with `TypeError` or similar because `MFNetwork(match_method=...)` isn't a valid kwarg yet. That's expected — the test exists to gate Task B3.

### Task B3: Wire `match_method` parameter into `MFNetwork`

**Files:**
- Modify: `stisim/networks/mf.py`

Add `match_method` to `MFPars` defaults; rewrite `MFNetwork.match_pairs` as a dispatcher.

- [ ] **Step 1: Add `match_method` to `MFPars.__init__`**

In `stisim/networks/mf.py`, inside `MFPars.__init__`, add the line just before `self.update(kwargs)`:

```python
        # Pair-formation algorithm: string key in matchers.MATCHERS, or a callable.
        self.match_method = 'kdtree_nn'
```

- [ ] **Step 2: Replace the body of `MFNetwork.match_pairs`**

In `stisim/networks/mf.py`, find `MFNetwork.match_pairs(self)` and replace its entire body with:

```python
    def match_pairs(self):
        """Dispatch to the matcher named by ``pars.match_method``.

        ``match_method`` is either a string key in ``matchers.MATCHERS``
        or a callable ``f(net) -> (p1, p2)``. See ``matchers.py``.
        """
        from .matchers import MATCHERS
        m = self.pars.match_method
        if callable(m):
            return m(self)
        return MATCHERS[m](self)
```

(Local import inside the method to avoid a circular import between `mf.py` and `matchers.py` — `matchers.py` doesn't import `mf` at module load, but if anyone later adds that dependency the local import keeps things robust.)

- [ ] **Step 3: Run the regression test**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/test_match_methods.py -v
```

Expected: all pass. Specifically:
- `test_sort_bisect_dispatch_matches_legacy` confirms `sort_bisect` reproduces the legacy output.
- The parametrized `test_each_matcher_produces_valid_pairs` confirms each of the seven methods runs.
- `test_callable_match_method_is_invoked` confirms the escape hatch.
- `test_unknown_string_raises_keyerror` confirms unknown methods fail cleanly.

- [ ] **Step 4: Confirm default sims now use `kdtree_nn`**

```bash
cd /Users/robynstuart/gf/stisim && python -c "
import stisim as sti
net = sti.MFNetwork()
print('default match_method =', net.pars.match_method)
assert net.pars.match_method == 'kdtree_nn'
print('OK')
"
```

Expected: `default match_method = kdtree_nn` and `OK`.

- [ ] **Step 5: Commit**

```bash
git add stisim/networks/mf.py stisim/tests/test_match_methods.py
git commit -m "feat: match_method parameter on MFNetwork (default kdtree_nn)"
```

### Task B4: Delete `stisim/pfa_variants.py` and its public exports

**Files:**
- Delete: `stisim/pfa_variants.py`
- Modify: `stisim/__init__.py`
- Delete: `stisim/tests/test_pfa_variants.py`

The seven subclasses are no longer needed — the registry covers them. Removing them shrinks the public API surface.

- [ ] **Step 1: Remove the `from .pfa_variants import *` line in `stisim/__init__.py`**

In `stisim/__init__.py`, delete the line:

```python
from .pfa_variants  import *
```

The line above (`from .networks import *`) brings in everything you need.

- [ ] **Step 2: Delete the legacy variants module and test file**

```bash
cd /Users/robynstuart/gf/stisim
git rm stisim/pfa_variants.py tests/test_pfa_variants.py
```

- [ ] **Step 3: Confirm no remaining references to `MFNetwork_*`**

```bash
grep -rn "MFNetwork_\|pfa_variants" stisim/ tests/ 2>/dev/null | grep -v "__pycache__\|\.pyc"
```

Expected: zero matches in `stisim/` and `tests/`. (Matches in `tests/devtests/` are fixed in Task B5.)

- [ ] **Step 4: Run the full test suite**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/ -x -q 2>&1 | tail -20
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add stisim/__init__.py
git commit -m "feat: remove MFNetwork_* subclasses (replaced by match_method)"
```

### Task B5: Update `devtest_pfa_comparison.py` to use string-based selection

**Files:**
- Modify: `stisim/tests/devtests/devtest_pfa_comparison.py`

Replace the `VARIANTS = [...]` list of `(name, cls)` tuples with method strings. Delete the `make_variant_class` shim — `StructuredSexual` accepts `match_method=...` directly through its `MFNetwork` parent.

- [ ] **Step 1: Edit `tests/devtests/devtest_pfa_comparison.py`**

Replace these lines:

```python
VARIANTS = [
    ('SortBisect',       sti.MFNetwork_SortBisect),
    ('SortPair',         sti.MFNetwork_SortPair),
    ('LSA',              sti.MFNetwork_LSA),
    ('DesiredAgeBucket', sti.MFNetwork_DesiredAgeBucket),
    ('GreedyOldEnough',  sti.MFNetwork_GreedyOldEnough),
    ('KDTreeNN',         sti.MFNetwork_KDTreeNN),
    ('BandMatch',        sti.MFNetwork_BandMatch),
]
```

with:

```python
# (display name, match_method string)
VARIANTS = [
    ('SortBisect',       'sort_bisect'),
    ('SortPair',         'sort_pair'),
    ('LSA',              'lsa'),
    ('DesiredAgeBucket', 'desired_age_bucket'),
    ('GreedyOldEnough',  'greedy_old_enough'),
    ('KDTreeNN',         'kdtree_nn'),
    ('BandMatch',        'band_match'),
]
```

In `run_generic`, replace `networks=cls()` with `networks=sti.MFNetwork(match_method=method)`:

```python
    for n in n_agents_list:
        for name, method in VARIANTS:
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
                    networks=sti.MFNetwork(match_method=method),
                    diseases=sti.HIV(),
                    ...
```

Delete the entire `make_variant_class` function — it's no longer needed.

In `run_zimbabwe`, replace the variant-class swap logic with a `match_method` parameter override on the existing `StructuredSexual`:

```python
def run_zimbabwe(n_reps, n_agents=10_000, checkpoint=None):
    """Calibrated benchmark on hivsim.demo('zimbabwe')."""
    import hivsim  # noqa: F401
    results = sc.objdict()
    for name, method in VARIANTS:
        if name == 'LSA':
            continue
        for rep in range(n_reps):
            key = (name, 'zimbabwe', rep)
            print(f'Running {key} ...')
            sim = hivsim.demo('zimbabwe', run=False, plot=False, n_agents=n_agents)
            orig_net = None
            for net in sim.pars['networks']:
                if isinstance(net, sti.StructuredSexual):
                    orig_net = net
                    break
            if orig_net is None:
                raise RuntimeError('No StructuredSexual found in hivsim.demo(zimbabwe) networks')
            orig_net.pars.match_method = method
            sim.pars['rand_seed'] = rep
            sim.pars['analyzers'] = list(sim.pars.get('analyzers') or []) + [
                PartnersLastYearAnalyzer(),
                PairFormationAgesAnalyzer(),
                PairPrevalenceAnalyzer(),
            ]
            with sc.timer() as t:
                sim.run()
            mf_net = next(nw for nw in sim.networks() if isinstance(nw, sti.MFNetwork))
            results[key] = sc.objdict(
                wall_time=t.elapsed,
                hiv_prevalence=sim.results.hiv.prevalence.values.copy(),
                partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                pair_formation_ages=sim.analyzers.pairformationagesanalyzer.records,
                pair_prevalence=sim.analyzers.pairprevalenceanalyzer.records,
                lifetime_partners=mf_net.lifetime_partners.values.copy(),
                sex_male=sim.people.male.values.copy(),
                age=sim.people.age.values.copy(),
            )
            print(f'    wall_time={t.elapsed:.2f}s')
            if checkpoint is not None:
                checkpoint(results)
    return results
```

- [ ] **Step 2: Run the devtest in quick mode to confirm it imports and runs end-to-end**

```bash
cd /Users/robynstuart/gf/stisim/tests/devtests && python devtest_pfa_comparison.py --quick 2>&1 | tail -20
```

Expected: it runs through all seven variants on 1k agents with 1 rep and produces a `pfa_comparison_results.obj` file. Wall times print without error.

- [ ] **Step 3: Commit**

```bash
git add tests/devtests/devtest_pfa_comparison.py
git commit -m "devtest: PFA comparison uses match_method string API"
```

### Task B6: Update the notebook builder's algorithm-description references

**Files:**
- Modify: `stisim/tests/devtests/_build_pfa_comparison_notebook.py`

The notebook still calls out `MFNetwork_*` subclasses in its prose. Update the references so it points at the new API.

- [ ] **Step 1: Find and update references**

```bash
cd /Users/robynstuart/gf/stisim && grep -n "MFNetwork_" tests/devtests/_build_pfa_comparison_notebook.py
```

For each occurrence (currently at lines 180 and 198), edit the markdown text so it says "the seven `match_method` choices" or similar — drop the `MFNetwork_*` phrasing.

Edit line 180:
- Old: `The generic block uses bare \`MFNetwork_*\` variants which have no SW layer, so`
- New: `The generic block uses bare \`MFNetwork\` (no SW layer), so`

Edit line 198:
- Old: `The generic block runs with bare \`MFNetwork_*\` variants. Zimbabwe uses the`
- New: `The generic block runs with bare \`MFNetwork\` (one method per run). Zimbabwe uses the`

- [ ] **Step 2: Rebuild the notebook from the script (to update the markdown cells)**

```bash
cd /Users/robynstuart/gf/stisim/tests/devtests && python _build_pfa_comparison_notebook.py
```

Expected: prints `Wrote .../pfa_comparison.ipynb`. The notebook is regenerated with updated prose; we do NOT re-execute it as part of this commit (the existing results file is fine).

- [ ] **Step 3: Commit**

```bash
git add tests/devtests/_build_pfa_comparison_notebook.py tests/devtests/pfa_comparison.ipynb
git commit -m "docs: notebook references match_method API instead of MFNetwork_* subclasses"
```

### Task B7: Final full-suite sanity check

**Files:** none

- [ ] **Step 1: Run the entire test suite**

```bash
cd /Users/robynstuart/gf/stisim && pytest tests/ -q 2>&1 | tail -10
```

Expected: all pass. Tally the test count and confirm `test_match_methods.py` contributed 11 new test cases (1 + 7 parametrized + 1 callable + 1 keyerror + 1 dispatch regression = 11; actually 4 functions, one parametrized over 7 → 10 reports).

- [ ] **Step 2: Confirm public API smoke**

```bash
cd /Users/robynstuart/gf/stisim && python -c "
import stisim as sti
# Symbols that must still exist after the refactor
for name in ['MFNetwork', 'SWNetwork', 'StructuredSexual', 'PriorPartners',
             'AgeMatchedMSM', 'AgeApproxMSM', 'NoPartnersFound', 'BaseNetwork',
             'BasePars', 'MFPars', 'SWPars', 'NetworkPars']:
    assert hasattr(sti, name), f'missing: {name}'

# Symbols that must NOT exist any more
for name in ['MFNetwork_LSA', 'MFNetwork_SortPair', 'MFNetwork_KDTreeNN',
             'MFNetwork_BandMatch', 'MFNetwork_DesiredAgeBucket',
             'MFNetwork_SortBisect', 'MFNetwork_GreedyOldEnough']:
    assert not hasattr(sti, name), f'should be removed: {name}'

# New API works
net = sti.MFNetwork()
assert net.pars.match_method == 'kdtree_nn'
print('OK')
"
```

Expected: prints `OK`.

- [ ] **Step 3: Show the final commit log on the branch**

```bash
git log --oneline rc1.5.6..HEAD | head -20
```

Confirm the commits land in the expected order: six refactor commits (A2–A7), then four feature commits (B1, B3, B4, B5, B6 — Task B2 doesn't commit on its own; its test ships with B3).

---

## Out of scope (followups, not in this PR)

These are noted in the spec and the project memory; do NOT do them in this PR:

- Recalibrating Zimbabwe / Eswatini / Kenya / Zambia models against the new default.
- Eswatini's `_build_pfa_comparison_notebook.py` update (separate repo).
- Resolving the duplicated relationship-acceptance Bernoulli in `desired_age_bucket` vs `MFNetwork.add_pairs`.
- Renaming `SWNetwork` → `FSWNetwork`.
- CHANGELOG entry (per the project memory `feedback_changelog_timing.md`, that happens at rc merge, only if asked).
