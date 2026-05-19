# `match_method` API + `networks/` package layout

**Date**: 2026-05-18
**Branch**: `feat/pfa-comparison` (off `rc1.5.6`)
**Predecessor**: [2026-05-15-pfa-comparison-design.md](./2026-05-15-pfa-comparison-design.md) and [tests/devtests/PFA_RESULTS.md](../../../tests/devtests/PFA_RESULTS.md)

## Goal

Expose the seven pair-formation algorithm (PFA) variants benchmarked on this
branch as a clean user-facing choice on `MFNetwork`, and reorganise the
network module so the new entry point is discoverable.

Two coupled changes:

1. **`networks/` package** — split `stisim/networks.py` (~1000 LoC, 9
   classes) into a sub-package mirroring `diseases/` and `interventions/`.
2. **`match_method` parameter** — replace the seven `MFNetwork_*`
   subclasses (currently in `stisim/pfa_variants.py`) with a string
   parameter on `MFNetwork`, dispatched through a registry of matcher
   functions. Default changes from current production (`sort_bisect`) to
   `kdtree_nn`.

## Motivation

`pfa_variants.py` makes PFA selection feel like a power-user feature:
users must import a different `MFNetwork_*` class. With the
`match_method` parameter, picking a PFA is the same shape as picking any
other modelling option (a string in `pars`).

The PFA comparison work showed the current production (`sort_bisect`)
does not honour `age_diff_pars` — it collapses the mean M-F age gap to
~0 instead of the intended +7yr. `kdtree_nn` honours `age_diff_pars`
with the cleanest implementation among the candidates that do. Making
it the default fixes a real bug in the calibrated behaviour of every
downstream STIsim model, at the cost of needing to recalibrate them.

The `networks/` reorganisation is opportunistic: we're already touching
this code, and `networks.py` has grown beyond comfortable single-file
size. The new layout makes `matchers.py` a natural sibling rather than
a top-level oddity.

## Non-goals

- Renaming `SWNetwork` → `FSWNetwork` (or anything else). The file moves
  into `fsw.py` for discoverability; the class name stays.
- Resolving `DesiredAgeBucket`'s duplicated relationship-acceptance
  logic with `MFNetwork.add_pairs`. Tracked as a follow-up; doesn't bite
  while `kdtree_nn` is the default.
- Changing `SWNetwork.match_pairs` or the MSM subclasses' `match_pairs`.
  The registry is MF-only.
- Recalibrating any downstream model. That happens after this lands, in
  the respective localisation repos.

## Package layout

```
stisim/networks/
├── __init__.py            # re-exports — public API stays flat
├── README.md              # describes the layout, especially layered networks
├── base.py                # NoPartnersFound, BasePars, NetworkPars, BaseNetwork
├── mf.py                  # MFPars, _mf_states, MFNetwork
├── fsw.py                 # SWPars, _sw_states, SWNetwork
├── msm.py                 # AgeMatchedMSM, AgeApproxMSM
├── layered_networks.py    # StructuredSexual, PriorPartners
└── matchers.py            # PFA functions + MATCHERS registry
```

**Source mapping from current `networks.py`:**

| Current location           | New location           |
|----------------------------|------------------------|
| `NoPartnersFound`          | `base.py`              |
| `_mf_states`               | `mf.py`                |
| `_sw_states`               | `fsw.py`               |
| `BasePars`                 | `base.py`              |
| `MFPars`                   | `mf.py`                |
| `SWPars`                   | `fsw.py`               |
| `NetworkPars`              | `base.py`              |
| `BaseNetwork`              | `base.py`              |
| `MFNetwork`                | `mf.py`                |
| `SWNetwork`                | `fsw.py`               |
| `StructuredSexual`         | `layered_networks.py`  |
| `PriorPartners`            | `layered_networks.py`  |
| `AgeMatchedMSM`            | `msm.py`               |
| `AgeApproxMSM`             | `msm.py`               |

**`pfa_variants.py` is deleted.** Its content moves to `matchers.py` as
module-level functions; the seven `MFNetwork_*` subclasses are gone.

**Public API stays flat.** `stisim/__init__.py` keeps `from .networks
import *`. So `sti.MFNetwork`, `sti.SWNetwork`, `sti.StructuredSexual`,
`sti.PriorPartners`, `sti.NoPartnersFound`, `sti.BaseNetwork` all keep
working unchanged.

`stisim/networks/__init__.py` re-exports the classes for the flat API
*and* exposes `matchers` as an attribute via `from . import matchers`,
so power users can reach matcher functions directly:
`sti.networks.matchers.kdtree_nn` or `from stisim.networks import matchers`.

## `match_method` API

### Matcher signature

Each matcher is a module-level function in `matchers.py`:

```python
def kdtree_nn(net):
    """Match each woman to her nearest-age available man via 1-NN."""
    ppl = net.sim.people
    f_looking, m_eligible = net._get_eligible()
    desired_ages = net._sample_desired_ages(f_looking)
    ...
    return p1, p2  # both ss.uids
```

The `net` parameter is the `MFNetwork` instance. Matchers use the
existing helpers (`_get_eligible`, `_sample_desired_ages`) which already
live on `MFNetwork`. They raise `NoPartnersFound` to signal an empty match.

### Registry

```python
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

### Parameter

Added to `MFPars` defaults:

```python
self.match_method = 'kdtree_nn'
```

### Dispatch

`MFNetwork.match_pairs` becomes a thin dispatcher:

```python
def match_pairs(self):
    m = self.pars.match_method
    if callable(m):
        return m(self)
    return MATCHERS[m](self)
```

String for the standard set; callable as the escape hatch for
user-defined matchers (no subclassing required). Unknown strings raise
`KeyError` at first call — naturally informative because the dict shows
the valid keys in the error.

### Method-specific constants

`band_match` uses a 5-year band width — currently the `band_width = 5`
class attribute on `MFNetwork_BandMatch`. After the refactor this lives
as a module-level constant in `matchers.py`:

```python
BAND_MATCH_WIDTH = 5
```

The matcher function reads it directly; no tuning knob exposed to
users. The `make_variant_class` shim in
`devtest_pfa_comparison.py` (which currently copies the `band_width`
class attribute when building `StructuredSexual_BandMatch`) goes away
entirely — devtest moves to string-based selection.

If someone needs a tunable band width later we lift it to a
`match_kwargs` dict on `MFPars`.

## Default change

`match_method = 'kdtree_nn'` is a behaviour change from current
production. Downstream calibrated models will produce different HIV/STI
trajectories:

- HIV prevalence shifts by ~18% in Zimbabwe, ~59% in Eswatini (from
  [PFA_RESULTS.md](../../../tests/devtests/PFA_RESULTS.md)).
- Stisim models pulling the next stisim release after this PR ships must
  recalibrate against `kdtree_nn` (or pin to `match_method='sort_bisect'`
  if they want to keep current behaviour).

Documented in:

- The CHANGELOG entry for the rc that includes this PR (touched at
  rc-merge time per the existing convention, not in this PR).
- The `MFNetwork` docstring lists the seven options.
- The `networks/README.md`.

## Internal call-sites that must update

After the refactor + feature, the following files reference the old
names and must move to the new API:

| File                                              | Change                                                                |
|---------------------------------------------------|-----------------------------------------------------------------------|
| `stisim/__init__.py`                              | Drop `from .pfa_variants import *`. `from .networks import *` stays. |
| `tests/test_pfa_variants.py`                      | Rename → `tests/test_match_methods.py`. Parametrize over `MATCHERS.keys()`. |
| `tests/devtests/devtest_pfa_comparison.py`        | Use `match_method='...'` strings. Drop `make_variant_class` shim.   |
| `tests/devtests/_build_pfa_comparison_notebook.py`| Update algorithm-description table to reference the new API.        |

External (not in this PR, but flagged for the user):

| File                                              | Change                                                              |
|---------------------------------------------------|---------------------------------------------------------------------|
| `hivsim_eswatini/run_pfa_comparison.py`           | Use `match_method='...'` strings.                                  |
| `hivsim_eswatini/_build_pfa_comparison_notebook.py`| Update algorithm-description table.                                |

The other localisations (Zimbabwe, Kenya, Zambia, syph_dx_zim, etc.)
don't reference `MFNetwork_*` directly — they just use `StructuredSexual`,
which gets the new default automatically.

## Testing

### `tests/test_match_methods.py` (new, replaces `test_pfa_variants.py`)

Parametrized over `matchers.MATCHERS.keys()`. For each method:

- 5-year sim, n=2k agents, `HIV` disease (per
  [feedback_diseases_for_mfnetwork](../../../../../.claude/projects/-Users-robynstuart-gf-stisim/memory/feedback_diseases_for_mfnetwork.md)
  — don't use bare SIS with MFNetwork).
- Assert edges are formed (non-empty `network.edges.p1` at sim end).
- Assert no error raised.

Plus two narrow tests:

- **Dispatch regression**: `match_method='sort_bisect'` reproduces the
  legacy `MFNetwork.match_pairs` output bit-for-bit on a fixed-seed run.
  Confirms the refactor didn't drift.
- **Callable escape hatch**: passing a callable `match_method` works
  and is invoked. Use a no-op matcher returning empty `ss.uids`.

### Existing tests

- `tests/devtests/devtest_pfa_comparison.py` continues to exercise all
  seven methods end-to-end via the new string API. No assertion changes
  needed.

## Sequencing

Two commits on `feat/pfa-comparison`:

1. **`refactor: split networks.py into networks/ package`** — pure file
   moves. `git diff -M` should show renames + a thin `__init__.py`.
   CI green after this commit. No behaviour change.
2. **`feat: match_method parameter on MFNetwork`** — adds
   `matchers.py`, `MFPars.match_method = 'kdtree_nn'`, dispatcher in
   `MFNetwork.match_pairs`. Deletes `pfa_variants.py` and the seven
   subclasses. Updates internal tests + devtest.

Each commit is independently revertable; the second is the
behaviour-changing one.

## Risks & open questions

- **CI uses calibrated downstream sims?** If any internal CI hits
  Zimbabwe/Eswatini calibrated values, those tests will fail when the
  default changes. Need to audit `tests/` for that pattern before the
  second commit.
- **`NetworkPars` consumer.** `sti.Sim` uses `NetworkPars` to filter
  kwargs into network buckets. Confirm the import path still resolves
  after the refactor (`from stisim.networks import NetworkPars`).
- **`from stisim.networks import X` users**. External code that does
  `from stisim.networks import MFNetwork` keeps working via the package
  `__init__.py`. Code that does `import stisim.networks as nets;
  nets.MFNetwork` also keeps working. Code that does `from
  stisim.networks.networks import MFNetwork` (unlikely but possible
  given the file used to be `networks.py`) would break — worth a quick
  grep across the wider repo set.
