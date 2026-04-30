# SW Participation Modeling — Design Bake-Off

Issue: [#307](https://github.com/starsimhub/stisim/issues/307)

---

## Slide 1 — Current overview: how transactional sex is modeled today

**State** ([networks.py:181-182](../stisim/networks.py#L181-L182)): two `BoolArr` on `StructuredSexual` — `fsw`, `client`.

**Assignment** ([networks.py:329-333](../stisim/networks.py#L329-L333)) via `set_sex_work`:
- `fsw_shares = bernoulli(p=0.05)` — flipped once for each female at intake.
- `client_shares = bernoulli(p=0.12)` — flipped once for each male at intake.
- **Once True, never flips back.**

**Matching** ([networks.py:542-586](../stisim/networks.py#L542-L586)):
- `active_fsw = over_debut & fsw` (filters by sexual debut only).
- `sw_intensity ~ Uniform(0,1)` sampled for each active FSW each step.
- Clients with `sw_seeking_rate` draw FSW partners, weighted by intensity.
- Edges typed as `sw`, condom use separately parameterized.

**Downstream** (~14 sites): PrEP eligibility, syphilis/STI testing, HIV stratified results, analyzers — all read `fsw` / `client` directly.

---

## Slide 2 — The problem

**Headline: FSW prevalence drifts over a sim.**
- `fsw_shares=0.05` is applied at **intake** only.
- With pop growth + no exits + differential aging, the FSW share of the adult population grows unboundedly.
- Example: a 40-year sim starting at 5% FSW can end >10% just from accumulation — and we have no lever to stop it short of lowering `fsw_shares` (breaks initial calibration).

**Secondary problems**
- "Lifelong" vs "currently" conflated across ~14 call sites. A 55-yr-old who sold sex at 22 is still PrEP-eligible.
- No re-entry, no exit programs, no response to economic shocks.
- Can only calibrate one of (lifetime prev, point prev) — not both.
- `sw_intensity` is the only dial for activity, and it's sampled fresh each step — no persistence.

*Supporting evidence*: `syph_dx_zim/archive/fsw_stats.py` dumps FSW % of all females vs 15-49, plus FSW age distribution at end-of-sim — shows the drift and the aged-out-but-still-flagged agents.

---

## Slide 3 — Planned architecture: MF + SW network split

Regardless of which SW mechanism is chosen, `StructuredSexual` will be refactored into two composable networks:

**`MFNetwork`** — heterosexual partnerships among 3 risk groups (low / medium / high). Contains all non-SW matching logic (`add_pairs_nonsw`, risk group assignment, duration/concurrency pars). No `fsw`/`client` state.

**`SWNetwork`** — sex work partnerships. Contains `fsw`, `client`, `sw_intensity` states and all SW matching logic (`set_sex_work`, `match_sex_workers`, `add_pairs_sw`). Standalone — no dependency on `MFNetwork.risk_group` (confirmed: `match_sex_workers` does not use risk group at all).

**`StructuredSexual`** — thin wrapper that instantiates both. Existing code continues to work unchanged.

**Two valid modeling modes:**
1. `MFNetwork` only — model SW implicitly via risk group 2 parameterization (see Option 3 below).
2. `MFNetwork` + `SWNetwork` — model SW explicitly as a separate network with its own matching logic (Options 1 or 2 below).

This split is orthogonal to the mechanism options below: whichever mechanism is chosen, it lives inside `SWNetwork`.

---

## Slide 5 — Option 1: Lifetime fate + time window

**Mechanism**
- At birth, draw:
  - `ever_fsw` (bool)
  - `age_sw_start` (float, years)
  - `dur_sw` (float, years)
- Computed: `current_fsw(t) = ever_fsw & (age_sw_start ≤ age < age_sw_start + dur_sw)`.
- `sw_intensity(t)` modulated within window (e.g., lower during pregnancy).

**Pros**
- Eliminates the lifelong-accumulation drift.
- `ever_fsw` gives a clean lifetime prevalence target; window width gives a point prevalence lever.
- Simple, deterministic, low implementation cost.

**Cons**
- ⚠️ **Residual drift**: `ever_fsw` rate is a lifetime-incidence parameter, not point prevalence. In a growing/aging population the age pyramid shifts, so `current_fsw %` still moves over the sim (smaller than today's drift, but real). Mitigation: calibrate `p_ever` per birth cohort, or accept drift as demographically driven.
- One contiguous window per person → no re-entry.
- Cycling only via intensity dial — intensity≈0 ≠ "not currently SW" semantically.
- Doesn't naturally respond to time-varying drivers (shocks, policy changes).

---

## Slide 6 — Option 2: Continuous propensity

**Mechanism**
- At birth: `sw_propensity ~ Bimodal(…)` (most near 0, long tail).
- Each step, determine `current_fsw` via either:
  - **Force-fit**: rank by `sw_propensity × modifier(t)`, take top-N to hit target point prevalence (same pattern as ART targeting).
  - **Model-based**: `p_current(t) = f(sw_propensity, age, pregnancy, partnered, …)`, sampled per step.

**Pros**
- ✅ **Force-fit variant fully eliminates drift** — point prevalence is pinned each step by construction.
- Continuous heterogeneity; natural cycling (in/out across steps).
- Easy to add time-varying drivers via `modifier(t)`.
- Logistic variant extensible and familiar.

**Cons**
- Without smoothing, FSW identity churns step-to-step → partnership/intervention continuity issues.
- More parameters; propensity shape needs data.
- "Ever FSW" is derived (ever selected), not set.

---

## Slide 7 — Option 3: Implicit behavioral definition (MF-only mode)

**Mechanism**
- No `SWNetwork` instantiated. Risk group 2 females are the effective "high-risk" population.
- SW-like behavior is defined by parameters: high concurrency, short relationship duration, high partner turnover.
- Optionally, add a derived property `is_keypop` based on observable behavioral thresholds (e.g., `partners_12 > N` or `concurrency > C`) rather than a categorical label.
- Downstream consumers (PrEP eligibility, testing) key off `is_keypop` or `risk_group == 2` rather than `fsw`.

**Pros**
- No separate network or matching logic — zero added complexity.
- Drift is a non-issue: risk group assignment is stable, and behavioral properties are naturally bounded by the MF network's dynamics.
- Natural fit for settings where the SW/non-SW distinction is analytically unimportant or poorly measured.
- `is_keypop` threshold can be tuned to match point prevalence of "high activity" without requiring a separate prevalence calibration target.

**Cons**
- No explicit SW edges — can't separately parameterize condom use or act frequency for transactional sex.
- No clean "ever FSW" lifetime prevalence target.
- Loses the semantic signal that a specific transmission route (transactional sex) is being modeled.
- Less appropriate when SW-specific interventions (targeted PrEP, outreach programs) need to be modeled distinctly.

---

## Slide 8 — Comparison

| Criterion | Opt 1: Window | Opt 2: Propensity | Opt 3: Implicit |
|---|---|---|---|
| Fixes prevalence drift | partial — eliminates accumulation drift, inherits demographic drift | ✅ full (force-fit) | ✅ not applicable |
| Point prevalence calibration | medium | **easy** (force-fit) | medium (via `is_keypop` threshold) |
| Lifetime prevalence calibration | **easy** | derived | not applicable |
| Cycling / re-entry | poor | **good** | natural (behavioral) |
| Separate SW transmission route | ✅ yes | ✅ yes | ✗ no |
| SW-specific interventions | ✅ yes | ✅ yes | partial (via `is_keypop`) |
| Intervention flexibility | medium | **high** | low |
| Implementation effort | **low** | medium | **lowest** |
| Semantic clarity | **high** | high | low |
| Requires SWNetwork | yes | yes | **no** |

**Open questions**
1. Data on SW duration distributions? (Defaults matter for Opt 1.)
2. Do we need re-entry in v1 or is one window enough?
3. Role of `sw_intensity` — modulation within window (Opt 1), or step-level weight (Opt 2)?
4. Calibration targets: point prevalence, lifetime prevalence, or both?
5. Is Opt 3 a first-class supported mode, or just "use MF and tune risk group 2 yourself"?
