# Calibration knobs — the levers STIsim actually exposes

Organized by what the analyst is trying to change, not by module. For
each use case, this document lists the parameters and data levers that
actually move the outcome, notes which internal pathway each acts on
(this matters — mis-targeting the pathway is a common calibration
failure mode), and points at the canonical source for the parameter's
definition.

Not a parameter dictionary. Deep semantics live in the source; use
`canonical-sources.md` to navigate.

---

## Framing — read this first

**STIsim's default parameter values are best-informed placeholders, not
values you are expected to accept as-is.** They are set from the
general literature so the model runs sensibly out of the box, but they
are almost never the values that a specific analysis should use. In
particular:

- **Discovering that a default doesn't fit your setting is not a bug.**
  It means you have found a parameter that varies across settings and
  that you need to set for your setting. That is the normal case, not
  a defect in STIsim.
- **Do not "calibrate" a parameter that should come from data.** Many
  parameters — especially anything network- or behavior-related (age
  of sexual debut, partnership durations, sex-work durations, risk-
  group proportions, condom use) — vary by country. These should be
  filled in from local surveys (DHS, PHIA, IBBS, etc.), not treated as
  free parameters for a calibrator to tune. Feeding them from data
  keeps the calibrator focused on the parameters that genuinely lack
  observational anchoring.
- **If a genuine global default looks wrong, fix it upstream.** Do not
  patch it in your local checkout. Other people use this tool. When
  you find a default that misses a known epidemiological signal, open
  a PR — the whole network of expert users benefits from the
  correction, and your own analysis benefits from having the fix in
  the canonical release rather than in a private branch.

Two different modes of intervention, then:

| Mode | What you're doing | Example |
|---|---|---|
| **Feed data** | Replace a default with a local value from a survey or data source. Not a calibration. | Age of sexual debut for your country from DHS. |
| **Calibrate** | Tune a parameter that genuinely lacks observational anchoring, to match a target output. | HIV `beta_m2f`, tuned against prevalence trajectory. |
| **Upstream fix** | Correct an STIsim default that misses a known signal. Belongs in a PR, not your workspace. | A default value contradicted by a well-established literature source. |

The sections below use these categories to flag which levers belong in
which mode.

---

## HIV mortality

There are three distinct classes of lever, and calibrating mortality
usually means using more than one.

### 1. Deduplicate the input data first (data step, not a calibration)

The all-cause background mortality rates STIsim consumes (via
`ss.Deaths`, from files like `<country>_deaths.csv`) include AIDS
deaths in high-HIV-burden countries. The HIV module *also* kills
agents via `p_hiv_death` and the `ti_zero` pathway. That is a
double-count and inflates apparent HIV mortality regardless of what
knobs are turned.

Fix: pre-process the deaths file with `stisim.data.dedup_deaths`
before running. This removes the AIDS component from the input rates
by interpolating between a pre-epidemic and post-epidemic anchor.

```python
import pandas as pd
from stisim.data import dedup_deaths, deleted_fraction

all_cause = pd.read_csv('data/<country>_deaths.csv')
hiv_deleted = dedup_deaths(all_cause, base_year=1985, end_year=2025)

all_cause.to_csv('data/<country>_deaths_all_cause.csv', index=False)
hiv_deleted.to_csv('data/<country>_deaths.csv', index=False)

# Sanity check without running the sim: AIDS share should peak in
# the mid-30s (~70–85%) and fall to ~0 at 80+.
diag = deleted_fraction(all_cause, hiv_deleted)
```

Do this **before** attempting to calibrate any HIV mortality parameter.
Otherwise the calibrator will be fighting a data-side overcount using
model-side knobs, and no parameter setting will reconcile them.

### 2. The rate knobs — act on `p_hiv_death` (calibration)

These scale the per-timestep CD4-stratified death hazards.

| Parameter | Effect | Notes |
|---|---|---|
| `rel_death` | Multiplies all HIV death probabilities (off- and on-ART) | The general "scale up or down HIV mortality" knob. |
| `rel_death_f` | Additional female multiplier | Folded into the off-ART rate; applies to both off- and on-ART pathways. |
| `art_death_age` | List of `(age_lo, age_hi, mult)` tuples for on-ART deaths by age | Analogous to `rel_sus_age`. Use to shape on-ART mortality's age gradient. |
| `rel_art_mortality_unsupp_m/f` | Sex-specific multiplier on unsuppressed on-ART deaths | Only affects the fraction not virally suppressed. |

Canonical source: `stisim/diseases/hiv.py`.

### 3. The duration knobs — act on `ti_zero` (calibration)

These do not appear on any rate table. They determine *when* the
`ti_zero` AIDS death fires, which in the current model accounts for
the majority of HIV deaths.

| Parameter | Effect | Notes |
|---|---|---|
| `dur_falling` | Duration of late-stage HIV (falling CD4) | Determines how much time between `ti_falling` and `ti_zero`. Shorter = earlier deaths. |
| `dur_post_art` | Post-ART-interruption survival distribution | Sets the mean / scale for the delay to death after an ART interruption. |
| `dur_post_art_scale_factor` | Variance parameter for `dur_post_art` | Controls the spread around the mean. |

If your AIDS-mortality deficit persists after tuning the rate knobs in
(2), check these. The rate knobs cannot fix a deficit whose source is
the `ti_zero` pathway.

Canonical source: `stisim/diseases/hiv.py`.

---

## HIV transmission

**All betas are constant over time.** There is no time-varying beta.
When time-varying transmission is needed (secular declines, epidemic
waves), it comes from adjusting *inputs that change over time* — not
from making beta a function of year.

The two levers that give time-varying effective transmission:

- **Condom use** — provided as a time-series input file (e.g.
  `condom_use.csv`). Fill this from DHS or similar; changes in condom
  use over time flow through `eff_condom` to change effective
  transmission.
- **ART coverage** — as coverage rises, viral suppression rises and
  effective onward transmission falls. Provide the ART scale-up
  schedule to the ART intervention; the intervention drives the
  time-varying transmission indirectly.

The static beta parameters (calibration knobs):

| Parameter | Effect |
|---|---|
| `beta_m2f` | Male-to-female per-act transmission probability |
| `rel_beta_f2m` | Female-to-male, relative to `beta_m2f` |
| `beta_m2c` | Mother-to-child (vertical) |
| `beta_breastfeed` | Breastfeeding transmission |
| `beta_m2m` | MSM transmission |
| `eff_condom` | Condom effectiveness |
| `rel_init_prev` | Scales the loaded initial prevalence |

Canonical source: `stisim/diseases/sti.py` (`BaseSTIPars`) and
`stisim/diseases/hiv.py` for HIV-specific overrides.

---

## Sexual networks (mostly data-driven, not calibration)

Network parameters are the largest cluster of "should come from data,
not the calibrator" levers in STIsim. They vary substantially by
country and by population, and treating them as calibration knobs
usually means fitting on a parameter that has an observed answer in
the local literature.

### Age of sexual debut

Base parameters (`stisim/networks/base.py`):

| Parameter | Default | Notes |
|---|---|---|
| `debut_f` | `ss.lognorm_ex(20, 3)` | Female age at sexual debut. Country-specific; fill from DHS. |
| `debut_m` | `ss.lognorm_ex(21, 3)` | Male age at sexual debut. Country-specific; fill from DHS. |
| `acts` | `ss.lognorm_ex(freqperyear(80), freqperyear(30))` | Annual coital frequency per partnership. Country-specific; fill from behavioral surveys where available. |

### Risk-group composition and partnership dynamics

MF-network parameters (`stisim/networks/mf.py`):

| Parameter | Default | Notes |
|---|---|---|
| `prop_f0` | 0.85 | Proportion of women in the low-risk group. Country-specific; DHS partnership counts inform. |
| `prop_m0` | 0.80 | Proportion of men in the low-risk group. |
| `prop_f2` | 0.01 | Proportion of women in the high-risk group. |
| `prop_m2` | 0.02 | Proportion of men in the high-risk group. |
| `concurrency_dist` | `ss.poisson(lam=1)` | Concurrent-partnership distribution. Country- and population-specific. |
| `dur_dist` | `ss.lognorm_ex()` | Partnership duration distribution. |
| `p_matched_stable`, `p_mismatched_casual` | see source | Assortative-mixing probabilities by risk group. |

### Sex-work-specific parameters

FSW-network parameters (`stisim/networks/fsw.py`):

| Parameter | Default | Notes |
|---|---|---|
| `dur_sw` | `ss.lognorm_ex(mean=5, std=3)` years | Duration in sex work. Highly country-specific; fill from IBBS. |
| `dur_client` | `ss.lognorm_ex(mean=10, std=5)` years | Duration as a client of sex work. |
| `age_sw_start` | derived | Age of entry into sex work. |

Canonical source: `stisim/networks/base.py`, `mf.py`, `fsw.py`, `msm.py`.

**How to use these.** For a country study, walk the network parameters
one by one and ask: *do I have local data for this?* If yes, feed it
in and take it out of the calibration set. If no, note it as a
calibration parameter and be explicit in the analysis-brief about
which network parameters lack observational anchoring.

---

## Care-seeking

STIsim's practical convention is to parameterize care-seeking on the
disease modules and testing interventions directly, not via the
standalone `CareSeeking` module.

**This is a natural data-entry point.** Care-seeking probabilities and
testing coverage over time are commonly observed (DHS "sought care
for STI symptoms" tables, ART coverage from Spectrum, HTS testing
volumes from national program data). Feed observed values in
rather than calibrating them from prevalence.

| Lever | Where |
|---|---|
| `p_symp_care` on SEIS diseases | Sets the probability that a symptomatic agent seeks care. Fill from DHS or national STI program data. |
| `rel_test` on testing interventions | Multiplies testing rates for scale-up scenarios; typically driven from a time-series of testing coverage. |
| Project-level multipliers (e.g. `care_seek_mult`) | Applied above the disease parameters when needed. |

Canonical source: `stisim/diseases/sti.py`, `stisim/interventions/*.py`.

---

## Demographics

Not typically calibration territory. The demographic inputs (age
structure, ASFR, deaths, migration) are data — pulled from UN WPP
via `stisim.data.downloaders` or downloaded manually. Once the input
files are correct (including HIV-deletion for high-burden settings),
demographic behavior usually needs only light tuning.

The main lever:

| Lever | Where |
|---|---|
| `dem_pars={'rel_migration': ...}` | Location-level scaling on net migration to match a specific dataset (e.g. `rel_migration=0.5` in the Zimbabwe project to match UN WPP after preserving the age structure). |

Canonical source: `stisim/demographics.py`, `stisim/data/downloaders.py`.

---

## How to extend this document

Add a section per new use case. Each section should:

1. Name the outcome the analyst is trying to change.
2. List the levers that actually move it, grouped by pathway or role.
3. **Classify each lever** as data-driven (fill from local sources),
   calibration (tune to fit output), or upstream-fix (if the default
   is wrong).
4. Flag any data-side pre-processing needed (like `dedup_deaths` for
   mortality).
5. Point at canonical sources rather than paraphrase them.

The value of this document is knowing **which knob matches which
role** — the risk of listing them out of context is that an analyst
calibrates the wrong thing and concludes the model is broken.
