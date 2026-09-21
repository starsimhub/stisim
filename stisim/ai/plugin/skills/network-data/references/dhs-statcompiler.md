# DHS Statcompiler — age at first sexual intercourse

Country-specific age-of-debut data for the `debut_f` and `debut_m`
network parameters. Source: DHS Statcompiler
(statcompiler.com), which is the browsable frontend for the DHS
Program indicator database.

## What we need from the data

STIsim's default age of debut is:

```python
debut_f = ss.lognorm_ex(20, 3)   # female mean 20, std 3
debut_m = ss.lognorm_ex(21, 3)   # male mean 21, std 3
```

For a country study we want:
- **Median age at first sexual intercourse** for women (20–49 aggregate)
  and for men (20–49 aggregate). Used to set the mean of the lognormal.
- **Cumulative first-sex fractions by exact age** (15, 18, 20, 22, 25).
  Used to characterise the spread — from the cumulative curve the
  standard deviation of the underlying age-of-debut distribution can
  be estimated, replacing stisim's default `std=3`.

For a lighter touch, the median alone gives a decent mean estimate;
default `std=3` can be kept until a spread-informed value earns its
place.

## Manual workflow (via the Statcompiler web UI)

1. Go to https://statcompiler.com/.
2. Click **Choose Indicator**.
3. Select **Complete list**.
4. Navigate to **Sexual behavior**.
5. Select:
   - *Median age at first sexual intercourse [Women]: 20–49*
   - *Median age at first sexual intercourse [Men]: 20–49*
6. Select the country of interest.
7. Record the values (per survey year); use the most recent survey (or
   average across recent surveys if the analysis is over a stretch of
   years).

For the spread-informed parameterisation, also retrieve:

- *First sexual intercourse by exact age 15 [Young women]*
- *First sexual intercourse by exact age 18 [Young women]*
- *First sexual intercourse by exact age 20 [Young women]*

(and 22, 25 if available — the "Young women" cohort is 15–24, so
values above ~24 are censored). These cumulative-by-exact-age values
give a shape estimate.

## Programmatic workflow (via the DHS API)

DHS operates a REST API at `https://api.dhsprogram.com/` that returns
the same indicator data Statcompiler surfaces. No authentication is
required for basic use; registering as a partner
(email api@dhsprogram.com) enlarges the per-page result cap.

**Base URL:** `https://api.dhsprogram.com/`
**Data endpoint:** `/rest/dhs/data`
**Indicator catalogue endpoint:** `/rest/dhs/indicators`
**Formats:** `f=json`, `f=xml`, `f=csv`, `f=html`

### Indicator IDs

The relevant indicator IDs for age-of-debut work:

| IndicatorId | Meaning |
|---|---|
| `SX_AAFS_W_M2A` | Median age at first sexual intercourse [Women]: 20–49 |
| `SX_AAFS_M_M2A` | Median age at first sexual intercourse [Men]: 20–49 (54, 59) |
| `SX_AAFS_W_M20` … `SX_AAFS_W_M45` | Median age, women, by 5-year age band |
| `SX_AAFS_M_M20` … `SX_AAFS_M_M60` | Median age, men, by 5-year age band |
| `SX_SBAG_W_B15` … `SX_SBAG_W_B25` | Cumulative first-sex fraction by exact age (women, all ages) |
| `SX_SBAY_W_B15` … `SX_SBAY_W_B20` | Cumulative first-sex fraction by exact age (young women 15–24) |
| `SX_SBAY_W_MSX` | Median age at first sexual intercourse [Young women] |

To discover other indicators, query the catalogue:

```bash
curl -s 'https://api.dhsprogram.com/rest/dhs/indicators?f=json&perpage=5000' \
  | python -c "import json,sys; \
    [print(f\"{i['IndicatorId']:20} {i['Label']}\") \
     for i in json.load(sys.stdin)['Data'] \
     if 'first sex' in i['Label'].lower()]"
```

### Example: retrieve medians for Zimbabwe

```python
import requests, pandas as pd

resp = requests.get(
    'https://api.dhsprogram.com/rest/dhs/data',
    params={
        'indicatorIds': 'SX_AAFS_W_M2A,SX_AAFS_M_M2A',
        'countryIds': 'ZW',
        'breakdown': 'national',
        'f': 'json',
        'perpage': 50,
    },
)
df = pd.DataFrame(resp.json()['Data'])
print(df[['SurveyYear', 'Indicator', 'Value']])
```

Country codes are two-letter (ZW for Zimbabwe, ZM for Zambia, MW for
Malawi, etc.). Full country list via `/rest/dhs/countries?f=json`.

Sample output (Zimbabwe, current as of the most recent DHS rounds):

| SurveyYear | Indicator | Value |
|---|---|---|
| 1994 | Median age at first sex [Women]: 20–49 | 18.4 |
| 1999 | Median age at first sex [Women]: 20–49 | 18.8 |
| 2005 | Median age at first sex [Women]: 20–49 | 18.7 |
| 2010 | Median age at first sex [Women]: 20–49 | 18.9 |
| 2015 | Median age at first sex [Women]: 20–49 | 18.7 |
| 1994 | Median age at first sex [Men]: 20–49 | 19.3 |
| 1999 | Median age at first sex [Men]: 20–49 | 19.6 |

Notice the Zimbabwe medians (~18.7 women, ~19.5 men) run meaningfully
lower than stisim's defaults (20 women, 21 men) — a concrete case of
why country-specific data matters.

## Translating to stisim parameters

**Median → mean of the lognormal.** For an `ss.lognorm_ex(mean, std)`
distribution, `mean` is the arithmetic mean and `std` is the standard
deviation of the distribution — not of the underlying normal. The DHS
median is close enough to the mean for practical use in stisim's
lognormal shape (the distributions are moderately skewed but not
extreme).

Simple substitution:

```python
debut_f = ss.lognorm_ex(18.7, 3)   # Zimbabwe women 20–49
debut_m = ss.lognorm_ex(19.5, 3)   # Zimbabwe men 20–49
```

**Cumulative curve → distribution parameters.** With `F(15)`, `F(18)`,
`F(20)`, `F(22)`, `F(25)` (fractions who have had first sex by that
exact age), two points identify a two-parameter distribution.
`stisim.data.logn_percentiles_to_pars` inverts a pair of quantiles to
scipy-style lognormal parameters; a short conversion then gives the
`(mean, std)` form `ss.lognorm_ex` wants.

```python
import numpy as np
from stisim.data import logn_percentiles_to_pars

# Two points from the cumulative curve for women in your country,
# e.g. Zimbabwe DHS: 18% had first sex by 15, 60% by 18.
s, scale = logn_percentiles_to_pars(x1=15, p1=0.18, x2=18, p2=0.60)

# Convert scipy (s, scale) to ss.lognorm_ex (mean, std).
mean = scale * np.exp(s ** 2 / 2)
std = mean * np.sqrt(np.exp(s ** 2) - 1)
debut_f = ss.lognorm_ex(mean, std)
```

Pick two points that flank the mass — early (below the median) and
late (above) — for a well-conditioned fit. Additional points can be
used to check the fit: predict `F(x)` at the held-out ages from the
fitted lognormal and compare to the observed cumulative fractions.

Log the source in `analysis-brief.md`:

> Age of debut set from DHS Zimbabwe (most recent round used: DHS 2015,
> indicator SX_AAFS_W_M2A = 18.7, SX_AAFS_M_M2A = 19.5). Standard
> deviation left at stisim default (3) — no spread fit performed for
> this pass.

## Cautions

- **Country and survey coverage.** DHS surveys are not annual and not
  universal. Some countries have one round, some have none. Check
  `/rest/dhs/countries?f=json` and the country's survey history before
  assuming an indicator is available.
- **Reporting bias.** Self-reported age at first sex is known to have
  systematic biases (recall error, social-desirability effects) that
  differ by sex and setting. Use as the best available number, not as
  ground truth.
- **Cohort effects.** Median age at first sex changes across cohorts.
  The 20–49 aggregate mixes cohorts born 1965–1995 in a 2015 survey;
  for a projection whose recent cohorts differ from older ones, the
  5-year band indicators (SX_AAFS_W_M20, SX_AAFS_W_M25, etc.) are
  more informative than the aggregate.
