# Pair matching

Each timestep, [`MFNetwork`](structured_sexual.md) forms new heterosexual partnerships by matching *looking* women to *eligible* men. The algorithm that does this is selected by the `match_method` parameter (default `'closest_age_tapered_seeking'`). All matchers share the same inputs and contract; they differ only in how they turn age preferences into pairs.

## The matching problem

```
each timestep
      │
      ▼
 eligible women ──► sample desired partner age ──► desired ages ┐
 (post-debut,        (own age + gap drawn from                  │
  under-partnered,    age_diff_pars, by age group               ├─► match_method ─► (p1, p2)
  "looking")          × risk group)                             │        pairs
                                                                │
 eligible men ───────────────────────────────────► actual ages ┘
 (post-debut, under-partnered)                                          │
                                                                        ▼
                                                              add_pairs assigns
                                                              type / duration / acts
```

- **Eligible women** are post-debut, under-partnered (`partners < concurrency`), and pass the `p_pair_form` "looking" filter. Each draws a **desired** male partner age — her own age plus an age gap sampled from `age_diff_pars` (varies by her age group and risk group; see [age mixing](structured_sexual.md#age-mixing)).
- **Eligible men** are post-debut and under-partnered, with their **actual** ages.
- A matcher returns `(p1, p2)` — equal-length `ss.uids` arrays of male and female partners — or raises `NoPartnersFound` if it can pair no one. Relationship type (stable/casual/one-time), duration, and coital acts are assigned afterwards by `add_pairs`, not by the matcher.

## Default: `closest_age_tapered_seeking`

The default matcher pairs each woman with the closest-aged available man at or above her target, and damps older women's search so realized age gaps stay faithful to `age_diff_pars` even when older men are scarce. It runs in four steps.

**1. Taper older women's search.** A woman's chance of looking this step falls with age, reaching zero at `f_partnership_taper_cut` (default 55):

```
P(woman looks)
  1 ┤●●●●●●●●●●●●●●●●
    │                ●●
    │                  ●●●        ∝ ((cut − age) / offset) ^ 1.5
  0 ┤────────────────────●●●●●──► age
    0                  ↑          ↑
              taper onset        cut  (f_partnership_taper_cut)
            (cut − offset)
```

`offset` is `mean_gap + 3·sd`, the widest plausible age gap across age/risk groups. Without this taper, older women who draw large gaps cannot find partners, biasing realized gaps downward.

**2. Sample desired ages** for the women who are looking (own age + gap), and collect eligible men's actual ages.

**3. Closest-age sweep.** Sort women by desired age and men by actual age, then sweep once with a two-pointer scan: each woman takes the nearest still-available man near her target (no man is reused). A woman whose nearest free man is more than `max_deviation` (1 year) off her target is skipped — later, higher-target women may still match.

```
F desired (sorted):  24    27    33        M actual (sorted): 22  25  28  31  34

   24 ─► nearest free man near 24 ─► 25     ✓  gap 1
   27 ─► nearest free man near 27 ─► 28     ✓  gap 1
   33 ─► nearest free man near 33 ─► 34     ✓  gap 1
                                            22, 31 left unmatched
```

**4. Trim extremes.** Drop any pair whose realized age gap exceeds `mean_gap + 3·sd`, removing unrealistic matches the sweep may have forced.

A small residual downward bias (<1 year) in the female–male age gap remains, driven mainly by debut-age limits truncating the low end of young women's preferences. Note also that in calibration `age_diff_pars` is not orthogonal to the concurrency and duration parameters, since those shape the pool of people looking each step.

## Other matchers

`match_method` accepts any of these registry keys (or a callable — see below):

| `match_method` | How it pairs | Age-gap fidelity |
|----------------|--------------|------------------|
| `closest_age_tapered_seeking` *(default)* | sorted closest-age sweep + older-woman taper + extreme-gap trim | highest |
| `kdtree_nn` | KD-tree nearest-neighbour on male age; contested men go to the closest woman | high |
| `lsa` | globally optimal assignment over the full age-distance matrix; `O(n³)`, reference only | optimal (slow) |
| `greedy_old_enough` | each woman (by ascending target) takes the youngest free man ≥ her target | one-sided |
| `sort_pair` | sort both groups by age, zip, truncate to the shorter | moderate |
| `sort_bisect` | sort + trim non-overlapping age tails + subsample (legacy production) | low (boundaries only) |
| `desired_age_bucket` | 1-year desired-age buckets, sample a man within | coarse ⁽¹⁾ |
| `band_match` | 5-year age bands, shuffle and zip within each band | coarse ⁽¹⁾ |

⁽¹⁾ `desired_age_bucket` and `band_match` draw from their own RNG seeded per `(rand_seed, ti)` — reproducible for a fixed run, but not branching-stable under Starsim's Common Random Numbers (adding an intervention can reshuffle pairings).

## Custom matchers

Pass a callable as `match_method`. It receives the network and returns `(p1, p2)`. Two helpers cover the common setup:

```python
import stisim as sti
from stisim.networks import NoPartnersFound

def my_matcher(net):
    f_looking, m_eligible = net._get_eligible()      # women: ss.uids; men: bool mask (use m_eligible.uids)
    desired = net._sample_desired_ages(f_looking)    # target male age per looking woman
    ...                                              # your pairing logic
    if nothing_matched:
        raise NoPartnersFound()
    return p1, p2                                     # equal-length male, female ss.uids

net = sti.MFNetwork(match_method=my_matcher)
```

`p1` must be male uids and `p2` female uids, of equal length. Raise `NoPartnersFound` when no pairing is possible for the timestep.
