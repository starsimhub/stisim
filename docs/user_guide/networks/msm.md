# MSM networks

STIsim provides three networks for modeling sexual contact between men. `AgeMatchedMSM` and `AgeApproxMSM` extend [`MFNetwork`](structured_sexual.md), inheriting risk groups, concurrency, and partnership classification, and differ only in how partners are matched. `MSMScaleFreeNetwork` extends `BaseNetwork` directly and uses a preferential-attachment kernel instead of risk groups.

| Class | Matching strategy | When to use |
|-------|------------------|-------------|
| `sti.AgeMatchedMSM` | Sort eligible males by age and pair sequentially | Lightweight; partners always have similar ages |
| `sti.AgeApproxMSM` | Reuses `MFNetwork` age-difference preferences to match split groups of eligible males | More flexible age mixing |
| `sti.MSMScaleFreeNetwork` | Preferential-attachment (rich-get-richer) degree kernel with Markovian edge deletion | Scale-free degree distribution; concentrated transmission |

## Participation

All three networks expose a `p_msm` parameter (default `ss.bernoulli(p=0.015)`) which sets the fraction of males that participate in MSM partnerships. Participation is assigned at sexual debut and is fixed for the lifetime of the agent.

```python
msm = sti.AgeMatchedMSM(p_msm=ss.bernoulli(p=0.05))  # 5% of males
```

## Combining with the heterosexual network

MSM networks are typically added alongside the default `StructuredSexual` network so that bisexual men contribute to transmission in both:

```python
sim = sti.Sim(
    diseases='hiv',
    networks=[sti.StructuredSexual(), sti.AgeMatchedMSM()],
)
```

Each network maintains its own edges and partnership classification. Disease modules see contacts from all networks each timestep.

## Scale-free network

`sti.MSMScaleFreeNetwork` (Whittles-2019 S2 kernel) forms edges via a rate proportional to a degree-based pair weight, producing a scale-free degree distribution where a minority of high-degree agents drive most transmission. Edges are deleted at random subject to a hard duration cap.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_msm` | `ss.bernoulli(p=0.015)` | Fraction of males in the MSM pool |
| `target_mean_degree` | 2.0 | Target mean concurrent partners per pool agent |
| `target_mean_dur` | `ss.years(2)` | Target mean edge duration |
| `max_edge_dur` | `ss.years(10)` | Hard cap on edge persistence |
| `phi` | 1.0 | Turnover parameter (sets initial network density) |

```python
msm = sti.MSMScaleFreeNetwork(target_mean_degree=3.0)
```

Note: unlike the other networks, this one is not branching-stable under Common Random Numbers — adding an intervention that perturbs the population reshuffles edges from that point on. Reproducibility holds per `(rand_seed, ti)`. It also does not record edge expirations, so analyzers that depend on them (e.g. `PartnershipFormationAnalyzer`) skip it.

## API

- `sti.AgeMatchedMSM` — exact-age matching; cheap and deterministic given inputs.
- `sti.AgeApproxMSM` — distribution-based matching; reuses age-preference machinery from `MFNetwork`.
- `sti.MSMScaleFreeNetwork` — preferential-attachment kernel; scale-free degree distribution.
