# MSM networks

STIsim provides two networks for modeling sexual contact between men. Both extend [`StructuredSexual`](structured_sexual.md), inheriting risk groups, concurrency, and partnership classification, but differ in how partners are matched.

| Class | Matching strategy | When to use |
|-------|------------------|-------------|
| `sti.AgeMatchedMSM` | Sort eligible males by age and pair sequentially | Lightweight; partners always have similar ages |
| `sti.AgeApproxMSM` | Reuses `StructuredSexual` age-difference preferences to match split groups of eligible males | More flexible age mixing |

## Participation

Both networks expose an `msm_share` parameter (default `ss.bernoulli(p=0.015)`) which sets the fraction of males that participate in MSM partnerships. Participation is assigned at sexual debut and is fixed for the lifetime of the agent.

```python
msm = sti.AgeMatchedMSM(msm_share=ss.bernoulli(p=0.05))  # 5% of males
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

## API

- `sti.AgeMatchedMSM` — exact-age matching; cheap and deterministic given inputs.
- `sti.AgeApproxMSM` — distribution-based matching; reuses age preference machinery from `StructuredSexual`.
