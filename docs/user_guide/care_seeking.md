# Care-seeking

`CareSeeking` assigns each agent a cross-disease care-seeking propensity — a latent
tendency to seek testing and treatment. Interventions can read this propensity so
that the same agents are consistently more (or less) likely to engage with care
across diseases and over time, rather than being re-randomised at every contact.

```python
sim = sti.Sim(
    diseases='hiv',
    interventions=[sti.HIVTest(test_prob_data=0.2)],
    demographics=[sti.CareSeeking()],   # shared propensity used by interventions
)
```

The module holds a per-agent propensity that downstream interventions (e.g. testing
and ART prioritisation) can use to break ties or weight eligibility, so high-engagers
and low-engagers behave coherently across the whole cascade.

> **Stub** — expand with the propensity distribution, how interventions consume it,
> and a worked example showing correlated testing/treatment uptake. See the API
> reference for [`care_seeking`](../api/care_seeking.qmd).
