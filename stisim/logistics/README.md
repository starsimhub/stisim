# Logistics modules

STIsim's `logistics` module models the supply chain behind distribution interventions: a **product**, a **supply** quantity of it, a **collection** of those supplies, and an abstract **intervention** that hands them out to eligible, care-seeking agents.

| File | Contents |
|------|----------|
| `product.py` | `Product` — a distributable item (e.g. a PrEP dose, an appointment slot) with per-time-index efficacy (`eff_by_ti`) and a per-unit `cost` |
| `product_category.py` | `ProductCategory` — enum grouping products by function (`PREP`, `ART`, `BARRIER`); passed as a member, not a string |
| `delivery_mode.py` | `DeliveryMode` — enum for how a product is delivered (`SHOT`, `PILL`, `TOPICAL`, `EDUCATIONAL`) |
| `supply.py` | `Supply` — a non-negative (or `np.inf`) quantity of one `Product`, tracking use and accrued cost |
| `supplies.py` | `Supplies` — an immutable, name-keyed collection of `Supply` objects (`len`/`in`/iteration/`[]`); may be **shared** across interventions, in which case quantity drawdown and `accrued_cost` aggregate across all sharers |
| `supplied_intervention.py` | `SuppliedIntervention` — abstract base (`ss.Intervention`); subclasses **must** implement `step()` to drive the per-step flow `eligibilities → determine_care_seeking → calc_supply_distribution → distribute → use` |

## Notes

- **`SuppliedIntervention` is abstract** and cannot be instantiated directly — subclass it and implement `step()`.
- **Enums are strict**: `Product` accepts `ProductCategory`/`DeliveryMode` *members* only; raw strings raise `ValueError`.
- **Cost granularity**: `SuppliedIntervention.accrued_cost` is one intervention's own spend; `Supplies.accrued_cost` is the pooled total across every intervention sharing that `Supplies`.
