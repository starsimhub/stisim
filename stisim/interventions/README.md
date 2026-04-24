# Intervention modules

STIsim separates **products** (diagnostics with sensitivity/specificity) from **interventions** (delivery strategies defining eligibility, timing, and coverage).

| File | Contents |
|------|----------|
| `base_interventions.py` | `STIDx` (diagnostic product), `STITest` (base testing), `SymptomaticTesting`, `STITreatment`, `PartnerNotification`, `ProductMix` |
| `hiv_interventions.py` | `HIVTest`, `ART`, `VMMC`, `PrEP` and coverage parsing utilities |
| `syphilis_interventions.py` | Syphilis-specific testing and treatment (ANC screening, newborn treatment) |
| `gonorrhea_interventions.py` | Gonorrhea-specific treatment |
| `bv_interventions.py` | BV-specific treatment |
