# Disease modules

Each file defines a disease model that can be added to an STIsim simulation.

| File | Disease | Model type |
|------|---------|------------|
| `sti.py` | The SEIS template (`SEIS`), STIsim's results (`STIResults`), and betas by network (`sexual_beta()`) | Template for SEIS-type STIs |
| `hiv.py` | HIV | CD4-based progression (acute/latent/falling) |
| `syphilis.py` | Syphilis | Staged (primary/secondary/latent/tertiary) |
| `chlamydia.py` | Chlamydia | SEIS with sex-stratified symptoms |
| `gonorrhea.py` | Gonorrhea | SEIS with sex-stratified symptoms |
| `trichomoniasis.py` | Trichomoniasis | SEIS with persistence |
| `bv.py` | Bacterial vaginosis | CST-based microbiome model |
| `gud.py` | Genital ulcer disease | Simple SIS |

The diseases are declared (`ss.Disease`: states and transitions), with what is not a state or transition as code in the same class; the bacterial STIs pass their defaults to `SEIS` from `sti.py`.
