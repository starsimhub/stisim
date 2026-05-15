# Disease modules

Each file defines a disease model that can be added to an STIsim simulation.

| File | Disease | Model type |
|------|---------|------------|
| `sti.py` | Base classes (`BaseSTI`, `SEIS`, `BaseSTIPars`, `STIPars`) | Template for SEIS-type STIs |
| `hiv.py` | HIV | CD4-based progression (acute/latent/falling) |
| `syphilis.py` | Syphilis | Staged (primary/secondary/latent/tertiary) |
| `chlamydia.py` | Chlamydia | SEIS with sex-stratified symptoms |
| `gonorrhea.py` | Gonorrhea | SEIS with sex-stratified symptoms |
| `trichomoniasis.py` | Trichomoniasis | SEIS with persistence |
| `bv.py` | Bacterial vaginosis | CST-based microbiome model |
| `gud.py` | Genital ulcer disease | Simple SIS |

All disease classes extend `BaseSTI` (or `SEIS` for bacterial STIs) from `sti.py`. Each disease has a companion `*Pars` class defining default parameters.
