# STIsim package

Source for the `stisim` Python package. STIsim extends [Starsim](https://docs.starsim.org) with sexually-transmitted-infection-specific modules.

| File / folder | Contents |
|---|---|
| `sim.py` | The `Sim` class — STIsim's entry point, subclassing `ss.Sim` to wire in default networks, demographics, and connectors. |
| `parameters.py` | `SimPars` and parameter-handling utilities. |
| `diseases/` | Disease modules (HIV, syphilis, chlamydia, gonorrhea, trichomoniasis, BV, GUD) and base classes (`BaseSTI`, `SEIS`). See `diseases/README.md`. |
| `networks.py` | Sexual contact networks: `StructuredSexual`, `PriorPartners`, `AgeMatchedMSM`, `AgeApproxMSM`. |
| `interventions/` | Testing, treatment, partner notification, and HIV-specific interventions (ART, VMMC, PrEP). See `interventions/README.md`. |
| `connectors/` | Coinfection connectors that adjust susceptibility/transmissibility between co-circulating diseases. See `connectors/README.md`. |
| `analyzers.py` | Result-tracking analyzers (coinfection statistics, sex-worker statistics, etc.). |
| `demographics.py` | Pregnancy and migration modules. |
| `care_seeking.py` | Care-seeking behavior used by symptomatic-testing interventions. |
| `calibration.py` | `Calibration` class plus build/evaluation helpers built on top of Starsim's calibration framework. |
| `data/` | Data loaders and downloaders for country-specific demographic and epidemiological inputs. |
| `utils.py` | Shared utilities. |
| `version.py` | Package version string. |
