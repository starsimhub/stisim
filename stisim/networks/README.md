# Network modules

Sexual contact networks for STIsim simulations. Each file defines one or more
network classes; the table below shows where to find each.

| File | Provides |
|------|----------|
| `base.py` | `BaseNetwork`, `BasePars`, `NetworkPars`, `NoPartnersFound` |
| `mf.py` | `MFNetwork`, `MFPars` — heterosexual partnerships |
| `fsw.py` | `SWNetwork`, `SWPars` — sex-work (FSW–client) partnerships |
| `msm.py` | `AgeMatchedMSM`, `AgeApproxMSM` — men who have sex with men |
| `layered_networks.py` | `StructuredSexual`, `PriorPartners` — composite networks |
| `matchers.py` | Pair-formation algorithms used by `MFNetwork.match_pairs` |

## Layered networks

A *layered network* bundles two or more network behaviours into a single class so a sim can model their combined dynamics without juggling multiple network objects.

`StructuredSexual` is the canonical example and the default choice in most STIsim models: it layers heterosexual partnerships (from `MFNetwork`) with sex-work partnerships (from `SWNetwork`) on a single edge list. Real populations have both regular and transactional partnerships, and STIs transmit across both — so most calibrated models reach for `StructuredSexual` rather than `MFNetwork` alone.

`PriorPartners` is a supplementary network you add *alongside* an MF network when `recall_prior=True`. It stores recently-dissolved edges so partner-notification interventions can trace them.

## Pair-formation algorithms

`MFNetwork.match_pairs` dispatches through a registry in `matchers.py`, selected by the `match_method` parameter (string or callable):

```python
import stisim as sti

# Default (kdtree_nn — nearest-neighbour by age)
net = sti.MFNetwork()

# Explicit choice
net = sti.MFNetwork(match_method='desired_age_bucket')

# Custom matcher
def my_matcher(net):
    ...
    return p1, p2  # both ss.uids
net = sti.MFNetwork(match_method=my_matcher)
```

See `matchers.py` for the seven built-in methods and the tradeoffs between them.
