# HIVsim

HIVsim is a thin convenience wrapper that ships sensible defaults for HIV simulations. It is included in the STIsim package and imported separately:

```python
import stisim as sti
import hivsim
```

The import alias convention is to use `import hivsim` (no alias), because aliases like `hiv` would collide with disease instances such as `hiv = hivsim.HIV()`.

## What HIVsim adds over `sti.Sim`

When you create `hivsim.Sim()` with no arguments, you get an HIV-focused simulation with these defaults pre-wired:

- **Disease**: `sti.HIV` (inserted as the first disease)
- **Demographics**: `ss.Pregnancy()`, `ss.Deaths()`
- **Networks**: `sti.StructuredSexual()`, `ss.MaternalNet()`, and `ss.BreastfeedingNet()` (only when pregnancy is present, to support breastfeeding-mediated MTCT)
- **Interventions**: `sti.HIVTest()`, `sti.ART()`, `sti.VMMC()`, `sti.Prep()`

Compared to `sti.Sim`, HIVsim:

- Always includes HIV first in the disease list.
- Splits parameter passing: any parameter belonging to `HIVPars` is routed to the HIV module, anything else goes to the sim.
- Provides a `demo()` helper for runnable examples.

## Quick start

```python
import hivsim

sim = hivsim.Sim(n_agents=1000, dur=40)
sim.run()
sim.plot('hiv', annualize=True)
```

## Routing parameters to HIV vs the sim

Any keyword argument that matches a key in `sti.HIVPars()` is routed to HIV; the rest are passed to `sti.Sim`. To pass them explicitly:

```python
sim = hivsim.Sim(
    sim_pars={'n_agents': 5000, 'dur': 50},
    hiv_pars={'init_prev': 0.1, 'beta_m2f': 0.06},
)
```

## Demos

`hivsim.demo()` runs a packaged example end-to-end:

```python
hivsim.demo()              # Minimal HIV sim
hivsim.demo('zimbabwe')    # Zimbabwe model with calibrated parameters and UNAIDS data
sim = hivsim.demo('zimbabwe', run=False, n_agents=500)  # Just construct it
```

Available examples are listed in `hivsim.sim.EXAMPLES`. Each example is implemented as a standalone module under `hivsim_examples/`.

## When to use HIVsim vs `sti.Sim`

Use **HIVsim** when:

- HIV is the primary disease of interest.
- The default network/intervention bundle is appropriate.
- You want a quick path from script to HIV epidemic projection.

Use **`sti.Sim`** directly when:

- You are modeling multiple co-circulating STIs.
- You need to override the default networks or interventions before the sim is built.
- You want full control over module ordering.
