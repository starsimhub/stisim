# Documentation assets

## `architecture.{dot,svg,png}`

Source-of-truth diagram of the Starsim ecosystem (engine → disease packages → localizations → analyses), with the HIVsim/STIsim relationship made explicit. Used on docs.stisim.org (see `docs/user_guide/index.md`).

Source: `architecture.dot` (Graphviz). Regenerate after editing:

```bash
cd docs/assets
dot -Tsvg architecture.dot > architecture.svg
dot -Tpng -Gdpi=150 architecture.dot > architecture.png
```

Install Graphviz with either:

```bash
brew install graphviz                            # macOS Homebrew
conda install -c conda-forge graphviz            # any conda env (works offline of brew)
```

The `stisim` cluster's module list (diseases/networks/interventions/logistics/connectors) is the part most likely to drift as features are added — check it against `stisim/__init__.py` and the CHANGELOG when regenerating for a release. The other `*sim` package boxes (covasim/hpvsim/fpsim/tbsim) describe sibling repos and aren't verified from this repo; update them only with input from those packages' maintainers.
