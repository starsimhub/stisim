# STIsim

[![tests](https://github.com/starsimhub/stisim/actions/workflows/tests.yaml/badge.svg)](https://github.com/starsimhub/stisim/actions/workflows/tests.yaml)
[![PyPI](https://img.shields.io/pypi/v/stisim?label=PyPI)](https://pypi.org/project/stisim/)


STIsim is an agent-based modeling framework in which users can design and configure simulations of co-circulating sexually-transmitted diseases. STIsim uses the [Starsim](https://starsim.org) architecture, and belongs to the Starsim model suite which also includes [Covasim](https://covasim.org), [HPVsim](https://hpvsim.org), and [FPsim](https://fpsim.org).

## Requirements

Python 3.9-3.14 or R.

We recommend, but do not require, installing STIsim in a virtual environment, such as [Miniconda](https://docs.anaconda.com/miniconda/).

## Installation

### Python

STIsim is most easily installed via [PyPI](https://pypi.org/project/stisim/):
```sh
pip install stisim
```

Or with [uv](https://github.com/astral-sh/uv):
```sh
uv init example
cd example
uv add stisim
```

STIsim can also be installed locally (including optional dependencies for testing and documentation). To do this, clone first this repository, then run:
```sh
pip install -e .[dev]
```

### R

STIsim can be used from R via [rstarsim](https://github.com/starsimhub/rstarsim), which calls the Python engine through reticulate. The Python packages are required regardless of whether you use the R or Python interface.

Install the R packages:
```r
install.packages(c("reticulate", "devtools"))
devtools::install_github("starsimhub/rstarsim")
```

On first use, `rstarsim` will set up a conda environment automatically if needed. To use an existing environment instead:
```r
library(starsim)
load_starsim("my_env_name")
```

See [r.starsim.org](https://r.starsim.org) for more information on using Starsim from R.

## Usage and documentation

Documentation, including tutorials and a user guide, is available at https://docs.stisim.org. Additional resources:
1. The [examples](#examples) below show how to set up, calibrate, and run country-level models
2. Read the articles that have been published about analyses using STIsim (see [references](#references) below)
3. Email us: [info@starsim.org](mailto:info@starsim.org)

## Development roadmap
The roadmap for future model development can be viewed [here](https://github.com/orgs/starsimhub/projects/26/views/6).

## References

Publications using STIsim include:

1. **Reduction in overtreatment of gonorrhoea and chlamydia through point-of-care testing compared with syndromic management for vaginal discharge: a modelling study for Zimbabwe** (2026). Stuart RM, Newman LM, Manguro G, Dziva Chikwari C, Marks M, Peters RPH, Klein D, Snyder L, Kerr C, Rao DW. *Sex Transm Infect* https://doi.org/10.1136/sextrans-2025-056646. Preprint: https://doi.org/10.21203/rs.3.rs-8843262/v1 

2. **Point-of-care testing to strengthen sexually transmitted infection case management in resource-constrained settings** (2026). Peters RPH, Manguro G, Ong'wen PA, Mdingi MM, Applegate TL, Stuart R, Harding-Esch EM, Manabe YC, Ndowa F, Van Der Pol B. *Sex Transm Infect*, https://doi.org/10.1136/sextrans-2025-056833.

3. **Estimating the value of novel tests for active syphilis in Zimbabwe: how much overtreatment can be avoided?** (2026). Stuart RM, Marks M, Dziva Chikwari C, Muellenmeister AM, Abeysuriya RG, Manguro G, Peters RPH, Rao DW, Newman LM. Forthcoming in *STD* (accepted May 3, 2026). Preprint: https://doi.org/10.21203/rs.3.rs-9559421/v1


## Examples

The following repositories contain end-to-end analyses built with STIsim, and are a good starting point for understanding how to set up, calibrate, and run country-level models:

- **[hiv_kenya](https://github.com/starsimhub/hiv_kenya)** -- HIV transmission model for Kenya with structured sexual networks, testing (FSW-targeted and general population), ART, and PrEP. Includes both Python and R interfaces and Optuna-based calibration.
- **[hiv_zambia](https://github.com/starsimhub/hiv_zambia)** -- HIV transmission model for Zambia. Similar structure to the Kenya model, useful for comparing how the same framework is adapted to a different country context.
- **[stisim_vddx_zim](https://github.com/starsimhub/stisim_vddx_zim)** -- Multi-STI model (gonorrhea, chlamydia, trichomoniasis, BV) for Zimbabwe evaluating point-of-care diagnostics vs. syndromic management for vaginal discharge. Demonstrates co-circulating STIs, intervention comparison, and scenario analysis.
- **[syph_dx_zim](https://github.com/starsimhub/syph_dx_zim)** -- Joint HIV-syphilis model for Zimbabwe evaluating syphilis diagnostic algorithms. Features coinfection connectors, stage-specific syphilis transmission, congenital syphilis, and treatment pathway analysis.
- **[hiv_vmb_southafrica](https://github.com/starsimhub/hiv_vmb_southafrica)** -- HIV-vaginal microbiome transmission model calibrated to South Africa. This model was developed to evaluate novel products to shift the vaginal microbiome and their impact on HIV and pre-term birth.

## Contributing

We welcome all contributions to STIsim! Please refer to our [code of conduct](https://github.com/starsimhub/stisim/blob/main/code_of_conduct.md) and [contributors's guide](https://github.com/starsimhub/stisim/blob/main/contributing.md). You'll find information there about our style guide, which is essential reading prior to contributing. Questions or comments can be directed to [info@starsim.org](mailto:info@starsim.org), or on this project's [GitHub](https://github.com/starsimhub/stisim) page. 

See `.github/workflows/README.md` for details on publishing new releases of STIsim.

## Disclaimer

The code in this repository was developed by IDM, the Burnet Institute, and other collaborators to support our joint research on flexible agent-based modeling. We've made it publicly available under the MIT License to provide others with a better understanding of our research and an opportunity to build upon it for their own work. We make no representations that the code works as intended or that we will provide support, address issues that are found, or accept pull requests. You are welcome to create your own fork and modify the code to suit your own modeling needs as permitted under the MIT License.
