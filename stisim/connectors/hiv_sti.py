"""
Syphilis-HIV connector for running coinfection analyses
"""

import starsim as ss
from stisim.diseases.syphilis import SyphilisPlaceholder


__all__ = ['hiv_syph', 'hiv_tv', 'hiv_ng', 'hiv_ct', 'hiv_bv']


class hiv_syph(ss.Connector):
    """Connector for bidirectional HIV-syphilis coinfection interactions.

    Models how syphilis increases susceptibility to HIV (default 2.67x) and
    HIV transmissibility (1.2x), and how HIV/AIDS status can modify syphilis
    susceptibility and transmissibility.

    Args:
        hiv_module: The HIV disease module instance.
        syphilis_module: The syphilis disease module instance.
        pars (dict): Optional parameter overrides for relative susceptibility
            and transmissibility multipliers.
    """

    def __init__(self, hiv_module, syphilis_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.syphilis = syphilis_module
        self.define_pars(
            dt='month',

            # Changes to HIV due to syphilis coinfection
            rel_sus_hiv_syph=2.67,  # Relative increase in susceptibility to HIV due to syphilis
            rel_trans_hiv_syph=1.2,  # Relative increase in transmission due to syphilis

            # Changes to syphilis due to HIV coinfection
            rel_sus_syph_hiv=1,         # People with HIV are x more likely to acquire syphilis
            rel_sus_syph_aids=1,        # People with AIDS are x more likely to acquire syphilis
            rel_trans_syph_hiv=1,       # People with HIV are x more likely to transmit syphilis
            rel_trans_syph_aids=1,      # People with AIDS are x more likely to transmit syphilis
        )
        self.update_pars(pars, **kwargs)

        self.define_states(
            ss.FloatArr('rel_sus_hiv_syph', default=self.pars.rel_sus_hiv_syph),
            ss.FloatArr('rel_trans_hiv_syph', default=self.pars.rel_trans_hiv_syph),
            ss.FloatArr('rel_sus_syph_hiv', default=self.pars.rel_sus_syph_hiv),
            ss.FloatArr('rel_sus_syph_aids', default=self.pars.rel_sus_syph_aids),
            ss.FloatArr('rel_trans_syph_hiv', default=self.pars.rel_trans_syph_hiv),
            ss.FloatArr('rel_trans_syph_aids', default=self.pars.rel_trans_syph_aids),
        )

        return

    def step(self):

        # HIV changes due to syphilis
        syphilis = self.syphilis.active
        self.hiv.rel_sus[syphilis] *= self.rel_sus_hiv_syph[syphilis]
        self.hiv.rel_trans[syphilis] *= self.rel_trans_hiv_syph[syphilis]

        # Syphilis changes due to HIV
        if isinstance(self.syphilis, SyphilisPlaceholder):
            return

        hiv = self.hiv.cd4 < 500
        self.syphilis.rel_sus[hiv] *= self.rel_sus_syph_hiv[hiv]
        self.syphilis.rel_trans[hiv] *= self.rel_trans_syph_hiv[hiv]

        aids = self.hiv.cd4 < 200
        self.syphilis.rel_sus[aids] *= self.rel_sus_syph_aids[aids]
        self.syphilis.rel_trans[aids] *= self.rel_trans_syph_aids[aids]

        return


class hiv_tv(ss.Connector):
    """Connector for HIV-trichomoniasis coinfection interaction.

    Models how trichomoniasis increases susceptibility to HIV acquisition
    (default 1.5x). Effect is unidirectional: trichomoniasis affects HIV
    susceptibility only.

    Args:
        hiv_module: The HIV disease module instance.
        tv_module: The trichomoniasis disease module instance.
        pars (dict): Optional parameter overrides.
    """

    def __init__(self, hiv_module, tv_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.tv = tv_module
        self.define_pars(
            # Changes to HIV due to trichonomiasis coinfection
            # Sources:
            #   - https://www.who.int/news-room/fact-sheets/detail/trichomoniasis
            #   - https://pubmed.ncbi.nlm.nih.gov/30341233/
            rel_sus_hiv_tv=1.5,  # Trich infections: 1.5x risk of HIV acquisition.
        )
        self.update_pars(pars, **kwargs)
        return

    def step(self):
        tv = self.tv.infected
        self.hiv.rel_sus[tv] *= self.pars.rel_sus_hiv_tv
        return


class hiv_ng(ss.Connector):
    """Connector for HIV-gonorrhea coinfection interaction.

    Models how gonorrhea increases both susceptibility to HIV acquisition
    (default 1.2x) and HIV transmissibility (default 1.2x).

    Args:
        hiv_module: The HIV disease module instance.
        ng_module: The gonorrhea disease module instance.
        pars (dict): Optional parameter overrides.
    """

    def __init__(self, hiv_module, ng_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.ng = ng_module
        self.define_pars(
            # Changes to HIV due to gonorrhea coinfection
            rel_sus_hiv_ng=1.2,    # Having NG increases risk of HIV acquisition
            rel_trans_hiv_ng=1.2,  # Having NG increases HIV transmissibility
        )
        self.update_pars(pars, **kwargs)

        return

    def step(self):
        ng = self.ng.infected
        self.hiv.rel_sus[ng] *= self.pars.rel_sus_hiv_ng
        self.hiv.rel_trans[ng] *= self.pars.rel_trans_hiv_ng
        return


class hiv_ct(ss.Connector):
    """Connector for HIV-chlamydia coinfection interaction.

    Models how chlamydia infection modifies susceptibility to HIV acquisition.
    The default relative susceptibility multiplier is 1 (no effect), serving
    as a placeholder for user-defined values.

    Args:
        hiv_module: The HIV disease module instance.
        ct_module: The chlamydia disease module instance.
        pars (dict): Optional parameter overrides.
    """

    def __init__(self, hiv_module, ct_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.ct = ct_module
        self.define_pars(
            # Changes to HIV due to chlamydia coinfection
            rel_sus_hiv_ct=1,  # Having an CT infection leads to XXx risk of HIV acquisition.
        )
        self.update_pars(pars, **kwargs)

        return

    def step(self):
        ct = self.ct.infected
        self.hiv.rel_sus[ct] *= self.pars.rel_sus_hiv_ct
        return


class hiv_simplebv(ss.Connector):

    def __init__(self, hiv_module, bv_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.bv = bv_module
        self.define_pars(
            # Changes to HIV due to BV
            rel_sus_hiv_bv=2,       # PLACEHOLDER: BV leads to 2x risk of HIV acquisition
            rel_trans_hiv_bv=2,     # PLACEHOLDER: BV leads to 2x risk of HIV transmission
        )
        self.update_pars(pars, **kwargs)

        return

    def step(self):
        bv = self.bv.infected
        self.hiv.rel_trans[bv] *= self.pars.rel_trans_hiv_bv
        self.hiv.rel_sus[bv] *= self.pars.rel_sus_hiv_bv
        return


class hiv_bv(hiv_simplebv):
    """Connector for HIV-bacterial vaginosis coinfection interaction.

    Extends ``hiv_simplebv`` to use the CST-IV vaginal microbiome state
    (rather than a simple infected flag) to determine which agents receive
    increased HIV susceptibility (default 2x) and transmissibility (default 2x).

    Args:
        hiv_module: The HIV disease module instance.
        bv_module: The BV disease module instance.
        pars (dict): Optional parameter overrides.
    """

    def step(self):
        cst4 = self.bv.cst4.uids
        self.hiv.rel_sus[cst4] *= self.pars.rel_sus_hiv_bv
        self.hiv.rel_trans[cst4] *= self.pars.rel_trans_hiv_bv
        return
