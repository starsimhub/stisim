"""
Syphilis-HIV connector for running coinfection analyses
"""

import starsim as ss
from stisim.diseases.syphilis import SyphilisPlaceholder


__all__ = ['hiv_syph', 'hiv_trich']


class hiv_syph(ss.Connector):

    def __init__(self, hiv_module, syphilis_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.syphilis = syphilis_module
        self.default_pars(
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

        self.add_states(
            ss.FloatArr('rel_sus_hiv_syph', default=self.pars.rel_sus_hiv_syph),
            ss.FloatArr('rel_trans_hiv_syph', default=self.pars.rel_trans_hiv_syph),
            ss.FloatArr('rel_sus_syph_hiv', default=self.pars.rel_sus_syph_hiv),
            ss.FloatArr('rel_sus_syph_aids', default=self.pars.rel_sus_syph_aids),
            ss.FloatArr('rel_trans_syph_hiv', default=self.pars.rel_trans_syph_hiv),
            ss.FloatArr('rel_trans_syph_aids', default=self.pars.rel_trans_syph_aids),
        )

        return

    def update(self):

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


class hiv_trich(ss.Connector):

    def __init__(self, hiv_module, trich_module, pars=None, **kwargs):
        super().__init__()

        self.hiv = hiv_module
        self.trich = trich_module
        self.default_pars(
            # Changes to HIV due to trichonomiasis coinfection
            # Sources:
            #   - https://www.who.int/news-room/fact-sheets/detail/trichomoniasis
            #   - https://pubmed.ncbi.nlm.nih.gov/30341233/
            rel_sus_hiv_trich=ss.normal(loc=1.5, scale=0.25),  # Trich infections: 1.5x risk of HIV acquisition.
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            ss.FloatArr('rel_sus_hiv_trich', default=self.pars.rel_sus_hiv_trich),
        )

        return

    def update(self):
        trich = self.trich.infected
        self.hiv.rel_sus[trich] *= self.rel_sus_hiv_trich[trich]
        return

