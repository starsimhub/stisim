import starsim as ss
import stisim as sti

class gud_syph(ss.Connector):

    def __init__(self, pars=None, **kwargs):
        super().__init__(label='GUD-Syphilis', requires=[sti.GUD, sti.Syphilis])
        self.define_pars(
            rel_sus_syph_gud=2,     # People with GUD are 2x more likely to acquire syphilis
            rel_trans_syph_gud=2,   # People with GUD are 2x more likely to transmit syphilis
            rel_sus_gud_syph=2,     # People with syphilis are 2x more likely to acquire GUD
            rel_trans_gud_syph=2,   # People with syphilis are 2x more likely to transmit GUD
        )
        self.update_pars(pars, **kwargs)

        return

    def update(self, sim):
        """
        GUD-syphilis interactions
        """
        # People with GUD are more likely to acquire syphilis
        sim.diseases.syphilis.rel_sus[sim.people.gud.infected] = self.pars.rel_sus_syph_gud

        # People with GUD are more likely to transmit syphilis
        sim.diseases.syphilis.rel_trans[sim.people.gud.infected] = self.pars.rel_trans_syph_gud

        # People with active syphilis are more likely to acquire GUD
        sim.diseases.gud.rel_sus[sim.diseases.syphilis.active] = self.pars.rel_sus_gud_syph

        # People with active syphilis are more likely to transmit GUD
        sim.diseases.gud.rel_trans[sim.diseases.syphilis.active] = self.pars.rel_trans_gud_syph

        return
