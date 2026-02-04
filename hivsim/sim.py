import starsim as ss
import stisim as sti

class Sim(sti.Sim):
    def __init__(self, pars=None, **kwargs):
        hiv_pars = {}
        default_hiv_keys = sti.diseases.hiv.HIVPars().keys()
        for key in kwargs.keys():
            if key in default_hiv_keys:
                hiv_pars[key] = kwargs.pop(key)
        hiv = sti.HIV(pars=hiv_pars)
        super().__init__(pars=pars, diseases=[hiv], **kwargs)
        return