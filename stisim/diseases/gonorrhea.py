"""
Gonorrhea disease module
"""

import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['Gonorrhea']


class Gonorrhea(SEIS):
    """
    Neisseria gonorrhoeae, as an SEIS infection (see `sti.SEIS`)

    With `amr=True`, infections carry a strain, `amr` (`'sens'` or `'res'`), passed on by
    whoever infected them; treatment can then fail more often for resistant strains, and
    resistance can be acquired on failure (see `sti.GonorrheaTreatment`).

    Args:
        amr (bool/list): stratify infections by drug resistance: True for `['sens', 'res']`, or the levels
        init_res (float): with `amr`, the share of initial infections that are resistant
        kwargs (dict): passed to `sti.SEIS`
    """
    defaults = dict(
        dur_exp = ss.constant(0),
        p_symp = [0.35, 0.65], # https://doi.org/10.1371/journal.pone.0143304
        dur_presymp = [[ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(1)]], # Women, men: symptoms appear within days
        p_symp_care = [0.66, 0.83], # See Table 2: https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
        dur_symp2care = [[ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]],
        dur_asymp2clear = [[ss.months(8), ss.months(2)], [ss.months(6), ss.months(3)]],
        dur_symp2clear = [[ss.months(9), ss.months(2)], [ss.months(6), ss.months(3)]], # Assumption
        p_pid = 0.0, # TODO
        dur_prepid = ss.lognorm_ex(ss.months(1.5), ss.months(3)),
        beta_m2f = 0.06,
        init_prev = 0.01,
        eff_condom = 0.9,
    )

    def __init__(self, name='ng', amr=False, init_res=0.0, **kwargs):
        super().__init__(name=name, **kwargs)
        if amr:
            levels = ['sens', 'res'] if amr is True else list(amr)
            self.stratify(amr=levels, kind='strain', split=[1 - init_res, init_res] + [0]*(len(levels) - 2))
        return
