"""
Trichomoniasis disease module
"""

import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['Trichomoniasis']


class Trichomoniasis(SEIS):
    """
    Trichomonas vaginalis, as an SEIS infection (see `sti.SEIS`)

    Most asymptomatic infections in women do not clear spontaneously (`p_clear`), and
    persist for `dur_persist`; men's do.

    Args:
        kwargs (dict): passed to `sti.SEIS`
    """
    defaults = dict(
        dur_exp = ss.constant(0),
        p_symp = [0.4, 0.5], # https://sti.bmj.com/content/76/4/248
        dur_presymp = [[ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(1)]], # Women, men: symptoms appear within days
        p_symp_care = [0.39, 0.27],
        dur_symp2care = [[ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]],
        dur_asymp2clear = [[ss.months(48), ss.months(6)], [ss.weeks(26), ss.weeks(4)]], # At least 3-5 years in women, about 4 months in men (https://sti.bmj.com/content/76/4/248)
        dur_symp2clear = [[ss.weeks(20), ss.weeks(4)], [ss.weeks(18), ss.weeks(4)]], # Assumptions
        p_clear = 0.1, # Most women do not clear spontaneously, men do (https://sti.bmj.com/content/76/4/248)
        dur_persist = ss.constant(ss.years(100)),
        p_pid = 0.0,
        dur_prepid = ss.lognorm_ex(ss.weeks(6), ss.weeks(4)),
        eff_condom = 0.0,
        beta_m2f = 0.1,
        init_prev = 0.01,
    )

    def __init__(self, name='tv', **kwargs):
        super().__init__(name=name, **kwargs)
        return
