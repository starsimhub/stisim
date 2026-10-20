"""
Chlamydia trachomatis disease module
"""

import starsim as ss
from stisim.diseases.sti import SEIS

__all__ = ['Chlamydia']


class Chlamydia(SEIS):
    """
    Chlamydia trachomatis, as an SEIS infection (see `sti.SEIS`)

    Args:
        kwargs (dict): passed to `sti.SEIS`
    """
    defaults = dict(
        dur_exp = ss.constant(ss.weeks(1)),
        p_symp = [0.20, 0.54], # https://doi.org/10.1016/j.epidem.2010.04.002
        dur_presymp = [[ss.weeks(1), ss.weeks(10)], [ss.weeks(0.25), ss.weeks(1)]], # Women, men: symptoms appear within days
        p_symp_care = [0.42, 0.83], # See Table 2 in https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
        dur_symp2care = [[ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]],
        dur_asymp2clear = [[ss.months(18), ss.months(1)], [ss.months(12), ss.months(1)]], # Women: 433 days (https://doi.org/10.1016/j.epidem.2010.04.002)
        dur_symp2clear = [[ss.months(18), ss.months(1)], [ss.months(12), ss.months(1)]], # Assumption
        p_pid = 0.0, # Assumption used in https://doi.org/10.1086/598983, based on https://doi.org/10.1016/s0029-7844(02)02118-x
        dur_prepid = ss.lognorm_ex(ss.months(1.5), ss.months(3)),
        beta_m2f = 0.06,
        init_prev = 0.01,
        eff_condom = 0.4, # doi:10.1001/archpedi.159.6.536
    )

    def __init__(self, name='ct', **kwargs):
        super().__init__(name=name, **kwargs)
        return
