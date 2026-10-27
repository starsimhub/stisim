"""
Chlamydia trachomatis disease module
"""

import sciris as sc
import starsim as ss
from stisim.diseases.sti import SEIS, by_sex, lognorm_by_sex

__all__ = ['Chlamydia']


class Chlamydia(SEIS):
    """
    Chlamydia trachomatis disease module.

    Models chlamydia as an SEIS infection (see `sti.SEIS`) with a latent period of a week,
    sex-specific symptoms and care seeking, and natural clearance.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'ct'``)
        kwargs (dict): parameter values (e.g. ``beta_m2f=0.08``), and the arguments of `sti.SEIS`
    """
    def __init__(self, pars=None, name='ct', **kwargs):
        defaults = dict(
            dur_exp = ss.constant(ss.weeks(1)),

            # Symptoms
            p_symp = by_sex(0.20, 0.54), # https://doi.org/10.1016/j.epidem.2010.04.002
            dur_presymp = lognorm_by_sex([ss.weeks(1), ss.weeks(10)], [ss.weeks(0.25), ss.weeks(1)]), # Men: symptoms should appear within days

            # Care seeking
            p_symp_care = by_sex(0.42, 0.83), # See Table 2 in https://docs.google.com/document/d/16t46nTL2qMHmA0C1gSPz8OhI6ccy6vVv3OCfkmYFUtw/edit?tab=t.0
            dur_symp2care = lognorm_by_sex([ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]),

            # Clearance
            dur_asymp2clear = lognorm_by_sex([ss.months(18), ss.months(1)], [ss.months(12), ss.months(1)]), # Women: 433 days (https://doi.org/10.1016/j.epidem.2010.04.002); men as above
            dur_symp2clear = lognorm_by_sex([ss.months(18), ss.months(1)], [ss.months(12), ss.months(1)]), # Assumption

            # PID
            p_pid = 0.0, # Assumption used in https://doi.org/10.1086/598983, based on https://doi.org/10.1016/s0029-7844(02)02118-x
            dur_prepid = ss.lognorm_ex(ss.months(1.5), ss.months(3)),

            # Transmission and initial conditions
            beta_m2f = 0.06,
            init_prev = 0.01,
            eff_condom = 0.4, # doi:10.1001/archpedi.159.6.536
        )
        super().__init__(name=name, **sc.mergedicts(defaults, pars, kwargs))
        return
