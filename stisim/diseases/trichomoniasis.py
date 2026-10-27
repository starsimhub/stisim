"""
Trichomoniasis disease module
"""

import sciris as sc
import starsim as ss
from stisim.diseases.sti import SEIS, by_sex, lognorm_by_sex

__all__ = ['Trichomoniasis']


class Trichomoniasis(SEIS):
    """
    Trichomonas vaginalis disease module.

    Models trichomoniasis as an SEIS infection (see `sti.SEIS`) with sex-specific symptoms and
    clearance. Most asymptomatic infections in women do not clear on their own (`p_clear`), and
    persist for `dur_persist`.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'tv'``)
        kwargs (dict): parameter values (e.g. ``beta_m2f=0.1``), and the arguments of `sti.SEIS`
    """
    def __init__(self, pars=None, name='tv', **kwargs):
        defaults = dict(
            dur_exp = ss.constant(0), # Infectious at once

            # Symptoms
            p_symp = by_sex(0.4, 0.5), # https://sti.bmj.com/content/76/4/248
            dur_presymp = lognorm_by_sex([ss.weeks(1), ss.weeks(12)], [ss.weeks(0.25), ss.weeks(1)]), # Men: symptoms should appear within days

            # Care seeking
            p_symp_care = by_sex(0.39, 0.27),
            dur_symp2care = lognorm_by_sex([ss.months(2), ss.months(1)], [ss.weeks(1), ss.weeks(2)]),

            # Clearance. The average duration of infection in women is at least 3–5 years, and about 4 months in men (https://sti.bmj.com/content/76/4/248)
            dur_asymp2clear = lognorm_by_sex([ss.months(48), ss.months(6)], [ss.weeks(26), ss.weeks(4)]),
            dur_symp2clear = lognorm_by_sex([ss.weeks(20), ss.weeks(4)], [ss.weeks(18), ss.weeks(4)]), # Assumptions
            p_clear = 0.1, # Most women do not spontaneously clear, men do (https://sti.bmj.com/content/76/4/248)
            dur_persist = ss.constant(ss.years(100)),

            # PID
            p_pid = 0.0,
            dur_prepid = ss.lognorm_ex(ss.weeks(6), ss.weeks(4)),

            # Transmission and initial conditions
            beta_m2f = 0.1,
            init_prev = 0.01,
            eff_condom = 0.0,
        )
        super().__init__(name=name, **sc.mergedicts(defaults, pars, kwargs))
        self.requires = 'structuredsexual'
        return
