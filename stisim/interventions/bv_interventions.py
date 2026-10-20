"""
BV treatment
"""

import starsim as ss
import sciris as sc

__all__ = ['treat_BV']


class treat_BV(ss.treat):
    """
    Treatment of bacterial vaginosis (metronidazole), for newly symptomatic women

    Each woman who develops symptoms (`bv.symptoms`) seeks care with probability
    `p_seek_care` and is treated: she is cured with an efficacy by her stable CST
    (`tx_effectiveness`) and moves along `bv.tx` into the post-treatment state, from which
    she returns to her CST (see `sti.BV`, `dur_post_tx`). Each course is `tx_doses` doses.

    Args:
        pars (dict): parameter values, as below
        p_seek_care (float): probability that a newly symptomatic woman is treated
        tx_effectiveness (dict/list): efficacy by stable CST: `[cst1, cst3, cst4]`, or Starsim 3's dict of `stable_cst1`... (numbers or bernoullis)
        tx_doses (int): doses per course
        start_year (float): first year of treatment (default: the start of the sim)
        stop_year (float): last year of treatment
        eligibility (Cond/func): who is treated, instead of the newly symptomatic (then at each step they are eligible)
        disease (str): the name of the BV module
        kwargs (dict): passed to `ss.treat` (e.g. `capacity=`, `name=`)

    Results: those of `ss.treat` (`delivered`, `cured`, `failed`, `overtreated`), `new_treated`, and `new_doses`.

    Examples:
        ```python
        sim = sti.Sim(diseases='bv', interventions=sti.treat_BV(p_seek_care=0.5, start_year=2005))
        ```
    """
    def __init__(self, pars=None, p_seek_care=0.9, tx_effectiveness=None, tx_doses=14, start_year=None, stop_year=None,
                 eligibility=None, disease='bv', durable_cst1_change=False, name='treat_bv', **kwargs):
        pars = sc.mergedicts(pars)
        p_seek_care = pars.pop('p_seek_care', p_seek_care)
        tx_effectiveness = pars.pop('tx_effectiveness', tx_effectiveness)
        tx_doses = pars.pop('tx_doses', tx_doses)
        start_year = pars.pop('start_year', start_year)
        stop_year = pars.pop('stop_year', stop_year)
        durable_cst1_change = pars.pop('durable_cst1_change', durable_cst1_change)
        for key in ['tx_delay', 'tx_effect_delay', 'dt']: # Unused in Starsim 3
            pars.pop(key, None)
        if durable_cst1_change or set(pars) & {'rr_douching', 'rr_poor_menstrual_hygiene', 'baseline_dur_stable_cst_change'}:
            errormsg = 'A durable change to CST 1 after treatment is not supported in Starsim 4; give sti.BV a post-treatment window (dur_post_tx=) instead'
            raise NotImplementedError(errormsg)
        if pars:
            errormsg = f'Unknown parameters for treat_BV: {sc.strjoin(pars.keys())}'
            raise ValueError(errormsg)
        p_seek_care = getattr(p_seek_care, 'pars', dict(p=p_seek_care))['p'] if isinstance(p_seek_care, ss.Dist) else p_seek_care
        eff = sc.ifelse(tx_effectiveness, [0.95, 0.90, 0.90])
        if isinstance(eff, dict):
            eff = [eff[k] for k in ['stable_cst1', 'stable_cst3', 'stable_cst4']]
        eff = [e.pars.p if isinstance(e, ss.Dist) else e for e in eff]
        self.tx_doses = tx_doses
        product = ss.Treatment(moves=f'{disease}.tx', efficacy=ss.by(stable=list(eff)), name=f'{name}_mtz')
        trigger = dict(where=eligibility, repeat=True) if eligibility is not None else dict(on=ss.state(f'{disease}.symptoms'))
        super().__init__(product, coverage=None if p_seek_care == 1 else p_seek_care, start=start_year, stop=stop_year,
                         disease=disease, name=name, **trigger, **kwargs)
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        for key,mult,label in [('new_treated', 1, 'New people treated'), ('new_doses', self.tx_doses, 'New doses administered')]:
            r = res.delivered
            self.results += ss.Result(key, values=mult*sc.dcp(r.values), scale=False, label=label, timevec=r.timevec, auto_plot=False, module=self.name)
        return
