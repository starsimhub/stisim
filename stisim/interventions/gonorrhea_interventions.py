"""
Gonorrhea treatment, with antimicrobial resistance
"""

import starsim as ss
from stisim.interventions.sti_interventions import STITreatment

__all__ = ['GonorrheaTreatment']


class GonorrheaTreatment(STITreatment):
    """
    Treatment for gonorrhea, with antimicrobial resistance as a strain

    With `amr=True` (which needs `sti.Gonorrhea(amr=True)`), the probability of cure is
    `base_treat_eff` for sensitive infections and `res_treat_eff` for resistant ones.
    Resistance acquired on failure is a step of its own, `acquire`: with probability
    `p_acquire`, an infection that is not cured becomes resistant. Pass it to the sim
    with the treatment. Resistance is passed on by whoever infects someone. Otherwise
    it is `sti.STITreatment` for NG.

    Args:
        base_treat_eff (float): the probability of cure (of a sensitive infection)
        res_treat_eff (float): with `amr`, the probability of cure of a resistant infection
        p_acquire (float): with `amr`, the probability that a treatment failure becomes resistant
        amr (bool): whether NG is stratified by resistance (`sti.Gonorrhea(amr=True)`)
        kwargs (dict): passed to `sti.STITreatment` (e.g. `who=`, `eligibility=`, `name=`)

    Examples:
        ```python
        ng = sti.Gonorrhea(amr=True, init_res=0.01)
        ng_tx = sti.GonorrheaTreatment(amr=True, res_treat_eff=0.3, p_acquire=0.1, who=panel.positive('ng'), name='ng_tx')
        sim = sti.Sim(diseases=ng, interventions=[panel, ng_tx, ng_tx.acquire])
        ```
    """
    def __init__(self, base_treat_eff=0.96, res_treat_eff=0.3, p_acquire=0.0, amr=False, pars=None, name='ng_tx', **kwargs):
        if pars: # As in STIsim 3
            base_treat_eff = pars.pop('base_treat_eff', base_treat_eff)
        eff = ss.by(amr=[base_treat_eff, res_treat_eff]) if amr else base_treat_eff
        super().__init__(diseases='ng', treat_eff=eff, pars=pars, name=name, **kwargs)
        self.acquire = None
        if amr and p_acquire:
            self.acquire = ss.Delivery(move=dict(amr='sens -> res'), who=self.failed, coverage=p_acquire, name=f'{name}_fail')
        return
