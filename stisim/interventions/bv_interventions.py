"""
Define BV interventions
"""

import starsim as ss
import numpy as np
import sciris as sc


# %% Treatment classes
# TODO: figure out whether this can be replaced by STITreatment
__all__ = ["treat_BV"]


class treat_BV(ss.Intervention):
    """
    Treat BV
    """

    def __init__(self, *args, pars=None, eligibility=None, **kwargs):
        super().__init__(*args)

        self.define_pars(
            dt='month',
            # Treatment effects
            start_year=0,  # Day to start treatment
            stop_year=None,  # Day to stop treatment
            p_seek_care=ss.bernoulli(p=0.9),  # Distribution of care-seeking behavior
            tx_delay=ss.weeks(1),  # Delay between symptoms and treatment
            tx_effect_delay=ss.weeks(1),  # Delay between treatment and effect
            tx_doses=14,  # Number of doses of treatment administered (2x/day for 7 days)
            tx_effectiveness=sc.objdict(
                stable_cst1=ss.bernoulli(
                    p=0.95  # Effectiveness of treatment for stably CST 1
                ),
                stable_cst3=ss.bernoulli(
                    p=0.9  # Effectiveness of treatment for stably CST 3
                ),
                stable_cst4=ss.bernoulli(
                    p=0.9  # Effectiveness of treatment for stably CST 4
                ),
            ),
            durable_cst1_change=False,
            rr_douching=0.5,  # Relative risk of stably transitioning to CST 1 due to douching
            rr_poor_menstrual_hygiene=0.5,  # Relative risk of stably transitioning to CST 1 due to poor menstrual hygiene
            baseline_dur_stable_cst_change=ss.months(
                6
            ),  # Duration of stable CST change - based on preliminary LACTIN V results
        ),
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility

        self.define_states(
            ss.FloatArr("ti_treated"),
            ss.FloatArr("orig_CST"),
            ss.BoolArr("on_tx"),
            ss.FloatArr("ti_clearance"),
            ss.FloatArr(
                "ti_return_stable_cst"
            ),  # Timestep to return to original stable CST when post-treatment effect wanes
            ss.BoolArr(
                "post_tx_effect"
            ),  # Flag for ongoing post-treatment effect (durable CST change)
        )

        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result("new_doses", dtype=int, label="New doses administered"),
            ss.Result("new_treated", dtype=int, label="New people treated"),
        ]
        self.define_results(*results)
        return

    def init_pre(self, sim):
        if self.pars.stop_year is None:
            self.pars.stop_year = sim.t.npts - 1
        else:
            ti = sc.findfirst(sim.t.yearvec, self.pars.stop_year)
            self.pars.stop_year = ti

        if self.pars.start_year:
            ti = sc.findfirst(sim.t.yearvec, self.pars.start_year)
            self.pars.start_year = ti
        super().init_pre(sim)
        return

    def init_post(self):
        super().init_post()
        sim = self.sim
        self.orig_CST[:] = sim.diseases.bv.stable_cst

    def check_eligibility(self):
        if self.eligibility is None:
            return (self.sim.diseases.bv.ti_symptomatic == self.ti).uids
        else:
            return self.eligibility(self.sim)

    def tx_cst1_change_duration(self, uids):
        """
        Calculate and set the duration of transitioning to a stable CST 1 based on
        douching and menstrual hygiene practices for individual agents.
        """
        # Baseline duration for  given CST state
        dur_stable_cst_change = np.full_like(
            uids, fill_value=self.pars.baseline_dur_stable_cst_change, dtype=float
        )
        # Adjust duration based on agent-specific characteristics
        douching = self.sim.diseases.bv.douching[uids]
        poor_hygiene = self.sim.diseases.bv.poor_menstrual_hygiene[uids]
        dur_stable_cst_change[
            douching
        ] *= self.pars.rr_douching  # halves the duration of stable CST change
        dur_stable_cst_change[
            poor_hygiene
        ] *= (
            self.pars.rr_poor_menstrual_hygiene
        )  # halves the duration of stable CST change
        return dur_stable_cst_change

    def step(self):
        bv = self.sim.diseases.bv
        if self.pars.stop_year >= self.ti >= self.pars.start_year:
            # Identify eligible agents for treatment
            eligible_inds = self.check_eligibility()
            seeks_care = self.pars.p_seek_care.filter(eligible_inds)
            if len(seeks_care):
                self.ti_treated[seeks_care] = self.ti
                bv.ti_treated[seeks_care] = self.ti
                self.on_tx[seeks_care] = True
                bv.on_tx[seeks_care] = True
                self.results["new_doses"][self.ti] += (
                    len(seeks_care) * self.pars.tx_doses
                )
                self.results["new_treated"][self.ti] += len(seeks_care)

        # Identify agents who have completed treatment
        done_treated = (self.on_tx & (self.ti_treated <= self.ti)).uids
        if len(done_treated):
            self.on_tx[done_treated] = False
            bv.on_tx[done_treated] = False
            self.ti_treated[done_treated] = np.nan

            # Process each CST state
            for stable_cst_val, stable_cst_attr in zip(
                ["stable_cst1", "stable_cst3", "stable_cst4"],
                [bv.stable_cst1, bv.stable_cst3, bv.stable_cst4],
            ):
                treated_cst = done_treated.intersect(stable_cst_attr)
                if len(treated_cst):
                    # Determine treatment effectiveness
                    tx_effective = self.pars.tx_effectiveness[stable_cst_val].filter(
                        treated_cst
                    )
                    if len(tx_effective):
                        # Transition CST states
                        if stable_cst_val == "stable_cst1":
                            bv.cst[tx_effective] = 1
                        else:
                            bv.cst[tx_effective] = 3
                            if self.pars.durable_cst1_change:
                                bv.cst[tx_effective] = 1
                                if len(tx_effective):
                                    self.orig_CST[tx_effective] = bv.stable_cst[
                                        tx_effective
                                    ]
                                    bv.stable_cst[tx_effective] = 1
                                    dur_stable_cst_change = (
                                        self.tx_cst1_change_duration(tx_effective)
                                    )
                                    self.ti_return_stable_cst[tx_effective] = np.round(
                                        self.ti + dur_stable_cst_change
                                    )
                                    self.post_tx_effect[tx_effective] = True

                        # Update symptomatic and clearance times
                        bv.symptomatic[tx_effective] = False
                        self.ti_clearance[tx_effective] = np.nan
                        bv.ti_clearance[tx_effective] = np.nan

                        # Update pregnancy module
                        if "pregnancy" in self.sim.demographics:
                            self.sim.demographics.pregnancy.rel_sus_ptb[
                                tx_effective
                            ] = 1
                            cst_pregnant_uids = tx_effective[
                                self.sim.demographics.pregnancy.pregnant[tx_effective]
                            ]
                            if len(cst_pregnant_uids):
                                self.sim.demographics.pregnancy.set_prognoses(
                                    cst_pregnant_uids, bv_update=True
                                )

            # Return agents to their original CST state
            treatment_waned = (
                self.post_tx_effect & (self.ti_return_stable_cst <= self.ti)
            ).uids
            if len(treatment_waned):
                bv.stable_cst[treatment_waned] = self.orig_CST[treatment_waned]
                self.post_tx_effect[treatment_waned] = False
        return


