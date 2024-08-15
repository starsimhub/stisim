"""
Template class for SEIS-type STIs
Used for chlamydia, gonorrhea, and trich
"""

import numpy as np
import starsim as ss
import sciris as sc

__all__ = ['SEIS']


class SEIS(ss.Infection):

    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.requires = 'structuredsexual'

        self.default_pars(

            # Durations
            dur_exp=ss.constant(1/52),
            dur_inf=[
                ss.lognorm_ex(52/52, 5/52),  # Women
                ss.lognorm_ex(52/52, 5/52),  # Men
            ],

            # Transmission
            beta=1.0,  # Placeholder
            beta_m2f=None,
            beta_f2m=None,
            beta_m2c=None,

            # Symptoms
            p_symp=[
                ss.bernoulli(p=0.375),  # Women
                ss.bernoulli(p=0.375),  # Men
            ],
            p_pid=ss.bernoulli(p=0.2),  # Probability depends on duration of infection
            dur_prepid=ss.lognorm_ex(3/52, 6/52),

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            # Natural history
            ss.BoolArr('exposed'),
            ss.BoolArr('asymptomatic'),
            ss.BoolArr('symptomatic'),
            ss.BoolArr('pid'),
            ss.FloatArr('dur_inf'),
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_symptomatic'),
            ss.FloatArr('ti_pid'),
            ss.FloatArr('ti_clearance'),
        )

        # Results by age and sex
        self.sex_results = sc.objdict()
        self.age_sex_results = sc.objdict()
        self.age_bins = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100])
        self.age_sex_result_keys = [
            'incidence',
            'prevalence',
            'symp_prevalence',
        ]

        return

    @property
    def treatable(self):
        """ Active bacterial presence -- includes exposed and infected, and responds to treatment """
        return self.exposed | self.infected

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.beta_m2f is not None:
            self.pars.beta['structuredsexual'][0] *= self.pars.beta_m2f
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        return

    def init_post(self):
        """ Make initial cases """
        super().init_post()
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        npts = self.sim.npts
        self.results += ss.Result(self.name, 'symp_prevalence', npts, dtype=float, scale=False, label="Symptomatic prevalence")
        self.results += ss.Result(self.name, 'incidence', npts, dtype=float, scale=False, label="Incidence")

        # Age/sex results
        for rkey in self.age_sex_result_keys:
            self.sex_results[rkey] = sc.objdict()
            self.age_sex_results[rkey] = sc.objdict()
            for skey in ['female', 'male', 'both']:
                self.sex_results[rkey][skey] = np.zeros(len(self.sim.yearvec))
                self.age_sex_results[rkey][skey] = np.zeros((len(self.age_bins)-1, len(self.sim.yearvec)))

        return
 
    def clear_infection(self, uids):
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.pid[uids] = False
        self.susceptible[uids] = True
        self.ti_clearance[uids] = self.sim.ti

    def update_pre(self):
        """ Updates prior to interventions """
        ti = self.sim.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Exposed -> infected
        new_infected = (self.exposed & (self.ti_infected <= ti)).uids
        if len(new_infected):
            self.exposed[new_infected] = False
            self.infected[new_infected] = True
            self.ti_infected[new_infected] = ti

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Progress symptoms
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        self.asymptomatic[new_symptomatic] = False
        self.symptomatic[new_symptomatic] = True
        self.ti_symptomatic[new_symptomatic] = ti

        # Progress PID
        new_pid = (self.infected & (self.ti_pid <= ti)).uids
        self.pid[new_pid] = True
        self.ti_pid[new_pid] = ti

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        bins = self.age_bins

        self.results['symp_prevalence'][ti] = self.results['n_symptomatic'][ti] / np.count_nonzero(ppl.alive)
        self.results['incidence'][ti] = self.results['new_infections'][ti] / self.results['n_susceptible'][ti]

        rmap = {'alive': 'both', 'female': 'female', 'male': 'male'}

        # Incidence and prevalence by age and sex
        for pkey, rkey in rmap.items():
            new_inf = ((self.ti_infected == ti) & ppl[pkey]).uids
            new_inf_ages = ppl.age[new_inf]
            n_sus = (self.susceptible & ppl[pkey]).uids
            n_sus_ages = ppl.age[n_sus]
            num, _ = np.histogram(new_inf_ages, bins=bins)
            denom, _ = np.histogram(n_sus_ages, bins=bins)
            self.age_sex_results['incidence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['incidence'][rkey][ti] = sc.safedivide(len(new_inf), len(n_sus))

        # Prevalence by age and sex
        for pkey, rkey in rmap.items():
            n_inf = (self.infected & ppl[pkey]).uids
            n_inf_ages = ppl.age[n_inf]
            num, _ = np.histogram(n_inf_ages, bins=bins)
            denom, _ = np.histogram(ppl.age[ppl[pkey]], bins=bins)
            self.age_sex_results['prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['prevalence'][rkey][ti] = sc.safedivide(len(n_inf), np.count_nonzero(ppl[pkey]))

            n_symp = (self.symptomatic & ppl[pkey]).uids
            n_symp_ages = ppl.age[n_symp]
            num, _ = np.histogram(n_symp_ages, bins=bins)
            self.age_sex_results['symp_prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['symp_prevalence'][rkey][ti] = sc.safedivide(len(n_symp), np.count_nonzero(ppl[pkey]))

        return

    def finalize_results(self):
        super().finalize_results()
        return

    def set_exposure(self, uids):
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.asymptomatic[uids] = True
        self.ti_exposed[uids] = self.sim.ti
        dur_exp = self.pars.dur_exp.rvs(uids)
        self.ti_infected[uids] = self.sim.ti + dur_exp/self.sim.dt
        return

    def set_symptoms(self, p, f_uids, m_uids):
        f_symp = p.p_symp[0].filter(f_uids)
        m_symp = p.p_symp[1].filter(m_uids)
        self.ti_symptomatic[f_symp] = self.ti_infected[f_symp]
        self.ti_symptomatic[m_symp] = self.ti_infected[m_symp]
        return

    def set_duration(self, p, f_uids, m_uids):
        dur_inf_f = p.dur_inf[0].rvs(f_uids)
        dur_inf_m = p.dur_inf[1].rvs(m_uids)
        self.dur_inf[f_uids] = dur_inf_f
        self.dur_inf[m_uids] = dur_inf_m
        return

    def set_pid(self, p, f_uids):
        pid = p.p_pid.filter(f_uids)
        dur_prepid = np.minimum(p.dur_prepid.rvs(pid), self.dur_inf[pid])
        self.ti_pid[pid] = self.ti_infected[pid] + dur_prepid/self.sim.dt
        return

    def wipe_dates(self, uids):
        """ Clear all previous dates """
        self.ti_exposed[uids] = np.nan
        self.ti_infected[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        self.ti_pid[uids] = np.nan
        self.ti_clearance[uids] = np.nan
        return

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected
        """
        super().set_prognoses(uids, source_uids)
        self.wipe_dates(uids)

        ppl = self.sim.people
        p = self.pars
        f_uids = ppl.female.uids.intersect(uids)
        m_uids = ppl.male.uids.intersect(uids)

        self.set_exposure(uids)
        self.set_symptoms(p, f_uids, m_uids)
        self.set_duration(p, f_uids, m_uids)
        self.set_pid(p, f_uids)

        # Determine when people recover
        self.ti_clearance[uids] = self.ti_infected[uids] + self.dur_inf[uids]/self.sim.dt

        return

    def plot(self, key=None, fig=None, style='fancy', fig_kw=None, plot_kw=None):
        """
        Plot all results in the Sim object

        Args:
            key (str): the results key to plot (by default, all)
            fig (Figure): if provided, plot results into an existing figure
            style (str): the plotting style to use (default "fancy"; other options are "simple", None, or any Matplotlib style)
            fig_kw (dict): passed to ``plt.subplots()``
            plot_kw (dict): passed to ``plt.plot()``

        """
        # Configuration
        flat = self.results.flatten()

        n_cols = np.ceil(np.sqrt(len(flat)))  # Number of columns of axes
        default_figsize = np.array([8, 6])
        figsize_factor = np.clip((n_cols-3)/6+1, 1, 1.5) # Scale the default figure size based on the number of rows and columns
        figsize = default_figsize*figsize_factor
        fig_kw = sc.mergedicts({'figsize':figsize}, fig_kw)
        plot_kw = sc.mergedicts({'lw':2}, plot_kw)

        # Do the plotting
        with sc.options.with_style(style):

            yearvec = self.sim.yearvec

            # Get the figure
            fig, axs = sc.getrowscols(len(flat), make=True, **fig_kw)
            if isinstance(axs, np.ndarray):
                axs = axs.flatten()
            else:
                axs = fig.axes
            if not sc.isiterable(axs):
                axs = [axs]

            # Do the plotting
            for ax, (key, res) in zip(axs, flat.items()):
                ax.plot(yearvec, res, **plot_kw, label=self.label)
                title = getattr(res, 'label', key)
                modtitle = res.module
                title = f'{modtitle}: {title}'
                ax.set_title(title)
                ax.set_xlabel('Year')

        sc.figlayout(fig=fig)

        return fig

    def plot(self):
        with sc.options.with_style('fancy'):
            flat = sc.flattendict(self.results, sep=': ')
            yearvec = self.sim.yearvec
            fig, axs = sc.getrowscols(len(flat), make=True)
            for ax, (k, v) in zip(axs.flatten(), flat.items()):
                ax.plot(yearvec, v)
                ax.set_title(k)
                ax.set_xlabel('Year')
        return fig
