"""
Common analyzers for STI analyses
"""
import warnings

# %% Imports and settings
import numpy as np
import sciris as sc
import starsim as ss
import pandas as pd

import stisim as sti
import pylab as pl
from collections import defaultdict, Counter

__all__ = ["result_grouper", "coinfection_stats", "sw_stats", "RelationshipDurations", "NetworkDegree", "DebutAge", "partner_age_diff", "TimeBetweenRelationships", "art_coverage", "PartnershipFormationAnalyzer"]

class result_grouper(ss.Analyzer):
    """Base analyzer providing conditional probability utilities for grouped results.

    Provides a ``cond_prob`` static method that computes the proportion of a
    numerator group within a denominator group. Intended as a base class for
    analyzers that compute stratified prevalence or infection statistics.
    """

    @staticmethod
    def cond_prob(numerator, denominator):
        numer = len((denominator & numerator).uids)
        denom = len(denominator.uids)
        out = sc.safedivide(numer, denom)
        return out


class coinfection_stats(result_grouper):
    """
    Generates stats for the coinfection of two diseases.
    This is useful for looking at the coinfection of HIV and syphilis, for example.

    Args:
        disease1 (str | ss.Disease): name of the first disease
        disease2 (str | ss.Disease): name of the second disease
        disease1_infected_state_name (str): name of the infected state for disease1 (default: 'infected')
        disease2_infected_state_name (str): name of the infected state for disease2 (default: 'infected')
        age_limits (list): list of two integers that define the age limits for the denominator.
        denom (function): function that returns a boolean array of the denominator, usually the relevant population.
            default: lambda self: (self.sim.people.age >= 15) & (self.sim.people.age < 50)
        *args, **kwargs : optional, passed to ss.Analyzer constructor
    """
    def __init__(self, disease1, disease2, disease1_infected_state_name='infected', disease2_infected_state_name='infected',
                 age_limits=None, denom=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if disease1 is None or disease2 is None:
            raise ValueError('Coinfection stats requires exactly 2 diseases')

        self.disease1 = disease1
        self.disease2 = disease2

        # if the diseases are objects, get their names and store them instead of the objects
        if isinstance(self.disease1, ss.Disease):
            self.disease1 = self.disease1.name
        if isinstance(self.disease2, ss.Disease):
            self.disease2 = self.disease2.name

        self.disease1_infected_state_name = disease1_infected_state_name
        self.disease2_infected_state_name = disease2_infected_state_name
        self.age_limits = age_limits or [15, 50]
        default_denom = lambda self: (self.sim.people.age >= self.age_limits[0]) & (self.sim.people.age < self.age_limits[1])
        self.denom = denom or default_denom

        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result(f'{self.disease1}_prev_no_{self.disease2}', dtype=float, scale=False),
            ss.Result(f'{self.disease1}_prev_has_{self.disease2}', dtype=float, scale=False),
            ss.Result(f'{self.disease1}_prev_no_{self.disease2}_f', dtype=float, scale=False),
            ss.Result(f'{self.disease1}_prev_has_{self.disease2}_f', dtype=float, scale=False),
            ss.Result(f'{self.disease1}_prev_no_{self.disease2}_m', dtype=float, scale=False),
            ss.Result(f'{self.disease1}_prev_has_{self.disease2}_m', dtype=float, scale=False),
        ]
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti = self.ti
        disease1name = self.disease1
        disease2name = self.disease2
        disease1obj = getattr(self.sim.diseases, self.disease1)
        disease2obj = getattr(self.sim.diseases, self.disease2)

        ppl = sim.people

        denom = self.denom(self)
        has_disease2 = getattr(disease2obj, self.disease2_infected_state_name)  # Adults with HIV
        has_disease1 = getattr(disease1obj, self.disease1_infected_state_name)  # Adults with syphilis

        has_disease1_f = denom & has_disease1 & ppl.female  # Women with dis1
        has_disease1_m = denom & has_disease1 & ppl.male  # Men with dis1
        has_disease2_f = denom & has_disease2 & ppl.female  # Women with dis2
        has_disease2_m = denom & has_disease2 & ppl.male  # Men with dis2
        no_disease2    = denom & (~has_disease2)  # Adults without dis2
        no_disease2_f  = no_disease2 & ppl.female  # Women without dis2
        no_disease2_m  = no_disease2 & ppl.male  # Men without dis2

        self.results[f'{disease1name}_prev_no_{disease2name}'][ti] = self.cond_prob(has_disease1, no_disease2)
        self.results[f'{disease1name}_prev_has_{disease2name}'][ti] = self.cond_prob(has_disease1, has_disease2)
        self.results[f'{disease1name}_prev_no_{disease2name}_f'][ti] = self.cond_prob(has_disease1_f, no_disease2_f)
        self.results[f'{disease1name}_prev_has_{disease2name}_f'][ti] = self.cond_prob(has_disease1_f, has_disease2_f)
        self.results[f'{disease1name}_prev_no_{disease2name}_m'][ti] = self.cond_prob(has_disease1_m, no_disease2_m)
        self.results[f'{disease1name}_prev_has_{disease2name}_m'][ti] = self.cond_prob(has_disease1_m, has_disease2_m)

        return


class sw_stats(result_grouper):
    """Track new infections and transmissions among sex workers and their clients.

    At each timestep, records the number and share of new infections and
    transmissions attributable to female sex workers (FSW), clients, and
    non-sex-worker populations, disaggregated by disease.

    Args:
        diseases (list): List of disease names (str) to track, e.g. ``['hiv', 'syphilis']``.
    """

    def __init__(self, diseases=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'sw_stats'
        self.diseases = diseases
        return

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        for d in self.diseases:
            results += [
                ss.Result('share_new_infections_fsw_'+d, scale=False, summarize_by='mean', auto_plot=False),
                ss.Result('share_new_infections_client_'+d,scale=False, summarize_by='mean', auto_plot=False),
                ss.Result('new_infections_fsw_'+d, dtype=int, auto_plot=False),
                ss.Result('new_infections_client_'+d, dtype=int, auto_plot=False),
                ss.Result('new_infections_non_fsw_'+d, dtype=int, auto_plot=False),
                ss.Result('new_infections_non_client_'+d, dtype=int, auto_plot=False),
                ss.Result('new_transmissions_fsw_'+d, dtype=int, auto_plot=False),
                ss.Result('new_transmissions_client_'+d, dtype=int, auto_plot=False),
                ss.Result('new_transmissions_non_fsw_'+d, dtype=int, auto_plot=False),
                ss.Result('new_transmissions_non_client_'+d, dtype=int, auto_plot=False),
            ]
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti = self.ti

        if ti > 0:

            for d in self.diseases:
                dis = sim.diseases[d]
                nw = sim.networks.structuredsexual
                adult = sim.people.age > 0
                fsw = nw.fsw & adult
                client = nw.client & adult
                non_fsw = sim.people.female & ~nw.fsw & adult
                non_client = sim.people.male & ~nw.client & adult
                newly_infected = (dis.ti_exposed == ti) & adult
                new_trans = dis.ti_transmitted_sex == ti
                total_acq = len(newly_infected.uids)

                newly_transmitting_fsw = (dis.ti_transmitted_sex == ti) & fsw
                newly_transmitting_clients = (dis.ti_transmitted_sex == ti) & client
                newly_transmitting_non_fsw = (dis.ti_transmitted_sex == ti) & non_fsw
                newly_transmitting_non_client = (dis.ti_transmitted_sex == ti) & non_client

                new_transmissions_fsw = dis.new_transmissions_sex[newly_transmitting_fsw]
                new_transmissions_client = dis.new_transmissions_sex[newly_transmitting_clients]
                new_transmissions_non_fsw = dis.new_transmissions_sex[newly_transmitting_non_fsw]
                new_transmissions_non_client = dis.new_transmissions_sex[newly_transmitting_non_client]

                self.results['share_new_infections_fsw_'+d][ti] = self.cond_prob(fsw, newly_infected)
                self.results['share_new_infections_client_'+d][ti] = self.cond_prob(client, newly_infected)

                self.results['new_infections_fsw_'+d][ti] = len((fsw & newly_infected).uids)
                self.results['new_infections_client_'+d][ti] = len((client & newly_infected).uids)
                self.results['new_infections_non_fsw_'+d][ti] = len((non_fsw & newly_infected).uids)
                self.results['new_infections_non_client_'+d][ti] = len((non_client & newly_infected).uids)

                self.results['new_transmissions_fsw_'+d][ti] = sum(new_transmissions_fsw)
                self.results['new_transmissions_client_'+d][ti] = sum(new_transmissions_client)
                self.results['new_transmissions_non_fsw_'+d][ti] = sum(new_transmissions_non_fsw)
                self.results['new_transmissions_non_client_'+d][ti] = sum(new_transmissions_non_client)

                # Note: total_trans only counts sexual transmission, not MTC (congenital).
                # So we don't validate against newly_infected which includes MTC.

        return


class RelationshipDurations(ss.Analyzer):
    """Analyze relationship durations in a StructuredSexual network.

    Records the mean and median duration of all relationships (stable, casual,
    etc.) at each timestep, disaggregated by sex. Durations are extracted from
    the network's ``relationship_durs`` tracking dict.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('mean_duration', dtype=float, scale=False, auto_plot=False),
            ss.Result('median_duration', dtype=float, scale=False, auto_plot=False),
        )
        return

    def step(self):
        pass

    def update_results(self):
        sim = self.sim
        ti = self.ti
        nw = sim.networks.structuredsexual
        rel_durations = self.get_relationship_durations()
        self.results['mean_duration'][ti] = np.mean(rel_durations)
        self.results['median_duration'][ti] = np.median(rel_durations)
        return

    def plot(self):
        sim = self.sim
        ti = self.ti
        female_relationship_durs, male_relationship_durs = self.get_relationship_durations()
        all_durations = female_relationship_durs + male_relationship_durs
        pl.figure(1)
        pl.hist(all_durations, bins=(max(all_durations) - min(all_durations)))
        pl.xlabel('Relationship Duration')
        pl.ylabel('Frequency')
        pl.title('Distribution of Relationship Durations')

        pl.figure(2)
        pl.hist(female_relationship_durs, bins=(max(all_durations) - min(all_durations)))
        pl.xlabel('Female Relationship Duration')
        pl.ylabel('Frequency')
        pl.title('Distribution of Female Relationship Durations')

        pl.figure(3)
        pl.hist(male_relationship_durs, bins=(max(all_durations) - min(all_durations)))
        pl.xlabel('Male Relationship Duration')
        pl.ylabel('Frequency')
        pl.title('Distribution of Male Relationship Durations')
        pl.show()


    def get_relationship_durations(self):
        """
        Returns the durations of all relationships, separated by sex.

        If include_current is False, return the duration of only relationships that have ended

        returns:
            female_durations: list of durations of relationships
            male_durations: list of durations of relationships
        """

        # Get the current duration of all relationships
        male_durations = []
        female_durations = []
        for pair, relationships in self.sim.networks.structuredsexual.relationship_durs.items():
            durs = [relationship['dur'] for relationship in relationships]

            # assign the durations to male and female lists
            for uid in pair:
                if self.sim.people.female[uid]:
                    female_durations.extend(durs)
                else:
                    male_durations.extend(durs)

        return female_durations, male_durations


class NetworkDegree(ss.Analyzer):
    """Analyze lifetime partner count distributions in a StructuredSexual network.

    At a specified year, records the number of lifetime partners per agent,
    disaggregated by sex and relationship type (stable, casual, one-time, sex
    work). Results are binned into a histogram for plotting and comparison
    with survey data.

    Args:
        year (float): Calendar year at which to record partner counts. Defaults
            to the last year of the simulation.
        bins (array): Bin edges for the partner count histogram. Defaults to
            ``[0, 1, ..., 20, 100]``.
        relationship_types (list): Relationship types to track, e.g.
            ``['stable', 'casual']``. Use ``'partners'`` for combined
            stable + casual counts.
    """

    def __init__(self, year=None, bins=None, relationship_types=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.year = year

        if bins is None:
            bins = np.concatenate([np.arange(21),[100]])
        self.bins = bins

        if relationship_types is None:
            relationship_types = ['stable', 'casual'] # Other options are 'partners' (stable+casual), 'onetime', 'sw'

        self.relationship_types = []
        if 'partners' in relationship_types:
            relationship_types.remove('partners')
            self.relationship_types.append('lifetime_partners')
        [self.relationship_types.append(f'lifetime_{relationship_type}_partners') for relationship_type in relationship_types]

        for relationship_type in self.relationship_types:
            setattr(self, f'{relationship_type}_f', [])
            setattr(self, f'{relationship_type}_m', [])

        return

    def init_results(self):
        """
        Add results for `n_rships`, separated for males and females
        Optionally disaggregate for risk level / age?
        """
        super().init_results()
        for relationship_type in self.relationship_types:
            self.results += [
                ss.Result(f'{relationship_type}_f', dtype=int, scale=False, shape=len(self.bins), auto_plot=False),
                ss.Result(f'{relationship_type}_m', dtype=int, scale=False, shape=len(self.bins), auto_plot=False),
            ]
        return

    def init_pre(self, sim, **kwargs):
        """
        Initialize the analyzer
        """
        super().init_pre(sim, **kwargs)
        self.year = sim.t.yearvec[-1]
        return


    def step(self):
        """
        record lifetime_partners for the user-specified year
        """
        if self.sim.t.yearvec[self.ti] == self.year:
            for relationship_type in self.relationship_types:
                # Get the number of partners, disaggregated by sex. We can't use a Result object for this because we
                # don't know how many agents there will be at any given time step. We can use Results for the binned
                # counts.

                female_partners = getattr(self.sim.networks.structuredsexual, relationship_type)[self.sim.people.female]
                male_partners = getattr(self.sim.networks.structuredsexual, relationship_type)[self.sim.people.male]

                getattr(self, f'{relationship_type}_f').extend(female_partners)
                getattr(self, f'{relationship_type}_m').extend(male_partners)

                # bin the data by number of partners
                female_counts, female_bins = np.histogram(female_partners, bins=self.bins)
                male_counts, male_bins = np.histogram(male_partners, bins=self.bins)

                for i, female_count, male_count in zip(range(len(self.bins)), female_counts, male_counts):
                    self.results[f'{relationship_type}_f'][i] = female_count
                    self.results[f'{relationship_type}_m'][i] = male_count
        return

    def plot(self):
        """
        Plot histograms and stats by sex and relationship type
        """

        for relationship_type in self.relationship_types:
            fig, axes = pl.subplots(1, 2, figsize=(9, 5), layout="tight")
            axes = axes.flatten()
            for ai, sex in enumerate(['f', 'm']):
                counts = self.results[f'{relationship_type}_{sex}'].values
                bins=self.bins

                total = sum(counts)
                counts = counts / total
                counts[-2] = counts[-2:].sum()
                counts = counts[:-1]

                axes[ai].bar(bins[:-1], counts)
                axes[ai].set_xlabel(f'Number of {relationship_type}')
                axes[ai].set_title(f'Distribution of partners, {sex}')
                axes[ai].set_ylim([0, 1])

                sex_counts = np.array(getattr(self, f'{relationship_type}_{sex}'))
                stats = f"Mean: {np.mean(sex_counts):.1f}\n"
                stats += f"Median: {np.median(sex_counts):.1f}\n"
                stats += f"Std: {np.std(sex_counts):.1f}\n"
                stats += f"%>20: {np.count_nonzero(sex_counts >= 20) / total * 100:.2f}\n"
                axes[ai].text(15, 0.5, stats)

            pl.show()

        return

class TimeBetweenRelationships(ss.Analyzer):
    """
    Analyzes the time between relationships in a structuredsexual network.
    Each timestep, for each debuted agent, check if they are in a relationship of the provided type.
    If not, increment the counter
    Otherwise, reset the counter to 0 and append the counter to the list of times between relationships for that agent.
    """
    def __init__(self, relationship_type='stable', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.relationship_type = relationship_type
        self.results['times_between_relationships'] = defaultdict(self.zero_list)
        return

    @staticmethod
    def zero_list():
        return [0]  # Initialize the list with a single zero


    def step(self):
        """
        For each debuted agent, check if they are in a relationship.
        If they are not, increment the time since last relationship by 1.
        If they are and time since last relationship is greater than 0, append the time to the list of times between relationships.
        """

        sim = self.sim
        ti = self.ti
        nw = sim.networks.structuredsexual

        debuted = (nw.debut <= sim.people.age).uids
        stable_relationships = nw.edges.edge_type == nw.edge_types[self.relationship_type]  # Get stable relationships
        p1_stable = nw.edges.p1[stable_relationships]  # Get p1s in stable relationships
        p2_stable = nw.edges.p2[stable_relationships]  # Get p2s in stable relationships

        stable_relationship_uids = set(np.concatenate([p1_stable, p2_stable]))  # Combine p1 and p2 uids in stable relationships
        not_stable_relationship_uids = set(debuted) - stable_relationship_uids  # Get uids not in stable relationships

        for uid in not_stable_relationship_uids:
            self.results['times_between_relationships'][uid][-1] += 1  # Increment time since last relationship for those not in stable relationships

        for uid in stable_relationship_uids:
            if self.results['times_between_relationships'][uid][-1] > 0:
                # If the agent is in a stable relationship, reset the time since last relationship and append to the list
                self.results['times_between_relationships'][uid].append(0)

        return

      
class partner_age_diff(ss.Analyzer):
    """Analyze age differences between sexual partners.

    Records the mean, median, and standard deviation of male-female age
    differences (male age minus female age) at each timestep. At a specified
    year, stores full age-difference distributions disaggregated by female
    age group for detailed plotting.

    Args:
        year (float): Calendar year at which to store detailed age-difference
            distributions. Defaults to 2000.
        age_bins (list): Female age group names matching the network's
            ``f_age_group_bins`` keys. Defaults to ``['teens', 'young', 'adult']``.
        network (str): Name of the network to analyze. Defaults to
            ``'structuredsexual'``.
    """

    def __init__(self, year=2000, age_bins=['teens', 'young', 'adult'], network='structuredsexual', *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.year = year
        self.network = network
        self.age_diffs = {}
        self.age_bins = age_bins
        return

    def init_results(self):
        """
        Initialize the results for the age differences.
        """
        super().init_results()
        self.define_results(
            ss.Result('age_diff_mean', dtype=float, scale=False, auto_plot=False),
            ss.Result('age_diff_median', dtype=float, scale=False, auto_plot=False),
            ss.Result('age_diff_std', dtype=float, scale=False, auto_plot=False),
        )
        return

    def step(self):
        """
        Record the age differences between partners in the specified year.
        """

        net = self.sim.networks[self.network]
        relationships = net.edges.dur > 1
        p1 = net.p1[relationships]
        p2 = net.p2[relationships]

        age_diffs = (self.sim.people.age[p1] - self.sim.people.age[p2])

        f_ages = self.sim.people.age[p2]

        # bin the female ages by the bins used in the structured sexual network
        # age_bins = sorted([bin[0] for bin in self.sim.networks.structuredsexual.pars.f_age_group_bins.values()])
        age_bin_limits = [net.pars.f_age_group_bins[bin][0] for bin in self.age_bins]
        age_bin_indices = np.digitize(f_ages, age_bin_limits) - 1

        self.results['age_diff_mean'][self.ti] = np.mean(age_diffs)
        self.results['age_diff_median'][self.ti] = np.median(age_diffs)
        self.results['age_diff_std'][self.ti] = np.std(age_diffs)

        if self.sim.t.yearvec[self.ti] == self.year:
            for bin in self.age_bins:
                self.age_diffs[bin] = age_diffs[age_bin_indices == self.age_bins.index(bin)]

    def plot(self):
        """
        Plot histograms of the age differences between partners.
        """
        if len(self.age_diffs) > 0:
            pl.figure(figsize=(8, 5))
            pl.hist(list(self.age_diffs.values()), label=list(self.age_diffs.keys()), bins=30, edgecolor='black', alpha=0.7)
            pl.legend()
            pl.xlabel('Age Difference (years)')
            pl.ylabel('Frequency')
            pl.title(f'Age Differences Between Partners in {self.year} (Male Age - Female Age)')
            pl.grid(True)
            pl.show()
        else:
            print("No age differences recorded for the specified year.")

        return


class DebutAge(ss.Analyzer):
    """Analyze the proportion of agents who have sexually debuted by age.

    Tracks the share of agents who are sexually active (past their debut age)
    at each single-year age bin, disaggregated by sex and birth cohort. Useful
    for validating debut age distributions against survey data such as DHS.

    Args:
        bins (array): Age bins to evaluate, e.g. ``np.arange(12, 31)``.
            Defaults to ages 12-30.
        cohort_starts (array): Birth-year cohort start years. Defaults to all
            cohorts that fit within the simulation timespan.
    """
    def __init__(self, bins=None, cohort_starts=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.bins = bins or np.arange(12, 31, 1)
        self.binspan = self.bins[-1] - self.bins[0]
        self.cohort_starts = cohort_starts

        return

    def init_pre(self, sim, force=False):
        if self.cohort_starts is None:
            first_cohort = sim.t.start.years
            last_cohort = sim.t.stop.years - self.binspan
            self.cohort_starts = sc.inclusiverange(first_cohort, last_cohort)
            self.cohort_ends = self.cohort_starts + self.binspan
            self.n_cohorts = len(self.cohort_starts)
            self.cohort_years = np.array([sc.inclusiverange(i, i + self.binspan) for i in self.cohort_starts])

        self.prop_active_f = np.zeros((self.n_cohorts, self.binspan + 1))
        self.prop_active_m = np.zeros((self.n_cohorts, self.binspan + 1))
        super().init_pre(sim, force=force)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('prop_active_f', dtype=float, scale=False, shape=len(self.cohort_starts), auto_plot=False),
            ss.Result('prop_active_m', dtype=float, scale=False, shape=len(self.cohort_starts), auto_plot=False),
        )
        return

    def step(self):

        sim = self.sim
        ppl = sim.people
        if sim.t.yearvec[sim.ti] in self.cohort_years:
            cohort_inds, bin_inds = sc.findinds(self.cohort_years, sim.t.yearvec[sim.ti])
            for ci, cohort_ind in enumerate(cohort_inds):
                bin_ind = bin_inds[ci]
                bin = self.bins[bin_ind]

                # all females cohort:
                conditions_f = ppl.female * ppl.alive * (ppl.age >= (bin - 1)) * (
                        ppl.age < bin)
                cohort_f_count = sum(conditions_f)
                # all active females in cohort:
                num_conditions_f = conditions_f * sim.networks.structuredsexual.over_debut
                debut_f_count = sum(num_conditions_f)

                self.prop_active_f[cohort_ind, bin_ind] = (debut_f_count) / (cohort_f_count) if cohort_f_count > 0 else 0

                # all males cohort:
                conditions_m = ~sim.people.female * sim.people.alive * (sim.people.age >= (bin - 1)) * (
                        sim.people.age < bin)
                cohort_m_count = sum(conditions_m)
                # all active males in cohort:
                num_conditions_m = conditions_m * sim.networks.structuredsexual.over_debut
                debut_m_count = sum(num_conditions_m)
                self.prop_active_m[cohort_ind, bin_ind] = (debut_m_count) / (cohort_m_count) if cohort_m_count > 0 else 0
        return

    def plot(self):
        """
        Plot the proportion of active agents by cohort and debut age
        """
        pl.figure(1)
        for row in self.prop_active_f:
            pl.plot(self.bins, row)
        pl.xlabel('Age')
        pl.ylabel('Share')
        pl.title('Proportion of females who are sexually active')
        pl.show()

        pl.figure(2)
        for row in self.prop_active_m:
            pl.plot(self.bins, row)
        pl.xlabel('Age')
        pl.ylabel('Share')
        pl.title('Proportion of males who are sexually active')
        pl.show()


class art_coverage(ss.Analyzer):
    """
    Track ART coverage (number and proportion) by sex and age bin.

    Results are stored as time series per stratum, accessible via:
        analyzer.results['n_art_f_15_25']    # Women 15-25 on ART (count)
        analyzer.results['p_art_m_25_35']    # Men 25-35 on ART (proportion of infected)

    Args:
        age_bins (list): age bin edges, e.g. [15, 25, 35, 45, 65]. Default: [15, 25, 35, 45, 65].
            Bins are half-open intervals: [lo, hi), i.e. lo <= age < hi.
    """

    def __init__(self, age_bins=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'art_coverage'
        self.age_bins = age_bins or [15, 25, 35, 45, 65]
        return

    def init_results(self):
        super().init_results()
        results = []

        # Aggregate results
        results.append(ss.Result('n_art',   dtype=int,   label='Total on ART'))
        results.append(ss.Result('p_art',   scale=False, label='ART coverage (proportion of infected)'))
        results.append(ss.Result('n_art_f', dtype=int,   label='Women on ART'))
        results.append(ss.Result('n_art_m', dtype=int,   label='Men on ART'))
        results.append(ss.Result('p_art_f', scale=False, label='ART coverage (women)'))
        results.append(ss.Result('p_art_m', scale=False, label='ART coverage (men)'))

        # Per age bin × sex
        for i in range(len(self.age_bins) - 1):
            lo, hi = self.age_bins[i], self.age_bins[i + 1]
            for sex in ['f', 'm']:
                results.append(ss.Result(f'n_art_{sex}_{lo}_{hi}', dtype=int, label=f'On ART {sex.upper()} {lo}-{hi}'))
                results.append(ss.Result(f'p_art_{sex}_{lo}_{hi}', scale=False, label=f'ART coverage {sex.upper()} {lo}-{hi}'))

        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti  = self.ti
        ppl = sim.people
        hiv = sim.diseases.hiv

        on_art   = hiv.on_art
        infected = hiv.infected
        female   = ppl.female
        male     = ppl.male
        age      = ppl.age

        # Aggregate
        n_art = len(on_art.uids)
        n_inf = len(infected.uids)
        self.results['n_art'][ti]   = n_art
        self.results['p_art'][ti]   = sc.safedivide(n_art, n_inf)
        self.results['n_art_f'][ti] = len((on_art & female).uids)
        self.results['n_art_m'][ti] = len((on_art & male).uids)
        self.results['p_art_f'][ti] = sc.safedivide(len((on_art & female).uids), len((infected & female).uids))
        self.results['p_art_m'][ti] = sc.safedivide(len((on_art & male).uids), len((infected & male).uids))

        # Per age bin × sex
        for i in range(len(self.age_bins) - 1):
            lo, hi = self.age_bins[i], self.age_bins[i + 1]
            age_mask = (age >= lo) & (age < hi)
            for sex, sex_mask in [('f', female), ('m', male)]:
                stratum     = age_mask & sex_mask
                n_on        = len((on_art & stratum).uids)
                n_inf_strat = len((infected & stratum).uids)
                self.results[f'n_art_{sex}_{lo}_{hi}'][ti] = n_on
                self.results[f'p_art_{sex}_{lo}_{hi}'][ti] = sc.safedivide(n_on, n_inf_strat)

        return

    def plot(self, by_age=True):
        """
        Plot ART coverage over time.

        Creates a 2-panel figure: aggregate coverage (left) and by age/sex (right).
        If by_age=False, only plots aggregate.

        Example::

            sim.run()
            sim.analyzers.art_coverage.plot()
        """
        yearvec = self.sim.t.yearvec
        sex_colors = {'f': '#d46e9c', 'm': '#4a90d9'}
        sex_labels = {'f': 'Women', 'm': 'Men'}

        if by_age and len(self.age_bins) > 2:
            fig, axes = pl.subplots(1, 3, figsize=(16, 5))
        else:
            fig, axes = pl.subplots(1, 1, figsize=(6, 5))
            axes = [axes]

        # Panel 1: aggregate
        ax = axes[0]
        ax.plot(yearvec, self.results['p_art'],   color='k',              linewidth=2, label='Overall')
        ax.plot(yearvec, self.results['p_art_f'], color=sex_colors['f'], linewidth=1.5, label='Women')
        ax.plot(yearvec, self.results['p_art_m'], color=sex_colors['m'], linewidth=1.5, label='Men')
        ax.set_xlabel('Year')
        ax.set_ylabel('ART coverage (proportion of infected)')
        ax.set_title('ART coverage')
        ax.set_ylim(0, 1)
        ax.legend(frameon=False)

        # Panels 2-3: by age (women, men)
        if by_age and len(self.age_bins) > 2 and len(axes) > 1:
            for sex, ax in zip(['f', 'm'], axes[1:]):
                for i in range(len(self.age_bins) - 1):
                    lo, hi = self.age_bins[i], self.age_bins[i + 1]
                    key = f'p_art_{sex}_{lo}_{hi}'
                    ax.plot(yearvec, self.results[key], label=f'{lo}-{hi}')
                ax.set_xlabel('Year')
                ax.set_ylabel('ART coverage')
                ax.set_title(f'{sex_labels[sex]} by age')
                ax.set_ylim(0, 1)
                ax.legend(frameon=False, fontsize=9)

        pl.tight_layout()
        return fig


class PartnershipFormationAnalyzer(ss.Analyzer):
    """
    Track partnership formation per network, gender, and age bin.

    Works with any analyzer that subclasses stisim BaseNetwork (that does not have class attribute
    records_all_expirations set to False on purpose. This attribute (when False) explicitly tells
    PartnershipFormationAnalyzer the network is not configured to work with it (yet)).

    Stores **one row per distinct edge** with its active interval, rather than
    re-snapshotting the active edge set every timestep. Formation is detected
    at the step where ``ti_formed == ti``; the end of the interval comes from
    the network's ``expired_this_loop`` (populated by ``_on_edge_dissolution`` when
    ``record_expired=True``), so ``step()`` only touches the per-step *flows*
    (newly-formed and newly-expired edges), never the active *stock*. Each
    partner's sex and age-at-formation are read from the edge/people directly, so
    the analyzer makes no p1=male/p2=female assumption and works for any gender
    combination (male-female, male-male).

    Recorded state lives in ``self.results['relationships'][nw]``: a list of one mutable
    row per distinct edge::

        [p1, p2, ti_formed, ti_expired, age_p1, age_p2, female_p1, female_p2]

    ``ti_expired`` is the timestep at which the edge was observed in
    ``expired_this_loop`` (so its active interval is the half-open
    ``[ti_formed, ti_expired)``), or ``None`` while the edge is still active
    at the end of the run. Getters map ``None -> final_ti + 1`` so an ongoing
    edge tests active through the last timestep.

    .. warning::
        **Single-edge-per-pair-per-timestep assumption.** Edges are keyed by
        ``(p1, p2, ti_formed)`` (see ``_indicies_of_rels_in_table_for_nw``). If a network ever forms
        **more than one edge between the same pair of agents at the same ti**,
        those edges collide on this key and **expiry data is corrupted**: only
        the last-formed colliding edge remains reachable for ``ti_expired``
        recording, so the others keep ``ti_expired = None`` (wrongly treated as
        active through run end) and their expiry events are dropped. Formation
        *counts* are unaffected (every edge is still appended as its own row), so
        only the active-interval results (``get_unique_partners_active_per_agent``
        and any duration computed from ``ti_expired``) are affected. The bundled
        networks do not do this; the assumption holds as long as a pair forms at
        most one edge per network per timestep (the no-concurrent-duplicate-edge
        invariant noted in ``BaseNetwork.add_pairs``). Forming a second edge
        between the same pair at a *different* ti is fine (distinct key).

    **Duration semantics.** ``ti_expired`` is defined as one past the last
    timestep the edge was present and transmission-capable during the disease
    transmission step (step 9 of the 16-step integration loop). Any duration
    reported or computed from this analyzer is therefore::

        duration = ti_expired - ti_formed
                 = number of active transmission steps the relationship existed

    where an "active transmission step" means the edge was present in the
    network with **both partners alive at step 9** (i.e. structurally available
    to transmit) -- not that a transmission event actually fired (that further
    depends on beta, acts, condom use, and infection status). This count is
    consistent regardless of when a relationship formed, when it dissolved, or
    why it dissolved:

      - duration expiry at R: active ``{F..R-1}`` -> ``duration = R - F``;
      - partner death/removal at T: the dying agent's ``alive`` flag flips at
        step 10 (*after* step 9), so the edge is still active at T -> active
        ``{F..T}``, ``ti_expired = T+1``, ``duration = (T+1) - F``;
      - formed and ended-by-death in the same timestep T: active ``{T}`` ->
        ``duration = 1`` (never a zero/empty interval);
      - still active at run end: active ``{F..final_ti}`` ->
        ``duration = (final_ti+1) - F`` (right-censored at the run boundary).

    Note this realized duration can be shorter than the network's drawn ``dur``
    parameter (the *intended* duration) when a partner dies or leaves early; the
    analyzer records the realized count, not the drawn one.

    Example (``self.results`` layout, network 'mfnetwork')::

        ana.results['relationships']['mfnetwork'][0]
        # -> [9, 17, 0, 4, 27.23, 19.62, False, True]
        #    p1=9, p2=17, ti_formed=0, ti_expired=4 (active ti 0..3),
        #    ages at formation, female flags. ti_expired is None if still active
        #    at run end.

    Reporting is via getters (all return ``{network_name: result}``; ``female``
    selects the subject sex):

      - ``get_n_partnerships_formed(female, age_bins=None, window_months=None)``
        -- ``{nw: {age_bin: array}}``. With ``window_months=None`` the inner
        array is the full per-timestep series (length ``n_ti``); with an int it
        is a single trailing-window sum (length 1).
      - ``get_n_partnerships_formed_per_agent(female, window_months=None)`` --
        ``{nw: array}`` of non-unique formation counts per agent over the window.
      - ``get_unique_partners_active_per_agent(female, window_months=None)`` --
        ``{nw: array}`` of unique partners per agent on edges active during the
        window (half-open interval test).

    Args:
        age_bins (list[str]): Age bin strings parsable by ``ss.parse_age_range``
            (e.g. ``['0-5', '5-10', ..., '70-75']``). Bins may overlap; a
            partner is counted in every bin containing its age. Defaults to
            five-year bins from ``'0-5'`` through ``'70-75'``.
        networks (list[str] | None): Network names to track. ``None`` (default)
            tracks every stisim ``BaseNetwork``-derived sexual network in the
            sim. Named missing networks raise ``ValueError``; named networks that
            are not ``BaseNetwork``-derived warn and are skipped. A network that
            cannot record all edge expirations (``records_all_expirations`` is
            False, e.g. ``MSMScaleFreeNetwork``) is **skipped with a warning**,
            since its active-interval results would be incomplete. The analyzer
            requires ``BaseNetwork`` because it reads BaseNetwork-only edge meta
            (``ti_formed``, ``age_p1``, ``age_p2``) and the ``expired_this_loop``
            buffer, and sets ``record_expired = True`` on each tracked network.
    """

    DEFAULT_AGE_BINS = [f'{lo}-{lo+5}' for lo in range(0, 75, 5)]

    # Edge-row column indices (rows stored in self.results['relationships'][nw]).
    _P1, _P2, _FORMATION_TI, _TI_EXPIRED, _AGE1, _AGE2, _FEM1, _FEM2 = range(8)

    def __init__(self, age_bins: list[str] = None, networks: list[str] = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.age_bins = list(age_bins) if age_bins is not None else list(self.DEFAULT_AGE_BINS)

        # Age-bin strings are used as result dict keys, so duplicates are not
        # allowed (they would collide). Detect via multiset-minus-set difference:
        # Counter(all occurrences) - Counter(one of each unique) leaves only the
        # surplus occurrences, whose keys are exactly the duplicated values.
        dups = sorted((Counter(self.age_bins) - Counter(set(self.age_bins))).keys())
        if len(dups) > 0:
            raise ValueError(f"PartnershipFormationAnalyzer: duplicate age_bins are not allowed: {dups}")

        self.age_ranges = [ss.parse_age_range(b) for b in self.age_bins]
        self.target_nw_names = networks
        self._tracked_nw_names = []                     # populated in init_results
        self._skipped_nw_names = []                  # networks skipped (can't record expirations)
        self._indicies_of_rels_in_table_for_nw = {}                        # nw -> {(p1,p2,ti_formed): row_idx}, pruned on expiry
        self._final_uids = {'f': None, 'm': None}    # surviving-agent uids by sex, set in step
        self._final_ti = None
        return

    def init_results(self):
        super().init_results()

        # Validate passed-in network names for tracking. The analyzer requires a
        # stisim BaseNetwork-derived sexual network (it reads BaseNetwork-only
        # edge meta -- ti_formed, age_p1, age_p2 -- and the expired_this_loop
        # buffer), so a plain ss.SexualNetwork is not enough.
        all_nw_names = self.sim.networks.keys()
        if self.target_nw_names is None:
            # default to all BaseNetwork-derived sexual networks
            candidates = [nw_name for nw_name in all_nw_names
                          if isinstance(self.sim.networks[nw_name], sti.BaseNetwork)]
        else:
            # ensure all requested networks exist
            missing = [nw_name for nw_name in self.target_nw_names if nw_name not in all_nw_names]
            if len(missing) > 0:
                missing_str = ', '.join(missing)
                raise ValueError(f"PartnershipFormationAnalyzer: network(s) '{missing_str}' not found in sim.networks. "
                                 f"Available: {sorted(all_nw_names)}.")

            # warn for any requested network that isn't BaseNetwork-derived; it will be ignored.
            unsupported = [nw_name for nw_name in self.target_nw_names
                           if not isinstance(self.sim.networks[nw_name], sti.BaseNetwork)]
            if len(unsupported) > 0:
                unsupported_str = ', '.join(unsupported)
                warnings.warn(f"PartnershipFormationAnalyzer: network(s) '{unsupported_str}' not a stisim "
                              f"BaseNetwork-derived sexual network, ignoring ...")
            candidates = [nw_name for nw_name in self.target_nw_names
                          if nw_name not in missing and nw_name not in unsupported]

        # Warn-and-skip networks that cannot record all edge expirations (e.g.
        # MSMScaleFreeNetwork removes edges outside _on_edge_dissolution), since
        # their active-interval results (ti_expired) would be incomplete.
        selected = []
        for nw_name in candidates:
            nw = self.sim.networks[nw_name]
            if not nw.records_all_expirations:
                self._skipped_nw_names.append(nw_name)
                warnings.warn(f"PartnershipFormationAnalyzer: network '{nw_name}' "
                              f"({type(nw).__name__}) does not record all edge expirations "
                              f"(it removes edges outside _on_edge_dissolution), so "
                              f"active-interval results would be incomplete; skipping.")
                continue
            selected.append(nw_name)
        self._tracked_nw_names = selected

        # Enable expiration recording on tracked networks (pure read-side-effect;
        # supplies ti_expired via nw.expired_this_loop each step).
        for nw_name in self._tracked_nw_names:
            self.sim.networks[nw_name].record_expired = True

        # One mutable row per distinct edge, plus a key->row index (pruned on
        # expiry) used only to patch ti_expired in step().
        self.results['relationships'] = {nw_name: [] for nw_name in self._tracked_nw_names}
        self._indicies_of_rels_in_table_for_nw = {nw_name: {} for nw_name in self._tracked_nw_names}
        return

    def _key_for_rel(self, p1, p2, ti_formed):
        return (int(p1), int(p2), int(ti_formed))

    def step(self):
        ti = self.ti
        ppl = self.sim.people

        for nw_name in self._tracked_nw_names:
            nw = self.sim.networks[nw_name]
            rel_table = self.results['relationships'][nw_name]
            rel_to_index = self._indicies_of_rels_in_table_for_nw[nw_name]

            # (1) append newly formed rels this ti to the rels table with ti_expired=None, and index by key for later
            # recording of actual ti_expired.
            new_rel_mask = np.asarray(nw.edges.ti_formed, dtype=np.int64) == ti
            if new_rel_mask.any():
                p1arr = np.asarray(nw.edges.p1, dtype=np.int64)[new_rel_mask]
                p2arr = np.asarray(nw.edges.p2, dtype=np.int64)[new_rel_mask]
                form_ti_arr = np.asarray(nw.edges.ti_formed, dtype=np.int64)[new_rel_mask]
                age1arr = np.asarray(nw.edges.age_p1, dtype=float)[new_rel_mask]
                age2arr = np.asarray(nw.edges.age_p2, dtype=float)[new_rel_mask]
                # Female sex (T/F) of agents in the relationships; not assuming p1/p2 sex for sexual network generality.
                fem1arr = np.asarray(ppl.female[ss.uids(p1arr)], dtype=bool)
                fem2arr = np.asarray(ppl.female[ss.uids(p2arr)], dtype=bool)
                # insert the new edges into the relationship table, recording its table index in the lookup dict
                # so we can find it later (to update relationship end time). None below denotes as-yet unknown end time.
                for (p1, p2, form_ti, age1, age2, fem1, fem2) in zip(p1arr, p2arr, form_ti_arr, age1arr, age2arr, fem1arr, fem2arr):
                    rel_key = self._key_for_rel(p1=p1, p2=p2, ti_formed=form_ti)
                    rel_to_index[rel_key] = len(rel_table)
                    # Ensure the order of data here matches the edge row column indicies noted at top of this analyzer
                    rel_table.append([int(p1), int(p2), int(form_ti), None, float(age1), float(age2), bool(fem1), bool(fem2)])

            # (2) Edges that expired this step (before disease transmission): record ti_expired = ti and prune
            # the index. Read-only -- the network accumulates into the buffer and clears it itself in finish_step
            # (step 14); see _record_expired.
            self._record_expired(nw, rel_table, rel_to_index, ti_expired=ti)

        # Cache surviving-agent uids and final ti for the getters.
        self._final_uids = {'f': np.asarray(ppl.female.uids, dtype=np.int64),
                            'm': np.asarray(ppl.male.uids, dtype=np.int64)}
        self._final_ti = ti
        return

    def _record_expired(self, nw, rel_table, rel_to_index, ti_expired):
        """Record ``ti_expired`` on each edge in ``nw.expired_this_loop``.

        Edges are matched to their row by the ``(p1, p2, ti_formed)`` key and
        the index entry is pruned. This is **read-only** with respect to the
        network: the network accumulates expirations across ``end_pairs`` and
        ``remove_uids`` and clears the buffer itself in ``finish_step`` (step 14),
        so several analyzers can each read the same buffer, and recording works
        with no analyzer present. Edges never placed in the buffer (i.e. still
        active) are left untouched, keeping ``ti_expired = None``.
        """
        # grab the expired relationship info
        expired_rels = nw.expired_this_loop
        # if at least one relationship expired, record known end dates.
        if len(expired_rels.get('p1', ())) > 0:
            p1_arr = np.asarray(expired_rels['p1'], dtype=np.int64)
            p2_arr = np.asarray(expired_rels['p2'], dtype=np.int64)
            ti_formed_arr = np.asarray(expired_rels['ti_formed'], dtype=np.int64)
            for p1, p2, form_ti in zip(p1_arr, p2_arr, ti_formed_arr):
                # Removing the relationship from our rel_to_index mapping, as we only need to record the end date once.
                rel_key = self._key_for_rel(p1=p1, p2=p2, ti_formed=form_ti)
                row_idx = rel_to_index.pop(rel_key, None)
                # set the end date to the specified ti
                if row_idx is not None:
                    rel_table[row_idx][self._TI_EXPIRED] = ti_expired
        return

    def finalize(self):
        """
        Resolve still-buffered relationships that expired at the end of the simulation.

        After the loop, the only edges left in a network's ``self.expired_this_loop`` are those removed by death at
        the **final** step 15 — i.e. after the last ``step()`` read. They were active through the disease phase of the
        final ti (step 9), so they get ``ti_expired = final_ti + 1``. Relationships/edges still active at sim end were
        never buffered for removal recording, so they keep ``ti_expired = None`` (the right-censored marker; getter
        methods treat ``None`` as still-active relationships at simulation end).
        """
        ti_expired = self._final_ti + 1
        for nw_name in self._tracked_nw_names:
            nw = self.sim.networks[nw_name]
            rel_table = self.results['relationships'][nw_name]
            rel_to_index = self._indicies_of_rels_in_table_for_nw[nw_name]
            self._record_expired(nw, rel_table, rel_to_index, ti_expired=ti_expired)

        super().finalize()
        return

    def _threshold(self, window_months):
        """First ti included by a trailing window (None => whole run)."""
        if window_months is None:
            return None
        return self._final_ti - window_months + 1

    def _rel_element_arrs(self, nw_name):
        """Columnar numpy arrays for one network's edge table.

        ``ti_expired`` is materialized with ``None -> final_ti + 1`` so the
        half-open active test ``ti_formed <= ti < ti_expired`` treats an
        ongoing edge as active through the final timestep.
        """
        rows = self.results['relationships'][nw_name]
        n = len(rows)
        sentinel = self._final_ti + 1
        if n == 0:
            z_int = np.empty(0, dtype=np.int64)
            z_flt = np.empty(0, dtype=float)
            z_bool = np.empty(0, dtype=bool)
            return dict(p1=z_int, p2=z_int.copy(), ti_formed=z_int.copy(),
                        ti_expired=z_int.copy(), age1=z_flt, age2=z_flt.copy(),
                        fem1=z_bool, fem2=z_bool.copy())
        cols = list(zip(*rows))  # tuple per column
        ti_expired = np.fromiter((sentinel if v is None else v for v in cols[self._TI_EXPIRED]),
                                 dtype=np.int64, count=n)
        return dict(
            p1=np.fromiter(cols[self._P1], dtype=np.int64, count=n),
            p2=np.fromiter(cols[self._P2], dtype=np.int64, count=n),
            ti_formed=np.fromiter(cols[self._FORMATION_TI], dtype=np.int64, count=n),
            ti_expired=ti_expired,
            age1=np.fromiter(cols[self._AGE1], dtype=float, count=n),
            age2=np.fromiter(cols[self._AGE2], dtype=float, count=n),
            fem1=np.fromiter(cols[self._FEM1], dtype=bool, count=n),
            fem2=np.fromiter(cols[self._FEM2], dtype=bool, count=n),
        )

    @staticmethod
    def _counts_for_surviving_agents(uids, surviving_agents):
        """Count occurrences of each uid, aligned to ``surviving_agents`` order."""
        counts = np.zeros(len(surviving_agents), dtype=np.int64)
        if len(uids) == 0 or len(surviving_agents) == 0:
            return counts
        tally = Counter(int(u) for u in uids)
        for i, u in enumerate(surviving_agents):
            counts[i] = tally.get(int(u), 0)
        return counts

    def get_n_partnerships_formed(self, female=False, age_bins=None, window_months=None):
        """
        Total count of partnerships formed (non-unique pairings), per network and age bin.

        When window_months is not used (None), one timeseries of partnership formation counts will be returned per
        age bin. When window_months is an int, N, it will return single-element lists of partnership formation counts
        during the final N months of the simulation, one single-element list per age bin.

        Args:
            female (bool): Subject sex (True=female, False=male).
            age_bins (list[str] | None): Bins to report; ``None`` (default) uses
                the analyzer's construction bins. A partner is counted in every
                bin containing its age at formation (overlapping bins each count
                it); a partner whose age is outside every bin is not counted.
            window_months (int | None): If ``None`` (default), the inner array is
                the full per-timestep series (length ``n_ti``, indexed by ti). If
                an int, the inner array is a single trailing-window sum
                (length 1).

        Returns:
            dict ``{network_name: {age_bin: ndarray}}``. The inner array is
            length ``n_ti`` when ``window_months is None``, else length 1.

        Example (per timestep, whole run)::

            ana.get_n_partnerships_formed(female=False, age_bins=['15-50', '40-50'])
            # -> {'mfnetwork': {'15-50': array([0, 2, 0, 0, 4]),   # by timestep
            #                   '40-50': array([0, 0, 1, 0, 0])}}

        Example (single trailing-window sum)::

            ana.get_n_partnerships_formed(female=False, age_bins=['15-50', '40-50'],
                                          window_months=12)
            # -> {'mfnetwork': {'15-50': array([6]),   # sum over last 12 steps
            #                   '40-50': array([1])}}
        """
        bins = self.age_bins if age_bins is None else age_bins
        ranges = self.age_ranges if age_bins is None else [ss.parse_age_range(b) for b in bins]
        n_ti = self._final_ti + 1
        threshold = self._threshold(window_months)
        lo_ti = 0 if threshold is None else max(0, threshold)

        out = {}
        for nw_name in self._tracked_nw_names:
            arr = self._rel_element_arrs(nw_name)
            # Each subject-sex endpoint of every relationship contributes its (ti_formed, age)
            # Note: a same-sex relationship contributes both partners (both "formed a relationship" (double counting, but they will be sorted out to potentially different age bins).
            subject1 = arr['fem1'] if female else ~arr['fem1']
            subject2 = arr['fem2'] if female else ~arr['fem2']
            ti_formeds = np.concatenate([arr['ti_formed'][subject1], arr['ti_formed'][subject2]])
            ages = np.concatenate([arr['age1'][subject1], arr['age2'][subject2]])
            nw_out = {}
            for age_bin, (lo, hi) in zip(bins, ranges):
                age_mask = (ages >= lo) & (ages < hi)
                # Per-timestep counts: a formation at ti contributes to index ti.
                per_ti = np.bincount(ti_formeds[age_mask], minlength=n_ti).astype(np.int64)
                if window_months is None:
                    # here we store the timeseries of relationship formation for return
                    nw_out[age_bin] = per_ti
                else:
                    # window_months was specified, so instead sum the counts in the final window_months timesteps.
                    nw_out[age_bin] = np.array([int(per_ti[lo_ti:].sum())])
            out[nw_name] = nw_out
        return out

    def get_n_partnerships_formed_per_agent(self, female=False, window_months=None):
        """
        Partnerships formed (non-unique) per agent over a trailing window, for agents alive at simulation end.

        Equivalent to asking living agents, "How many partnerships have you formed over the last X months?"

        Args:
            female (bool): Subject sex.
            window_months (int | None): Trailing ti window length for partnership counting, None for whole simulation.

        Returns:
            dict ``{network_name: array}`` indexed like the subject-sex
            surviving-agent uids (``ana._final_uids['f'|'m']``).

        Example (whole run)::

            ana.get_n_partnerships_formed_per_agent(female=False)
            # -> {'mfnetwork': array([0, 0, 1, 0, 4, ...])}
            #    aligned to ana._final_uids['m']; e.g. that male formed 4
            #    partnerships over the run (non-unique: a re-formed pair counts twice).
        """
        threshold = self._threshold(window_months)
        surviving_agent_uids = self._final_uids['f' if female else 'm']
        out = {}
        for nw_name in self._tracked_nw_names:
            arrs = self._rel_element_arrs(nw_name)
            # we do not presume sex of p1 and p2, so we grab the selected sex agents from both (p1 and p2 agents,
            # independently, will be 0-matching or all-matching the female sex bool passed in).
            sex1 = arrs['fem1'] if female else ~arrs['fem1']
            sex2 = arrs['fem2'] if female else ~arrs['fem2']
            uids = np.concatenate([arrs['p1'][sex1], arrs['p2'][sex2]])
            ti_formeds = np.concatenate([arrs['ti_formed'][sex1], arrs['ti_formed'][sex2]])
            if threshold is not None:
                keep = ti_formeds >= threshold
                uids = uids[keep]
            out[nw_name] = self._counts_for_surviving_agents(uids, surviving_agent_uids)
        return out

    def get_unique_partners_active_per_agent(self, female=False, window_months=None):
        """
        Unique partner counts per agent over a trailing window, for agents alive at simulation end.

        Equivalent to asking living agents, "How many unique partners have you had over the last X months?"

        An edge counts as active in the trailing window iff it had not expired
        before the window opened, i.e. ``ti_expired > win_start`` (its half-open
        active interval ``[ti_formed, ti_expired)`` overlaps the window). The
        window always ends at the final timestep, so no upper bound on
        ``ti_formed`` is needed; an ongoing edge carries
        ``ti_expired = final_ti + 1`` and therefore always counts. For
        ``window_months=None`` the window is the whole run (``win_start = 0``).
        Sex-general: the partner is the other endpoint regardless of its sex, so
        this works for same-sex networks too. Dedups by partner uid (a network is
        assumed not to hold concurrent duplicate edges between the same pair, so
        unique partners == unique active partnerships).

        Args:
            female (bool): Subject sex.
            window_months (int | None): Trailing ti window length for partner counting, None for whole simulation.

        Returns:
            dict ``{network_name: array}`` indexed like the subject-sex
            surviving-agent uids.

        Example (trailing 12 timesteps)::

            ana.get_unique_partners_active_per_agent(female=False, window_months=12)
            # -> {'mfnetwork': array([1, 0, 2, ...])}
            #    aligned to ana._final_uids['m']; e.g. the male at index 2 had 2
            #    distinct partners across edges active during the last 12 steps.
        """
        threshold = self._threshold(window_months)
        win_start = 0 if threshold is None else threshold
        surviving_agent_uids = self._final_uids['f' if female else 'm']
        out = {}
        for nw_name in self._tracked_nw_names:
            arr = self._rel_element_arrs(nw_name)
            # ti_expired > win_start drops edges that ended before the window opened;
            # ongoing edges carry ti_expired = final_ti + 1 (see _rel_element_arrs) so they pass.
            active_in_win = arr['ti_expired'] > win_start

            # subject_mask1/subject_mask2 select active edges whose p1/p2 (respectively)
            # is the requested sex (female if female=True, else male).
            subject_mask1 = active_in_win & (arr['fem1'] if female else ~arr['fem1'])
            subject_mask2 = active_in_win & (arr['fem2'] if female else ~arr['fem2'])

            # female (bool) identifies the subject -- subjects can be in p1 and/or p2 as we do not assume sex of p1/p2
            subjects = np.concatenate([arr['p1'][subject_mask1], arr['p2'][subject_mask2]])
            # swapping p1/p2 identifies the partners of the subjects (in subject-order of arrays)
            partners = np.concatenate([arr['p2'][subject_mask1], arr['p1'][subject_mask2]])

            if len(subjects) > 0:
                # Distinct (subject, partner) pairs -> count partners per subject.
                # Create 2-element pairing arrays from aligned subjects and partners arrays ([[sub0, part0], ..., [subN, partN]]
                pairings = np.stack([subjects, partners], axis=1)

                # treats each [subI, partI] element as items to be uniquified
                unique_pairs = np.unique(pairings, axis=0)
                out[nw_name] = self._counts_for_surviving_agents(unique_pairs[:, 0], surviving_agent_uids)
            else:
                out[nw_name] = np.zeros(len(surviving_agent_uids), dtype=np.int64)
        return out

