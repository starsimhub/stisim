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
    """Track partnership formation per network, gender, age bin, and timestep.

    Each timestep, for every tracked network, records the active edges and which
    of them were newly formed (``formation_ti == ti``). Each partner's sex and
    age-at-formation are taken from the edge/people directly, so the analyzer
    makes no p1=male/p2=female assumption and works for any gender combination
    (male-female, male-male, female-female).

    All recorded data lives in ``self.results``:

      - ``self.results['n_formed'][nw][sex][age_bin]`` -> ``array(n_ti)``:
        per-timestep count of partnerships formed in which a partner of that
        ``sex`` (``'f'``/``'m'``) fell in ``age_bin`` at formation. Both
        partners of every new edge are counted (a new partnership updates two
        cells).
      - ``self.results['edges'][nw]`` -> list of per-timestep records
        ``(ti, p1, p2, formation_ti, age_p1, age_p2, female_p1, female_p2)``,
        used by the per-agent / windowed getters.

    Example (``self.results`` layout, default bins, network 'structuredsexual')::

        # per-timestep male-side formations in age bin '20-25':
        ana.results['n_formed']['structuredsexual']['m']['20-25']
        # -> array([0, 1, 0, 2, 0, ...])   # 1 at ti=1, 2 at ti=3, ... (len == n_ti)

        # the raw edge record captured at one timestep (here a single edge):
        ana.results['edges']['structuredsexual'][0]
        # -> (0,                       # ti
        #     array([9]), array([17]),   # p1, p2  (edge endpoints' uids)
        #     array([0]),               # formation_ti (formed at ti=0)
        #     array([27.23]),           # age_p1 at formation
        #     array([19.62]),           # age_p2 at formation
        #     array([False]),           # female_p1 (this p1 is male)
        #     array([ True]))           # female_p2 (this p2 is female)

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
        window.

    Args:
        age_bins (list[str]): Age bin strings parsable by ``ss.parse_age_range``
            (e.g. ``['0-5', '5-10', ..., '70-75']``). Bins may overlap; a
            partner is counted in every bin containing its age. Defaults to
            five-year bins from ``'0-5'`` through ``'70-75'``.
        networks (list[str] | None): Network names to track. ``None`` (default)
            tracks every ``ss.SexualNetwork`` in the sim. Named networks that
            are missing raise ``ValueError``; named non-sexual networks warn and
            are skipped. Tracked networks' edges must carry ``formation_ti``,
            ``age_p1``, ``age_p2`` (BaseNetwork provides these).
    """

    DEFAULT_AGE_BINS = [f'{lo}-{lo+5}' for lo in range(0, 75, 5)]

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
        self.target_network_names = networks
        self._tracked_names = []                     # populated in init_results
        self._final_uids = {'f': None, 'm': None}    # subject universes, set in step
        self._final_ti = None
        return

    def init_results(self):
        super().init_results()

        # Validate passed-in network names for tracking.
        all_nw_names = self.sim.networks.keys()
        if self.target_network_names is None:
            # default to all sexual networks
            selected = [nw_name for nw_name in all_nw_names if isinstance(self.sim.networks[nw_name], ss.SexualNetwork)]
        else:
            # ensure all requested networks exist
            missing = [nw_name for nw_name in self.target_network_names if nw_name not in all_nw_names]
            if len(missing) > 0:
                missing_str = ', '.join(missing)
                raise ValueError(f"PartnershipFormationAnalyzer: network(s) '{missing_str}' not found in sim.networks. "
                                 f"Available: {sorted(all_nw_names)}.")

            # warn for any requested network that is non sexual; it will be ignored.
            non_sexual = [nw_name for nw_name in self.target_network_names
                          if not isinstance(self.sim.networks[nw_name], ss.SexualNetwork)]
            if len(non_sexual) > 0:
                non_sexual_str = ', '.join(non_sexual)
                warnings.warn(f"PartnershipFormationAnalyzer: network(s) '{non_sexual_str}' not a sexual network, ignoring ...")
            selected = [nw_name for nw_name in self.target_network_names
                        if nw_name not in missing and nw_name not in non_sexual]

        self._tracked_names = selected

        # Canonical per-timestep formation time series (the stisim-norm storage):
        #   self.results['n_formed'][nw_name][sex][age_bin] -> array(n_ti)
        n_ti = len(self.t.tvec)
        n_formed = {}
        for nw_name in self._tracked_names:
            n_formed[nw_name] = {
                'f': {age_bin: np.zeros(n_ti, dtype=int) for age_bin in self.age_bins},
                'm': {age_bin: np.zeros(n_ti, dtype=int) for age_bin in self.age_bins},
            }
        self.results['n_formed'] = n_formed

        # Raw per-timestep edge records, for the per-agent / windowed getters:
        #   self.results['edges'][nw_name] -> list of
        #     (ti, p1, p2, formation_ti, age_p1, age_p2, female_p1, female_p2)
        self.results['edges'] = {nw_name: [] for nw_name in self._tracked_names}
        return

    def step(self):
        ti = self.ti
        ppl = self.sim.people

        for nw_name in self._tracked_names:
            nw = self.sim.networks[nw_name]
            p1 = np.asarray(nw.edges.p1, dtype=np.int64)
            p2 = np.asarray(nw.edges.p2, dtype=np.int64)
            if len(p1) == 0:
                continue
            formation_ti = np.asarray(nw.edges.formation_ti, dtype=np.int64)
            # age_p1 / age_p2 are stored on the edge at formation, so they give
            # each partner's age at the moment the partnership formed.
            age_p1 = np.asarray(nw.edges.age_p1, dtype=float)
            age_p2 = np.asarray(nw.edges.age_p2, dtype=float)
            # Determine each partner's sex from people (no p1/p2 assumption), so
            # same-sex networks are handled correctly.
            female_p1 = np.asarray(ppl.female[ss.uids(p1)], dtype=bool)
            female_p2 = np.asarray(ppl.female[ss.uids(p2)], dtype=bool)

            self.results['edges'][nw_name].append(
                (ti, p1, p2, formation_ti, age_p1, age_p2, female_p1, female_p2))

            # Per-timestep formation counts: both partners of each new edge
            # contribute to their own (sex, age-bin) cell.
            new_mask = formation_ti == ti
            if new_mask.any():
                nfm = self.results['n_formed'][nw_name]
                for ages, fems in ((age_p1[new_mask], female_p1[new_mask]),
                                   (age_p2[new_mask], female_p2[new_mask])):
                    for sex_key, sex_mask in (('f', fems), ('m', ~fems)):
                        if not sex_mask.any():
                            continue
                        bin_counts = self._bin_counts(ages[sex_mask], self.age_ranges)
                        for i, age_bin in enumerate(self.age_bins):
                            if bin_counts[i]:
                                nfm[sex_key][age_bin][ti] += bin_counts[i]

        # Cache subject universes and final ti for the per-agent getters.
        self._final_uids = {'f': np.asarray(ppl.female.uids, dtype=np.int64),
                            'm': np.asarray(ppl.male.uids, dtype=np.int64)}
        self._final_ti = ti
        return

    @staticmethod
    def _bin_counts(ages, ranges):
        """Count ages into (possibly overlapping) bins; returns array(len(ranges))."""
        counts = np.zeros(len(ranges), dtype=np.int64)
        for i, (lo, hi) in enumerate(ranges):
            counts[i] = int(np.count_nonzero((ages >= lo) & (ages < hi)))
        return counts

    def _threshold(self, window_months):
        """First ti included by a trailing window (None => whole run)."""
        if window_months is None:
            return None
        return self._final_ti - window_months + 1

    def _iter_formations(self, nw_name, female):
        """Yield (ti, subject_uids, subject_ages) for newly-formed edges whose
        endpoint matches the subject sex (ages are ages at formation)."""
        for ti, p1, p2, fti, age1, age2, fem1, fem2 in self.results['edges'][nw_name]:
            new_mask = fti == ti
            if not new_mask.any():
                continue
            m1 = new_mask & (fem1 if female else ~fem1)
            m2 = new_mask & (fem2 if female else ~fem2)
            uids = np.concatenate([p1[m1], p2[m2]])
            if len(uids) == 0:
                continue
            ages = np.concatenate([age1[m1], age2[m2]])
            yield ti, uids, ages

    def _formed_per_ti(self, nw_name, female, ranges, age_bins):
        """Per-(age-bin, ti) formation counts for one network: array(n_bins, n_ti).

        Uses the stored ``n_formed`` time series when ``age_bins`` are the
        construction bins; recomputes from the raw edge records otherwise.
        """
        if age_bins is self.age_bins:
            nfm = self.results['n_formed'][nw_name]['f' if female else 'm']
            return np.vstack([nfm[age_bin] for age_bin in self.age_bins])
        # final_ti is the last stepped timestep, so n_ti == final_ti + 1 (matches
        # len(self.t.tvec) after a run, and avoids needing self.t in unit tests).
        n_ti = self._final_ti + 1
        arr = np.zeros((len(ranges), n_ti), dtype=np.int64)
        for ti, uids, ages in self._iter_formations(nw_name, female):
            arr[:, ti] += self._bin_counts(ages, ranges)
        return arr

    def get_n_partnerships_formed(self, female=False, age_bins=None, window_months=None):
        """Partnerships formed (non-unique), per network and age bin.

        Args:
            female (bool): Subject sex (True=female, False=male).
            age_bins (list[str] | None): Bins to report; ``None`` (default) uses
                the analyzer's construction bins. A partner is counted in every
                bin containing its age at formation (overlapping bins each count
                it).
            window_months (int | None): If ``None`` (default), the inner array is
                the full per-timestep series (length ``n_ti``, indexed by ti). If
                an int, the inner array is a single trailing-window sum
                (length 1).

        Returns:
            dict ``{network_name: {age_bin: ndarray}}``. The inner array is
            length ``n_ti`` when ``window_months is None``, else length 1.

        Example (per timestep, whole run)::

            ana.get_n_partnerships_formed(female=False, age_bins=['15-50', '40-50'])
            # -> {'structuredsexual': {'15-50': array([0, 2, 0, 0, 4]),   # by timestep
            #                          '40-50': array([0, 0, 1, 0, 0])}}

        Example (single trailing-window sum)::

            ana.get_n_partnerships_formed(female=False, age_bins=['15-50', '40-50'],
                                          window_months=12)
            # -> {'structuredsexual': {'15-50': array([6]),   # sum over last 12 steps
            #                          '40-50': array([1])}}
        """
        bins = self.age_bins if age_bins is None else age_bins
        ranges = self.age_ranges if age_bins is None else [ss.parse_age_range(b) for b in bins]
        threshold = self._threshold(window_months)

        out = {}
        for nw_name in self._tracked_names:
            per_ti = self._formed_per_ti(nw_name, female, ranges, bins)  # (n_bins, n_ti)
            if window_months is None:
                out[nw_name] = {age_bin: per_ti[i] for i, age_bin in enumerate(bins)}
            else:
                lo = max(0, threshold)
                out[nw_name] = {age_bin: np.array([int(per_ti[i, lo:].sum())])
                                for i, age_bin in enumerate(bins)}
        return out

    def get_n_partnerships_formed_per_agent(self, female=False, window_months=None):
        """Partnerships formed (non-unique) per agent over a trailing window.

        Args:
            female (bool): Subject sex.
            window_months (int | None): Trailing window length in timesteps;
                ``None`` aggregates over the whole run.

        Returns:
            dict ``{network_name: array}`` indexed like the subject-sex uid
            universe (``ana._final_uids['f'|'m']``).

        Example (whole run)::

            ana.get_n_partnerships_formed_per_agent(female=False)
            # -> {'structuredsexual': array([0, 0, 1, 0, 4, ...])}
            #    aligned to ana._final_uids['m']; e.g. that male formed 4
            #    partnerships over the run (non-unique: a re-formed pair counts twice).
        """
        threshold = self._threshold(window_months)
        uid_universe = self._final_uids['f' if female else 'm']
        out = {}
        for nw_name in self._tracked_names:
            counts_by_uid = defaultdict(int)
            for ti, uids, ages in self._iter_formations(nw_name, female):
                if threshold is not None and ti < threshold:
                    continue
                for u in uids:
                    counts_by_uid[int(u)] += 1
            out[nw_name] = np.array([counts_by_uid.get(int(u), 0) for u in uid_universe],
                                    dtype=np.int64)
        return out

    def get_unique_partners_active_per_agent(self, female=False, window_months=None):
        """Per-agent count of unique partners over a trailing window.

        Counts partners from edges *active* during the window (not just newly
        formed). Sex-general: the partner is the other endpoint regardless of
        its sex, so this works for same-sex networks too. Dedups by partner
        uid (a network is assumed not to hold concurrent duplicate edges
        between the same pair, so unique partners == unique active partnerships).

        Args:
            female (bool): Subject sex.
            window_months (int | None): Trailing window length; ``None`` uses
                the whole run.

        Returns:
            dict ``{network_name: array}`` indexed like the subject-sex uid
            universe.

        Example (trailing 12 timesteps)::

            ana.get_unique_partners_active_per_agent(female=False, window_months=12)
            # -> {'structuredsexual': array([1, 0, 2, ...])}
            #    aligned to ana._final_uids['m']; e.g. the male at index 2 had 2
            #    distinct partners across edges active during the last 12 steps
            #    (re-forming with the same partner still counts as 1 unique).
        """
        threshold = self._threshold(window_months)
        uid_universe = self._final_uids['f' if female else 'm']
        out = {}
        for nw_name in self._tracked_names:
            partner_sets = defaultdict(set)
            for ti, p1, p2, fti, age1, age2, fem1, fem2 in self.results['edges'][nw_name]:
                if threshold is not None and ti < threshold:
                    continue
                # Subject = endpoints matching the requested sex; the partner is
                # the opposite endpoint (whatever its sex).
                sub1 = fem1 if female else ~fem1
                sub2 = fem2 if female else ~fem2
                for s, partner in zip(p1[sub1], p2[sub1]):
                    partner_sets[int(s)].add(int(partner))
                for s, partner in zip(p2[sub2], p1[sub2]):
                    partner_sets[int(s)].add(int(partner))
            out[nw_name] = np.array([len(partner_sets.get(int(u), ())) for u in uid_universe],
                                    dtype=np.int64)
        return out

