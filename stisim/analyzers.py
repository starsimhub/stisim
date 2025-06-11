"""
Common analyzers for STI analyses
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import starsim as ss
import pandas as pd

import stisim as sti
import pylab as pl
from collections import defaultdict

__all__ = ["result_grouper", "coinfection_stats", "sw_stats", "RelationshipDurations", "NetworkDegree", "TimeBetweenRelationships"]


class result_grouper(ss.Analyzer):
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
        self.name = 'coinfection_stats'
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
        default_denom = lambda self: (self.sim.people.age >= self.age_limits[0]) & (self.sim.people.age < self.age_limits[0])
        self.denom = denom or default_denom

        return

    def init_results(self):
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
        has_disease2 = getattr(disease2obj, self.disease2_infected_state_name) # Adults with HIV
        has_disease1 = getattr(disease1obj, self.disease1_infected_state_name)  # Adults with syphilis

        has_disease1_f = denom & has_disease1 & ppl.female  # Women with dis1
        has_disease2_m = denom & has_disease1 & ppl.male  # Men with dis1
        has_disease2_f = denom & has_disease2 & ppl.female  # Women with dis2
        has_disease2_m = denom & has_disease2 & ppl.male  # Men with dis2
        no_disease2    = denom & ~has_disease2  # Adults without dis2
        no_disease2_f  = no_disease2 & ppl.female  # Women without dis2
        no_disease2_m  = no_disease2 & ppl.male  # Men without dis2

        self.results[f'{disease1name}_prev_no_{disease2name}'][ti] = self.cond_prob(has_disease1, no_disease2)
        self.results[f'{disease1name}_prev_has_{disease2name}'][ti] = self.cond_prob(has_disease1, has_disease2)
        self.results[f'{disease1name}_prev_no_{disease2name}_f'][ti] = self.cond_prob(has_disease1_f, no_disease2_f)
        self.results[f'{disease1name}_prev_has_{disease2name}_f'][ti] = self.cond_prob(has_disease1_f, has_disease2_f)
        self.results[f'{disease1name}_prev_no_{disease2name}_m'][ti] = self.cond_prob(has_disease2_m, no_disease2_m)
        self.results[f'{disease1name}_prev_has_{disease2name}_m'][ti] = self.cond_prob(has_disease2_m, has_disease2_m)

        return


class sw_stats(result_grouper):
    def __init__(self, diseases=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'sw_stats'
        self.diseases = diseases
        return

    def init_results(self):
        results = sc.autolist()
        for d in self.diseases:
            results += [
                ss.Result('share_new_infections_fsw_'+d, scale=False, summarize_by='mean'),
                ss.Result('share_new_infections_client_'+d,scale=False, summarize_by='mean'),
                ss.Result('new_infections_fsw_'+d, dtype=int),
                ss.Result('new_infections_client_'+d, dtype=int),
                ss.Result('new_infections_non_fsw_'+d, dtype=int),
                ss.Result('new_infections_non_client_'+d, dtype=int),
                ss.Result('new_transmissions_fsw_'+d, dtype=int),
                ss.Result('new_transmissions_client_'+d, dtype=int),
                ss.Result('new_transmissions_non_fsw_'+d, dtype=int),
                ss.Result('new_transmissions_non_client_'+d, dtype=int),
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

                total_trans = sum(new_transmissions_fsw) + sum(new_transmissions_client) + sum(new_transmissions_non_fsw) + sum(new_transmissions_non_client)
                if total_trans != len(newly_infected.uids):
                    errormsg = f'Infections acquired should equal number transmitted: {total_acq} vs {total_trans}'
                    raise ValueError(errormsg)

        return


class RelationshipDurations(ss.Analyzer):
    """
    Analyzes the durations of relationships in a structuredsexual network.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return

    def init_results(self):
        self.define_results(
            ss.Result('mean_duration', dtype=float, scale=False),
            ss.Result('median_duration', dtype=float, scale=False),
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
                ss.Result(f'{relationship_type}_f', dtype=int, scale=False, shape=len(self.bins)),
                ss.Result(f'{relationship_type}_m', dtype=int, scale=False, shape=len(self.bins)),
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
        self.times_between_relationships = defaultdict(self.zero_list)
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
            self.times_between_relationships[uid][-1] += 1  # Increment time since last relationship for those not in stable relationships

        for uid in stable_relationship_uids:
            if self.times_between_relationships[uid][-1] > 0:
                # If the agent is in a stable relationship, reset the time since last relationship and append to the list
                self.times_between_relationships[uid].append(0)

        return
