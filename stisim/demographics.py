"""
Define demographic modules used within STIsim
"""
import numpy as np
import starsim as ss
import sciris as sc
import pandas as pd
import stisim as sti


__all__ = ['Pregnancy', 'Migration']


class Pregnancy(ss.Pregnancy):
    """
    Pregnancy module for STIsim - includes preterm birth and gestational age
    """
    @staticmethod
    def p_ptb(self, sim, uids):
        ptb_probs = self.rel_sus_ptb[uids] * self.pars["p_ptb"]
        return ptb_probs

    def __init__(self, pars=None, metadata=None, **kwargs):
        super().__init__(pars=pars, metadata=metadata, **kwargs)
        self.define_states(
            ss.FloatArr("rel_sus_ptb", default=1),
            ss.BoolState("ptb", default=False, label="Preterm birth"),
            ss.FloatArr("trimester", default=0),
        )

        self.define_pars(
            p_ptb=0.057,  # Probability of preterm birth
        )

        self.p_ptb = ss.bernoulli(p=self.p_ptb)
        return

    @property
    def gestational_age(self):
        return self.ti - self.ti_pregnant

    def init_results(self):
        """
        Results could include a range of birth outcomes e.g. LGA, stillbirths, etc.
        Still unclear whether this logic should live in the pregnancy module, the
        individual disease modules, the connectors, or the sim.
        """
        super().init_results()
        self.define_results(
            ss.Result("ptbs", dtype=int, scale=True, label="New preterm births"),
        )
        return

    def set_prognoses(self, uids, bv_update=False):
        if bv_update is False:
            super().set_prognoses(uids)

            # Set trimester of pregnancy
            self.trimester[uids] = 1

        # Determine pre-term birth
        ptb, no_ptb = self.p_ptb.split(uids)
        self.ptb[ptb] = True
        self.ptb[no_ptb] = False
        return

    def update_states(self):
        ptb_deliveries = self.pregnant & self.ptb & (self.ti_delivery <= self.ti)
        self.results["ptbs"][self.ti] = np.count_nonzero(ptb_deliveries)
        super().update_states()
        self.ptb[ptb_deliveries] = False

        # Find pregnant women and update trimester
        pregnant_uids = self.pregnant.uids
        time_pregnant = self.ti - self.ti_pregnant[pregnant_uids]
        trimester = np.select(
            [time_pregnant <= 3, time_pregnant <= 6], [1, 2], default=3
        )
        self.trimester[pregnant_uids] = trimester

        postpartum_uids = self.postpartum.uids
        self.trimester[postpartum_uids] = 0

        return


class Migration(ss.Demographics):
    """
    Remove / add migrants
    Assumes a datafile with the number of migrants each year
    """
    def __init__(self, pars=None, migration_data=None, **kwargs):
        super().__init__()
        self.define_pars(
            # Migration parameters
            dt = ss.years(1),  # Unit of time for migration data
            migration_propensity = ss.normal(loc=1, scale=0.1),  # Propensity to emigrate
            slot_scale = 5, # Random slots will be assigned to newborn agents between min=n_agents and max=slot_scale*n_agents
            min_slots  = 100, # Minimum number of slots, useful if the population size is very small
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('immigrant', label='Immigrant'),
            ss.FloatArr('migration_propensity', label='Propensity to migrate'),
            ss.FloatArr('ti_emigrate', label='Time of emigration'),
        )
        self.migration_data = migration_data
        self.choose_slots = None # Distribution for choosing slots; set in self.init()
        return

    def init_migration_propensity(self):
        """
        Set individual's propensity to migrate
        This is currently a random variable but could be defined as a function of age/sex/other properties
        """
        uids = ss.uids(self.migration_propensity.isnan)
        self.migration_propensity[uids] = self.pars.migration_propensity.rvs(uids)
        return

    def init_pre(self, sim):
        super().init_pre(sim)

        # Slots
        low = sim.pars.n_agents + 1
        high = int(self.pars.slot_scale*sim.pars.n_agents)
        high = np.maximum(high, self.pars.min_slots) # Make sure there are at least min_slots slots to avoid artifacts related to small populations
        self.choose_slots = ss.randint(low=low, high=high, sim=sim, module=self)

        # Validate and process data
        data = self.migration_data
        if data is None:
            errormsg = 'Migration data not provided. Please provide dataframe with "year" as the index'
            raise ValueError(errormsg)
        else:
            # Check that the columns contain 'Time' and 'Value'
            if 'Time' not in data.columns or 'Value' not in data.columns:
                errormsg = 'Migration data must contain columns labeled "Time" and "Value"'
                raise ValueError(errormsg)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_migrants', dtype=int,   scale=True,  label='New migrants'),
        )
        return

    def get_migrants(self):
        """ Get the number of migrants for this timestep """
        mrd = self.migration_data
        if isinstance(mrd, (pd.Series, pd.DataFrame)):
            year_ind = sc.findnearest(mrd.Time, self.t.now('year'))
            nearest_year = mrd.index[year_ind]
            n_migrants = mrd.loc[nearest_year].Value
        elif sc.isnumber(mrd):
            n_migrants = mrd
        scaled_migrants = int(n_migrants*self.t.dt_year/self.sim.pars.pop_scale)
        return scaled_migrants

    def make_immigrants(self, twin_uids):
        """ Make immigrants by making copies of existing agents """
        people = self.sim.people
        n_new = len(twin_uids)
        if n_new == 0:
            new_uids = ss.uids()
        else:
            # Choose slots for the unborn agents
            new_slots = self.choose_slots.rvs(twin_uids)

            # Grow the arrays and set properties for the unborn agents
            new_uids = people.grow(len(new_slots), new_slots)
            people.age[new_uids] = people.age[twin_uids]
            people.slot[new_uids] = new_slots  # Before sampling female_dist
            people.female[new_uids] = people.female[twin_uids]

        return new_uids

    def step(self):
        """ Perform all updates """
        self.init_migration_propensity()  # Set propensity to migrate
        new_migrants = self.step_migration()  # Step through the actual migration logic
        self.remove_emigrants()  # Remove emigrants
        self.results.new_migrants[self.ti] = new_migrants
        return

    def step_migration(self):
        """ Select people to migrate """

        # Get the number of people who'll be coming in / going out this timestep
        new_migrants = self.get_migrants()
        uids = self.sim.people.alive
        weights = self.migration_propensity[uids]

        # new_migrants is positive -> add people to the population
        if new_migrants > 0:
            # Choose people to twin
            choices = np.argsort(-weights)[:new_migrants]
            twin_uids = uids.uids[choices]
            new_people = self.make_immigrants(twin_uids)

        # new_migrants is negative -> remove people from the population
        elif new_migrants < 0:
            choices = np.argsort(-weights)[:-new_migrants]
            migrant_uids = uids.uids[choices]
            self.ti_emigrate[migrant_uids] = self.ti

        return new_migrants

    def remove_emigrants(self):
        """ Remove people who've decided to emigrate this timestep """
        ti = self.ti
        emigrants = (self.ti_emigrate == ti).uids
        self.sim.people.request_removal(emigrants)
        return

