import starsim as ss


class IReporter(ss.Analyzer):
    """
    Core class for creating a consistent and flexible reporter interface. Designed to be subclassed, not used directly.
    """

    def __init__(self, select, when=None, subpop=None):
        """
        Args:
            select: (callable) a function that returns desired data to be tracked by uid, given a sim object
                For example, the CD4 counts of agents:
                    lambda sim: sim.diseases.hiv.cd4
            when: (optional) (callable) a function that returns True/False if the current timestep is to be recorded, given
                a sim object
                For example, timesteps from [10, 20):
                    lambda sim: (sim.ti >= 10) and (sim.ti < 20)
                Default: True for all time indices
            subpop: (optional) (callable, iterable) an iterable of uids OR a function that returns a an iterable of uids
                given a sim object.
                For example, HIV infected agents that are not on ART
                    lambda sim: ( sim.diseases.hiv.infected & (~sim.diseases.hiv.on_art) ).uids
                Default: all sim uids
        """

        super().__init__()
        self.select = select
        self.subpop = self.all_people if subpop is None else subpop
        self.when = self.all_time if when is None else when

    def update_results(self):
        sim = self.sim
        if self.when(sim=sim):
            data = self.select(sim=sim)
            uids = self.subpop(sim=sim)
            self.set_result(sim=sim, data=data, uids=uids)
        else:
            # record nothing
            pass

    def set_result(self, sim, data, uids):
        raise NotImplementedError("Subclasses of IAnalyzer must implement set_result.")

    @property
    def all_people(self):
        return lambda sim: sim.people.uids

    @property
    def all_time(self):
        return lambda sim: True

    def step(self):
        pass
