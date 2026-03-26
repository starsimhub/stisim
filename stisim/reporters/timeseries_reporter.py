import starsim as ss

from stisim.itimeseries_reporter import ITimeSeriesReporter


class TimeSeriesReporter(ITimeSeriesReporter):
    """
    An analyzer that generates and reports custom time-index-series data.

    Results are stored as a series in self.results dict accessed by key 'timeseries' .

    For example, to record the number of infected agents per timestep:

    selector = lambda sim: len(sim.diseases.hiv.infected.uids)

     In a three-timestep sim infecting agents on step 1, this could store the following in self.results['ts']

    [0.0, 14.0, 17.0]
    """
    def __init__(self, agg=None, *args, **kwargs):
        """
        Args:
            agg: (optional) (callable) a function that aggregates incoming data (data parameter) at a timestep,
                returning a result to be recorded.
                Default: no aggregation, selected data recorded as-is.
        """
        super().__init__(*args, **kwargs)
        self.agg = self.no_agg if agg is None else agg

    def set_result(self, sim, data, uids):
        uids = ss.arrays.uids(uids)  # this type casting is required for the next line to work right
        uid_filtered_data = data[uids]
        agg_data = self.agg(data=uid_filtered_data)

        self._set_result(result=self.ts_results, value=agg_data, ti=self.ti)

    @property
    def no_agg(self):
        return lambda data: data
