from stisim.reporters.ireporter import IReporter


class ITimeSeriesReporter(IReporter):
    """
    Core class for creating a consistent and flexible timeseries-based reporter interface. Designed to be subclassed,
    not used directly.
    """

    value_key = "val"
    ti_key = "ti"
    result_name = "ts"

    def init_results(self):
        super().init_results()
        self.results[self.result_name] = {}

    @classmethod
    def _init_result(cls):
        # standard creation of a blank item for recording
        return {cls.value_key: [], cls.ti_key: []}

    @classmethod
    def _set_result(cls, result, value, ti):
        # standard addition of new data item to results
        result[cls.value_key].append(value)
        result[cls.ti_key].append(ti)

    @property
    def ts_results(self):
        # root where all timeseries reporters record their result
        return self.results[self.result_name]
