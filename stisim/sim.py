import starsim as ss
import stisim as sti
import sciris as sc
from itertools import combinations
from .data import loaders as stidata
from .data import downloaders as stidl


class Sim(ss.Sim):
    """
    A subclass of starsim.Sim that is specifically designed for STI simulations.
    It initializes with a structured sexual network and includes STI-specific analyzers.

    Args:
        pars (dict): Parameters for the simulation
        label (str): A label for the simulation.
        people (ss.People): People object containing the agents in the simulation.
        demographics (ss.Demographics, list): Demographic modules to include.
        diseases (ss.Disease, list): Disease modules to include.
        networks (ss.Network, list): Network modules to include.
        interventions (ss.Intervention, list): Intervention modules to include.
        analyzers (ss.Analyzer, list): Analyzer modules to include.
        connectors (bool or list): If True, use default connectors; otherwise, provide a list of connectors.
        copy_inputs (bool): Whether to copy input parameters or not.
        data: Additional data to be used in the simulation.

        This class provides flexibility to initialize starsim modules in various ways. Default values are provided,
        so sti.Sim() can be called and generate a reasonable simulation without any inputs. Alternatively, modules can
        be passed in as strings by name, or as starsim modules. If a module is passed in as a string, it will be loaded
        from the stisim or starsim modules, and initialized with any parameters provided in the pars dictionary.

        Pars for each module can be provided via kwargs labeled with the format <module_type>_pars, e.g. disease_pars
        or network_pars, and should be a dictionary with keys corresponding to the module names.
        e.g.:

        ```python
        disease_pars = {
            'hiv': {'init_prev': 0.1,},
            'ng': {'init_prev': 0.05,}
        }
        ```

        Example usage:
        ```python
        sim = sti.Sim(diseases=['hiv', 'syphilis'])

        # The above will initialize a simulation for HIV and syphilis, using default parameters, including networks,
        # demographics, and connectors. If you want to specify parameters for the diseases, you can do so like this:

        sim = sti.Sim(diseases=['hiv', 'syphilis'],
                    disease_pars={
                        'hiv': {'init_prev': 0.1},
                        'syphilis': {'init_prev': 0.05}})
        ```

        You can also pass in custom parameters for the default modules, such as networks or demographics, like this:
        ```python
        sim = sti.Sim(network_pars={
                    'structuredsexual' = {(sw_seeking_rate=0.5)},
                    'maternalnet': {'maternal_age': 30}
        })
        ```

        These are all in addition to the standard starsim input style using module instances directly, e.g.:
        ```python
        sim = sti.Sim(diseases=[sti.HIV(init_prev=0.1), sti.Syphilis(init_prev=0.05)],
                      networks=[sti.StructuredSexual(sw_seeking_rate=0.5)])
        ```

        If parameters are provided for the same module type in more than one location (e.g. both in the pars dictionary
        and in the kwargs), the order of precedence is as follows:
        1. Parameters provided in the kwargs (e.g. disease_pars, network_pars)
        2. Parameters provided in the input args (e.g. diseases, networks)
        3. Parameters provided in the pars dictionary
        4. Default parameters defined in the sti.Sim.default_pars() method.
        5. Default parameters defined in each module's own class.
    """

    def __init__(self, pars=None, sim_pars=None, sti_pars=None, nw_pars=None,
                 label=None, people=None, demographics=None, diseases=None, networks=None,
                 interventions=None, analyzers=None, connectors=None, datafolder=None, **kwargs):

        # Inputs and defaults
        self.nw_pars = None     # Parameters for the networks - processed later
        self.sti_pars = None    # Parameters for the STIs - processed later
        self.pars = None        # Parameters for the simulation - processed later
        self.datafolder = datafolder
        # self.stis = ss.ndict()   # Set during init, after processing the STIs

        # Call the constructor of the parent class WITHOUT pars or module args
        super().__init__(pars=None, label=label)
        self.pars = sti.SimPars()  # Make default parameters

        # Separate the parameters, storing sim pars now and saving module pars to process in init
        sim_kwargs = dict(label=label, people=people, demographics=demographics, diseases=diseases, networks=networks,
                    interventions=interventions, analyzers=analyzers, connectors=connectors)
        sim_kwargs = {key: val for key, val in sim_kwargs.items() if val is not None}
        updated_pars = self.separate_pars(pars, sim_pars, sti_pars, nw_pars, sim_kwargs, **kwargs)
        self.pars.update(updated_pars)
        return

    def separate_pars(self, pars=None, sim_pars=None, sti_pars=None, nw_pars=None, sim_kwargs=None, **kwargs):
        """
        Create a nested dict of parameters that get passed to Sim constructor and the component modules
        Prioritization:
            - If any key appears in both pars and *_pars, the value from *_pars will be used.
            - If any key appears in both pars and kwargs, the value from kwargs will be used.
        """
        # Merge in pars and kwargs
        all_pars = sc.mergedicts(pars, sim_pars, sti_pars, nw_pars, sim_kwargs, kwargs)
        all_pars = self.remap_pars(all_pars)  # Remap any parameter names

        # Deal with sim pars
        if 'dur' in all_pars.keys():
            self.pars['stop'] = None
        user_sim_pars = {k: v for k, v in all_pars.items() if k in self.pars.keys()}
        for k in user_sim_pars: all_pars.pop(k)
        sim_pars = sc.mergedicts(user_sim_pars, sim_pars, _copy=True)  # Don't merge with defaults, those are set above

        # Deal with STI pars
        all_sti_pars = sti.merged_sti_pars()  # All STI parameters, ignoring duplicates
        user_sti_pars = {}
        for k, v in all_pars.items():
            if k in all_sti_pars.keys(): user_sti_pars[k] = v  # Just set
            if sc.checktype(v, dict):  # See whether it contains STI pars
                user_sti_pars[k] = {gk: gv for gk, gv in v.items() if gk in all_sti_pars}
        for k in user_sti_pars: all_pars.pop(k)
        sti_pars = sc.mergedicts(user_sti_pars, sti_pars, _copy=True)

        # Deal with network pars
        default_nw_pars = sti.NetworkPars()
        user_nw_pars = {k: v for k, v in all_pars.items() if k in default_nw_pars.keys()}
        for k in user_nw_pars: all_pars.pop(k)
        nw_pars = sc.mergedicts(default_nw_pars, user_nw_pars, nw_pars, _copy=True)

        # Raise an exception if there are any leftover pars
        if all_pars:
            raise ValueError(f'Unrecognized parameters: {all_pars.keys()}. Refer to parameters.py for parameters.')

        # Store the parameters for the modules - thse will be fed into the modules during init
        self.sti_pars = sti_pars    # Parameters for the STI modules
        self.nw_pars = nw_pars      # Parameters for the networks

        return sim_pars

    @staticmethod
    def remap_pars(pars):
        """
        Remap any parameter names to match the expected format for STIs and networks.
        This is useful for ensuring that parameters are correctly interpreted by the modules.
        """
        if 'beta' in pars and sc.isnumber(pars['beta']):
            pars['beta_m2f'] = pars.pop('beta')
        if 'location' in pars:
            pars['demographics'] = pars.pop('location')
        return pars

    def init(self, force=False, **kwargs):
        """
        Perform all initializations for the sim
        """
        # Process the STIs
        self.pars['diseases'] = self.process_stis()
        connectors = self.process_connectors()
        self.pars['connectors'] += connectors

        # Process the network
        self.pars['networks'] = self.process_networks()

        # Process the demographics
        demographics, people, total_pop = self.process_demographics()
        self.pars['demographics'] = demographics
        self.pars['people'] = people
        self.pars['total_pop'] = total_pop
        if self.pars['demographics'] is not None:
            self.pars.use_aging = True

        # Reset n_agents
        if self.pars['people'] is not None: self.pars['n_agents'] = len(self.pars['people'])

        super().init(force=force, **kwargs)  # Call the parent init method

        return self

    def process_networks(self):
        """
        Process the network parameters to create network module.
        If networks are provided, they will be used; otherwise, use default networks (usual case)
        """
        if len(self.pars['networks']):
            # If networks are provided, use them directly
            networks = self.pars['networks']

        else:
            # If no networks are provided, create them based on the network parameters
            networks = ss.ndict(
                sti.StructuredSexual(pars=self.nw_pars),
                ss.MaternalNet(),
            )
        return networks

    def process_demographics(self):
        """ Process the location to create people and demographics if not provided. """

        # If it's a string, do lots of work
        if sc.checktype(self.pars['demographics'], str):
            location = self.pars.pop('demographics')
            self.pars['demographics'] = ss.ndict()
            demographics = sc.autolist()

            if self.datafolder is None:
                # Check that the necessary data files are available
                indicators = ['age', 'deaths']
                if self.pars['use_pregnancy']: indicators.append('asfr')
                else: indicators.append('births')
                if self.pars['use_migration']: indicators.append('migration')
                start_year = ss.date(self.pars['start']).year
                ok, missing = stidl.check_downloaded(location, indicators, year=start_year)

                # If they aren't available, try to download them
                if not ok:
                    printmsg = (f'Could not find demographic data files for "{location}", attempting to download. '
                                f'Note that this requires an internet connection.')
                    print(printmsg, end='')
                    stidl.download_data(location=location, indicators=missing, start=start_year)

            # Load birth or fertility rates and turn into module
            if self.pars['use_pregnancy']:
                fertility_rates = stidata.get_rates(location, 'asfr', self.datafolder)
                pregnancy = sti.Pregnancy(fertility_rate=fertility_rates, metadata=dict(data_cols=dict(year='Time', age='AgeStart', value='Value')),)
                demographics += pregnancy
            else:
                birth_rates = stidata.get_rates(location, 'births', self.datafolder)
                births = ss.Births(birth_rate=birth_rates, metadata=dict(data_cols=dict(year='year', value='cbr')))
                demographics += births

            # Load death rates and turn into a module
            death_rates = stidata.get_rates(location, 'death', self.datafolder)
            deaths = ss.Deaths(death_rate=death_rates, rate_units=1, metadata=dict(data_cols=dict(year='Time', sex='Sex', age='AgeStart', value='Value')))
            demographics += deaths

            # Optionally add migration
            if self.pars['use_migration']:
                migration_data = stidata.get_rates(location, 'migration', self.datafolder)
                migration = sti.Migration(migration_data=migration_data)
                demographics += migration

            # Load age data and create people
            age_data = stidata.get_age_distribution(location, year=self.pars.start, datafolder=self.datafolder)
            total_pop = int(age_data.value.sum())
            age_data['value'] /= sum(age_data['value'])  # Normalize the age distribution
            people = ss.People(self.pars.n_agents, age_data=age_data)

        else:
            demographics = self.pars['demographics']
            people = self.pars['people']
            total_pop = self.pars['total_pop']

        return demographics, people, total_pop

    def process_stis(self):
        """
        Look up a disease by its name and return the corresponding module.
        """
        disease_pars = self.pars['diseases']
        if not disease_pars:  # Handles e.g. ss.ndict()
            disease_pars = []
        else:
            disease_pars = sc.tolist(disease_pars)  # Ensure it's a list
        stis = sc.autolist()
        if len(disease_pars) == 0:
            return stis

        all_sti_pars = sti.merged_sti_pars()  # All STI parameters, ignoring duplicates
        sti_main_pars = {k: v for k, v in self.sti_pars.items() if k in all_sti_pars}
        remaining_pars = {k: v for k, v in self.sti_pars.items() if k not in all_sti_pars}

        # Check that the remaining parameters are keyed by STI, remapping them if needed
        sti_options, sti_mapping = sti.sti_aliases()  # Get the options and mapping
        sti_pars = {}
        for sparname, spardict in remaining_pars.items():
            if sparname not in sti_mapping.keys():
                raise ValueError(f'Parameters for STI {sparname} were provided, but this is not an inbuilt STI')
            else:
                sti_pars[sti_mapping[sparname]] = spardict

        # Construct or interpret the STIs from the pars
        for stidis in disease_pars:

            # If it's a string, convert to a module
            if sc.checktype(stidis, str):
                if stidis not in sti_options.keys():
                    errormsg = f'STI {stidis} is not one of the inbuilt options.'
                    raise ValueError(errormsg)

                # See if any parameters have been provided for this STI
                this_stitype_pars = {}
                if stidis in sti_pars.keys():
                    this_stitype_pars = sti_pars[stidis]
                final_pars = sc.mergedicts(sti_main_pars, this_stitype_pars)
                stis += sti.make_sti(stidis, pars=final_pars)

            elif isinstance(stidis, ss.Disease):
                stis += stidis
            else:
                raise ValueError(f"Invalid STI type: {type(stidis)}. Must be str or sti.BaseSTI.")

        return stis

    def process_connectors(self):
        """
        Get the default connectors for the diseases in the simulation.
        Connectors are loaded based on the disease names or modules provided in the format <d1>_<d2>.
        """
        connectors = []
        parsed_diseases = []
        add_connectors = False

        if isinstance(self.pars.connectors, bool) and self.pars.connectors:
            errormsg = ('STIsim does not currently support automatically adding connectors. This feature'
                        ' will be added in a future release. For the time being, please add connectors '
                        'manually, e.g. by passing sti.hiv_syph(hiv, syphilis) to the sim.')
            raise NotImplementedError(errormsg)

            # TODO: The remaining code will be re-enabled once debugged
        #     self.pars['connectors'] = ss.ndict()  # Reset
        #     add_connectors = True
        #
        # if add_connectors:
        #     for disease in self.pars['diseases']:
        #         if isinstance(disease, str):
        #             parsed_diseases.append(disease.lower())
        #         if isinstance(disease, ss.Module):
        #             parsed_diseases.append(disease.name.lower())
        #
        #     # sort the diseases and then get all combinations of their pairs
        #     disease_pairs = combinations(parsed_diseases, 2)
        #
        #     # TODO: this does not quite work as intended because the ordering matters...
        #     for (d1, d2) in disease_pairs:
        #         try:
        #             connector = getattr(sti, f'{d1}_{d2}')
        #         except:
        #             try:
        #                 connector = getattr(sti, f'{d2}_{d1}')
        #             except:
        #                 continue
        #         connectors.append(connector(d1, d2))
        return connectors

    def case_insensitive_getattr(self, searchspace, attrname):
        """
        Find a class in the given package that matches the name, ignoring case.

        Args:
            searchspace (list): A list of classes or modules to search through.
            attrname (str): The name of the attribute to find, case-insensitive.
        """
        if not isinstance(searchspace, list):
            searchspace = [searchspace]

        for classname in searchspace:
            for attr in dir(classname):
                if attr.lower() == attrname.lower():
                    return getattr(classname, attr)
        return None