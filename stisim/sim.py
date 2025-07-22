import pandas as pd
import starsim as ss

import stisim as sti
import sciris as sc
from itertools import combinations


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

    def __init__(self, pars=None, sim_pars=None, sti_pars=None, nw_pars=None, label=None, people=None, demographics=None, diseases=None, networks=None,
                 interventions=None, analyzers=None, connectors=None, copy_inputs=True, data=None, **kwargs):

        # Inputs and defaults
        self.nw_pars = None     # Parameters for the networks - processed later
        self.sti_pars = None    # Parameters for the STIs - processed later
        self.pars = None        # Parameters for the simulation - processed later
        # Call the constructor of the parent class WITHOUT pars or module args
        super().__init__(pars=None, label=label)
        self.pars = sti.make_sim_pars()  # Make default parameters using values from parameters.py

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
        # Marge in pars and kwargs
        all_pars = sc.mergedicts(pars, sim_pars, sti_pars, nw_pars, sim_kwargs, kwargs)

        # Deal with sim pars
        user_sim_pars = {k: v for k, v in all_pars.items() if k in self.pars.keys()}
        for k in user_sim_pars: all_pars.pop(k)
        sim_pars = sc.mergedicts(user_sim_pars, sim_pars, _copy=True)  # Don't merge with defaults, those are set above

        # Deal with STI pars
        # TODO

        # Deal with network pars
        default_nw_pars = sti.make_network_pars()
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

    def init(self, force=False, **kwargs):
        """
        Perform all initializations for the sim
        """
        # Process the STIs
        stis = self.process_stis()
        self.pars['diseases'] += stis

        # Process the network
        self.pars['networks'] = self.process_networks()

        # Process the demographics
        demographics, people, total_pop = self.process_demographics()
        self.pars['demographics'] += demographics
        self.pars['people'] = people
        self.pars['total_pop'] = total_pop

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
                ss.Maternal(pars=self.nw_pars),
            )
        return networks

    def process_demographics(self):
        """ Process the location to create people and demographics if not provided. """

        # If it's a string, do lots of work
        if sc.checktype(self.pars['demographics'], str):
            location = self.pars.pop('demographics')
            self.pars['demographics'] = ss.ndict()
            birth_rates, death_rates = self.get_births_deaths(location)

            # Create modules for births and deaths
            births = ss.Births(birth_rate=birth_rates, metadata=dict(data_cols=dict(year='year', value='cbr')))
            deaths = ss.Deaths(death_rate=death_rates)

            try:
                age_data = hpdata.get_age_distribution(location, year=self.pars.start)
                pop_trend = hpdata.get_total_pop(location)
                pop_age_trend = hpdata.get_age_distribution_over_time(location)
                total_pop = int(age_data.value.sum())*1e3
            except ValueError as E:
                warnmsg = f'Could not load age data for requested location "{location}" ({str(E)}); using default'
                raise ValueError(warnmsg) from E

            people = ss.People(self.pars.n_agents, age_data=age_data)
            demographics = [births, deaths]

        else:
            demographics = self.pars['demographics']
            people = self.pars['people']
            total_pop = self.pars['total_pop']

        return demographics, people, total_pop

    @staticmethod
    def get_births_deaths(location, verbose=1, by_sex=True, overall=False):
        '''
        Get mortality and fertility data by location if provided

        Args:
            location (str):  location
            verbose (bool):  whether to print progress
            by_sex   (bool): whether to get sex-specific death rates (default true)
            overall  (bool): whether to get overall values ie not disaggregated by sex (default false)

        Returns:
            lx (dict): dictionary keyed by sex, storing arrays of lx - the number of people who survive to age x
            birth_rates (arr): array of crude birth rates by year
        '''

        if verbose:
            print(f'Loading location-specific demographic data for "{location}"')
        try:
            death_rates = hpdata.get_death_rates(location=location, by_sex=by_sex, overall=overall)
            birth_rates = hpdata.get_birth_rates(location=location)
            return birth_rates, death_rates
        except ValueError as E:
            warnmsg = f'Could not load demographic data for requested location "{location}" ({str(E)})'
            print(warnmsg)

    @staticmethod
    def init_modules(mod_type):
        # init diseases
        initialized_modules = []
        mod_type_notplural = mod_type[:-1]  # Remove the 's' at the end of the type
        pars_name = f'{mod_type_notplural}_pars'  # e.g. 'disease_pars' or 'network_pars'

        if mod_type in all_input_args:
            # if the inputs are not a list, e.g. a single disease or network, convert to a list
            if not isinstance(all_input_args[mod_type], list):
                all_input_args[mod_type] = [all_input_args[mod_type]]

            for mod in all_input_args[mod_type]:
                mod_getter = getattr(self, f'get_{mod_type_notplural}')
                mod_name = mod.name if isinstance(mod, ss.Module) else mod.lower()  # Get the name of the module

                pars = None
                if pars_name in all_input_args:
                    # If there are parameters for the disease or network, pop them
                    pars = all_input_args[pars_name].pop(mod_name, None)

                initialized_mod = mod_getter(mod, pars)  # Look up the disease module
                if isinstance(initialized_mod, list):
                    initialized_modules.extend(initialized_mod)
                else:
                    initialized_modules.append(initialized_mod)

        elif mod_type in self.pars:
            # If diseases are defined in the pars, use those (but only if they are not already in the input args!)
            initialized_modules = self.pars.pop(mod_type)

        all_input_args[mod_type] = initialized_modules

        if pars_name in all_input_args:
            for mod in all_input_args[mod_type]: #and initialized_mod.name in all_input_args[pars_name]:
                if mod.name in all_input_args[mod_type]:
                    # if the disease or network is already initialized, update its parameters
                    mod.update_pars(all_input_args[pars_name][mod.name])  # Update disease parameters if provided
            all_input_args.pop(pars_name)

        init_modules('diseases')
        init_modules('networks')
        init_modules('demographics')
        init_modules('interventions')
        init_modules('analyzers')

        # init connectors
        if all_input_args['connectors'] is not None and all_input_args['connectors'] == True:
            all_input_args['connectors'] = self.get_connectors(all_input_args['diseases'])

        input_pars = sc.mergedicts(pars, all_input_args, _copy=copy_inputs)

        super().__init__(pars=input_pars, data=data)
        return

    def process_stis(self):
        """
        Look up a disease by its name and return the corresponding module.
        """

        self.pars['genotypes'] = sc.tolist(self.pars['genotypes'])  # Ensure it's a list
        genotypes = sc.autolist()

        dis_dict = {
            'bv': sti.BV,
            'ct': sti.Chlamydia,
            'ctbl': sti.ChlamydiaBL,
            'ng': sti.Gonorrhea,
            'hiv': sti.HIV,
            'syphilis': sti.Syphilis,
            'tv': sti.Trichomoniasis,
        }

        if isinstance(disease, str):
            disease = disease.lower()
            if disease in dis_dict:
                return dis_dict[disease](pars=disease_pars)
            else:
                raise ValueError(f"Disease '{disease}' not found in STIsim diseases.")
        elif isinstance(disease, ss.Module):
            disease.update_pars(pars=disease_pars)
            return disease

        else:
            raise TypeError("Disease name must be a string or a starsim Module instance.")

    def get_demographic(self, demog, demog_pars=None):
        """
        Get the demographics modules for the simulation.
        """
        initialized_demographics = []

        # if the input is a string, assume it is a location name and load the corresponding demographics
        if isinstance(demog, str):
            # try to load the files
            try:
                location = demog
                path = ''
                if demog_pars is not None and 'data' in demog_pars:
                    path = demog_pars['data']
                fertility_rates = pd.read_csv(f'{path}{location}_asfr.csv')
                death_rates = pd.read_csv(f'{path}{location}_deaths.csv')

                initialized_demographics.append(sti.Pregnancy(pars={'fertility_rate': fertility_rates}, ))
                initialized_demographics.append(ss.Deaths(pars={'death_rate': death_rates, 'rate_units': 1}))
            except:
                print("Warning: Location demographic data files not found. Assuming demographic module name instead.")

                # if the files are not found, try to load the demographic from the stisim or starsim modules
                match = self.case_insensitive_getattr([sti, ss], demog)
                if match is not None:
                    demog = match(pars=demog_pars)
                    initialized_demographics.append(demog)
                else:
                    raise ValueError(f"Demographic module '{demog}' not found in STIsim demographics or data files.")

            return initialized_demographics
        else:
            demog.update_pars(demog_pars)  # Update demographic parameters if provided
            return demog


    def get_network(self, network, network_pars=None):
        """
        Look up a network by its name and return the corresponding module.
        """
        if isinstance(network, str):
            network = self.case_insensitive_getattr([sti, ss], network)
            if network is None:
                raise ValueError(f"Network '{network}' not found in STIsim networks.")
            return network(network_pars)  # Initialize the network class

        elif isinstance(network, ss.Network):
            network.update_pars(network_pars)
            return network
        else:
            raise TypeError("Network must be a string or a starsim Network instance.")



    def get_connectors(self, diseases):
        """
        Get the default connectors for the diseases in the simulation.
        Connectors are loaded based on the disease names or modules provided in the format <d1>_<d2>.

        Args:
            diseases (list): A list of disease names or disease modules to connect.
        """
        connectors = []
        parsed_diseases = []
        for disease in diseases:
            if isinstance(disease, str):
                parsed_diseases.append(disease.lower())
            if isinstance(disease, ss.Module):
                parsed_diseases.append(disease.name.lower())

        # sort the diseases and then get all combinations of their pairs
        disease_pairs = combinations(parsed_diseases, 2)

        for (d1, d2) in disease_pairs:
            try:
                connector = getattr(sti, f'{d1}_{d2}')
            except:
                try:
                    connector = getattr(sti, f'{d2}_{d1}')
                except:
                    continue
            connectors.append(connector(d1, d2))
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