"""
Download data needed for STIsim.

To use this, you will need to get an auth_key from the UN Data Portal API.

"""

import os
import numpy as np
import pandas as pd
import sciris as sc
import requests


# Set up
base_url = "https://population.un.org/dataportalapi/api/v1/"
thisdir = sc.thispath(__file__)
filesdir = thisdir / 'files'
files = sc.objdict()
files.country_codes = 'locations.csv'
auth_key_file = thisdir / 'auth_key.txt'


def _authkeyerror(location):
    """ Base error message when no authorization key is present """
    return f'''
To make demographics for {location}, there are two options.

1. STIsim can pull demographic data from UN WPP if you register for an authorization
    token with them. You can do this by sending an email to population@un.org with the
    Subject 'Data Portal Token Request (https://population.un.org/dataportalapi/index.html).
    Once you receive the authorization token, you need to save it in a file named 
    "auth_key.txt", which you will need to save in the folder "stisim/data/files". 
    
2. You can manually download datafiles and save them in the folder "stisim/data/files".
    These should be named as "{location}_deaths.csv" for background death rates, 
    "{location}_asfr.csv" for age-specific fertility, and "{location}_births.csv"
    if you want to use births instead of fertility. Refer to the files within the 
    folder "tests/test_data/", including "zimbabwe_asfr.csv, "zimbabwe_births.csv", 
    and "zimbabwe_deaths.csv" to see how these files should be formatted. You can 
    obtain global data from "https://population.un.org/dataportal/home".
'''


def get_auth_key(location):
    # Read in auth_key
    have_auth_key = os.path.exists(auth_key_file)
    if not have_auth_key:
        errormsg = _authkeyerror(location)
        raise ValueError(errormsg)
    with open(auth_key_file, 'r') as f:
        for line in f:
            line = line.strip()
            auth_key = line
    return auth_key


def get_filename(location, indicator, year=None):
    """
    Get the filename for a given location, indicator, and year.
    """
    if year is None:
        return filesdir / f'{location}_{indicator}.csv'
    else:
        return filesdir / f'{location}_{indicator}_{year}.csv'


indicator_codes = dict(
    deaths=80,  # Mortality
    asfr=68,  # Fertility
    age=47,  # Population size by age
    pop=49,  # Total population
    cbr=55,  # Crude birth rate
    cmr=59,  # Crude mortality rate
    migration=65,  # Net migration
)


def check_downloaded(location, indicators=None, year=None, verbose=1):
    """
    Check if data is downloaded for this location
    Args:
        location (str): the location to check
        indicators (list): the indicators to check for; if None, checks all
        verbose (int): detail to print (0 = none, 1 = reason for failure, 2 = everything)
    """

    # Do file checks
    exists = dict()
    if indicators is None:
        indicators = indicator_codes.keys()

    for key in indicators:
        if key == 'age':
            fn = get_filename(location, key, year)
        else:
            fn = get_filename(location, key)
        exists[key] = os.path.exists(fn)
        if verbose > 1:
            print(f'STIsim data: checking {fn}: {exists[key]}')
    ok = all(list(exists.values()))
    if not ok and verbose:
        print(f'STIsim data: at least one file missing: {exists}')
    elif ok and verbose > 1:
        print('STIsim data: all files exist')

    missing = [key for key, val in exists.items() if not val]

    return ok, missing


def get_country_codes(auth_key):
    country_code_file = filesdir / files.country_codes
    have_file = os.path.exists(country_code_file)
    if not have_file:
        try:
            get_locations(auth_key)
        except:
            raise ValueError('Could not get country codes')
    df = pd.read_csv(country_code_file)
    return df


def get_country_aliases():
    """ Define aliases for countries """
    country_mappings = {
       'Bolivia':        'Bolivia (Plurinational State of)',
       'Burkina':        'Burkina Faso',
       'Cape Verde':     'Cabo Verdeo',
       'Hong Kong':      'China, Hong Kong Special Administrative Region',
       'Macao':          'China, Macao Special Administrative Region',
       "Cote d'Ivoire":  "Côte d'Ivoire",
       "Cote dIvoire":   "Côte d'Ivoire",
       "Ivory Coast":    "Côte d'Ivoire",
       'DRC':            'Democratic Republic of the Congo',
       'Congo':          'Congo, Rep.',
       'Iran':           'Iran (Islamic Republic of)',
       'Laos':           "Lao People's Democratic Republic",
       'Micronesia':     'Micronesia (Federated States of)',
       'Korea':          'Republic of Korea',
       'South Korea':    'Republic of Korea',
       'Moldova':        'Republic of Moldova',
       'Russia':         'Russian Federation',
       'Palestine':      'State of Palestine',
       'Syria':          'Syrian Arab Republic',
       'Taiwan':         'Taiwan Province of China',
       'Macedonia':      'The former Yugoslav Republic of Macedonia',
       'UK':             'United Kingdom of Great Britain and Northern Ireland',
       'United Kingdom': 'United Kingdom of Great Britain and Northern Ireland',
       'Tanzania':       'United Republic of Tanzania',
       'USA':            'United States of America',
       'United States':  'United States of America',
       'Venezuela':      'Venezuela (Bolivarian Republic of)',
       'Vietnam':        'Viet Nam',
        }

    return country_mappings


def get_country_code(df, location):
    """
    Find a match between the data file and the provided location.
    """

    if sc.checktype(df, pd.DataFrame):
        countries = [l.lower() for l in df.name.values]
    else:
        errormsg = 'Type not supported'
        raise ValueError(errormsg)

    mapping = get_country_aliases()
    mapping = {key.lower(): val.lower() for key, val in mapping.items()}

    lloc = location.lower()
    if lloc not in countries and lloc in mapping:
        lloc = mapping[lloc]
    try:
        ind = countries.index(lloc)
        country_code = list(df['id'].values)[ind]
    except ValueError as E:
        suggestions = sc.suggest(location, countries, n=4)
        if suggestions:
            errormsg = f'Location "{location}" not recognized, did you mean {suggestions}? ({str(E)})'
        else:
            errormsg = f'Location "{location}" not found in data ({str(E)})'
        raise ValueError(errormsg)

    return country_code


def download_data(location, indicators=None, start=1950, stop=2100, step=10):
    """ Download data """

    # Get auth_key
    auth_key = get_auth_key(location)

    # Process indicators
    if indicators is None:
        indicators = indicator_codes.keys()

    # Get country code
    df = get_country_codes(auth_key)
    country_code = get_country_code(df, location)

    for indicator in indicators:
        icode = indicator_codes[indicator]
        dfs = sc.autolist()
        if indicator == 'age': years = [start]
        else: years = np.arange(start, stop, step)
        for year in years:
            target = f"data/indicators/{icode}/locations/{country_code}/start/{year}/end/{year}"
            dfs += get_data(auth_key, target, do_save=False)
        df = pd.concat(dfs)
        df = df.loc[df.variant == 'Median']
        label_dict = {'timeLabel': 'Time', 'sex': 'Sex', 'ageStart': 'AgeStart', 'value': 'Value'}
        for ol, nl in label_dict.items():
            df = df.rename(columns={ol: nl})
        df = df[label_dict.values()]
        if indicator == 'deaths':
            df = df.loc[df.Sex != 'Both sexes']
        if indicator == 'age':
            df = df.loc[df.Sex == 'Both sexes']
            fn = get_filename(location, indicator, year)
        else:
            fn = get_filename(location, indicator)
        df.to_csv(fn)

        print(f'Successfully downloaded {indicator} for {location}!')

    return


def get_indicators(base_url, auth_key, do_save=True):
    return get_available(base_url, auth_key, 'indicators/', do_save=do_save)


def get_locations(auth_key, do_save=True):
    return get_available(base_url, auth_key, 'locations/', do_save=do_save)


def get_data(auth_key, which, do_save=False):
    df = get_available(base_url, auth_key, which, filename='mortality', do_save=do_save)
    return df


def get_available(base_url, auth_key, which, filename=None, do_save=True):
    """ Save a list of indicators or locations to a csv file """
    first_url = base_url + which
    payload = {}
    headers = {'Authorization': auth_key}
    response = requests.request("GET", first_url, headers=headers, data=payload)
    j = response.json()

    # Set up list of dataframes
    dfs = sc.autolist()
    dfs += pd.json_normalize(j['data'])

    while j['nextPage'] != None:
        pageno = j['pageNumber']
        print(f'Processing page {pageno}')
        new_target = j['nextPage'].split('?')[1]
        new_url = base_url + which + '/?' + new_target
        response = requests.request("GET", new_url, headers=headers, data=payload)
        j = response.json()
        df_temp = pd.json_normalize(j['data'])
        dfs += df_temp

    df = pd.concat(dfs)
    # if do_save:
    #     if filename is None: filename = which.strip('/')
    #     df.to_csv(filename+'.csv')
    #
    return df


