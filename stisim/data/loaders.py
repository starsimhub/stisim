'''
Load data
'''

#%% Housekeeping
import pandas as pd
import sciris as sc


__all__ = ['get_age_distribution', 'get_rates']


thisdir = sc.thispath(__file__)
filesdir = thisdir / 'files'
files = sc.objdict()
files.age_dist_sex = 'populations_by_sex.obj'
files.birth = 'births.csv'
files.asfr = 'asfr.csv'
files.death = 'deaths.csv'
files.migration = 'migration.csv'


def get_age_distribution(location=None, year=None, datafolder=None):
    """ Load age distribution for a given location & year"""
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename= f'{location}_age_{year}.csv', folder=datafolder)
    try:
        raw_df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for age distribution.'
        raise ValueError(errormsg) from E

    # Remap column names
    raw_df = raw_df.rename(columns={'Time': 'year', 'AgeStart': 'age', 'Value': 'value'})
    raw_df = raw_df[['year', 'age', 'value']]

    return raw_df


def get_rates(location, which, datafolder=None):
    """ Load birth/death/fertility/migration rates for a given location """
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename=f'{location}_{files[which]}', folder=datafolder)
    try:
        df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for {which}.'
        raise ValueError(errormsg) from E
    return df

