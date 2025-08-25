'''
Load data
'''

#%% Housekeeping
import pandas as pd
import sciris as sc
import unicodedata
import re
import os


__all__ = ['get_age_distribution', 'get_rates']


thisdir = sc.thispath(__file__)
filesdir = thisdir / 'files'
files = sc.objdict()
files.age_dist_sex = 'populations_by_sex.obj'
files.birth = 'births.csv'
files.asfr = 'asfr.csv'
files.death = 'deaths.csv'
files.migration = 'migration.csv'


# Cache data as a dict
cache = dict()


def sanitizestr(string=None, alphanumeric=True, nospaces=True, asciify=True, lower=True, spacechar='_', symchar='_'):
    ''' Remove all non-printable characters from a string -- to be moved to Sciris eventually '''
    string = str(string)
    if asciify:
        string = unicodedata.normalize('NFKD', string).encode('ascii', 'ignore').decode()
    if nospaces:
        string = string.replace(' ', spacechar)
    if lower:
        string = string.lower()
    if alphanumeric:
        string = re.sub('[^0-9a-zA-Z ]', symchar, string)
    return string


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

