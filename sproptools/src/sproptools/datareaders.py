#! /usr/bin/env python3

import collections
import glob
import os
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from .bartelsrotation import br_to_date

"""
Tools for reading cosmic-ray data files in different formats:

- Simple text files with two-column data: `read_xy`

- Sunspot number file in SIDC format: `read_ssn_file`

- Files downloaded from ASI database at https://tools.ssdc.asi.it/CosmicRays/:
  - XML files: `read_asi_file`
  - A time series of XML files: `read_asi_timeseries`
  - CSV files for time-dependent data: `read_asi_csv_timeseries`

- Files downloaded from the CRDB database at https://lpsc.in2p3.fr/crdb:
  - USINE files: `read_crdb_file`

Data from ASI and CRDB files is returned as a pandas dataframe with metadata attached.

"""

def read_xy(filename, encoding='utf-8', dtype=float):
    '''
    Read a simple text file with two columns (e.g., energy vs flux).
    '''

    data = np.genfromtxt(filename, encoding=encoding, dtype=dtype)
    x = data[:, 0]
    y = data[:, 1]

    return x, y


def read_ssn_file(filename, encoding='utf-8'):
    '''
    Extract smoothed sunspot number data from file in SIDC format
    (e.g., the file in solarprop's data directory)

    Returns array of np.datetim64 objects and array of corresponding SSN values.
    '''

    data = np.genfromtxt(filename, encoding=encoding, usecols=(0,1,2,3), dtype='i,i,f,f', unpack=True)
    year = data[0]
    month = data[1]
    ssn = data[3]

    datetimes = np.array([np.datetime64(f'{y}-{m:02}') for y,m in zip(year, month)])

    return datetimes, ssn


def read_asi_csv_timeseries(csvfile):
    """
    Convert data from ASI time-dependent csv format to pandas dataframe.

    csvfile: Filename of input csv file.
    """

    attr = {}

    # parse particle and mission from filename
    filenametokens = os.path.splitext(os.path.split(csvfile)[-1])[0].split('_')
    if len(filenametokens) > 1:
        attr['particle'] = filenametokens[0]
        attr['mission'] = filenametokens[1]

    df = pd.read_csv(csvfile)

    # extract units and rename columns
    # assumption: each original header has the form "quantity unit"
    headers = [h for h in df.keys() if h not in ['bartels_rotation_number']]
    units = {x: y for x, y in [h.strip().split(' ', 1) for h in headers]}
    df.rename(columns={h: q for h, q in zip(headers, units.keys())}, inplace=True)

    # convert date string to datetime objects
    if 'date' in df.keys():
        df['date'] = pd.to_datetime(df['date'])

    # BR -> time range
    if 'bartels_rotation_number' in df.keys():
        df = df.assign(time_min=lambda x: br_to_date(x['bartels_rotation_number']),
                       time_max=lambda x: br_to_date(x['bartels_rotation_number'] + 1))

    attr['units'] = units
    df.attrs = attr
    return df


def read_asi_timeseries(filepattern):
    """
    Read a series of ASI xml files. Returns a single pandas dataframe.

    filepattern: pattern for file names to read (may contain wildcards)
    """

    results = []
    asifiles = sorted(glob.glob(filepattern))
    if not asifiles:
        raise ValueError('No files to read!')

    for asifile in asifiles:
        results.append(read_asi_file(asifile))

    missions = set(x.attrs['mission'] for x in results)
    if len(missions) != 1:
        raise ValueError(f'Mission not unique: {missions}')
    particles = set(x.attrs['particle'] for x in results)
    if len(particles) != 1:
        raise ValueError(f'Particle not unique: {particles}')
    # FIXME make sure units are identical

    return pd.concat(results)


def read_asi_file(inputfile):
    '''
    Parse XML file downloaded from ASI cosmic-ray database.

    Return data as a pandas dataframe, with metadata stored in the `attrs` attribute of the dataframe.

    The following simplified quantities are calculated automatically if possible:
       - statistical and systematical errors (arithmetic mean of the upper and lower error bars),
       - total error (quadratic sum of statistical and systematical errors).

    Usage example:

        # read XML file downloaded from ASI database
        data = read_asi_file(os.path.expandvars('$OPENDATA/data/asifile.xml'))
        # pretty-print data
        print(data)
        # example for data access
        plt.errorbar(data['rigidity'], data['flux'], yerr=data['flux_total_error'])

    '''

    # initialize XML parser
    root = ET.parse(inputfile).getroot()

    attr = {}

    # read header information
    header = root.find('HEADER')
    attr['mission'] = header.findtext('MISSION')
    attr['particle'] = header.findtext('PARTICLE')
    attr['units'] = {}

    units = header.find('UNITS')
    for x in units.iter():
        if x.tag and x.tag != 'UNITS':
            attr['units'][x.tag] = x.text

    # extract individual data points
    datapoints = collections.defaultdict(list)

    for d in root.findall('DATA'):

        for x in d.iter():
            if x.tag and x.tag != 'DATA':
                try:
                    # try interpreting data as float ...
                    datapoints[x.tag].append(float(x.text))
                except ValueError:
                    try:
                        # ... or as datetime64 ...
                        datapoints[x.tag].append(np.datetime64(x.text))
                    except ValueError:
                        # ... or as string
                        datapoints[x.tag].append(x.text)

    # create dataframe
    df = pd.DataFrame(datapoints)

    # calculate simplified errors (average of "low" and "high" errors)
    quantities = list(df.keys())
    for q in quantities:
        if 'error' in q and '_low' in q:
            errorname = q.replace('_low','')
            if f'{errorname}_high' in datapoints.keys():
                df = df.assign(**{errorname: lambda x: 0.5 * (np.abs(x[q]) + np.abs(x[f'{errorname}_high']))})

    # calculate total error if possible
    quantities = list(df.keys())
    for q in quantities:
        if q.endswith('_statistical_error'):
            systname = q.replace('_statistical_', '_systematical_')
            if systname in quantities:
                errorname = q.replace('_statistical_', '_total_')
                df = df.assign(**{errorname: lambda x: np.sqrt(x[q]**2 + x[systname]**2)})

    df.attrs = attr
    return df



# Format: USINE code
#   Col.1  -  QUANTITY NAME (case insensitive)
#   Col.2  -  SUB-EXP NAME (case insensitive, no space)
#   Col.3  -  EAXIS TYPE: EKN, EK, R, or ETOT
#   Col.4  -  <E>: mean value bin [GeV/n, GeV, GV, or GeV]
#   Col.5  -  EBIN_LOW
#   Col.6  -  EBIN_HIGH
#   Col.7  -  QUANTITY VALUE: [#/sr/s/m2/EAxis] if flux , no unit if ratio
#   Col.8  -  ERR_STAT-
#   Col.9  -  ERR_STAT+
#   Col.10 -  ERR_SYST-
#   Col.11 -  ERR_SYST+
#   Col.12 -  ADS URL FOR PAPER REF (no space)
#   Col.13 -  phi [MV]
#   Col.14 -  DISTANCE EXP IN SOLAR SYSTEM [AU]
#   Col.15 -  DATIMES: format = yyyy/mm/dd-hhmmss:yyyy/mm/dd-hhmmss;...
#   Col.16 -  IS UPPER LIMIT: format = 0 or 1

def read_crdb_file(inputfile):
    '''
    Parse text file downloaded from CRDB cosmic-ray database.

    Return data and metadata as an instance of CrdbData.
    '''

    d = np.genfromtxt(inputfile, dtype='U20,U50,U5,f,f,f,f,f,f,f,f,U40,f,f,U40,i', unpack=True)

    attr = {}
    attr['mission'] = str(d[1][0])
    attr['particle'] = str(d[0][0])
    attr['ads_url'] = str(d[11][0])

    eaxis_type = str(d[2][0])

    default_units = { 'EKN': 'GeV/n', 'EK': 'GeV', 'R': 'GV', 'ETOT': 'GeV' }
    default_units['value'] = f'#/sr/s/m2/{default_units[eaxis_type]}'
    attr['units'] = { eaxis_type: default_units[eaxis_type], 'value': default_units['value'] }
    
    mean_ebin = d[3]
    ebin_low = d[4]
    ebin_high = d[5]
    value = d[6]
    err_stat_minus = d[7]
    err_stat_plus = d[8]
    err_syst_minus = d[9]
    err_syst_plus = d[10]
    is_upper_limit = d[15]

    err_stat = 0.5 * (np.abs(err_stat_minus) + np.abs(err_stat_plus))
    err_syst = 0.5 * (np.abs(err_syst_minus) + np.abs(err_syst_plus))
    err_total = np.sqrt(err_stat**2 + err_syst**2)

    datimes = [str(x).split(':') for x in d[14]]
    time_min = pd.to_datetime([x[0] for x in datimes], format='%Y/%m/%d-%H%M%S').astype('datetime64[s]')
    time_max = pd.to_datetime([x[1] for x in datimes], format='%Y/%m/%d-%H%M%S').astype('datetime64[s]')
    
    df = pd.DataFrame({f'{eaxis_type}_min': ebin_low,
                       f'{eaxis_type}_max': ebin_high,
                       f'{eaxis_type}_mean': mean_ebin,
                       'value': value,
                       'err_stat_minus': np.abs(err_stat_minus),
                       'err_stat_plus': err_stat_plus,
                       'err_stat': err_stat,
                       'err_syst_minus': np.abs(err_syst_minus),
                       'err_syst_plus': err_syst_plus,
                       'err_syst': err_syst,
                       'err_total': err_total,
                       'is_UL': is_upper_limit,
                       'time_min': time_min,
                       'time_max': time_max,
                       })
    
    # check consistency: quantity, experiment, eaxis, times are the same in all rows
    if list(d[0]).count(attr['particle']) != len(df):
        raise ValueError('Inconsistency in particle names!')
    if list(d[1]).count(attr['mission']) != len(df):
        raise ValueError('Inconsistency in mission names!')
    if list(d[2]).count(eaxis_type) != len(df):
        raise ValueError('Inconsistency in eaxis type!')
    if list(d[11]).count(attr['ads_url']) != len(df):
        raise ValueError('Inconsistency in ADS URLs!')

    df.attrs = attr
    return df
