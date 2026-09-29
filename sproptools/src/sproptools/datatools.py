#! /usr/bin/env python3

import numpy as np
import pandas as pd

"""
Tools for operations on pandas dataframes related to cosmic-ray data.
"""

def for_timestamp(data, timestamp):
    """
    Select data valid for given timestamp (e.g., a np.datetime64 value).

    Returns a pandas dataframe.
    """

    if 'time_min' in data.keys() and 'time_max' in data.keys():

        selection = (timestamp >= data['time_min']) & (timestamp <= data['time_max'])
        df = data[selection]

    elif 'date' in data.keys():
        selection = (timestamp >= data['date']) & (timestamp <= data['date'])
        df = data[selection]

    else:
        raise ValueError('Cannot determine time ranges due to invalid column names!')

    if df.empty:
            raise ValueError(f'Cannot find data for timestamp {timestamp}!')

    df.attrs = data.attrs
    return df


def for_timerange(data, timerange):
    """
    Select data valid for given time range (e.g., a pair of np.datetime64 values).

    Returns a pandas dataframe.
    """

    if 'time_min' in data.keys() and 'time_max' in data.keys():
        # time_min and time_max are assumed to be exact time stamps, so we use '<=' here ...
        selection = (data['time_min'] >= timerange[0]) & (data['time_max'] <= timerange[1])
        df = data[selection]

    elif 'date' in data.keys():
        # ... but 'date' will only be exact to 1 day, to we use '<'
        selection = (data['date'] >= timerange[0]) & (data['date'] < timerange[1])
        df = data[selection]

    else:
        raise ValueError('Cannot determine time ranges due to invalid column names!')

    if df.empty:
            raise ValueError(f'Cannot find data for timestamp {timerange}!')

    df.attrs = data.attrs
    return df


def timerange(data):
    """
     Return overall time range covered by the dataframe as pair of numpy datetime64 values.
    """
    if 'time_min' in data.keys() and 'time_max' in data.keys():
        time_min = sorted(data['time_min'].unique())[0]
        time_max = sorted(data['time_max'].unique())[-1]
        return time_min, time_max

    if 'date' in data.keys():
        dates = sorted(data['date'].unique())
        return dates[0], dates[-1]

    raise ValueError('Cannot determine time range!')


def find_unit(data, quantity):
    """
    Determine unit for given quantity from the metadata of the dataframe.
    """

    if not 'units' in data.attrs:
        raise ValueError('No units stored in metadata of dataframe!')

    if quantity in data.attrs['units']:
        return data.attrs['units'][quantity]
    if '_min' in quantity:
        q = quantity.replace('_min','')
        if q in data.attrs['units']:
            return data.attrs['units'][q]
    if '_max' in quantity:
        q = quantity.replace('_max','')
        if q in data.attrs['units']:
            return data.attrs['units'][q]
    return None


def from_min_max(data, quantity, loc='geometric_mean', lw_spectral_index=3.0):
    """
    Calculate new series of data points, where x values are averaged based on x_min and x_max.

    For example, if the dataset only contains bin edges "rigidity_min" and "rigidity_max",
    a series of rigidity values for plotting can be calculated as `from_min_max('rigidity')`.

    The location of the new data points is determined by the `loc` argument:
      - geometric_mean: geometric mean of left and right bin edges
      - arithmetic_mean: arithmetic mean of left and right bin edges
      - lafferty_wyatt: find characteristic bin rigidities using method of Lafferty & Wyatt (1995), eq.(6),
        assuming spectral index (x^(-gamma)) specified by the argument `lw_spectral_index`.

    """
    xmin = f'{quantity}_min'
    xmax = f'{quantity}_max'
    if xmin not in data.keys() or xmax not in data.keys():
        raise KeyError(f'Data file does not contain {xmin} and {xmax} values!')
    qmin = data[xmin]
    qmax = data[xmax]

    if loc == 'geometric_mean':
        return np.sqrt(qmin * qmax)
    if loc == 'arithmetic_mean':
        return 0.5 * (qmin + qmax)
    if loc == 'lafferty_wyatt':

        gamma = lw_spectral_index

        if gamma == 1.0:
            return (qmax - qmin) / np.log(qmax / qmin)
        else:
            denom = np.power(qmax, 1.0 - gamma) - np.power(qmin, 1.0 - gamma)
            num = (qmax - qmin) * (1.0 - gamma)
            return np.power(num / denom, 1.0 / gamma)

    raise ValueError(f'Illegal location algorithm "{loc}"!')


def time_averaged(data, by, columns, weights_from):
    """
    Calculate weighted time averages of data quantities, grouped by selected variable.

    Typical usage:

        avg_data = time_averaged(data, by='rigidity_min', columns=['flux', 'syst_error'], weights_from='stat_error')

    or you can even average over Bartels rotations for the full data set like this:

        data['bartels_rotation'] = date_to_br(data['date'])
        avg_data = time_averaged(data, by=['bartels_rotation', 'rigidity_min'], columns=['flux', 'syst_error'], weights_from='stat_error')

    """

    df = data.copy()

    # set new time range
    data_timerange = timerange(data)
    df['time_min'] = data_timerange[0]
    df['time_max'] = data_timerange[1]

    # everything not used for the averaging will be copied over, except date, which will be obsolete after averaging
    copy_over = [k for k in df.keys() if k != by and k != weights_from and k not in columns and k != 'date']
    agg_dict = {k: 'first' for k in copy_over}

    # weights for weighted means
    agg_dict['w'] = 'sum'
    df['w'] = 1. / df[weights_from]**2

    # temporary columns: weight x quantity
    for col in columns:
        df[f'w_x_{col}'] = df['w'] * df[col]
        agg_dict[f'w_x_{col}'] = 'sum'

    # calculate weighted sums and sum of weights
    df = df.groupby(by).aggregate(agg_dict)

    # get uncertainty on main quantity (usually flux)
    df[weights_from] = 1. / np.sqrt(df['w'])

    # calculate final weighted means and remove temporary columns
    for col in columns:
        df[col] = df[f'w_x_{col}'] / df['w']
        del df[f'w_x_{col}']
    del df['w']

    # delete now-obsolete date
    if 'units' in df.attrs:
        if 'date' in df.attrs['units']:
            del df.attrs['units']['date']

    return df
