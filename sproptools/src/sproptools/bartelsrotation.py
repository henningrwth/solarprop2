#! /usr/bin/env python3

"""
Functions for converting unix time stamps to Bartels Rotation numbers, and vice versa.
"""

import numpy as np

def date_to_br(date):
    """
    Calculate Bartels rotation number from given date.

    Parameters
    ----------

    date : datetime
           Input date(s).


    Returns
    -------
    float
        Bartels rotation number (floating-point).

    """

    return (np.array(date).astype('datetime64[D]') - np.datetime64('1832-02-08')) / np.timedelta64(27, 'D') + 1.0


def br_to_date(br):
    """
    Calculate start date for given Bartels rotation number.
    """

    return np.datetime64('1832-02-08') + (br - 1) * np.timedelta64(27, 'D')


def br_to_timerange(br):
    """
    Calculate time range for given Bartels rotation number.

    Returns a tuple of start and end date.
    """

    return br_to_date(br), br_to_date(br + 1)
