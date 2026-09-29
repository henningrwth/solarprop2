#! /usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

class SolarpropInterpolator:
    """
    Handle interpolation of a dataframe containing
    the results of a time series of Solarprop runs.

    Logarithmic interpolation is used for fluxes.

    Example:

        interp = SolarpropInterpolator(sprop[sprop['Z'] == 1], 'Ekin', 'date', 'flux')
        ekin = np.geomspace(0.1, 100.0, num=50)
        interpolated_fluxes = interp(ekin, np.datetime64('2012-05-25'))

    """

    def __init__(self, dataframe, eaxis, taxis, vaxis):
        """
        Construct interpolator.

        dataframe: Dataframe of solarprop results for a single species.
        eaxis: Name of kinematic variable (e.g., kinetic energy or rigidity) in dataframe.
        taxis: Name of time axis (or dates) in dataframe.
        vaxis: Name of flux value axis in dataframe.
        """

        df = dataframe.copy()

        # prepare for log-interpolation
        log_eaxis = f'log_{eaxis}'
        log_vaxis = f'log_{vaxis}'
        df[log_eaxis] = np.log(df[eaxis])
        df[log_vaxis] = np.log(df[vaxis])

        pivottable = df.pivot(index=log_eaxis, columns=taxis, values=log_vaxis)
        self.interpolator = RegularGridInterpolator((pivottable.index, pivottable.columns), pivottable.to_numpy(), method='linear')
        self.datetype = pivottable.columns.dtype


    def __call__(self, x, t):
        """
        Evaluate flux(es) at given kinematic variable(s) x and date/time t.
        """
        try:
            v = self.interpolator((np.log(x), t.astype(self.datetype)))
        except ValueError as e:
            print(e)
            xgrid = np.exp(self.interpolator.grid[0])
            tgrid = self.interpolator.grid[1].astype(self.datetype)
            print(f'Grid range in x-direction from {xgrid[0]} to {xgrid[-1]}, called with x = {x}')
            print(f'Grid range in t-direction from {tgrid[0]} to {tgrid[-1]}, called with t = {t}')
            raise e
        return np.exp(v)
