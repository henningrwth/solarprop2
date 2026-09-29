#! /usr/bin/env python3

"""
Plot the corrected version of Fig. 2 in Bobik et al. (2012).

The plot depicts the magnetic field amplitudes of the modified Parker spiral field
used in Helmod at different heliospheric radii vs colatitude.

The curves for 5 AU and 10 AU in the paper are wrong though. This script shows the
correct curves.

"""

import os
import importlib.resources as resources

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species

from solarprop.units import AU, deg, nanotesla
import solarprop

plt.style.use('solarprop.mplstyle')

plt.rcParams['xaxis.labellocation'] = 'center'
plt.rcParams['yaxis.labellocation'] = 'center'
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['lines.linewidth'] = 2.5

proton = Species.from_name('proton')

parameters = { 'model': 'Helmod2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'angleFile': str(resources.files(solarprop).joinpath('solarprop-data/angle.dat')),
               'omniFile': str(resources.files(solarprop).joinpath('solarprop-data/omniweb.dat')),
               'ssnFile': str(resources.files(solarprop).joinpath('solarprop-data/ssn.dat')),
               'kappaScaling': 1.0,
               'BfieldScaling': 1.0,
               'year': 1998,
               'month': 6,
               'day': 7,
               'tiltModel': 'L',
               'vswMin': 400.,
               'Bearth': 5.0,
               'highSolarActivityThreshold': 1000.,
               'angle': 0.,
               'forceNumericalDerivatives': False,
              }

distances = np.array([1., 5., 10.])
colatitudes = np.linspace(0., 180., 360)

s = solarprop.Solarprop(parameters)
B_mod = [np.array([s.magnetic_field(r*AU, theta*deg, 0.) / nanotesla for theta in colatitudes]) for r in distances]

s.update_parameters({'delta0': 0.})
B_parker = [np.array([s.magnetic_field(r*AU, theta*deg, 0.) / nanotesla for theta in colatitudes]) for r in distances]

plt.figure()
plt.plot(colatitudes, B_parker[0], color='blue', ls=':')
plt.plot(colatitudes, B_mod[0], color='blue')
plt.plot(colatitudes, B_parker[1], color='red', ls=':')
plt.plot(colatitudes, B_mod[1], color='red')
plt.plot(colatitudes, B_parker[2], color='green', ls=':')
plt.plot(colatitudes, B_mod[2], color='green')
plt.annotate(f'{distances[0]:g} AU', [90., np.max(B_mod[0])], xytext=(0., 5.), textcoords='offset points', color='blue', ha='center')
plt.annotate(f'{distances[1]:g} AU', [90., np.max(B_mod[1])], xytext=(0., 5.), textcoords='offset points', color='red', ha='center')
plt.annotate(f'{distances[2]:g} AU', [90., np.max(B_mod[2])], xytext=(0., -15.), textcoords='offset points', color='green', ha='center', va='top')

plt.xlabel('Colatitude')
plt.gca().xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}°"))
plt.ylabel('$|B|$ (nT)')
plt.xlim(0., 180.)
plt.ylim(1.e-3, 1.e1)
plt.yscale('log')
plt.minorticks_on()
plt.grid(which='both')
plt.grid(which='minor', color='lightgrey', linestyle=':', linewidth=0.5)

plt.tight_layout()
plt.show()
