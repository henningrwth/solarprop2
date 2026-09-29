#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions to Figure 11 of

  R. D. Strauss et al., Astrophys. Space Sci. (2012) 339:223-236:

propagation times and energy losses for 100 MeV protons in the A<0 cycle.

"""

import argparse
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species
from sproptools.datareaders import read_xy

import solarprop
from solarprop.units import GeV, days

plt.style.use('solarprop.mplstyle')
plt.rcParams['font.size'] = 24

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--number', type=int, default=1000, help='Number of pseudo-particles per tilt angle.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use result files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

parameters = { 'model': 'Strauss2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'polarity': -1,
               'simpleHcsDrift': False,
               'forceNumericalDerivatives': False,
               'hcsFactor': 1.0,
              }

fname = 'test_strauss_2012_fig11_results.npz'

ekin = 0.1  # GeV

ref_angles, ref_proptime = read_xy('input/Strauss2012_Fig11_proptime.txt')
_, ref_deltaT = read_xy('input/Strauss2012_Fig11_deltaT.txt')
ref_deltaT /= 1000.  # MeV -> GeV


if args.reuse:
    print(f'Reading results from {fname}')
    npzfile = np.load(fname)
    tilt_angles= npzfile['tilt_angles']
    v_proptime= npzfile['v_proptime']
    v_err_proptime= npzfile['v_err_proptime']
    v_deltaT= npzfile['v_deltaT']
    v_err_deltaT= npzfile['v_err_deltaT']
else:
    tilt_angles = np.arange(10., 91., 10.)
    v_deltaT = np.zeros_like(tilt_angles)
    v_err_deltaT = np.zeros_like(tilt_angles)
    v_proptime = np.zeros_like(tilt_angles)
    v_err_proptime = np.zeros_like(tilt_angles)

    s = solarprop.Solarprop(parameters)

    for iangle, tilt_angle in enumerate(tilt_angles):

        s.update_parameters({'angle': tilt_angle})

        endpoints = s.endpoints(args.number, ekin)
        if not endpoints:
            sys.exit(1)

        deltaT = np.array([t.T() for t in endpoints if not t.isInsideHeliosphere()]) - ekin*GeV
        proptime = np.array([t.time() for t in endpoints if not t.isInsideHeliosphere()])

        mean_deltaT = np.mean(deltaT)
        med_deltaT = np.median(deltaT)
        err_deltaT = np.std(deltaT, ddof=1) / np.sqrt(len(deltaT))
        mean_proptime = np.mean(proptime)
        med_proptime = np.median(proptime)
        err_proptime = np.std(proptime, ddof=1) / np.sqrt(len(proptime))
        print(f'alpha= {tilt_angle} deg: mean delta T= {mean_deltaT/GeV:.3f} GeV, median= {med_deltaT/GeV:.3f} GeV, mean propagation time= {mean_proptime/days:.2f} days, median= {med_proptime/days:.2f} days')
        v_deltaT[iangle] = mean_deltaT
        v_err_deltaT[iangle] = err_deltaT
        v_proptime[iangle] = mean_proptime
        v_err_proptime[iangle] = err_proptime

    # save results
    np.savez(fname, tilt_angles=tilt_angles, v_proptime=v_proptime, v_err_proptime=v_err_proptime, v_deltaT=v_deltaT, v_err_deltaT=v_err_deltaT)
    print(f'Stored results in {fname}')

fig, ax = plt.subplots(1, 2, figsize=(18, 12))
ax[0].plot(ref_angles, ref_proptime, 'o', fillstyle='full', markerfacecolor='white', markeredgecolor='blue', markersize=12, label='Reference')
ax[0].errorbar(tilt_angles, v_proptime/days, yerr=v_err_proptime/days, fmt='s', color='blue', label='Solarprop 2.0')
ax[0].set_xlabel('HCS tilt angle')
ax[0].xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}°"))
ax[0].set_ylabel('Propagation time (days)')
ax[0].set_xlim(0., 99.)
ax[0].set_ylim(80., 180.)
ax[0].minorticks_on()
ax[0].grid()
ax[0].legend(loc='upper left')

ax[1].plot(ref_angles, ref_deltaT, 'o', fillstyle='full', markerfacecolor='white', markeredgecolor='red', markersize=12, label='Reference')
ax[1].errorbar(tilt_angles, v_deltaT/GeV, yerr=v_err_deltaT/GeV, fmt='s', color='red', label='Solarprop 2.0')
ax[1].set_xlabel('HCS tilt angle')
ax[1].xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:g}°"))
ax[1].yaxis.set_major_locator(mticker.MaxNLocator(nbins='auto', steps=[1, 2, 5, 10]))
ax[1].set_ylabel('Energy loss (GeV)')
ax[1].set_xlim(0., 99.)
ax[1].set_ylim(0.35, 2.5)
ax[1].minorticks_on()
ax[1].grid()
ax[1].legend(loc='upper left')
plt.tight_layout()

if args.batch:
    plt.savefig('test_strauss_2012_fig11.pdf')
else:
    plt.show()
