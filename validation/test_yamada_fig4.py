#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  Y. Yamada et al., Geophys. Res. Lett. Vol. 25 No. 13 (1998) 2353-2356, https://doi.org/10.1029/98GL51869

and compare results to their Figure 4.

"""

import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy

import solarprop
from solarprop.units import second, AU, MV, GV, cm2, km

plt.style.use('solarprop.mplstyle')

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--particles', type=int, default=1000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

def yamada_lis(p, m):
    return np.pow(m**2+p**2, -1.85) * p / 4. / np.pi  # flux  = p^2*f/4pi

ekin = np.geomspace(0.01, 50.0, num=32)  # Solarprop: 30 bins from 0.01 to 1278 (0.01*1.5**29) GeV
# apparently, in the Yamada paper, the transformation from
# differential flux in momentum to differential flux in kinetic energy
# is done ignoring the flux derivative ...
lis_T = yamada_lis(proton.convert(ekin, KinematicQuantity.KineticEnergy, KinematicQuantity.Momentum), proton.mass)

# we either save the Green function matrix file or re-use it, depending on whether the -r flag is set
rw_opt = 'green' if args.reuse else 'save'

parameters = { 'model': 'Yamada1998',
               'mass': proton.mass,
               'charge': proton.charge,
               'kappa0': 5.e22,
               'particlesPerBin': args.particles,
               rw_opt: 'yamada98.green',
              }

s = solarprop.Solarprop(parameters)
flux_mod = s.modulate(ekin, lis_T)

# phi = V(R-1)/3kappa
phi = 400.*km/second * (100.-1.)*AU / 3. / (parameters['kappa0'] * cm2/second/GV)
print(f'FF modulation potential phi = {phi / MV:.0f} MV')

def p2(T, m):
    return T**2 + 2*m*T

def forcefield(T, flux_lis, phi, m):
    T0 = T + phi
    flux0 = np.exp(np.interp(np.log(T0), np.log(T), np.log(flux_lis)))
    flux_mod = p2(T, m) / p2(T0, m) * flux0
    return T, flux_mod

ekin_ff, flux_ff = forcefield(ekin, lis_T, phi / GV, proton.mass)

yamada_flux_ekin, yamada_flux_vals = read_xy('input/Yamada1998_Fig4_flux.txt')
yamada_mod_ekin, yamada_mod_flux = read_xy('input/Yamada1998_Fig4_5e22.txt')

plt.figure(figsize=(20,10))
#plt.plot(yamada_flux_ekin, yamada_flux_vals, 'k-', label='LIS (paper)')  # digitized LIS
plt.plot(ekin, lis_T, 'b--', label='LIS')
plt.plot(yamada_mod_ekin, yamada_mod_flux, 'ko', label='Reference')
plt.plot(ekin_ff, flux_ff, 'g:', label='Force-field')
plt.plot(ekin, flux_mod, 'r-', label='Solarprop 2.0')
plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.0095, 40.)
plt.legend()
mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")
plt.gca().xaxis.set_major_locator(mylocator)
plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_yamada_fig4.pdf')
else:
    plt.show()
