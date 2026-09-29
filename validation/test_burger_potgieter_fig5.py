#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  R.A. Burger and M.S. Potgieter, ApJ 339 (1989) 501-511, https://adsabs.harvard.edu/full/1989apj...339..501b

and compare results. To be more precise, this model aims at
reproducing Fig. 5 in the Burger&Potgieter paper. The figure was done
with the same parameters as used in Kota&Jokipii, ApJ 265 (1983) 573,
but with the new implementation of HCS drift by Burger&Potgieter,
see also Burger&Hattingh, Astrophys. Space Sci. 230 (1995) 375-382.

"""

import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy

import solarprop

plt.style.use('solarprop.mplstyle')

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--particles', type=int, default=10000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

sf = 12.  # normalization from comparison to original paper

def bp_lis(T, m):
    return np.sqrt(T**2 + 2*T*m) * np.pow(T + m, -3.6) * sf

ekin = np.geomspace(0.01, 25.0, num=30)
v_lis_T = bp_lis(ekin, proton.mass)


parameters = { 'model': 'BP1989',
               'mass': proton.mass,
               'charge': proton.charge,
               'polarity': -1,
               'particlesPerBin': args.particles,
              }

# we either save the Green function matrix file or re-use it, depending on whether the -r flag is set
rw_opt = 'green' if args.reuse else 'save'

p0deg = { 'angle': 0, rw_opt: 'bp_0deg.green'}
p30deg = { 'angle': 30,  rw_opt: 'bp_30deg.green'}

s = solarprop.Solarprop(parameters | p0deg)
flux_mod_0deg = s.modulate(ekin, v_lis_T)

s.update_parameters(p30deg)
flux_mod_30deg = s.modulate(ekin, v_lis_T)

bp_flux_ekin, bp_flux_vals = read_xy('input/BurgerPotgieter1989_Fig5_flux.txt')
bp_0deg_ekin, bp_0deg_flux = read_xy('input/BurgerPotgieter1989_Fig5_0deg.txt')
bp_30deg_ekin, bp_30deg_flux = read_xy('input/BurgerPotgieter1989_Fig5_30deg.txt')

vanilla_lis_ekin, vanilla_lis_flux = read_xy('input/vanilla_ref4_LIS.dat')
vanilla_lis_flux = sf * vanilla_lis_flux
vanilla_0deg_ekin, vanilla_0deg_flux = read_xy('input/vanilla_ref4_0deg.dat')
vanilla_0deg_flux = sf * vanilla_0deg_flux
vanilla_30deg_ekin, vanilla_30deg_flux = read_xy('input/vanilla_ref4_30deg.dat')
vanilla_30deg_flux = sf * vanilla_30deg_flux

plt.figure(figsize=(20,10))
#plt.plot(bp_flux_ekin, bp_flux_vals, 'k-', label='Flux (paper)')
plt.plot(bp_0deg_ekin, bp_0deg_flux, 'go', label=r'Reference: $\mathregular{\alpha}$ = 0°')
plt.plot(bp_30deg_ekin, bp_30deg_flux, 'ro', label=r'Reference: $\mathregular{\alpha}$ = 30°')
plt.plot(ekin, v_lis_T, 'b--', label='LIS')

#plt.plot(vanilla_lis_ekin, vanilla_lis_flux, 'b:', label='vanilla LIS')
plt.plot(vanilla_0deg_ekin, vanilla_0deg_flux, 'g:', label='Solarprop: 0°')
plt.plot(vanilla_30deg_ekin, vanilla_30deg_flux, 'r:', label='Solarprop: 30°')

plt.plot(ekin, flux_mod_0deg, 'g-', label='Solarprop 2.0: 0°')
plt.plot(ekin, flux_mod_30deg, 'r-', label='Solarprop 2.0: 30°')
plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.01, 10.5)
plt.ylim(0.08, 6.)
plt.legend()
mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")
plt.gca().xaxis.set_major_locator(mylocator)
plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_burger_potgieter_fig5.pdf')
else:
    plt.show()
