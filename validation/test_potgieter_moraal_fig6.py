#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  M.S. Potgieter and H. Moraal, ApJ 294 (1985) 425-440, https://adsabs.harvard.edu/abs/1985ApJ...294..425P

and compare results to their Figure 6.

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
argparser.add_argument('-n', '--particles', type=int, default=1000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

sf = 12.  # normalization from comparison to original paper

def pm_lis(T):
    beta = proton.convert(T, KinematicQuantity.KineticEnergy, KinematicQuantity.Beta)
    m = proton.mass
    return beta * np.pow(T + m, -2.6) * sf

ekin = np.geomspace(0.01, 25.0, num=30)
lis_T = pm_lis(ekin)

parameters = { 'model': 'PM1985',
               'mass': proton.mass,
               'charge': proton.charge,
               'particlesPerBin': args.particles,
              }

# we either save the Green function matrix file or re-use it, depending on whether the -r flag is set
rw_opt = 'green' if args.reuse else 'save'

pneg = { 'polarity': -1, rw_opt: 'pm_neg.green' }
ppos = { 'polarity': 1,  rw_opt: 'pm_pos.green' }

s = solarprop.Solarprop(parameters | pneg)
flux_mod_neg = s.modulate(ekin, lis_T)

s.update_parameters(ppos)
flux_mod_pos = s.modulate(ekin, lis_T)

# unfortunately, Figure 6 in the scan of the original paper is heavily distorted above ~1 GeV,
# so digitizing values of Potgieter and Moraal is difficult
pm_flux_ekin, pm_flux_vals = read_xy('input/PotgieterMoraal1985_Fig6_flux.txt')
pm_neg_ekin, pm_neg_flux = read_xy('input/PotgieterMoraal1985_Fig6_neg.txt')
pm_pos_ekin, pm_pos_flux = read_xy('input/PotgieterMoraal1985_Fig6_pos.txt')

vanilla_lis_ekin, vanilla_lis_flux = read_xy('input/vanilla_ref3_LIS.dat')
vanilla_lis_flux = sf * vanilla_lis_flux
vanilla_neg_ekin, vanilla_neg_flux = read_xy('input/vanilla_ref3_neg.dat')
vanilla_neg_flux = sf * vanilla_neg_flux
vanilla_pos_ekin, vanilla_pos_flux = read_xy('input/vanilla_ref3_pos.dat')
vanilla_pos_flux = sf * vanilla_pos_flux

plt.figure(figsize=(20,10))
# comparison of the digitized and analytic LIS can be used to gauge effects of the distortion of the paper scan
#plt.plot(pm_flux_ekin, pm_flux_vals, 'k-', label='LIS (paper)')
plt.plot(ekin, lis_T, 'b--', label='LIS')
plt.plot(pm_neg_ekin, pm_neg_flux, 'go', label='Reference: D(-) (A<0)')
plt.plot(pm_pos_ekin, pm_pos_flux, 'ro', label='Reference: D(+) (A>0)')

#plt.plot(vanilla_lis_ekin, vanilla_lis_flux, 'b:', label='vanilla LIS')
plt.plot(vanilla_neg_ekin, vanilla_neg_flux, 'g:', label='Solarprop: A<0')
plt.plot(vanilla_pos_ekin, vanilla_pos_flux, 'r:', label='Solarprop: A>0')

plt.plot(ekin, flux_mod_neg, 'g-', label='Solarprop 2.0: A<0')
plt.plot(ekin, flux_mod_pos, 'r-', label='Solarprop 2.0: A>0')
plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.01, 10.5)
plt.ylim(0.01, 10.)
plt.legend()
mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")
plt.gca().xaxis.set_major_locator(mylocator)
plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_potgieter_moraal_fig6.pdf')
else:
    plt.show()
