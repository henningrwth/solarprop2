#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  J.R. Jokipii and D.A. Kopriva, ApJ 234 (1979) 384-392, https://adsabs.harvard.edu/pdf/1979ApJ...234..384J

and compare results to their Figure 2.

"""

import argparse
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy

import solarprop

plt.style.use('solarprop.mplstyle')

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-d', '--num-diff', action='store_true', help='Include results from numerical differentiation.')
argparser.add_argument('-n', '--particles', type=int, default=1000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

sf = 8.0  # LIS normalization from comparison to original paper
def jk_lis(p, m):
   return np.pow(m**2+p**2, -1.8) * p * sf  # flux ~ p^2*f

ekin = np.geomspace(0.01442, 10.0, num=21)
# apparently, in the model paper, the transformation from
# differential flux in momentum to differential flux in kinetic energy
# is done ignoring the flux derivative ...
lis_T = jk_lis(proton.convert(ekin, KinematicQuantity.KineticEnergy, KinematicQuantity.Momentum), proton.mass)

parameters = { 'model': 'JK1979',
               'mass': proton.mass,
               'charge': proton.charge,
               'particlesPerBin': args.particles,
              }

negative = { 'polarity': -1 }
positive = { 'polarity': 1 }

# we either save the Green function matrix file or re-use it, depending on whether the -r flag is set
rw_opt = 'green' if args.reuse else 'save'

fileneg = { rw_opt: 'jk_neg.green' }
filepos = { rw_opt: 'jk_pos.green' }
fileneg2 = { rw_opt: 'jk_neg2.green' }
filepos2 = { rw_opt: 'jk_pos2.green' }

s = solarprop.Solarprop(parameters)

s.update_parameters(negative | fileneg)
flux_mod_neg = s.modulate(ekin, lis_T)

s.update_parameters(positive | filepos)
flux_mod_pos = s.modulate(ekin, lis_T)

if args.num_diff:
   s.update_parameters(negative | fileneg2 | {'forceNumericalDerivatives': True})
   flux_mod_neg2 = s.modulate(ekin, lis_T)

   s.update_parameters(positive | filepos2 | {'forceNumericalDerivatives': True})
   flux_mod_pos2 = s.modulate(ekin, lis_T)

jk_flux_ekin, jk_flux_vals = read_xy('input/JokipiiKopriva1979_Fig2_flux.txt')
jk_neg_ekin, jk_neg_flux = read_xy('input/JokipiiKopriva1979_Fig2_qAneg.txt')
jk_pos_ekin, jk_pos_flux = read_xy('input/JokipiiKopriva1979_Fig2_qApos.txt')

vanilla_lis_ekin, vanilla_lis_flux = read_xy('input/vanilla_ref2_LIS.dat')
vanilla_lis_flux = sf * vanilla_lis_flux
vanilla_neg_ekin, vanilla_neg_flux = read_xy('input/vanilla_ref2_neg.dat')
vanilla_neg_flux = sf * vanilla_neg_flux
vanilla_pos_ekin, vanilla_pos_flux = read_xy('input/vanilla_ref2_pos.dat')
vanilla_pos_flux = sf * vanilla_pos_flux

plt.figure(figsize=(20,10))

# digitized results from reference
#plt.plot(jk_flux_ekin, jk_flux_vals, 'k-', label='LIS (paper)')
plt.plot(jk_neg_ekin, jk_neg_flux, 'go', label='Reference: qA<0')
plt.plot(jk_pos_ekin, jk_pos_flux, 'ro', label='Reference: qA>0')
plt.plot(ekin, lis_T, 'b--', label='LIS')

# results obtained with vanilla Solarprop
#plt.plot(vanilla_lis_ekin, vanilla_lis_flux, 'b:', label='vanilla LIS')
plt.plot(vanilla_neg_ekin, vanilla_neg_flux, 'g:', label='Solarprop: qA<0')
plt.plot(vanilla_pos_ekin, vanilla_pos_flux, 'r:', label='Solarprop: qA>0')

# results from Solarprop 2.0
plt.plot(ekin, flux_mod_neg, 'g-', label='Solarprop 2.0: qA<0')
if args.num_diff:
   plt.plot(ekin, flux_mod_neg2, '-', color='c', label='Num.diff. qA<0')
plt.plot(ekin, flux_mod_pos, 'r-', label='Solarprop 2.0: qA>0')
if args.num_diff:
   plt.plot(ekin, flux_mod_pos2, '-', color='m', label='Num.diff. qA>0')
plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.0145, 10.)
plt.ylim(0.02, 6.)
plt.legend()
mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")
plt.gca().xaxis.set_major_locator(mylocator)
plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
   plt.savefig('test_jokipii_kopriva_fig2.pdf')
else:
   plt.show()
