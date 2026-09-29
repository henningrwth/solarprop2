#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  R.A. Burger, ApJ 760:60 (2012)

and compare results to Figure 2 in the reference publication.

"""

import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species
from sproptools.datareaders import read_xy

import solarprop

plt.style.use('solarprop.mplstyle')

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--particles', type=int, default=10000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

def lis(T):
    m = proton.mass
    return np.sqrt(T**2 + 2*T*m) * np.pow(T + m, -3.6) * 10.

ekin = np.geomspace(0.01, 25.0, num=30)
v_lis_T = lis(ekin)

parameters = { 'model': 'Burger2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'polarity': -1,
               'particlesPerBin': args.particles,
              }

# we either save the Green function matrix file or re-use it, depending on whether the -r flag is set
rw_opt = 'green' if args.reuse else 'save'

p0deg = { 'angle': 0, rw_opt: 'burger2012_0deg.green' }
p30deg = { 'angle': 30, rw_opt: 'burger2012_30deg.green' }

s = solarprop.Solarprop(parameters | p0deg)
flux_mod_0deg = s.modulate(ekin, v_lis_T)
s.update_parameters(p30deg)
flux_mod_30deg = s.modulate(ekin, v_lis_T)

paper_flux_ekin, paper_flux_vals = read_xy('input/Burger2012_Fig2_flux.txt')
paper_0deg_ekin, paper_0deg_flux = read_xy('input/Burger2012_Fig2_0deg.txt')
paper_30deg_ekin, paper_30deg_flux = read_xy('input/Burger2012_Fig2_30deg.txt')

plt.figure(figsize=(20,10))
#plt.plot(paper_flux_ekin, paper_flux_vals, 'k-', label='Flux (paper)')
plt.plot(paper_0deg_ekin, paper_0deg_flux, 'go', label=r'Reference: $\mathregular{\alpha}$ = 0°')
plt.plot(paper_30deg_ekin, paper_30deg_flux, 'ro', label=r'Reference: $\mathregular{\alpha}$ = 30°')
plt.plot(ekin, v_lis_T, 'b--', label='LIS')

plt.plot(ekin, flux_mod_0deg, 'g-', label='Solarprop 2.0: 0°')
plt.plot(ekin, flux_mod_30deg, 'r-', label='Solarprop 2.0: 30°')
plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.01, 10.)
plt.ylim(0.01, 10.)
plt.legend()
mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")
plt.gca().xaxis.set_major_locator(mylocator)
plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_burger_2012_fig2.pdf')
else:
    plt.show()
