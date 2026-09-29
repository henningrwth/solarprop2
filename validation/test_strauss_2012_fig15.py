#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  R. D. Strauss et al., Astrophys. Space Sci. (2012) 339:223-236

and compare results to Figure 15 in that publication. 

"""

import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy

import solarprop

plt.style.use('solarprop.mplstyle')
plt.rcParams['font.size'] = 24

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

parameters = { 'model': 'Strauss2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'polarity': -1,
               'particlesPerBin': args.particles,
               'forceNumericalDerivatives': False,
               'kappaScaling': 0.66892,
               'heliosphereBoundary': 10.,
               'rho': 0.05,
               'kappaExponent': 0.5,
               'hcsFactor': 1.0,
              }

parameters_burger = { 'model': 'Burger2012',
                      'mass': proton.mass,
                      'charge': proton.charge,
                      'polarity': -1,
                      'particlesPerBin': args.particles,
                     }

rw_opt = 'green' if args.reuse else 'save'
p0deg_full = { 'angle': 0, 'simpleHcsDrift': False, rw_opt: 'strauss2012_fig15_0deg_full.green'}
p30deg_full = { 'angle': 30,  'simpleHcsDrift': False, rw_opt: 'strauss2012_fig15_30deg_full.green'}
p0deg_simple = { 'angle': 0, 'simpleHcsDrift': True, rw_opt: 'strauss2012_fig15_0deg_simple.green'}
p30deg_simple = { 'angle': 30,  'simpleHcsDrift': True, rw_opt: 'strauss2012_fig15_30deg_simple.green'}
p0deg_burger = { 'angle': 0, rw_opt: 'strauss2012_fig15_0deg_burger2012.green'}
p30deg_burger = { 'angle': 30, rw_opt: 'strauss2012_fig15_30deg_burger2012.green'}

s = solarprop.Solarprop(parameters)

s.update_parameters(parameters | p0deg_full)
flux_0deg_full = s.modulate(ekin, v_lis_T)
s.update_parameters(parameters | p30deg_full)
flux_30deg_full = s.modulate(ekin, v_lis_T)

s.update_parameters(parameters | p0deg_simple)
flux_0deg_simple = s.modulate(ekin, v_lis_T)
s.update_parameters(parameters | p30deg_simple)
flux_30deg_simple = s.modulate(ekin, v_lis_T)

s_burger = solarprop.Solarprop(parameters_burger)
s_burger.update_parameters(p0deg_burger)
flux_0deg_burger = s_burger.modulate(ekin, v_lis_T)
s_burger.update_parameters(p30deg_burger)
flux_30deg_burger = s_burger.modulate(ekin, v_lis_T)


# data points digitized from Fig. 15 of the Strauss+ (2012) paper
str_0deg_ekin, str_0deg_flux = read_xy('input/Strauss2012_Fig15_0deg.txt')
str_30deg_ekin, str_30deg_flux = read_xy('input/Strauss2012_Fig15_30deg.txt')
# data points digitized from Fig. 2 of Burger (2012)
burger_0deg_ekin, burger_0deg_flux = read_xy('input/Burger2012_Fig2_0deg.txt')
burger_30deg_ekin, burger_30deg_flux = read_xy('input/Burger2012_Fig2_30deg.txt')
# data points digitized from Fig. 5 of Pei+ (2012)
pei_lis_ekin, pei_lis_flux = read_xy('input/Pei2012_Fig5_lis.txt')
pei_0deg_ekin, pei_0deg_flux = read_xy('input/Pei2012_Fig5_0deg.txt')
pei_30deg_ekin, pei_30deg_flux = read_xy('input/Pei2012_Fig5_30deg.txt')


plt.figure(figsize=(9,12))
plt.plot(ekin, v_lis_T, 'k-', label='LIS')
#plt.plot(pei_lis_ekin, pei_lis_flux, 'r:', label='LIS (Pei)')
plt.plot(str_0deg_ekin, str_0deg_flux, 'ks', label=r'Strauss (2012), $\mathregular{\alpha}$ = 0°')
plt.plot(str_30deg_ekin, str_30deg_flux, 'ko', markerfacecolor='white', label=r'Strauss (2012), $\mathregular{\alpha}$ = 30°')
plt.plot(burger_0deg_ekin, burger_0deg_flux, 'go', label=r'Burger "New" (2012), $\mathregular{\alpha}$ = 0°')
plt.plot(burger_30deg_ekin, burger_30deg_flux, 'rs', label=r'Burger "New" (2012), $\mathregular{\alpha}$ = 30°')
plt.plot(pei_0deg_ekin, pei_0deg_flux, 'gs', label=r'Pei (2012), $\mathregular{\alpha}$ = 0°')
plt.plot(pei_30deg_ekin, pei_30deg_flux, 'ro', label=r'Pei (2012), $\mathregular{\alpha}$ = 30°')

plt.plot(ekin, flux_0deg_full, 'g-', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 0°')
plt.plot(ekin, flux_30deg_full, 'r-', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 30°')
plt.plot(ekin, flux_0deg_simple, 'g:', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 0° (simple)')
plt.plot(ekin, flux_30deg_simple, 'r:', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 30° (simple)')
plt.plot(ekin, flux_0deg_burger, c='palegreen', ls='--', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 0° (Burger2012)')
plt.plot(ekin, flux_30deg_burger, c='lightcoral', ls='--', label=r'Solarprop 2.0: $\mathregular{\alpha}$ = 30° (Burger2012)')

plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.01, 10.)
plt.ylim(0.01, 10.)
plt.legend(fontsize=15)
mylocator = mticker.LogLocator(subs=[1., 3., ])
myformatter = mticker.StrMethodFormatter("{x:g}")
#plt.gca().xaxis.set_major_locator(mylocator)
#plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_strauss_2012_fig15.pdf')
else:
    plt.show()
