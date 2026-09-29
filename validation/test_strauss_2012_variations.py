#! /usr/bin/env python3

"""
Test consistency of results under varying settings in the Strauss2012 model.

"""

import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species
from sproptools.datareaders import read_xy

import solarprop
from solarprop.units import days, second

plt.style.use('solarprop.mplstyle')
plt.rcParams['font.size'] = 24

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument('-a', '--angle', type=float, default=30., help='Tilt angle (in degrees).')
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--particles', type=int, default=5000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files from previous run.')
args = argparser.parse_args()

proton = Species.from_name('proton')

def langner_lis(T):
    # LIS parameterized in Langner and Potgieter, Journal Of Geophysical Research 109 (2004) A01103, doi:10.1029/2003JA010158
    t1 = np.exp(4.64 - 0.08*np.log(T)**2 - 2.91*np.sqrt(T))
    t2 = np.exp(3.22 - 2.86*np.log(T) - 1.50 / T)
    return t1 * (T<1.0) + t2 * (T>=1.0)

ekin = np.geomspace(0.001, 100.0, num=35)
lis_langner_T = langner_lis(ekin)

heliosphere_boundary = 30.  # AU

par = {
    'mass': proton.mass,
    'charge': proton.charge,
    'particlesPerBin': args.particles,
    'polarity': -1,
    'angle': args.angle,
    'heliosphereBoundary': heliosphere_boundary,
    'rho': 0.02,
    'kappaExponent': 1.0,
}

par_strauss = par | \
{ 'model': 'Strauss2012',
  'forceNumericalDerivatives': False,
  'hcsFactor': 1.0,
  'dt': 0.004*days/second,
  'simpleHcsDrift': False,
  'force2D': False,
  'kappaScaling': 1.0,
 }

par_burger = par | \
{ 'model': 'Burger2012',
  'kappaScaling':  1.0/0.66892,
}

rw_opt = 'green' if args.reuse else 'save'

s = solarprop.Solarprop(par_strauss)

s.update_parameters({rw_opt: 'var_strauss_std.green'})
flux_std = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_strauss | {'forceNumericalDerivatives': True, rw_opt: 'var_strauss_numdiff.green'})
flux_numdiff = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_strauss | {'simpleHcsDrift': True, rw_opt: 'var_strauss_simple.green'})
flux_simple = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_strauss | {'dt': 1200., rw_opt: 'var_strauss_dt1200.green'})
flux_dt_1200 = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_strauss | {'force2D': True, rw_opt: 'var_strauss_2d.green'})
flux_2d = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_burger | {rw_opt: 'var_strauss_burger.green'})
flux_burger = s.modulate(ekin, lis_langner_T)

# data points digitized from Fig. 5 of the Strauss (2012) paper
str_flux_ekin, str_flux_vals = read_xy('input/Strauss2012_Fig5_LIS.txt')
str_neg_10deg_ekin, str_neg_10deg_flux = read_xy('input/Strauss2012_Fig5_Aneg_10deg.txt')
str_neg_75deg_ekin, str_neg_75deg_flux = read_xy('input/Strauss2012_Fig5_Aneg_75deg.txt')


plt.figure(figsize=(9,12))
plt.plot(ekin, lis_langner_T, 'k-', label='LIS')
plt.plot(str_neg_10deg_ekin, str_neg_10deg_flux, 'D', color='darkgrey', label=r'Strauss (2012), $\mathregular{\alpha}$ = 10°')
plt.plot(str_neg_75deg_ekin, str_neg_75deg_flux, '^', color='darkgrey', label=r'Strauss (2012), $\mathregular{\alpha}$ = 75°')

plt.plot(ekin, flux_std, 'r-', label=rf'standard, $\mathregular{{\alpha}}$ = {args.angle:g}°')
plt.plot(ekin, flux_numdiff, 'g:', label=r'num. diff.')
plt.plot(ekin, flux_simple, 'b--', label=r'simple HCS drift')
plt.plot(ekin, flux_dt_1200, 'k', ls=(0, (1, 1)), label=r'dt = 1200 s')
plt.plot(ekin, flux_2d, 'm-.', label=r'2D')
plt.plot(ekin, flux_burger, 'c', ls=(0, (3, 1, 1, 1)), label=r'Burger2012 model')

plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.001, 100.)
plt.ylim(1.e-5, 1.e2)
plt.legend(fontsize='small')
mylocator = mticker.LogLocator(subs=[1., 3., ])
myformatter = mticker.StrMethodFormatter("{x:g}")
#plt.gca().xaxis.set_major_locator(mylocator)
#plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()

if args.batch:
    plt.savefig('test_strauss_2012_variations.pdf')
else:
    plt.show()
