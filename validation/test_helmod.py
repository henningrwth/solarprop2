#! /usr/bin/env -S python3 -u

"""
Calculate SOLARPROP predictions for the models

  - Helmod2012, aiming at an implementation of the Helmod model published in
    Bobik et al. ApJ 745 (2012) 132, https://doi.org/10.1088/0004-637X/745/2/132 and
  - Helmod2018, aiming at an implementation of the Helmod model published in
    Boschini et al., Adv. Space Res. 62 (2018) 2859-2879, https://doi.org/10.1016/j.asr.2017.04.017

With default arguments, this program aims at reproducing Figures 5, 6, 8, and 10 of
the Bobik et al. (2012) paper.

"""

import argparse
from copy import copy
from importlib import resources
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sunpy.coordinates.sun import carrington_rotation_number as cr

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy, read_crdb_file

import solarprop

plt.style.use('solarprop.mplstyle')
plt.rcParams['figure.figsize'] = (18, 10)
plt.rcParams['savefig.dpi'] = 140

mylocator = mticker.LogLocator(subs=[1., 2., 3., 5.])
myformatter = mticker.StrMethodFormatter("{x:g}")

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Do not display plots in the end.')
argparser.add_argument('-m', '--model', default='Helmod2012', choices=['Helmod2012', 'Helmod2018'], help='Model name of model to test.')
argparser.add_argument('-n', '--particles', type=int, default=1000, help='Number of pseudo-particles per bin.')
argparser.add_argument('-r', '--reuse', action='store_true', help='Re-use Green matrix files stored during previous run. Only use if nothing changed in the settings or model.')
argparser.add_argument('-q', '--quick', action='store_true', help='Skip additional calculations like numerical derivatives, varying timescales, etc?')
argparser.add_argument('-y', '--with-yamada', action='store_true', help='Also plot predictions by Yamada model?')
args = argparser.parse_args()

# choose correct option for saving/loading Green matrix files
rw_opt = 'green' if args.reuse else 'save'
# include predictions by Yamada model?
yamada = args.with_yamada

proton = Species.from_name('proton')

parameters = { 'model': args.model,
               'mass': proton.mass,
               'charge': proton.charge,
               'angleFile': str(resources.files(solarprop).joinpath('solarprop-data/angle.dat')),
               'omniFile': str(resources.files(solarprop).joinpath('solarprop-data/omniweb.dat')),
               'ssnFile': str(resources.files(solarprop).joinpath('solarprop-data/ssn.dat')),
               'kappaScaling': 0.65,
               'BfieldScaling': 1.0,
               'modelSpecificTimestep': False,
               'highSolarActivityThreshold': 150.,
               'particlesPerBin': args.particles,
              }

par_yamada = { 'model': 'Yamada1998',
               'mass': proton.mass,
               'charge': proton.charge,
               'nmFile': str(resources.files(solarprop).joinpath('solarprop-data/nm.dat')),
               'kappaScaling': 1.0,
               'particlesPerBin': args.particles,
              }

def lis(T, J0, gamma):
    R = proton.convert(T, KinematicQuantity.KineticEnergy, KinematicQuantity.Rigidity)

    flux_highR = J0*R**(-gamma)
    lnR = np.log(R)
    x = 9.472 - 1.999*lnR - 0.6938*lnR**2 + 0.2988*lnR**3 - 0.04714*lnR**4
    flux_lowR  = J0/(1.9e4) * np.exp(x)

    highR = (R>=7.)
    lowR = ~highR
    return lowR*flux_lowR + highR*flux_highR

ekin = np.geomspace(0.2, 200.0, num=23)

# values for LIS normalization and spectral index extracted from Fig. 3 of Bobik+ (2012)
#

# high solar activity:
#
# section 7.2 of Bobik+ (2012): using K_perp_theta=K_perp_r
# and solar latitudes lower than 30° (which seems to refer
# to the forward-in-time integration method and is not applicable here)

# Fig. 5: comparison to BESS-2000 data
print ('  ====================  BESS 2000 ====================  ')
date_bess00 = {'year': 2000,
               'month': 8,
               'day': 10,
               }

bess00 = parameters | date_bess00 | { 'tiltModel': 'L',
                                      'iotaPolarRegion': 1.0,
                                      rw_opt: f'{args.model}_bess00.green',
                                     }

lis_T_bess00 = lis(ekin, J0=1.74e4, gamma=2.78)
s_bess00 = solarprop.Solarprop(bess00)
flux_mod_bess00 = s_bess00.modulate(ekin, lis_T_bess00)

# Simplified HCS drift leads to artifacts at higher energies (need to debug this...),
# so we restrict the energy range for this case.
ekin_simplehcs = ekin[ekin<50.]

if not args.quick:
    s_bess00_nd = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_nd.green', 'forceNumericalDerivatives': True})
    flux_mod_bess00_nd = s_bess00_nd.modulate(ekin, lis_T_bess00)

    #s_bess00_dt500 = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_dt500.green', 'dt': 500.})
    #flux_mod_bess00_dt500 = s_bess00_dt500.modulate(ekin, lis_T_bess00)

    s_bess00_dt2k = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_dt2k.green', 'dt': 2000.})
    flux_mod_bess00_dt2k = s_bess00_dt2k.modulate(ekin, lis_T_bess00)

    #s_bess00_dt3k = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_dt3k.green', 'dt': 3000.})
    #flux_mod_bess00_dt3k = s_bess00_dt3k.modulate(ekin, lis_T_bess00)

    s_bess00_dt = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_dt.green', 'modelSpecificTimestep': True})
    flux_mod_bess00_dt = s_bess00_dt.modulate(ekin, lis_T_bess00)

    s_bess00_simplehcs = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_simplehcs.green', 'simpleHcsDrift': True})
    flux_mod_bess00_simplehcs = s_bess00_simplehcs.modulate(ekin_simplehcs, lis_T_bess00)

    s_bess00_kappa1 = solarprop.Solarprop(bess00 | {rw_opt: f'{args.model}_bess00_kappa1.green', 'kappaScaling': 1.0})
    flux_mod_bess00_kappa1 = s_bess00_kappa1.modulate(ekin, lis_T_bess00)

if yamada:
    s_bess00_yam = solarprop.Solarprop(par_yamada | date_bess00 | { rw_opt: 'yamada_bess00.green' })
    flux_mod_bess00_yam = s_bess00_yam.modulate(ekin, lis_T_bess00)


# Fig. 6: comparison to BESS-2002 data
print ('  ====================  BESS 2002 ====================  ')
date_bess02 = {'year': 2002,
               'month': 8,
               'day': 7,
               }
bess02 = parameters | date_bess02 | { 'tiltModel': 'L',
                                      'iotaPolarRegion': 1.0,
                                      rw_opt: f'{args.model}_bess02.green',
                                     }

lis_T_bess02 = lis(ekin, J0=1.74e4, gamma=2.75)
s_bess02 = solarprop.Solarprop(bess02)
flux_mod_bess02 = s_bess02.modulate(ekin, lis_T_bess02)

if not args.quick:
    s_bess02_nd = solarprop.Solarprop(bess02 | {rw_opt: f'{args.model}_bess02_nd.green', 'forceNumericalDerivatives': True})
    flux_mod_bess02_nd = s_bess02_nd.modulate(ekin, lis_T_bess02)

    s_bess02_dt2k = solarprop.Solarprop(bess02 | {rw_opt: f'{args.model}_bess02_dt2k.green', 'dt': 2000.})
    flux_mod_bess02_dt2k = s_bess02_dt2k.modulate(ekin, lis_T_bess02)

    s_bess02_dt = solarprop.Solarprop(bess02 | {rw_opt: f'{args.model}_bess02_dt.green', 'modelSpecificTimestep': True})
    flux_mod_bess02_dt = s_bess02_dt.modulate(ekin, lis_T_bess02)

    s_bess02_simplehcs = solarprop.Solarprop(bess02 | {rw_opt: f'{args.model}_bess02_simplehcs.green', 'simpleHcsDrift': True})
    flux_mod_bess02_simplehcs = s_bess02_simplehcs.modulate(ekin_simplehcs, lis_T_bess02)

    s_bess02_kappa1 = solarprop.Solarprop(bess02 | {rw_opt: f'{args.model}_bess02_kappa1.green', 'kappaScaling': 1.0})
    flux_mod_bess02_kappa1 = s_bess02_kappa1.modulate(ekin, lis_T_bess02)

if yamada:
    s_bess02_yam = solarprop.Solarprop(par_yamada | date_bess02 | { rw_opt: 'yamada_bess02.green' })
    flux_mod_bess02_yam = s_bess02_yam.modulate(ekin, lis_T_bess02)


#
# not dominated by high solar activity:
#

# Fig. 8: comparison to AMS-1998 data
print ('  ====================  AMS 1998 ====================  ')
date_ams98 = {'year': 1998,
              'month': 6,
              'day': 7,
              }
ams98 = parameters | date_ams98 | { 'tiltModel': 'L',
                                    rw_opt: f'{args.model}_ams98.green',
                                   }
if args.model == 'Helmod2018':
    ams98 |= {'kappaScaling': 1.5}

lis_T_ams98 = lis(ekin, J0=1.74e4, gamma=2.78)
s_ams98 = solarprop.Solarprop(ams98)
flux_mod_ams98 = s_ams98.modulate(ekin, lis_T_ams98)

if not args.quick:
    s_ams98_nd = solarprop.Solarprop(ams98 | {rw_opt: f'{args.model}_ams98_nd.green', 'forceNumericalDerivatives': True})
    flux_mod_ams98_nd = s_ams98_nd.modulate(ekin, lis_T_ams98)

    s_ams98_dt2k = solarprop.Solarprop(ams98 | {rw_opt: f'{args.model}_ams98_dt2k.green', 'dt': 2000.})
    flux_mod_ams98_dt2k = s_ams98_dt2k.modulate(ekin, lis_T_ams98)

    s_ams98_dt = solarprop.Solarprop(ams98 | {rw_opt: f'{args.model}_ams98_dt.green', 'modelSpecificTimestep': True})
    flux_mod_ams98_dt = s_ams98_dt.modulate(ekin, lis_T_ams98)

    s_ams98_simplehcs = solarprop.Solarprop(ams98 | {rw_opt: f'{args.model}_ams98_simplehcs.green', 'simpleHcsDrift': True})
    flux_mod_ams98_simplehcs = s_ams98_simplehcs.modulate(ekin_simplehcs, lis_T_ams98)

    s_ams98_kappa1 = solarprop.Solarprop(ams98 | {rw_opt: f'{args.model}_ams98_kappa1.green', 'kappaScaling': 1.0})
    flux_mod_ams98_kappa1 = s_ams98_kappa1.modulate(ekin, lis_T_ams98)

if yamada:
    s_ams98_yam = solarprop.Solarprop(par_yamada | date_ams98 | { rw_opt: 'yamada_ams98.green' })
    flux_mod_ams98_yam = s_ams98_yam.modulate(ekin, lis_T_ams98)


# Fig. 10: comparison to PAMELA-2006/08 data
print ('  ====================  PAMELA 2006/08 ====================  ')
date_pamela =  {'year': 2007,
                'month': 6,
                'day': 1,
                }
pamela = parameters | date_pamela | { 'tiltModel': 'L',
                                      rw_opt: f'{args.model}_pamela.green',
                                     }

lis_T_pamela = lis(ekin, J0=1.77e4, gamma=2.79)
s_pamela = solarprop.Solarprop(pamela)
flux_mod_pamela = s_pamela.modulate(ekin, lis_T_pamela)

if not args.quick:
    s_pamela_nd = solarprop.Solarprop(pamela | {rw_opt: f'{args.model}_pamela_nd.green', 'forceNumericalDerivatives': True})
    flux_mod_pamela_nd = s_pamela_nd.modulate(ekin, lis_T_pamela)

    s_pamela_dt2k = solarprop.Solarprop(pamela | {rw_opt: f'{args.model}_pamela_dt2k.green', 'dt': 2000.})
    flux_mod_pamela_dt2k = s_pamela_dt2k.modulate(ekin, lis_T_pamela)

    s_pamela_dt = solarprop.Solarprop(pamela | {rw_opt: f'{args.model}_pamela_dt.green', 'modelSpecificTimestep': True})
    flux_mod_pamela_dt = s_pamela_dt.modulate(ekin, lis_T_pamela)

    s_pamela_simplehcs = solarprop.Solarprop(pamela | {rw_opt: f'{args.model}_pamela_simplehcs.green', 'simpleHcsDrift': True})
    flux_mod_pamela_simplehcs = s_pamela_simplehcs.modulate(ekin_simplehcs, lis_T_pamela)

    s_pamela_kappa1 = solarprop.Solarprop(pamela | {rw_opt: f'{args.model}_pamela_kappa1.green', 'kappaScaling': 1.0})
    flux_mod_pamela_kappa1 = s_pamela_kappa1.modulate(ekin, lis_T_pamela)

if yamada:
    s_pamela_yam = solarprop.Solarprop(par_yamada | date_pamela | { rw_opt: 'yamada_pamela.green' })
    flux_mod_pamela_yam = s_pamela_yam.modulate(ekin, lis_T_pamela)


datestring = f'{ams98["year"]}-{ams98["month"]}-{ams98["day"]}'
cr_number = cr(datestring)
print(f'CR for {datestring} = {cr_number:.2f}')

helmod_lis_ams98_ekin, helmod_lis_ams98_flux = read_xy('input/Bobik2012_Fig8_flux.txt')
helmod_mod_ams98_ekin, helmod_mod_ams98_flux = read_xy('input/Bobik2012_Fig8_modulated.txt')
helmod_lis_bess00_ekin, helmod_lis_bess00_flux = read_xy('input/Bobik2012_Fig5_flux.txt')
helmod_mod_bess00_ekin, helmod_mod_bess00_flux = read_xy('input/Bobik2012_Fig5_modulated.txt')

pams01 = read_crdb_file('input/CRDB/AMS01_ProtonFlux_2000.dat')
pbess00 = read_crdb_file('input/CRDB/BESS00_ProtonFlux_2007.dat')
pbess02 = read_crdb_file('input/CRDB/BESS02_ProtonFlux_2007.dat')
ppamela = read_crdb_file('input/CRDB/PAMELA_ProtonFlux_2011.dat')

xlabel = r'$\mathregular{E_{kin}}$ (GeV)'
ylabel = r'Flux (${\mathregular{m}^\mathregular{-2}}\,\mathregular{sr}^\mathregular{-1}\,\mathregular{s}^\mathregular{-1}\,\mathregular{GeV}^\mathregular{-1}$)'

xmin = 0.14
xmax = 280.

thinline = 2.5

fig, ax = plt.subplots(num='AMS-1998')
ax.plot(ekin, lis_T_ams98, 'k-', label='LIS')
ax.plot(ekin, flux_mod_ams98, 'r-', label='Solarprop 2.0')
#ax.plot(helmod_lis_ams98_ekin, helmod_lis_ams98_flux, 'g:', label='Reference LIS')
ax.plot(helmod_mod_ams98_ekin, helmod_mod_ams98_flux, 'c--', lw=thinline, label='Reference')
if not args.quick:
    ax.plot(ekin, flux_mod_ams98_nd, 'g:', label='Num.diff.')
    ax.plot(ekin, flux_mod_ams98_dt2k, ls=(0, (1, 3)), color='coral', label='dt = 2000 s')
    ax.plot(ekin, flux_mod_ams98_dt, ls=(0, (4, 6)), lw=thinline, color='darkred', label=r'$\mathregular{d}t=r^2/\kappa_{rr}$')
    ax.plot(ekin_simplehcs, flux_mod_ams98_simplehcs, '-.', color='darkviolet', label='simple drifts')
    ax.plot(ekin, flux_mod_ams98_kappa1, ls=(0, (3, 2, 1, 2, 1, 2)), lw=thinline, color='darkgrey', label=r'no $\kappa$ scaling')
if yamada:
    ax.plot(ekin, flux_mod_ams98_yam, '--', color='grey', label='Yamada model')
ax.errorbar(pams01['EKN_mean'], pams01['value'], yerr=pams01['err_total'], fmt='o', color='red', label='AMS-1998')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(1.e-3, 1.e5)
ax.legend()
ax.xaxis.set_major_locator(copy(mylocator))
ax.xaxis.set_major_formatter(copy(myformatter))
ax.minorticks_on()
plt.tight_layout()
if args.batch:
    plt.savefig(f'test_helmod_{args.model}_ams98.pdf')

fig, ax = plt.subplots(num='BESS-2000', figsize=(14, 10))
ax.plot(ekin, lis_T_bess00, 'k-', label='LIS')
ax.plot(ekin, flux_mod_bess00, 'r-', label='Solarprop 2.0')
#ax.plot(helmod_lis_bess00_ekin, helmod_lis_bess00_flux, 'b--', label='Reference LIS')
ax.plot(helmod_mod_bess00_ekin, helmod_mod_bess00_flux, 'c--', lw=thinline, label='Reference')
if not args.quick:
    ax.plot(ekin, flux_mod_bess00_nd, 'g:', label='Num.diff.')
    #ax.plot(ekin, flux_mod_bess00_dt500, ':', color='skyblue', label='dt = 500 s')
    ax.plot(ekin, flux_mod_bess00_dt2k, ls=(0, (1, 3)), color='coral', label='dt = 2000 s')
    #ax.plot(ekin, flux_mod_bess00_dt3k, ':', color='thistle', label='dt = 3000 s')
    ax.plot(ekin, flux_mod_bess00_dt, ls=(0, (4, 6)), lw=thinline, color='darkred', label=r'$\mathregular{d}t=r^2/\kappa_{rr}$')
    ax.plot(ekin_simplehcs, flux_mod_bess00_simplehcs, '-.', color='darkviolet', label='simple drifts')
    ax.plot(ekin, flux_mod_bess00_kappa1, ls=(0, (3, 2, 1, 2, 1, 2)), lw=thinline, color='darkgrey', label=r'no $\kappa$ scaling')
if yamada:
    ax.plot(ekin, flux_mod_bess00_yam, '--', color='grey', label='Yamada model')
ax.errorbar(pbess00['EKN_mean'], pbess00['value'], yerr=pbess00['err_total'], fmt='o', color='red', label='BESS-2000')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, 26.)
ax.set_ylim(1., 1.e5)
ax.legend()
ax.xaxis.set_major_locator(copy(mylocator))
ax.xaxis.set_major_formatter(copy(myformatter))
ax.minorticks_on()
plt.tight_layout()
if args.batch:
    plt.savefig(f'test_helmod_{args.model}_bess00.pdf')

fig, ax = plt.subplots(num='BESS-2002')
ax.plot(ekin, lis_T_bess02, 'k-', label='LIS')
ax.plot(ekin, flux_mod_bess02, 'r-', label='Solarprop 2.0')
if not args.quick:
    ax.plot(ekin, flux_mod_bess02_nd, 'g:', label='Num.diff.')
    ax.plot(ekin, flux_mod_bess02_dt2k, ls=(0, (1, 3)), color='coral', label='dt = 2000 s')
    ax.plot(ekin, flux_mod_bess02_dt, ls=(0, (4, 6)), lw=thinline, color='darkred', label=r'$\mathregular{d}t=r^2/\kappa_{rr}$')
    ax.plot(ekin_simplehcs, flux_mod_bess02_simplehcs, '-.', color='darkviolet', label='simple drifts')
    ax.plot(ekin, flux_mod_bess02_kappa1, ls=(0, (3, 2, 1, 2, 1, 2)), lw=thinline, color='darkgrey', label=r'no $\kappa$ scaling')
if yamada:
    ax.plot(ekin, flux_mod_bess02_yam, '--', color='grey', label='Yamada model')
ax.errorbar(pbess02['EKN_mean'], pbess02['value'], yerr=pbess02['err_total'], fmt='o', color='red', label='BESS-2002')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(1.e-3, 1.e5)
ax.legend()
ax.xaxis.set_major_locator(copy(mylocator))
ax.xaxis.set_major_formatter(copy(myformatter))
ax.minorticks_on()
plt.tight_layout()
if args.batch:
    plt.savefig(f'test_helmod_{args.model}_bess02.pdf')

fig, ax = plt.subplots(num='PAMELA-2006/08')
ax.plot(ekin, lis_T_pamela, 'k-', label='LIS')
ax.plot(ekin, flux_mod_pamela, 'r-', label='Solarprop 2.0')
if not args.quick:
    ax.plot(ekin, flux_mod_pamela_nd, 'g:', label='Num.diff.')
    ax.plot(ekin, flux_mod_pamela_dt2k, ls=(0, (1, 3)), color='coral', label='dt = 2000 s')
    ax.plot(ekin, flux_mod_pamela_dt, ls=(0, (4, 6)), lw=thinline, color='darkred', label=r'$\mathregular{d}t=r^2/\kappa_{rr}$')
    ax.plot(ekin_simplehcs, flux_mod_pamela_simplehcs, '-.', color='darkviolet', label='simple drifts')
    ax.plot(ekin, flux_mod_pamela_kappa1, ls=(0, (3, 2, 1, 2, 1, 2)), lw=thinline, color='darkgrey', label=r'no $\kappa$ scaling')
if yamada:
    ax.plot(ekin, flux_mod_pamela_yam, '--', color='grey', label='Yamada model')
ax.errorbar(ppamela['EKN_mean'], ppamela['value'], yerr=ppamela['err_total'], fmt='o', color='red', label='PAMELA-2006/08')
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(1.e-3, 1.e5)
ax.legend()
ax.xaxis.set_major_locator(copy(mylocator))
ax.xaxis.set_major_formatter(copy(myformatter))
ax.minorticks_on()
plt.tight_layout()
if args.batch:
    plt.savefig(f'test_helmod_{args.model}_pamela.pdf')

if not args.batch:
    plt.show()
