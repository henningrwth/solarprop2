#! /usr/bin/env python3

"""
Calculate SOLARPROP predictions for the model published in

  R. D. Strauss et al., Astrophys. Space Sci. (2012) 339:223-236

and compare results to Figure 5 of this paper.

"""

import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from sproptools.kinematics import Species, KinematicQuantity
from sproptools.datareaders import read_xy

import solarprop

plt.style.use('solarprop.mplstyle')
plt.rcParams['font.size'] = 26

argparser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-b', '--batch', action='store_true', help='Save plot to file and do not draw.')
argparser.add_argument('-n', '--particles', type=int, default=1000, help='Number of pseudo-particles per bin.')
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


parameters = { 'model': 'Strauss2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'particlesPerBin': args.particles,
               'forceNumericalDerivatives': False,
               'hcsFactor': 1.0,
              }

rw_opt = 'green' if args.reuse else 'save'
par_pos_10deg_simple = { 'polarity': +1, 'angle': 10, 'simpleHcsDrift': True, rw_opt: 'strauss_pos_10deg_simple.green'}
par_pos_75deg_simple = { 'polarity': +1, 'angle': 75, 'simpleHcsDrift': True, rw_opt: 'strauss_pos_75deg_simple.green'}
par_neg_10deg_simple = { 'polarity': -1, 'angle': 10, 'simpleHcsDrift': True, rw_opt: 'strauss_neg_10deg_simple.green'}
par_neg_75deg_simple = { 'polarity': -1, 'angle': 75, 'simpleHcsDrift': True, rw_opt: 'strauss_neg_75deg_simple.green'}
par_pos_10deg_full = { 'polarity': +1, 'angle': 10, 'simpleHcsDrift': False, rw_opt: 'strauss_pos_10deg_full.green'}
par_pos_75deg_full = { 'polarity': +1, 'angle': 75, 'simpleHcsDrift': False, rw_opt: 'strauss_pos_75deg_full.green'}
par_neg_10deg_full = { 'polarity': -1, 'angle': 10, 'simpleHcsDrift': False, rw_opt: 'strauss_neg_10deg_full.green'}
par_neg_75deg_full = { 'polarity': -1, 'angle': 75, 'simpleHcsDrift': False, rw_opt: 'strauss_neg_75deg_full.green'}
par_pos_10deg_slow = { 'polarity': +1, 'angle': 10, 'simpleHcsDrift': False, 'hcsFactor': 0.5, rw_opt: 'strauss_pos_10deg_slow.green'}
par_pos_75deg_slow = { 'polarity': +1, 'angle': 75, 'simpleHcsDrift': False, 'hcsFactor': 0.5, rw_opt: 'strauss_pos_75deg_slow.green'}
par_neg_10deg_slow = { 'polarity': -1, 'angle': 10, 'simpleHcsDrift': False, 'hcsFactor': 0.5, rw_opt: 'strauss_neg_10deg_slow.green'}
par_neg_75deg_slow = { 'polarity': -1, 'angle': 75, 'simpleHcsDrift': False, 'hcsFactor': 0.5, rw_opt: 'strauss_neg_75deg_slow.green'}

s = solarprop.Solarprop(parameters)

s.update_parameters(par_pos_10deg_simple)
flux_mod_pos_10deg_simple = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_pos_75deg_simple)
flux_mod_pos_75deg_simple = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_10deg_simple)
flux_mod_neg_10deg_simple = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_75deg_simple)
flux_mod_neg_75deg_simple = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_pos_10deg_full)
flux_mod_pos_10deg_full = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_pos_75deg_full)
flux_mod_pos_75deg_full = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_10deg_full)
flux_mod_neg_10deg_full = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_75deg_full)
flux_mod_neg_75deg_full = s.modulate(ekin, lis_langner_T)

s.update_parameters(par_pos_10deg_slow)
flux_mod_pos_10deg_slow = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_pos_75deg_slow)
flux_mod_pos_75deg_slow = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_10deg_slow)
flux_mod_neg_10deg_slow = s.modulate(ekin, lis_langner_T)
s.update_parameters(par_neg_75deg_slow)
flux_mod_neg_75deg_slow = s.modulate(ekin, lis_langner_T)

# data points digitized from Fig. 5 of the Strauss (2012) paper
str_flux_ekin, str_flux_vals = read_xy('input/Strauss2012_Fig5_LIS.txt')
str_pos_10deg_ekin, str_pos_10deg_flux = read_xy('input/Strauss2012_Fig5_Apos_10deg.txt')
str_pos_75deg_ekin, str_pos_75deg_flux = read_xy('input/Strauss2012_Fig5_Apos_75deg.txt')
str_neg_10deg_ekin, str_neg_10deg_flux = read_xy('input/Strauss2012_Fig5_Aneg_10deg.txt')
str_neg_75deg_ekin, str_neg_75deg_flux = read_xy('input/Strauss2012_Fig5_Aneg_75deg.txt')

# LIS digitized from Moskalenko et al., The Astrophysical Journal 565 (2002) 280-296
mos_lis_ekin, mos_lis_data = read_xy('input/moskalenko2002_lis.txt')
mos_lis_vals = 1.e-3 * mos_lis_data * np.power(mos_lis_ekin, -2.)

# Simplified HCS drift leads to artifacts at higher energies for A<0, alpha=75deg (need to debug this...),
# so we restrict the energy range for this case.
sel_simple = ekin<30.

plt.figure(figsize=(12,12))
#plt.plot(str_flux_ekin, str_flux_vals, 'k-', label='Reference: LIS')
#plt.plot(mos_lis_ekin, mos_lis_vals, 'b:', label='LIS Moskalenko 2002')
plt_lis, = plt.plot(ekin, lis_langner_T, 'k-', label='LIS')  # LIS: Langner+ (2004)
plt_pos_10deg_ref, = plt.plot(str_pos_10deg_ekin, str_pos_10deg_flux, 'o', c='blue', label=r'Reference')
plt_pos_10deg_ful, = plt.plot(ekin, flux_mod_pos_10deg_full, '-', c='blue', label=r'Solarprop 2.0')
plt_pos_10deg_sim, = plt.plot(ekin, flux_mod_pos_10deg_simple, ls='-', color='deepskyblue', lw=2, label=r'simple HCS drift')
plt_pos_10deg_slo, = plt.plot(ekin, flux_mod_pos_10deg_slow, c='aqua', ls=':', lw=2.5, label=r'$v\to v/2$')

plt_pos_75deg_ref, = plt.plot(str_pos_75deg_ekin, str_pos_75deg_flux, 'gs', label=r'Reference')
plt_pos_75deg_ful, = plt.plot(ekin, flux_mod_pos_75deg_full, 'g--', label=r'Solarprop 2.0')
plt_pos_75deg_sim, = plt.plot(ekin, flux_mod_pos_75deg_simple, ls='--', color='mediumseagreen', lw=2, label=r'simple HCS drift')
plt_pos_75deg_slo, = plt.plot(ekin, flux_mod_pos_75deg_slow, c='lime', ls='-.', lw=2.5, label=r'$v\to v/2$')

plt_neg_10deg_ref, = plt.plot(str_neg_10deg_ekin, str_neg_10deg_flux, 'D', c='mediumvioletred', label=r'Reference')
plt_neg_10deg_ful, = plt.plot(ekin, flux_mod_neg_10deg_full, '-', c='mediumvioletred', label=r'Solarprop 2.0')
plt_neg_10deg_sim, = plt.plot(ekin, flux_mod_neg_10deg_simple, c='fuchsia', ls='-', lw=2, label=r'simple HCS drift')
plt_neg_10deg_slo, = plt.plot(ekin, flux_mod_neg_10deg_slow, c='orchid', ls=':', lw=2.5, label=r'$v\to v/2$')

plt_neg_75deg_ref, = plt.plot(str_neg_75deg_ekin, str_neg_75deg_flux, 'r^', label=r'Reference')
plt_neg_75deg_ful, = plt.plot(ekin, flux_mod_neg_75deg_full, 'r--', label=r'Solarprop 2.0')
plt_neg_75deg_sim, = plt.plot(ekin[sel_simple], flux_mod_neg_75deg_simple[sel_simple], c='chocolate', ls='--', lw=2, label=r'simple HCS drift')
plt_neg_75deg_slo, = plt.plot(ekin, flux_mod_neg_75deg_slow, c='tomato', ls='-.', lw=2.5, label=r'$v\to v/2$')

plt.xlabel(r'$\mathregular{E_{kin}}$ (GeV)')
plt.ylabel(r'Flux (arbitrary units)')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.001, 100.)
plt.ylim(1.e-5, 1.e2)

# legends
title_fontsize = 16
legend_fontsize = 14
lis_legend = plt.gca().legend(handles=[plt_lis], loc='upper right', edgecolor='black')
plt.gca().add_artist(lis_legend)
leg_pos_10deg = plt.gca().legend(handles=[plt_pos_10deg_ref, plt_pos_10deg_ful, plt_pos_10deg_sim, plt_pos_10deg_slo],
                                 bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., linewidth=0,
                                 title=r'A>0  $\mathregular{\alpha}$ = 10°', fontsize=legend_fontsize, title_fontsize=title_fontsize, alignment='left')
plt.gca().add_artist(leg_pos_10deg)
leg_pos_75deg = plt.gca().legend(handles=[plt_pos_75deg_ref, plt_pos_75deg_ful, plt_pos_75deg_sim, plt_pos_75deg_slo],
                                 bbox_to_anchor=(1.02, 0.80), loc='upper left', borderaxespad=0., linewidth=0,
                                 title=r'A>0  $\mathregular{\alpha}$ = 75°', fontsize=legend_fontsize, title_fontsize=title_fontsize, alignment='left')
plt.gca().add_artist(leg_pos_75deg)
leg_neg_10deg = plt.gca().legend(handles=[plt_neg_10deg_ref, plt_neg_10deg_ful, plt_neg_10deg_sim, plt_neg_10deg_slo],
                                 bbox_to_anchor=(1.02, 0.6), loc='upper left', borderaxespad=0., linewidth=0,
                                 title=r'A<0  $\mathregular{\alpha}$ = 10°', fontsize=legend_fontsize, title_fontsize=title_fontsize, alignment='left')
plt.gca().add_artist(leg_neg_10deg)
leg_neg_75deg = plt.gca().legend(handles=[plt_neg_75deg_ref, plt_neg_75deg_ful, plt_neg_75deg_sim, plt_neg_75deg_slo],
                                 bbox_to_anchor=(1.02, 0.4), loc='upper left', borderaxespad=0., linewidth=0,
                                 title=r'A<0  $\mathregular{\alpha}$ = 75°', fontsize=legend_fontsize, title_fontsize=title_fontsize, alignment='left')
plt.gca().add_artist(leg_neg_75deg)

mylocator = mticker.LogLocator(subs=[1., 3., ])
myformatter = mticker.StrMethodFormatter("{x:g}")
#plt.gca().xaxis.set_major_locator(mylocator)
#plt.gca().xaxis.set_major_formatter(myformatter)
plt.grid(which='both')
plt.tight_layout()
plt.subplots_adjust(right=0.8)

if args.batch:
    plt.savefig('test_strauss_2012_fig5.pdf')
else:
    plt.show()
