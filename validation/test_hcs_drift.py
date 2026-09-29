#! /usr/bin/env python3

"""
Plot the heliospheric current sheet and test tracking particles along it.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.optimize as opt
import numdifftools as nd

from sproptools.kinematics import Species, KinematicQuantity
from solarprop.units import AU, nanotesla
import solarprop

plt.style.use('solarprop.mplstyle')

argparser = argparse.ArgumentParser(description=__doc__)
argparser.add_argument('-a', '--tilt-angle', type=float, metavar='DEG', default=15., help='Tilt angle (in degrees).')
argparser.add_argument('-m', '--minimization-method', choices=['Nelder-Mead', 'Powell'], default='Nelder-Mead', help='Method for finding point of closest distance on HCS.')
argparser.add_argument('-r', '--r-max', type=float, default=10., metavar='AU', help='Radial extension of heliosphere (in AU).')
argparser.add_argument('-s', '--start-coords', nargs=3, type=float, default=[1.0, 90., 120.], help='Start coordinates: r in AU, theta in degrees, phi in degrees.')
argparser.add_argument('-v', '--vmax', type=float, default=5., help='Maximum for distance color scale (in AU).')
argparser.add_argument('--restricted-view', action='store_true', help='Zoom into 2D distance plot?')
args = argparser.parse_args()

proton = Species.from_name('proton')

omega_sun = 2.* np.pi / (25.4 * 24. * 3600.)  # rad s^(-1)  (solarprop value): 2.863073e-6 rad/s
##omega_sun = 2.66e-6  # # rad s^(-1) (Value from Raath PhD (2016))
#omega_sun = 2.9e-6  # rad s^(-1)  (value from Jokipii and Thomas (1981))

v_sw = 400. * 1.e3 / sc.au  # AU / s
tilt_angle = np.deg2rad(args.tilt_angle)
phi_0 = 0.

r_max = args.r_max
ndots = int(r_max*10)
clipping_factor = 1.0  # 0.75

parameters = { 'model': 'Strauss2012',
               'mass': proton.mass,
               'charge': proton.charge,
               'polarity': 1,
               'angle': args.tilt_angle,
               'forceNumericalDerivatives': False,
              }

def theta_prime_hcs(alpha, r, phi):
    return np.pi/2. + np.arcsin(np.sin(alpha) * np.sin(phi - phi_0 + omega_sun/v_sw * r))

def theta_prime_hcs_tan(alpha, r, phi):
    return np.pi/2. - np.arctan(np.tan(alpha) * np.sin(phi - phi_0 + omega_sun/v_sw * r))

def theta_prime_hcs_simple(alpha, r, phi):
    return np.pi/2. + alpha * np.sin(phi - phi_0 + omega_sun/v_sw * r)

def rho_phi(x, y):
    return np.hypot(x, y), np.atan2(y, x)

def x_y_z(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def beta_sign(r, phi):
    return np.sign(np.cos(phi + phi_0 + omega_sun/v_sw * r))

def tan_psi(r, theta):
    return omega_sun * r * np.sin(theta) / v_sw

def sin_psi(r, theta):
    x = tan_psi(r, theta)
    return x / np.sqrt(x**2 + 1.)

def cos_psi(r, theta):
    x = tan_psi(r, theta)
    return 1. / np.sqrt(x**2 + 1.)

def sin2_psi(r, theta):
    x = tan_psi(r, theta)
    x2 = x**2
    return x2 / (x2 + 1.)

def tanbeta2_func(r, theta_p):
    sin_thetap_sqr = np.sin(theta_p)**2
    cos_thetap_sqr = 1. - sin_thetap_sqr
    numer = np.clip(np.sin(tilt_angle)**2 - cos_thetap_sqr, 0., None)
    return (omega_sun*r/v_sw)**2 / sin2_psi(r, theta_p) * numer / sin_thetap_sqr

def cosbeta_func(r, theta_p):
    x2 = tanbeta2_func(r, theta_p)
    return 1./np.sqrt(x2+1.)

def sinbeta_func(r, theta_p, phi):
    sign = beta_sign(r, phi)
    x2 = tanbeta2_func(r, theta_p)
    return sign * np.sqrt(x2/(x2+1.))

def sincosbeta_func(r, theta_p, phi):
    sign = beta_sign(r, phi)
    x2 = tanbeta2_func(r, theta_p)
    cosbeta = 1./np.sqrt(x2+1.)
    sinbeta = sign * np.sqrt(x2) * cosbeta
    return sinbeta, cosbeta

def beta_func(r, theta_p, phi):
    sign = beta_sign(r, phi)
    return sign * np.atan(np.sqrt(tanbeta2_func(r, theta_p)))

def magfield_amp(r, theta):
    return 3.4e-9/r**2 * np.sqrt(1. + tan_psi(r, theta)**2)

def larmor_radius(species, Ekn_in_GeV, r, theta):
    rigidity = species.convert(Ekn_in_GeV, KinematicQuantity.KineticEnergyPerNucleon, KinematicQuantity.Rigidity)
    magfield = magfield_amp(r, theta)
    return (rigidity / sc.c) / magfield

def D2(r1, r2, theta1, theta2, phi1, phi2):
    return r1**2 + r2**2 - 2*r1*r2*(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2) + np.cos(theta1)*np.cos(theta2))

def D2_hcs(r1, r2, theta1, phi1, phi2):
    return D2(r1, r2, theta1, theta_prime_hcs(tilt_angle, r2, phi2), phi1, phi2)

def L2_hcs(r1, r2, theta1, phi):
    return r1**2 * (tilt_angle * np.sin(omega_sun/v_sw*r2 + phi) - theta1 + np.pi/2)**2 + (r2-r1)**2

def dL2_hcs(r1, r2, theta1, phi):
    k = omega_sun / v_sw
    return k*tilt_angle*r1**2 * (2*tilt_angle*np.sin(k*r2 + phi) - 2*theta1 + np.pi) * np.cos(k*r2 + phi) + 2*(r2-r1)

def distance_to_hcs(r, theta, phi, alternative_start=None, verbose=False):
    rphi0_default = [r, phi]
    res = opt.minimize(lambda x: D2_hcs(r, x[0], theta, phi, x[1]), x0=rphi0_default,
                   bounds=[(0., r_max), (0., 2.*np.pi)], method=args.minimization_method, options={'maxiter': 100})
    distance = np.sqrt(res.fun)
    better_start = '(point)'

    if alternative_start is not None:
        res_alt = opt.minimize(lambda x: D2_hcs(r, x[0], theta, phi, x[1]), x0=alternative_start,
                               bounds=[(0., r_max), (0., 2.*np.pi)], method=args.minimization_method, options={'maxiter': 100})
        distance_alt = np.sqrt(res_alt.fun)
        if distance_alt < distance:
            res = res_alt
            better_start = '(alt)'

    if not res.success:
        print(res.message)
    coords_hcs = (res.x[0], theta_prime_hcs(tilt_angle, res.x[0], res.x[1]), res.x[1])  # (r, theta, phi)
    if verbose:
        print(f'{r:g} AU {np.rad2deg(theta):g} deg {np.rad2deg(phi):g} deg: D= {distance:g} AU, {res.nit} iterations {better_start}, point at HCS: ({coords_hcs[0]:g} AU {np.rad2deg(coords_hcs[1]):g} deg {np.rad2deg(coords_hcs[2]):g} deg)')
    return distance, coords_hcs


ekn = 0.1  # GeV/n
r0 = float(args.start_coords[0])
phi0 = np.deg2rad(float(args.start_coords[2]))
theta0 = np.deg2rad(float(args.start_coords[1]))
start_point = np.array(x_y_z(r0, theta0, phi0))
speed_of_light = sc.c/sc.au
beta = proton.convert(ekn, KinematicQuantity.KineticEnergyPerNucleon, KinematicQuantity.Beta)
rL = larmor_radius(proton, ekn, r0, theta0)
print(f'({args.start_coords}): rL (Ekn={ekn:g} GeV)= {rL:g} AU, beta={beta:g}')
rL_b = larmor_radius(proton, ekn, r_max, np.pi/2)
print(f'rL at r={r_max:g} AU on the ecliptic: {rL_b:g} AU')

dt = 0.004 * 24. * 3600.

if not args.restricted_view:
    r_lin = np.linspace(0., r_max, ndots)
    phi_lin = np.linspace(0., clipping_factor*2.*np.pi, ndots)
else:
    r_lin = np.linspace(r0-2., r0+2., ndots)
    phi_lin = np.linspace(phi0-0.1, phi0+0.1, ndots)
R, Phi = np.meshgrid(r_lin, phi_lin)

theta = theta_prime_hcs(tilt_angle, R, Phi)
X, Y, Z = x_y_z(R, theta, Phi)

r_proj = np.linspace(0.01, r_max, 5*ndots)
theta_proj = theta_prime_hcs(tilt_angle, r_proj, 0.)
theta_proj_simple = theta_prime_hcs_simple(tilt_angle, r_proj, 0.)
theta_proj_tan = theta_prime_hcs_tan(tilt_angle, r_proj, 0.)

x_proj, _ , z_proj = x_y_z(r_proj, theta_proj, 0.)
x_proj_simple, _ , z_proj_simple= x_y_z(r_proj, theta_proj_simple, 0.)
x_proj_tan, _ , z_proj_tan= x_y_z(r_proj, theta_proj_tan, 0.)

dtheta_prime_dr_func = lambda r, phi: nd.Derivative(lambda r: theta_prime_hcs(tilt_angle, r, phi))(r)
dtheta_prime_dr_proj = dtheta_prime_dr_func(r_proj, 0.)


beta_proj = beta_func(r_proj, theta_proj, 0.)
sinbeta_proj, cosbeta_proj = sincosbeta_func(r_proj, theta_proj, 0.)

# track one particle along HCS
nsteps = 300
# particle coordinates after each step
r_p = np.zeros(nsteps)
theta_p = np.zeros(nsteps)
phi_p = np.zeros(nsteps)

r_p[0] = r0
phi_p[0] = phi0
theta_p[0] = theta0

# HCS footpoints for each step
r_foot = np.zeros(nsteps-1)
theta_foot = np.zeros(nsteps-1)
phi_foot = np.zeros(nsteps-1)

for istep in range(nsteps-1):
    r = r_p[istep]
    theta = theta_p[istep]
    phi = phi_p[istep]
    if r < r_max:
        alternative_start = None
        if istep > 0:
            alternative_start = [r_foot[istep-1], phi_foot[istep-1]]
        d, coords_hcs = distance_to_hcs(r, theta, phi, alternative_start)
        r_foot[istep], theta_foot[istep], phi_foot[istep] = coords_hcs

        sinbeta, cosbeta = sincosbeta_func(*coords_hcs)
        sinpsi = sin_psi(coords_hcs[0], coords_hcs[1])
        cospsi = cos_psi(coords_hcs[0], coords_hcs[1])
        vec_r = cosbeta * sinpsi
        vec_theta = sinbeta
        vec_phi = cosbeta * cospsi

        r_larmor = larmor_radius(proton, ekn, r, theta)
        x = np.clip(d / r_larmor, 0., 2.)
        #print('x=', x)
        velocity = (0.457 - 0.412*x + 0.0915*x**2) * beta * speed_of_light
        dr = velocity * dt * vec_r
        dtheta = velocity * dt * vec_theta / r
        dphi = velocity * dt * vec_phi / r / np.sin(theta)
    else:
        dr = dtheta = dphi = 0.
        r_foot[istep], theta_foot[istep], phi_foot[istep] = r_foot[istep-1], theta_foot[istep-1], phi_foot[istep-1]

    r_p[istep+1] = r + dr
    theta_p[istep+1] = theta + dtheta
    # reset particle to HCS instead:
    #theta_p[istep+1] = theta_prime_hcs(tilt_angle, r+dr, phi+dphi)
    phi_p[istep+1] = phi + dphi

xyz_particle = np.array([x_y_z(r_p, theta_p, phi_p)])[0]
xyz_foot = np.array([x_y_z(r_foot, theta_foot, phi_foot)])[0]

# calculate rho=sqrt(x^2+y^2) vs z and theta of particle track and corresponding values for HCS
theta_hcs_p = theta_prime_hcs(tilt_angle, r_p, phi_p)
z_hcs_p = x_y_z(r_p, theta_hcs_p, phi_p)[2]
rho_p = np.hypot(xyz_particle[0], xyz_particle[1])
rho_foot = np.hypot(xyz_foot[0], xyz_foot[1])

# more detailed analysis for start point
distance_0 = np.sqrt(D2_hcs(r0, R, theta0, phi0, Phi))
D_0, coords_hcs = distance_to_hcs(r0, theta0, phi0, alternative_start=None, verbose=True)
point_0 = np.array(x_y_z(r0, theta0, phi0))
point_hcs = np.array(x_y_z(*coords_hcs))

L_0 = np.sqrt(L2_hcs(r0, R, theta0, phi0))
L_0_lin = np.sqrt(L2_hcs(r0, r_lin, theta0, phi0))
dL2_0_lin = dL2_hcs(r0, r_lin, theta0, phi0)

#
# solarprop calculation for comparison
#
s = solarprop.Solarprop(parameters)
B_sprop = s.magnetic_field(r0 * AU, theta0, phi0) / nanotesla
print(B_sprop, 'nT')


# plotting
plt.rcParams['font.size'] = 14.0

# distance from start point to HCS
fig, ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'}, figsize=(12,10))
im = ax.pcolormesh(Phi, R, distance_0, cmap='magma_r', vmin=0., vmax=args.vmax)
if not args.restricted_view:
    ax.set_rlim(0., r_max)
else:
    ax.set_rlim(r0-2., r0+2.)
    ax.set_thetalim(phi0-0.1, phi0+0.1)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('$D$ (AU)')
ax.plot(phi0, r0, 'go')
ax.plot(coords_hcs[2], coords_hcs[0], 'r+')

# distance in simplified, 1d radial version
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(r_lin, dL2_0_lin, 'r--')
ax.plot(r_lin, L_0_lin, 'r-')
ax.set_xlabel('$r$ (AU)')
ax.set_ylabel('$L$ (AU)')
#ax.set_ylim(bottom=0.)
ax.axvline(r0)
ax.axvline(coords_hcs[0], color='g', ls=':')

# particle track and corresponding values for HCS: plot rho vs z
fig, ax = plt.subplots(figsize=(18,12))
ax.plot(rho_p, xyz_particle[2], 'g.-')
ax.plot(rho_foot, xyz_foot[2], 'r+')
ax.plot(rho_p, z_hcs_p, 'k:')
ax.set_xlabel(r'$\rho$ (AU)')
ax.set_ylabel('$z$ (AU)')
fig.subplots_adjust(hspace=0.01)
fig.tight_layout()
    
# reproduce Figure 4 of Strauss+ (2012) paper
fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(5, 12))
ax[0].plot(x_proj, z_proj, 'b-')
ax[0].plot(x_proj_simple, z_proj_simple, 'g:')
ax[0].plot(x_proj_tan, z_proj_tan, 'r-.')
ax[0].plot(r_proj,  r_proj * np.sin(tilt_angle), '--', color='lightgrey')
ax[0].plot(r_proj, -r_proj * np.sin(tilt_angle), '--', color='lightgrey')
ax[-1].set_xlabel('$x$ (AU)')
ax[-1].set_xlim(left=0.)
ax[0].set_ylabel('$z$ (AU)')
ax[0].grid()

ax[1].plot(x_proj, dtheta_prime_dr_proj, 'b-')
ax[1].set_ylabel(r'$\mathrm{d}\theta^\prime/\mathrm{d}r$')
ax[1].grid()

ax[2].plot(x_proj, beta_proj, 'b-')
ax[2].set_ylabel(r'$\beta$ (rad)')
ax[2].grid()

ax[3].plot(x_proj, sinbeta_proj, 'b-')
ax[3].plot(x_proj, cosbeta_proj, 'b-.')
ax[3].set_ylabel(r'$\sin\beta$ and $\cos\beta$')
ax[3].grid()

fig.subplots_adjust(hspace=0.01)
fig.tight_layout()

# 3D plot of HCS, start point, and particle track
fig = plt.figure(figsize=(16, 14))
ax = fig.add_subplot(projection='3d')
ax.set_zlim(-r_max, r_max)
ax.set_xlabel('$x$ (AU)')
ax.set_ylabel('$y$ (AU)')
ax.set_zlabel('$z$ (AU)')

ax.plot_surface(X, Y, Z, cmap=plt.cm.RdYlGn, alpha=0.1)
ax.plot(*start_point, 'bo')
ax.plot(*point_0, 'ko')
ax.plot(*point_hcs, 'r+')
ax.plot([point_0[0], point_hcs[0]], [point_0[1], point_hcs[1]], [point_0[2], point_hcs[2]], 'r-')
ax.plot(xyz_particle[0], xyz_particle[1], xyz_particle[2], 'g.--')
ax.plot(xyz_foot[0], xyz_foot[1], xyz_foot[2], 'r+')
plt.show()
