#include "modifiedspiral.h"

#include "particle.h"
#include "units.h"

#include <iostream>

using namespace Units;

/**
 * Constructor.
 *
 * @param modelInfo %Model information.
 * @param _isHighSolarActivity Flag: period of high solar activity? Determines the latitudinal dependence of the solar wind speed.
 * @param _polarity Heliospheric polarity (\c +1 or \c -1).
 * @param _phase Heliospheric phase. Possible values: \c +1 ascending, \c -1 descending, \c 0 unknown.
 * @param _tiltAngle Tilt angle of the heliospheric current sheet.
 * @param _Bearth Amplitude \f$ B_\oplus \f$ of the heliospheric magnetic field at Earth.
 * @param _vswMin \f$ v_\mathrm{omni} \f$, minimum solar wind speed.
 * @param _delta0 Amplitude \f$ \delta_0 \f$ for the polar field correction.
 *
 */
ModifiedSpiral::ModifiedSpiral(const ModelInformation& modelInfo, bool _isHighSolarActivity, int _polarity, int _phase, double _tiltAngle, double _Bearth,
			       double _vswMin, double _delta0) :
	Bearth(_Bearth),
	r0_sqr(std::pow(1.0 * AU, 2)),
	rSun(6.96342e8 * meter),
	solarWindSpeedMin(_vswMin),
	solarWindSpeedMax(760. * km / s),
	isHighSolarActivity(_isHighSolarActivity),
	OmegaSun(2. * pi * rad / (25.4 * days)),
	polarity(_polarity),
	phase(_phase),
	tiltAngle(_tiltAngle),
	delta0(_delta0),
	thetaThresholdForCorrection(1.7 * deg)
{
	double Gamma_e = OmegaSun / solarWindSpeed(Particle()) * 1. * AU;
	double t_e = delta0 * (1. * AU / rSun);
	double N_e = 1.0 + Gamma_e * Gamma_e + t_e * t_e;
	BfieldNormalization = Bearth / std::sqrt(N_e);
}

/** Calculate solar wind speed. */
double ModifiedSpiral::solarWindSpeed(const Particle& p) const
{
	if (isHighSolarActivity)
	{
		return solarWindSpeedMin;
	}

	double absCosTheta = std::abs(p.cosTheta);

	if (absCosTheta > 0.8660254)
	{
		// cos(30deg) = sqrt(3)/2 = 0.8660254
		return solarWindSpeedMax;
	}

	double V = solarWindSpeedMin * (1.0 + absCosTheta);

	return std::min(V, solarWindSpeedMax);
}

/** Calculate \f$ \mathrm{d}V_{sw}/\mathrm{d}\cos\theta \f$. */
double ModifiedSpiral::dvsw_dcostheta(const Particle& p) const
{
	if (isHighSolarActivity || p.cosTheta == 0.0)
	{
		return 0.0;
	}
	double absCosTheta = std::abs(p.cosTheta);

	if (absCosTheta > 0.8660254)
	{
		return 0.0;
	}

	if (vsw >= solarWindSpeedMax)
	{
		return 0.0;
	}

	return solarWindSpeedMin * p.cosTheta / absCosTheta;
}

/** Calculate cached values at the beginning of each pseudo-particle step.
 *
 * \attention Make sure this function is called before any of the getters in this class.
 */
void ModifiedSpiral::calculate(const Particle& p)
{
	vsw = solarWindSpeed(p);
	Gamma = OmegaSun / vsw * p.r * p.sinTheta;
	onePlusGamma2 = 1.0 + Gamma * Gamma;
	double rOverRsun = p.r / rSun;
	heavi = getHeaviside(p);
	sign = polarity * heavi;
	delta = deltaFunction(p);
	t = delta * rOverRsun;
	dtdmu = -rOverRsun * deltaPrimeFunction(p) / p.sinTheta;
	N = onePlusGamma2 + t * t;
	tau = 1.0 + t * t;
	Gammaprime = dGamma_dcostheta(p);

	driftConstant = 2. * polarity / 3. * p.rigidity * p.velocity / BfieldNormalization / r0_sqr * p.r * Gamma / pow(N, 2);
}

/** Component in radial direction of \f$ \vec{B}/B^2 \f$. */
double ModifiedSpiral::getBOverB2_r(const Particle& p) const
{
	return sign / BfieldNormalization / r0_sqr * std::pow(p.r, 2) / N;
}

/** Component in polar direction of \f$ \vec{B}/B^2 \f$. */
double ModifiedSpiral::getBOverB2_theta(const Particle& p) const
{
	return sign / BfieldNormalization / r0_sqr * std::pow(p.r, 3) / rSun * delta / N;
}

/** Component in azimuthal direction of \f$ \vec{B}/B^2 \f$. */
double ModifiedSpiral::getBOverB2_phi(const Particle& p) const
{
	return -sign / BfieldNormalization / r0_sqr * std::pow(p.r, 2) * Gamma / N;
}

/** Radial component of standard gradient and curvature drift velocity field outside HCS. */
double ModifiedSpiral::getDriftVelocity_r(const Particle& p) const
{
	double deltaPrime = deltaPrimeFunction(p);

	double term1 = delta != 0.0 ? t * t * deltaPrime / delta : 0.0;

	double vprime = -dvsw_dcostheta(p) * p.sinTheta;
	double term2 = 0.5 * vprime / vsw * (tau - Gamma * Gamma);

	double term3 = -p.getCotTheta() * tau;

	return heavi * driftConstant * (term1 + term2 + term3);
}

/** Polar component of standard gradient and curvature drift velocity field outside HCS. */
double ModifiedSpiral::getDriftVelocity_theta(const Particle& p) const
{
	return heavi * driftConstant * (1.0 + N);
}

/** Modification factor \f$ \delta(\theta) \f$ for the magnetic field. */
double ModifiedSpiral::deltaFunction(const Particle& p) const
{
	double absCosTheta = std::abs(p.cosTheta);

	// polar region?
	if (absCosTheta < 0.8660254) // cos(30deg) = sqrt(3)/2 = 0.8660254
		return 0.0;

	double cosThetaMax = std::cos(thetaThresholdForCorrection);

	if (std::abs(p.cosTheta) > cosThetaMax)
	{
		double deltaMax = delta0 / std::sin(thetaThresholdForCorrection);
		return heavi * deltaMax;
	}

	return delta0 * heavi / p.sinTheta;
}

/** Derivative of \f$ \delta(\theta) \f$ with respect to the polar angle \f$ \theta \f$. */
double ModifiedSpiral::deltaPrimeFunction(const Particle& p) const
{
	double absCosTheta = std::abs(p.cosTheta);

	// polar region?
	if (absCosTheta < 0.8660254)
		return 0.0;

	double cosThetaMax = std::cos(thetaThresholdForCorrection);

	if (std::abs(p.cosTheta) > cosThetaMax)
		return 0.0;

	return -delta0 * heavi * p.cosTheta / std::pow(p.sinTheta, 2);
}

/** Calculate \f$ \mathrm{d}\Gamma/\mathrm{d}\cos\theta \f$. */
double ModifiedSpiral::dGamma_dcostheta(const Particle& p) const
{
	double mu = p.getCosTheta();
	return -Gamma * (mu / (1.0 - mu * mu) + dvsw_dcostheta(p) / vsw);
}

/** Calculate \f$ B_\oplus / B \f$. */
double ModifiedSpiral::getBearthOverB(const Particle& p) const
{
	return Bearth / getB(p);
}

/** Larmor radius \f$ R_L = |R/B| \f$. */
double ModifiedSpiral::getLarmorRadius(const Particle& p) const
{
	return p.absRigidity / getB(p);
}

/** Absolute value of magnetic field. */
double ModifiedSpiral::getB(const Particle& p) const
{
	return BfieldNormalization * r0_sqr / std::pow(p.r, 2) * std::sqrt(N);
}

/** Heaviside function \f$ H(\theta - \theta^\prime) \f$ where \f$ \theta^\prime \f$ is the location of the HCS. */
int ModifiedSpiral::getHeaviside(const Particle& p) const
{
	double thetaPrime = getThetaPrime(p);
	return 1 - 2 * (p.getTheta() > thetaPrime);
}

/**
 * Location of heliospheric current sheet (simplified formula).
 *
 * The simplified formula is used to calculate the location of the current sheet:
 * \f[
 *   \theta^\prime = \frac{\pi}{2} + \alpha\sin\left(\phi+\frac{\Omega r}{V_{sw}}\right).
 * \f]
 */
double ModifiedSpiral::getThetaPrime(const Particle& p) const
{
	return pihalf + tiltAngle * sin(p.phi + OmegaSun / vsw * p.r);
}
