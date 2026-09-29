#include "parkerspiral.h"

#include "particle.h"
#include "units.h"

#include <iostream>

using namespace Units;

/**
 * Construct the magnetic field.
 *
 * The numerical constant \f$ B_0 \f$ for the amplitude of the heliospheric
 * magnetic field is chosen such that \f$ |B_\oplus|=5\,\mathrm{nT} \f$ at Earth.
 *
 * @param modelInfo %Model information.
 * @param _wavyCurrentSheet Use a wavy current sheet instead of a flat one?
 */
ParkerSpiral::ParkerSpiral(const ModelInformation& modelInfo, bool _wavyCurrentSheet) :
	wavyCurrentSheet(_wavyCurrentSheet)
{
	// amplitude such that B_earth = 5 nT
	static constexpr double B0 = 3.4 * nanotesla;

	double BfieldScaling = modelInfo.getFloatOptionWithDefault("BfieldScaling", 1.0, false);
	BfieldNormalization = BfieldScaling * B0;
	r0_sqr = std::pow(1.0 * AU, 2);

	OmegaSun = 2. * pi * rad / (25.4 * days);
	polarity = modelInfo.getPolarity();
	tiltAngle = modelInfo.getTiltAngle();

	double Gamma_e = OmegaSun / solarWindSpeed(Particle()) * 1. * AU;
	Bearth = BfieldNormalization * std::sqrt(1.0 + Gamma_e * Gamma_e);
	std::cout << "[ParkerSpiral] B_earth = " << Bearth / nanotesla << " nT" << std::endl;
}

double ParkerSpiral::solarWindSpeed(const Particle&) const
{
	static const double V = 400 * km / s;
	return V;
}

/** Calculate cached values at the beginning of each pseudo-particle step.
 *
 * \attention Make sure this function is called before any of the getters in this class.
 */
void ParkerSpiral::calculate(const Particle& p)
{
	Gamma = OmegaSun / solarWindSpeed(p) * p.r * p.sinTheta;
	onePlusGamma2 = 1.0 + Gamma * Gamma;

	driftConstant = polarity * 2. / 3. * p.rigidity * p.velocity / BfieldNormalization / r0_sqr * p.r * Gamma / pow(onePlusGamma2, 2);
	driftConstantTimesHeaviside = getHeaviside(p) * driftConstant;
}

double ParkerSpiral::getBOverB2_r(const Particle& p) const
{
	double r = p.r;
	int sign = polarity * getHeaviside(p);
	return sign * 1.0 / BfieldNormalization / r0_sqr * (r * r) / onePlusGamma2;
}

double ParkerSpiral::getBOverB2_phi(const Particle& p) const
{
	double r = p.r;
	int sign = polarity * getHeaviside(p);
	return -sign * 1.0 / BfieldNormalization / r0_sqr * (r * r) * Gamma / onePlusGamma2;
}

double ParkerSpiral::getDriftVelocity_r(const Particle& p) const
{
	return -driftConstantTimesHeaviside * p.cosTheta / p.sinTheta;
}

double ParkerSpiral::getDriftVelocity_theta(const Particle& p) const
{
	return driftConstantTimesHeaviside * (1.0 + onePlusGamma2);
}

double ParkerSpiral::getDriftVelocity_phi(const Particle& p) const
{
	return driftConstantTimesHeaviside * Gamma * p.cosTheta / p.sinTheta;
}

double ParkerSpiral::getBearthOverB(const Particle& p) const
{
	return Bearth / getB(p);
}

double ParkerSpiral::getLarmorRadius(const Particle& p) const
{
	return p.absRigidity / getB(p);
}

double ParkerSpiral::getB(const Particle& p) const
{
	return BfieldNormalization * r0_sqr / (p.r * p.r) * std::sqrt(onePlusGamma2);
}

int ParkerSpiral::getHeaviside(const Particle& p) const
{
	if (wavyCurrentSheet)
	{
		double thetaPrime = getThetaPrime(p);
		return 1 - 2 * (p.getTheta() > thetaPrime);
	}
	else
	{
		if (p.cosTheta < 0.0)
		{
			return -1;
		}
		else
		{
			return 1;
		}
	}
}

double ParkerSpiral::getTanPsi(const Particle& p) const
{
	return OmegaSun / solarWindSpeed(p) * p.r * p.sinTheta;
}

double ParkerSpiral::getSinPsi(const Particle& p) const
{
	double x = getTanPsi(p);
	return x / std::sqrt(x * x + 1.);
}

double ParkerSpiral::getCosPsi(const Particle& p) const
{
	double x = getTanPsi(p);
	return 1. / std::sqrt(x * x + 1.);
}

/**
 * The simplified formula is used to calculate the location of the current sheet:
 * \f[
 *   \theta^\prime = \frac{\pi}{2} + \alpha\sin\left(\phi+\frac{\Omega r}{V_{sw}}\right).
 * \f]
 */
double ParkerSpiral::getThetaPrime(const Particle& p) const
{
	return pi / 2.0 + tiltAngle * std::sin(p.phi + OmegaSun / solarWindSpeed(p) * p.r);
}
