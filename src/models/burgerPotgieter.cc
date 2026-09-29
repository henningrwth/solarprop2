#include "burgerPotgieter.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo %Model information
 */
BurgerPotgieter::BurgerPotgieter(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo, false)
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 3000., false) * s;

	kappa0 = 5.e21 * cm2/s;

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();
	if (modelInfoKappa0 > 0.0)
	{
		kappa0 = modelInfoKappa0;
	}

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0, false);
	kappa0 *= kappaScaling;

	perPart = 0.05;

	heliosphereBoundary = 10.*AU;
}


/**
 * Calculate the quantities needed several times in one particle step.
 */
void BurgerPotgieter::calculate(const Particle& p) {

	B.calculate(p);

	kappa_parallel = kappa0*sqrt(p.absRigidity/(1*GV/c))*p.beta*B.getBearthOverB(p);
	kappa_perpendicular = perPart*kappa_parallel;

	kappa_rr = 1.0/B.onePlusGamma2 * (kappa_parallel + kappa_perpendicular*B.Gamma*B.Gamma);
	kappa_thetatheta = kappa_perpendicular;
}


double BurgerPotgieter::getDriftR(const Particle& p) const
{
	return B.getDriftVelocity_r(p);
}


double BurgerPotgieter::getDriftTheta(const Particle& p) const
{
	return B.getDriftVelocity_theta(p);
}


double BurgerPotgieter::getDriftSheetR(const Particle& p) const {
	double alpha = B.tiltAngle/rad;
	double deltaTheta = 2.*p.absRigidity*B.solarWindSpeed(p) / (B.BfieldNormalization * B.r0_sqr * B.OmegaSun * cos(alpha));

	// Look if the particle is inside the "cone" of the wavy heliospheric current sheet
	if (std::abs(p.getTheta() - pi/2.) < alpha+deltaTheta)
	{
		return B.polarity*p.charge * p.velocity/6.0 * cos(alpha) * deltaTheta/sin(alpha+deltaTheta);
	}

	return 0.;
}


double BurgerPotgieter::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	double Gamma2 = B.onePlusGamma2 - 1.0;
	return kappa_rr/p.r * (4.+3.*Gamma2) / B.onePlusGamma2
			+ kappa_parallel / p.r *2.*Gamma2 / pow(B.onePlusGamma2,2) * (perPart-1.);
}


double BurgerPotgieter::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	return -1./(p.r*p.r) * kappa_perpendicular * p.cosTheta / B.onePlusGamma2 * (1.0 + B.onePlusGamma2);
}
