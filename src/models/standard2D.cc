#include "standard2D.h"
#include "units.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo %Model information
 */
Standard2D::Standard2D(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo, false)
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;

	if (modelInfo.hasOption("dt"))
	{
		dt_max = modelInfo.getFloatOption("dt") * second;
		dt_max_autoset = false;
	}

	kappa0 = (137.0*MV/modelInfo.getModulationPotential() - 0.061) * 0.5 * 7.5e-4 * AU*AU/s;
	perPart = 0.02;
	heliosphereBoundary = 100.*AU;

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();
	if (modelInfoKappa0 > 0.0)
	{
		kappa0 = modelInfoKappa0;
	}

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0);
	kappa0 *= kappaScaling;

	std::cout << "Diffusion tensor rescaling: " << kappaScaling << std::endl;
	std::cout << "Determined diffusion constant (cm^2 s^-1 GV^-1): " << kappa0 / (cm2/s/GV) << std::endl;

	alphaH = 0.;
	tanAlphaH = 0.;
	f = 0.;
	fprime = 0.;
}


/**
 * Calculate the quantities needed several times in one particle step.
 */
void Standard2D::calculate(const Particle& p)
{
	B.calculate(p);

	if (dt_max_autoset && B.polarity*p.charge > 0)
	{
		dt_max = 1000*s;
	}
	else
	{
		dt_max = 5000*s;
	}

	double rigidityValue = p.absRigidity / (GV/c);

	if (rigidityValue < 0.1)
		rigidityValue = 0.1;

	kappa_parallel = kappa0*rigidityValue*p.beta*B.getBearthOverB(p);

	if (B.polarity*p.charge > 0)
	{
		kappa_parallel *= 0.07;
	}

	kappa_perpendicular = perPart*kappa_parallel;

	// This is a bug in vanilla solarprop: 2r_L/r is in rad, but angle is given in deg.
	alphaH = acos(pow(pi/sin((B.tiltAngle/deg + 2.*B.getLarmorRadius(p)/p.r)*deg)-1.,-1));
	tanAlphaH = tan(alphaH);

	kappa_rr = (kappa_parallel + kappa_perpendicular*pow(B.Gamma,2)) / B.onePlusGamma2;
	kappa_thetatheta = kappa_perpendicular;

	double arg = (1. - 2. * p.getTheta() / pi)*tanAlphaH;

	f = atan(arg)/alphaH;
	fprime = -2./(alphaH*pi)*tanAlphaH/(1.0 + arg*arg);
}


double Standard2D::getDriftR(const Particle& p) const
{
	return -B.driftConstant * p.cosTheta/p.sinTheta * f;
}


double Standard2D::getDriftTheta(const Particle& p) const
{
	return B.driftConstant * (1.0 + B.onePlusGamma2) * f;
}

double Standard2D::getDriftSheetR(const Particle& p) const
{
	// suppress enormous HCS drift close to Earth, which leads to unphysical results
	if (p.r < 1.5*AU)
		return 0.0;

	return -B.driftConstant/2. * B.onePlusGamma2 * fprime;
}

double Standard2D::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	return 2./p.r * (kappa_rr +
			 (kappa_perpendicular - kappa_parallel) * pow(B.Gamma,2) / pow(B.onePlusGamma2,2) +
			 kappa_rr * (1.0 + 0.5*pow(B.Gamma,2)) / B.onePlusGamma2);
}


double Standard2D::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	return -1./(p.r*p.r)*kappa_perpendicular*p.cosTheta/B.onePlusGamma2*(1.0+B.onePlusGamma2);
}
