#include "jokipiiKopriva.h"

#include "numdiff.h"
#include "particle.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo %Model information
 */
JokipiiKopriva::JokipiiKopriva(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo, false)
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 1000.) * s;

	kappa0 = 5.0e+21 * cm2/s;
	perPart = 0.1;
	heliosphereBoundary = 10.*AU;

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();
	if (modelInfoKappa0 > 0.0)
	{
		kappa0 = modelInfoKappa0;
	}

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0);
	kappa0 *= kappaScaling;

	forceNumericalDerivatives = modelInfo.getBooleanOptionWithDefault("forceNumericalDerivatives", false);

	if (forceNumericalDerivatives)
	{
		std::cout << "[JokipiiKopriva] Forcing numerical calculation of derivatives in diffusion and drift terms..." << std::endl;
	}
}

/**
 * Calculate the quantities needed several times in one particle step.
 */
void JokipiiKopriva::calculate(const Particle& p)
{
	if (!nd && forceNumericalDerivatives)
	{
		nd = std::make_shared<NumericalDerivatives>(this);
	}

	B.calculate(p);

	kappa_parallel = kappa0*sqrt(p.absRigidity/(1*GV/c))*p.beta;
	kappa_perpendicular = perPart*kappa_parallel;

	kappa_rr = kappa_parallel/B.onePlusGamma2 + kappa_perpendicular*B.Gamma*B.Gamma/B.onePlusGamma2;
	kappa_thetatheta = kappa_perpendicular;
}


double JokipiiKopriva::getDriftR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateRadialGradientCurvatureDriftVelocity(p, false);

	return B.getDriftVelocity_r(p);
}

double JokipiiKopriva::getDriftTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculatePolarGradientCurvatureDriftVelocity(p);

	return B.getDriftVelocity_theta(p);
}

double JokipiiKopriva::getDriftSheetR(const Particle& p) const
{
	// Distance to the current sheet
	double d = fabs(p.r*p.cosTheta);
	double rL = B.getLarmorRadius(p);

	// There is only a drift effect close to the current sheet
	if (d < 2.*rL)
	{
		double x = d / rL;
		return B.polarity*p.charge*B.Gamma / sqrt(B.onePlusGamma2) *
				(0.457 - 0.412*x + 0.0915*x*x) * p.velocity;
	}


	return 0.;
}

double JokipiiKopriva::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	return 2./p.r*(kappa_rr + kappa_parallel*B.Gamma*B.Gamma/pow(B.onePlusGamma2,2)*(perPart-1.));
}

double JokipiiKopriva::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	return (-2.0*p.cosTheta)/(p.r*p.r)*kappa_perpendicular;
}
