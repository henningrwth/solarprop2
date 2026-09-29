#include "potgieterMoraal.h"
#include "numdiff.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo %Model information
 */
PotgieterMoraal::PotgieterMoraal(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo, false)
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 2000.) * s;

	kappa0 = 6.0e21 * cm2/s;

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();
	if (modelInfoKappa0 > 0.0)
	{
		kappa0 = modelInfoKappa0;
	}

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0);
	kappa0 *= kappaScaling;

	perPart = 0.0333333333333;
	heliosphereBoundary = 50.*AU;

	double thetaonehalf = 85./90.*pi/2.;
	alphaH = acos(pi/(2.*thetaonehalf)-1.);
	tanAlphaH = tan(alphaH);

	f = 0.0;
	fprime = 0.0;

	forceNumericalDerivatives = modelInfo.getBooleanOptionWithDefault("forceNumericalDerivatives", false);

	if (forceNumericalDerivatives)
	{
		std::cout << "[PotgieterMoraal] Forcing numerical calculation of derivatives in diffusion terms..." << std::endl;
	}
}

/**
 * Calculate the quantities needed several times in one particle step.
 */
void PotgieterMoraal::calculate(const Particle& p)
{
	if (!nd && forceNumericalDerivatives)
	{
		nd = std::make_shared<NumericalDerivatives>(this);
	}

	B.calculate(p);

	double rigidityValue = p.absRigidity / (GV/c);
	kappa_perpendicular = perPart*kappa0*rigidityValue*p.beta*B.getBearthOverB(p);

	if (rigidityValue < 0.4)
		rigidityValue = 0.4;

	double r_AU = p.r/AU;
	double theta = p.getTheta();
	kappa_parallel = kappa0*rigidityValue*p.beta*(1. + r_AU*r_AU);

	kappa_rr = 1.0/B.onePlusGamma2 * (kappa_parallel + kappa_perpendicular*B.Gamma*B.Gamma);
	kappa_thetatheta = kappa_perpendicular;

	f =  1./alphaH*atan((1.-2.*theta/pi)*tanAlphaH);
	fprime = -2.*pi/alphaH*tanAlphaH/(pi*pi+pow((pi-2.*theta),2)*pow(tanAlphaH,2));
}

double PotgieterMoraal::getDriftR(const Particle& p) const
{
	return -B.driftConstant * p.cosTheta/p.sinTheta * f;

}

double PotgieterMoraal::getDriftTheta(const Particle& p) const
{
	return B.driftConstant * (1.0 + B.onePlusGamma2) * f;
}

double PotgieterMoraal::getDriftSheetR(const Particle& p) const
{
	return -B.driftConstant/2. * B.onePlusGamma2 * fprime;
}

double PotgieterMoraal::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	double r = p.r;
	return 2./r*(kappa_rr + (kappa_perpendicular - kappa_parallel)*B.Gamma*B.Gamma/pow(B.onePlusGamma2,2)
		     + 1.*kappa_parallel*r*r/(1.+r*r)*1./B.onePlusGamma2
		     + 1./2.*kappa_perpendicular*(1+B.onePlusGamma2)*B.Gamma*B.Gamma/pow(B.onePlusGamma2,2));
}


double PotgieterMoraal::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	return -1./(p.r*p.r)*kappa_perpendicular*p.cosTheta/B.onePlusGamma2*(1.0+B.onePlusGamma2);
}
