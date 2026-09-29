#include "burger2012.h"
#include "numdiff.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor.
 *
 * @param modelInfo %Model information
 */
Burger2012::Burger2012(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo, false) // we take care of the waviness of the HCS ourselves
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;
	kappa_rphi = 0.0;
	kappa_phiphi = 0.0;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 1000., false) * s;

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

	kappa_exponent = modelInfo.getFloatOptionWithDefault("kappaExponent", 0.5, false);

	perPart = modelInfo.getFloatOptionWithDefault("rho", 0.05, false);

	heliosphereBoundary = modelInfo.getFloatOptionWithDefault("heliosphereBoundary", 10., false) * AU;

	forceNumericalDerivatives = modelInfo.getBooleanOptionWithDefault("forceNumericalDerivatives", false, false);

	std::cout << "[Burger2012] kappa0= " << kappa0/(cm2/s) << " cm^2/s, "
		  << "a= " << kappa_exponent << ", "
		  << "lambda_perp/lambda_par= " << perPart << ", R= " << heliosphereBoundary/AU << " AU"
		  << std::endl;

}

/**
 * Calculate the quantities needed several times in one particle step.
 */
void Burger2012::calculate(const Particle& p)
{
	if (!nd && forceNumericalDerivatives)
	{
		nd = std::make_shared<NumericalDerivatives>(this);
	}

	B.calculate(p);

	double rigidityValue = p.absRigidity / (GV/c);
	double kappa_R = std::pow(rigidityValue, kappa_exponent);

	kappa_parallel = kappa0 * kappa_R * p.beta * B.getBearthOverB(p);
	kappa_perpendicular = perPart * kappa_parallel;

	double Gamma2 = B.onePlusGamma2 - 1.0;
	kappa_rr = (kappa_parallel + kappa_perpendicular*Gamma2) / B.onePlusGamma2;
	kappa_thetatheta = kappa_perpendicular;

	kappa_rphi = (kappa_perpendicular - kappa_parallel) * B.Gamma / B.onePlusGamma2;
	kappa_phiphi = (kappa_parallel * Gamma2 + kappa_perpendicular) / B.onePlusGamma2;

	fdiffusion = kappa0 * std::sqrt(rigidityValue) * p.beta * B.Bearth / (B.BfieldNormalization * B.r0_sqr);

	double kappa_a = p.absRigidity * p.velocity / 3.0 / B.getB(p);

	parameter_k = rigidityValue > 3.5 ? 27.52 * std::pow(rigidityValue, -0.25) : 20.12;

	double tanalpha = std::tan(B.tiltAngle);
	double t1 = p.getR() * B.OmegaSun/B.solarWindSpeed(p);
	double phi_star = p.getPhi() + t1;
	double theta_ns = pihalf - std::atan(tanalpha * std::sin(phi_star));

	double n1 = -tanalpha * std::cos(phi_star);
	double d1 = 1.0 + std::pow(tanalpha * std::sin(phi_star), 2);
	double tanbeta = t1 * n1 / d1;
	double tangamma = n1 / d1 / p.getSinTheta();

	cos_nu = 1.0 / std::sqrt(1. + tanbeta*tanbeta + tangamma*tangamma);
	double argument = parameter_k * (theta_ns - p.getTheta()) * cos_nu;
	drift_coeff = std::tanh(argument);

	// cancel the Heaviside factor included in the ParkerSpiral class here,
	// since the tanh(...) term is supposed to take care of the sign change across the HCS
	drift_coeff *= B.getHeaviside(p);

	v_dns = B.polarity * p.charge * parameter_k * kappa_a / p.r / std::pow(std::cosh(argument), 2);
}

double Burger2012::getThetaPrimeHCS(const Particle& p) const
{
	double t1 = p.getR() * B.OmegaSun/B.solarWindSpeed(p);
	double phi_star = p.getPhi() + t1;

	return pihalf - std::atan(std::tan(B.tiltAngle) * std::sin(phi_star));
}

double Burger2012::getDriftR(const Particle& p) const
{
	return drift_coeff * B.getDriftVelocity_r(p);
}

double Burger2012::getDriftTheta(const Particle& p) const
{
	return drift_coeff * B.getDriftVelocity_theta(p);
}

double Burger2012::getDriftPhi(const Particle& p) const
{
	return drift_coeff * B.getDriftVelocity_phi(p);
}

double Burger2012::getDriftSheetR(const Particle& p) const
{
	double sinPsi = B.getSinPsi(p);
	return v_dns * sinPsi * cos_nu;
}

double Burger2012::getDriftSheetTheta(const Particle& p) const
{
	double sign = -std::copysign(1.0, std::cos(p.getPhi() + B.OmegaSun/B.solarWindSpeed(p) * p.r));
	double sin_nu = std::sqrt(1.0 - std::pow(cos_nu, 2));
	return sign * v_dns * sin_nu;
}

double Burger2012::getDriftSheetPhi(const Particle& p) const
{
	double cosPsi = B.getCosPsi(p);
	return v_dns * cosPsi * cos_nu;
}

double Burger2012::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	double Gamma2 = B.onePlusGamma2 - 1.0;

	return fdiffusion*p.r*(6*B.onePlusGamma2*Gamma2*perPart + 4*B.onePlusGamma2 - 3*Gamma2*Gamma2*perPart - 3*Gamma2) / pow(B.onePlusGamma2, 2.5);
}

double Burger2012::getLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRPhi(p);

	return 0.0;
}

double Burger2012::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	return -fdiffusion*p.cosTheta*perPart*(B.onePlusGamma2 + 1) / pow(B.onePlusGamma2, 1.5);
}

double Burger2012::getLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearAzimuthalDiffusionTermFromKappaRPhi(p);

	return fdiffusion*B.Gamma*(B.onePlusGamma2 + 3)*(perPart - 1) / (pow(B.onePlusGamma2, 2.5)*p.sinTheta);
}

double Burger2012::getLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearAzimuthalDiffusionTermFromKappaPhiPhi(p);

	return 0.0;
}
