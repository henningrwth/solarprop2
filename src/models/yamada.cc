#include "yamada.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo Model information
 */
Yamada::Yamada(const ModelInformation& modelInfo) :
	Model(modelInfo)
{
	FACTORY_INIT;

	kappa_rr = 0.0;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 15000.) * s;

	heliosphereBoundary = 100.*AU;

	// Force-field approximation: phi = V(R-1)/3kappa
	kappa0 = std::abs(Yamada::solarWindSpeed(Particle())) * (heliosphereBoundary - 1.0*AU) / 3.0 / modelInfo.getModulationPotential();

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();
	if (modelInfoKappa0 > 0.0)
	{
		kappa0 = modelInfoKappa0;
	}

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0);
	kappa0 *= kappaScaling;

	std::cout << "[Yamada] Diffusion constant (cm^2 s^-1 GV^-1): " << kappa0 / (cm2/s/GV) << std::endl;
}

/** Solar wind speed. */
double Yamada::solarWindSpeed(const Particle&) const
{
	static const double V = 400 * km/s;
	return V;
}

void Yamada::calculate(const Particle& p)
{
	kappa_rr = kappa0*p.velocity*p.absRigidity;
}

double Yamada::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	return 2.*kappa_rr/p.r;
}
