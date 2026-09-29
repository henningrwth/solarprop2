#include "custom.h"

#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/// \file custom.cc
/// \brief Member functions of the Custom model should be implemented here.

/**
 * Constructor for custom model implementation.
 *
 * All option extraction should be done here.
 *
 * @param modelInfo %Model information.
 */
Custom::Custom(const ModelInformation& modelInfo) :
	Model(modelInfo)
{
	// Keep this to make sure model is registered in model factory:
	FACTORY_INIT;

	// extract relevant parameter values from modelInfo and store them in class members...
	// ...

	// set maximum time step and extent of the heliosphere
	dt_max = modelInfo.getFloatOptionWithDefault("dt", 500.) * second;
	heliosphereBoundary = 100. * AU;
}

double Custom::solarWindSpeed(const Particle&) const
{
	// return solar wind speed at particle position

	static const double V = 400 * km / s;
	return V;
}

int Custom::getHeaviside(const Particle& p) const
{
	// return value of Heaviside function H for given particle (H=+1 above HCS, H=-1 below HCS)
	return 0;
}

void Custom::calculate(const Particle& p)
{
	// calculate cached values that are needed several times in a particle step

	// kappa_parallel = kappa0 * ... ;
	// kappa_perpendicular_r = ... ;
	// kappa_perpendicular_theta = ... ;

	// kappa_rr = ... ;
	// kappa_thetatheta = ... ;
}
