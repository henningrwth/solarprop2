#include "model.h"
#include "particle.h"
#include "sde.h"
#include "units.h"

using namespace Units;

/**
 * Base constructor that sets variables common to all model implementations.
 *
 * @param modelInfo %Model information.
 */
Model::Model(const ModelInformation& modelInfo)
{
	dt_max = 0;
	heliosphereBoundary = 0;

	dynamicStep = modelInfo.getBooleanOptionWithDefault("dynamicStep", false);
}

/** Calculate \f$ \frac{1}{r^2} \frac{\partial}{\partial r} (r^2 V_{sw}) \f$. */
double Model::getEnergyGainTerm(const Particle& p) const
{
	double r = p.getR();
	double V_sw = solarWindSpeed(p);

	return 2. * V_sw / r;
}
