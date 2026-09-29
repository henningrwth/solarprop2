#include "sde.h"

#include "model.h"
#include "random.h"
#include "units.h"

#include <cassert>
#include <iomanip>
#include <iostream>

using namespace Units;

/**
 * Constructor.
 *
 * @param m The model defining the terms in the %SDE.
 * @param _globalSeed Global seed for random numbers. The unique particle number will be added to this to obtain a unique random generator for each particle.
 *
 */
SDE::SDE(Model* m, unsigned int _globalSeed) :
	model(m),
	globalSeed(_globalSeed)
{
	dt = 0.0;

	assert(model);
}

/** Destructor. */
SDE::~SDE()
{
	// do nothing else
	// model and random objects are not owned by this class
}

/**
 *  Simulate trajectory from the earth to the heliospheric boundary.
 *
 *  \param p The pseudo-particle to be tracked.
 *  \param dumpIntermediateSteps Print all intermediate steps?
 *  \param dumpFinal Print the final step?
 *  \param historyLevel History level:
 *   - 0: do not store particle history
 *   - 1: store only endpoints
 *   - 2: store all steps
 */
void SDE::simulate(Particle& p, bool dumpIntermediateSteps, bool dumpFinal, int historyLevel)
{
	// Time limit after which the propagation of a pseudo-particle should be stopped
	static const double sTimeLimit = 50000000 * s;

	// Abort pseudo-particles that come too close to the Sun
	static const double sInnerRadiusLimit = 0.01 * AU;

	// random generator for particle
	// link seed to particle number to ensure reproducibility of results across runs
	Random r(globalSeed + p.getParticleNumber());

	bool storeAllSteps = (historyLevel > 1);
	if (storeAllSteps)
		history.emplace_back(Particle(p));

	// Simulate only inside the heliosphere and below the time limit
	while (p.getR() > sInnerRadiusLimit && p.getR() < model->getHeliosphereBoundary() && p.getTime() < sTimeLimit && std::abs(p.getCosTheta()) < 1.0)
	{
		if (dumpIntermediateSteps)
		{
			p.dump();
		}

		step(p, storeAllSteps, r);
	}

	if (historyLevel == 1)
		history.emplace_back(Particle(p));

	// also dump particle after final step?
	if (dumpFinal)
	{
		p.dump();
	}
}

/**
 * Perform a single pseudo-particle step.
 *
 *  \param p The pseudo-particle being tracked.
 *  \param storeHistory Record step in \c history vector?
 *  \param rndm The random generator used for the Wiener process.
 */
void SDE::step(Particle& p, bool storeHistory, Random& rndm)
{
	// Compute cached values.
	model->calculate(p);

	// Generate random numbers.
	generateRandomNumbers(rndm);

	// temporarily store variables
	double rTemp = p.getR();
	double timeTemp = p.getTime();
	double pTemp = p.getP();
	double cosThetaTemp = p.getCosTheta();
	double phiTemp = p.getPhi();

	// the first call may set dt such that numerical issues near r=0 or near the poles are avoided
	double deltaR = getDeltaR(p);
	double deltaTime = getDeltaTime(p);
	double deltaP = getDeltaP(p);
	double deltaCosTheta = getDeltaCosTheta(p);
	double deltaPhi = getDeltaPhi(p);
	double Rh = model->getHeliosphereBoundary();

	if (deltaR == 0. && deltaCosTheta == 0. && deltaPhi == 0.)
		throw std::runtime_error("Particle is stuck!");

	// Make sure that the particle is not carried beyond the heliosphere (where solar wind velocity changes!)
	if (rTemp + deltaR > Rh)
	{
		double fraction = (Rh - rTemp) / deltaR;
		fraction *= 1.0001; // avoid numerical problems: make sure R>Rh.
		deltaR *= fraction;
		deltaTime *= fraction;
		deltaP *= fraction;
		deltaCosTheta *= fraction;
		deltaPhi *= fraction;
	}

	// Change pseudo-particle values
	rTemp += deltaR;
	timeTemp += deltaTime;
	pTemp += deltaP;
	cosThetaTemp += deltaCosTheta;
	phiTemp += deltaPhi;

	// Take care that theta is always between 0 and pi (see Strauss et al., Kopp2012)
	if (cosThetaTemp > 1.0)
	{
		cosThetaTemp = 2.0 - cosThetaTemp;
		phiTemp += pi;
	}
	else if (cosThetaTemp < -1.0)
	{
		cosThetaTemp = -2.0 - cosThetaTemp;
		phiTemp += pi;
	}

	// Take care that phi is always between 0 and 2*pi (see Strauss et al.)
	if (phiTemp < 0)
	{
		phiTemp += 2 * pi;
	}
	else if (phiTemp > 2 * pi)
	{
		phiTemp -= 2 * pi;
	}

	// Define radial boundary condition for r = 0
	if (rTemp < 0)
	{
		rTemp = -1. * rTemp;
	}

	p.setR(rTemp);
	p.setCosTheta(cosThetaTemp);
	p.setPhi(phiTemp);
	p.setMomentum(pTemp);
	p.setTime(timeTemp);

	if (p.getR() >= model->getHeliosphereBoundary())
	{
		p.setInsideHeliosphere(false);
	}

	if (storeHistory)
		history.emplace_back(Particle(p));
}

/**
 * Calculate energy gain for current step.
 *
 * \note This function is not used, in order to avoid the extra loss term that appears in the %SDEs when
 * calculating \f$ \mathrm{d}T \f$ instead of \f$ \mathrm{d}p \f$.
 */
double SDE::getDeltaT(const Particle& p) const
{
	double T = p.getT();
	double T0 = p.getMass() * c2;

	return 1.0 / 3.0 * model->getEnergyGainTerm(p) * T * (T + 2. * T0) / (T + T0) * dt;
}

/**
 * Calculate momentum gain for current step.
 */
double SDE::getDeltaP(const Particle& p) const
{
	double mom = p.getP();

	return 1.0 / 3.0 * model->getEnergyGainTerm(p) * mom * dt;
}
