#ifndef DENSITY_H_
#define DENSITY_H_

#include "properties.h"

#include <string>
#include <vector>

class Model;
class Particle;
class Random;

/// \file density.h

/// \typedef std::vector<std::vector<Particle> > ParticleHistories
/// \brief Store intermediate steps of pseudo-particles.
/// To access the information:
/// - First index: Pseudo-particle.
/// - Second index: Step of given pseudo-particle, recorded during tracking.
using ParticleHistories = std::vector<std::vector<Particle>>;

/**
 * Perform %SDE integration for collection of test particles.
 */
class Density
{
	public:
	Density(const Model*, const ParticleProperties&, const std::vector<double>& kineticEnergies, unsigned int _particlesPerEnergyBin, bool verbose);
	Density(const Model*, const ParticleProperties&, unsigned int _particlesPerEnergyBin,
		const std::vector<EnergySmearingProperties>& kineticEnergySmearingProps, bool verbose);
	~Density();

	ParticleHistories simulate(bool dumpIntermediateSteps, bool dumpFinal, int historyLevel, int globalSeed);

	const Particle& getParticle(unsigned int energyBin, unsigned int particleIndexInEnergyBin) const;

	/** Get cosmic-ray particle properties. */
	ParticleProperties getProperties() const { return crProps; }

	/** Get number of pseudo-particles used for each bin of initial kinetic energies. */
	unsigned int numberOfParticlesPerEnergyBin() const { return particlesPerEnergyBin; }

	private:
	double smearEnergy(const EnergySmearingProperties& p, Random& rndm);

	private:

	/// List of pseudo-particles.
	std::vector<Particle> particles;

	/// Cosmic-ray particle for which this Density object is valid.
	const ParticleProperties crProps;

	/// The model governing the modulation process. Each thread will operate on a copy of this template.
	const Model* modelTemplate;

	/// Number of pseudo-particles used for each bin of initial kinetic energies.
	unsigned int particlesPerEnergyBin;
};

#endif
