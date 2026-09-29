#include "density.h"

#include "model.h"
#include "particle.h"
#include "random.h"
#include "sde.h"
#include "units.h"

#include <iostream>

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include <omp.h>

using namespace Units;
using namespace std::string_literals;

/**
 * Constructor without energy smearing.
 *
 * The inital pseudo-particles will have discrete energy values as dictated by the \p kineticEnergies vector.
 *
 * @param _modelTemplate The model governing the modulation process. Each thread will operate on a copy of this template.
 * @param prop Mass and charge of the cosmic-ray particle.
 * @param kineticEnergies Vector of initial kinetic energies.
 * @param _particlesPerEnergyBin Number of test particles to be tracked per energy bin.
 * @param verbose Print the number of particles per bin and the overall range of the energy binning?
 */
Density::Density(const Model* _modelTemplate, const ParticleProperties& prop, const std::vector<double>& kineticEnergies, unsigned int _particlesPerEnergyBin,
		 bool verbose) :
	crProps(prop),
	modelTemplate(_modelTemplate),
	particlesPerEnergyBin(_particlesPerEnergyBin)
{
	if (verbose)
	{
#ifdef _OPENMP
		const int nThreadsMax = omp_get_max_threads();
#else
		const int nThreadsMax = 1;
#endif
		if (kineticEnergies.size() > 1)
			std::cout << "Simulate " << particlesPerEnergyBin << " particles per energy bin in " << kineticEnergies.size() << " bins from "
				  << kineticEnergies.front() / GeV << " to " << kineticEnergies.back() / GeV << " GeV using "
				  << nThreadsMax << " threads..." << std::endl;
		else if (kineticEnergies.size() == 1)
			std::cout << "Simulate " << particlesPerEnergyBin << " particles at " << kineticEnergies.front() / GeV << " GeV using "
				  << nThreadsMax << " threads..." << std::endl;
	}

	particles.reserve(kineticEnergies.size() * particlesPerEnergyBin);

	for (unsigned int ebin = 0; ebin < kineticEnergies.size(); ++ebin)
	{
		for (unsigned int i = 0; i < particlesPerEnergyBin; ++i)
		{
			int particleNumber = ebin * particlesPerEnergyBin + i;
			particles.emplace_back(prop, kineticEnergies[ebin], particleNumber);
		}
	}
}

/**
 * @brief Constructor with energy smearing.
 *
 * The initial energies of the pseudo-particles will be smeared assuming a power-law distribution spanning the respective energy bin.
 *
 * @param _modelTemplate The model governing the modulation process. Each thread will operate on a copy of this template.
 * @param prop Mass and charge of the cosmic-ray particle.
 * @param _particlesPerEnergyBin Number of test particles to be tracked per energy bin.
 * @param kineticEnergySmearingProps Vector with range and spectral index for energy smearing, one entry per energy bin.
 * @param verbose Print the number of particles per bin and the overall range of the energy binning?
 */
Density::Density(const Model* _modelTemplate, const ParticleProperties& prop, unsigned int _particlesPerEnergyBin,
		 const std::vector<EnergySmearingProperties>& kineticEnergySmearingProps, bool verbose) :
	crProps(prop),
	modelTemplate(_modelTemplate),
	particlesPerEnergyBin(_particlesPerEnergyBin)
{
	if (verbose)
	{
#ifdef _OPENMP
		const int nThreadsMax = omp_get_max_threads();
#else
		const int nThreadsMax = 1;
#endif
		if (kineticEnergySmearingProps.size())
			std::cout << "Simulate " << particlesPerEnergyBin << " particles per energy bin in " << kineticEnergySmearingProps.size()
				  << " bins from " << kineticEnergySmearingProps.front().T_low / GeV << " to " << kineticEnergySmearingProps.back().T_up / GeV
				  << " GeV using " << nThreadsMax << " threads ..." << std::endl;
	}

	particles.reserve(kineticEnergySmearingProps.size() * particlesPerEnergyBin);

	for (unsigned int ebin = 0; ebin < kineticEnergySmearingProps.size(); ++ebin)
	{
		const EnergySmearingProperties& sm = kineticEnergySmearingProps[ebin];
		Random rndmEnergy(sm.randomSeed);

		for (unsigned int i = 0; i < particlesPerEnergyBin; ++i)
		{
			double T = smearEnergy(sm, rndmEnergy);
			int particleNumber = ebin * particlesPerEnergyBin + i;
			particles.emplace_back(prop, T, particleNumber);
		}
	}
}

/**
 * Destructor.
 *
 * Clear the particle vector.
 */
Density::~Density()
{
	particles.clear();
}

/**
 * Getter for a specific particle.
 *
 * \param energyBin Energy bin (initial energies) of the particle.
 * \param particleIndexInEnergyBin Index of the particle inside the energy bin specified by \p energyBin.
 *
 * \return A \c const reference to the Particle object.
 */
const Particle& Density::getParticle(unsigned int energyBin, unsigned int particleIndexInEnergyBin) const
{

	unsigned int index = energyBin * particlesPerEnergyBin + particleIndexInEnergyBin;

	if (particleIndexInEnergyBin >= particlesPerEnergyBin || index >= particles.size())
		throw std::runtime_error("Illegal particle index ("s + std::to_string(energyBin) + ","s + std::to_string(particleIndexInEnergyBin) + ")!"s);

	return particles[index];
}

/**
 *  Simulate trajectory for all particles from the earth to the heliospheric boundary.
 *
 *  If \p historyLevel is true, returns the particle histories during %SDE propagation
 *  (a vector of a vector of particles:
 *  first index is particle number, second index is step number).
 *  For a \p historyLevel of 1, only the endpoints are stored, for a \p historyLevel of 2,
 *  all steps are stored.
 *
 *  \param dumpIntermediateSteps Print all intermediate steps?
 *  \param dumpFinal Print the final step?
 *  \param historyLevel See the description above.
 *  \param globalSeed The global seed. The (unique) particle number will be added to this seed to create a unique random generator for each particle.
 *
 *  \returns A ParticleHistories object.
 */
ParticleHistories Density::simulate(bool dumpIntermediateSteps, bool dumpFinal, int historyLevel, int globalSeed)
{
	unsigned int size = particles.size();
#ifdef _OPENMP
	const int nThreadsMax = omp_get_max_threads();
#else
	const int nThreadsMax = 1;
#endif
	std::vector<std::unique_ptr<Model>> models;
	models.reserve(nThreadsMax);
	std::vector<std::unique_ptr<SDE>> sdes;
	sdes.reserve(nThreadsMax);

	for (int iThread = 0; iThread < nThreadsMax; ++iThread)
	{
		models.push_back(modelTemplate->clone());
		sdes.push_back(models[iThread]->createSDE(globalSeed));
	}

	#pragma omp parallel for schedule(dynamic, 2)
	for (unsigned int i = 0; i < size; ++i)
	{
#ifdef _OPENMP
		auto t = omp_get_thread_num();
		const std::unique_ptr<SDE>& sde = sdes[t];
#else
		const std::unique_ptr<SDE>& sde = sdes[0];
#endif
		sde->simulate(particles[i], dumpIntermediateSteps, dumpFinal, historyLevel);

#ifdef PYTHON_BINDINGS
		// let Ctrl+C interrupt modulation run
		if (t == 0)
		{
			if (PyErr_CheckSignals() != 0)
				throw py::error_already_set();
		}
#endif
	}

	ParticleHistories histories(historyLevel ? size : 0);

	for (int iThread = 0; iThread < nThreadsMax; ++iThread)
	{
		if (historyLevel)
		{
			const std::vector<Particle>& history = sdes[iThread]->getHistory();
			for (const Particle& p : history)
				histories.at(p.getParticleNumber()).emplace_back(p);
		}
	}

	return histories;
}

/**
 * Calculate a smeared initial particle energy.
 *
 * @param p EnergySmearingProperties object containing the relevant bin ranges and spectral index.
 * @param rndm Random generator object to be used for smearing.
 *
 * @return Smeared energy.
 */
double Density::smearEnergy(const EnergySmearingProperties& p, Random& rndm)
{
	double r = rndm.getUniformRandomNumber();

	if (std::abs(p.spectralIndex + 1.0) > 0.001)
	{
		// standard case
		double a = std::pow(p.T_low, p.spectralIndex + 1.0);
		double b = std::pow(p.T_up, p.spectralIndex + 1.0);
		return std::pow((b - a) * r + a, 1.0 / (p.spectralIndex + 1.0));
	}
	else
	{
		// E^(-1) special case
		return p.T_low * std::pow(p.T_up / p.T_low, r);
	}
}
