#ifndef SDE_H_
#define SDE_H_

#include "arrows.h"
#include "particle.h"

#include <vector>

class Model;
class Random;

/**
 *  Abstract base class for SDEs.
 *
 *  There are two derived classes:
 *  - SDE2d for two-dimensional integration in \f$ (r, \theta) \f$, including cases where \f$ \kappa_{r\theta}\neq 0 \f$.
 *  - SDE3d for three-dimensional integration, assuming no polar diffusion, \f$ \kappa_{r\theta}=\kappa_{\theta{}r}=\kappa_{\theta\phi}=\kappa_{\phi\theta}=0 \f$.
 */
class SDE
{
	public:
	SDE(Model* m, unsigned int _globalSeed);
	virtual ~SDE();

	virtual void simulate(Particle& p, bool dumpIntermediateSteps, bool dumpFinal, int historyLevel);
	virtual Arrows calculateArrows(const Particle& p) const = 0;

	/**
	 *  Getter for particle history (if stored).
	 *
	 *  \note This will contain steps for different particles (which can be disentangled by their particle number).
	 */
	std::vector<Particle> getHistory() const { return history; }

	protected:

	void step(Particle& p, bool storeHistory, Random& rndm);

	/** Generate and store random numbers needed for the SDE step. */
	virtual void generateRandomNumbers(Random& rndm) = 0;

	virtual double getDeltaR(const Particle& p) { return 0.0; }
	virtual double getDeltaCosTheta(const Particle& p) const { return 0.0; }
	virtual double getDeltaPhi(const Particle&) const { return 0.0; }

	virtual double getDeltaT(const Particle& p) const;
	virtual double getDeltaP(const Particle& p) const;

	/** Get the time step. */
	virtual double getDeltaTime(const Particle&) const { return dt; }

	protected:

	/// Adjusted dt to avoid numerical instability around \f$ r = 0 \f$.
	double dt;

	/// The model defining the terms in the %SDE.
	Model* model = nullptr;

	/// Global seed for random numbers. The unique particle number will be added to this to obtain a unique random generator for each particle.
	unsigned int globalSeed;

	/// Pseudo-particle steps are stored here if the relevant options are used.
	std::vector<Particle> history;
};

#endif
