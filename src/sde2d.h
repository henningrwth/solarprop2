#ifndef SDE_2D_H_
#define SDE_2D_H_

#include "sde.h"

/**
 * Set of SDEs for two-dimensional models.
 *
 * This assumes \f$ \phi=0 \f$, and \f$ \kappa_{r\phi} = \kappa_{\theta\phi} = \kappa_{\phi\phi} = 0 \f$.
 *
 * Only the \f$ r  \f$ and \f$ \mu=\cos\theta \f$ coordinates are considered. The full set of equations
 * for \f$ \kappa_{r\theta} \neq 0 \f$ is implemented.
 *
 * See the \c jupyter notebook \c calculations/sde_derivation_2d.ipynb for a derivation of the SDEs.
 *
 */
class SDE2d : public SDE
{
	public:
	SDE2d(Model* m, unsigned int _globalSeed);
	~SDE2d() override{};

	Arrows calculateArrows(const Particle& p) const override;

	private:

	void generateRandomNumbers(Random& rndm) override;

	double adjusted_dt(const Particle& p, double v, double driftR, double driftsheetR, double linearDiffusionTermRR, double kappaRR,
			   double linearCoefficient) const;

	double getDeltaR(const Particle& p) override;
	double getDeltaCosTheta(const Particle& p) const override;

	private:

	/// %Random number for Wiener process in radial direction.
	double wr;

	/// %Random number for Wiener process in polar direction.
	double wx;
};

#endif
