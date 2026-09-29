#ifndef SDE_3D_H_
#define SDE_3D_H_

#include "sde.h"

/**
 * Set of SDEs for simplified three-dimensional models.
 *
 * This assumes \f$ \kappa_{r\theta} = \kappa_{\theta\phi} = 0 \f$.
 *
 * See the \c jupyter notebook \c calculations/sde_derivation_3d_mu.ipynb for a derivation of the SDEs.
 *
 */
class SDE3d : public SDE
{
	public:
	SDE3d(Model* m, unsigned int _globalSeed);
	~SDE3d() override{};

	Arrows calculateArrows(const Particle& p) const override;

	private:

	void generateRandomNumbers(Random& rndm) override;

	double adjusted_dt(const Particle& p) const;

	double getDeltaR(const Particle& p) override;
	double getDeltaCosTheta(const Particle& p) const override;
	double getDeltaPhi(const Particle& p) const override;

	private:

	/// %Random number for Wiener process in radial direction.
	double w_r;

	/// %Random number for Wiener process in polar direction.
	double w_mu;

	/// %Random number for Wiener process in azimuthal direction.
	double w_phi;
};

#endif
