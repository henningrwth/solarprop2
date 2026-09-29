#ifndef HELMOD2012_H_
#define HELMOD2012_H_

#include <vector>

#include "helmod.h"

class NumericalDerivatives;
class Particle;

/**
 * Implementation the `%Helmod2012` model.
 *
 * Reference: Bobik et al., ApJ 745 (2012) 132, https://doi.org/10.1088/0004-637X/745/2/132
 *
 * See the \c jupyter notebook \c calculations/helmod2012_diffusion.ipynb for a derivation
 * of the terms involving the derivatives of the components of the diffusion coefficient in
 * heliospheric coordinates. The corresponding C++ source code is auto-generated in this file.
 * For a visualization and numerical investigation of the equations for heliospheric
 * current sheet drift, see \c calculations/helmod_hcs_drift.ipynb.
 *
 */
class Helmod2012 : public Helmod, RegisteredInFactory<Model, Helmod2012, const ModelInformation&>
{
	public:
	Helmod2012() = delete;
	Helmod2012(const ModelInformation&);
	Helmod2012(const Helmod2012&) = default;
	~Helmod2012() override {};
	std::unique_ptr<Model> clone() const override { return std::make_unique<Helmod2012>(*this); }

	void calculate(const Particle& particle) override;

	/** 1/r^2 d/dr (r^2 kappa_rr) */
	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;

	/** 1/(r^2) d/dcostheta (kappa_thetatheta*sin^2(theta)) */
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override;

	/** -d/dcostheta (kappa_r_theta*sin(theta)/r) */
	double getLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const override;

	/** -1/r^2 d/dr (r*kappa_r_theta*sin(theta)) */
	double getLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const override;

	/** Factory function for ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Helmod2012>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "Helmod2012"; }

	private:

	static double defaultKappaZeroFromSolarPolarityAndPhase(int polarity, int phase, double ssn);

	/// factor in front of diffusion coefficient: \f$ f = \beta \cdot K_0 \cdot K_P(P,t) \cdot B_\oplus / 3 \f$.
	double fdiffusion;
};

#endif
