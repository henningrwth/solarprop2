#ifndef HELMOD2018_H_
#define HELMOD2018_H_

#include <vector>

#include "helmod.h"

class NumericalDerivatives;
class Particle;

/**
 * Implementation the `%Helmod2018` model.
 *
 * This should correspond to their code version 3.
 *
 * Reference: Boschini et al., Adv. Space Res. 62 (2018) 2859-2879, https://doi.org/10.1016/j.asr.2017.04.017
 *
 * See the \c jupyter notebook \c calculations/helmod2018_diffusion.ipynb for a derivation
 * of the terms involving the derivatives of the components of the diffusion coefficient in
 * heliospheric coordinates. The corresponding C++ source code is auto-generated in this file.
 * For a visualization and numerical investigation of the equations for heliospheric
 * current sheet drift, see \c calculations/helmod_hcs_drift.ipynb.
 *
 */
class Helmod2018 : public Helmod, RegisteredInFactory<Model, Helmod2018, const ModelInformation&>
{
	public:
	Helmod2018() = delete;
	Helmod2018(const ModelInformation&);
	Helmod2018(const Helmod2018&) = default;
	~Helmod2018() override {};
	std::unique_ptr<Model> clone() const override { return std::make_unique<Helmod2018>(*this); }

	/** Calculate all quantities which are needed several times in one step. */
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
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Helmod2018>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "Helmod2018"; }

	private:

	static double defaultKappaZeroFromSolarPolarityAndPhase(int polarity, int phase, double ssn);

	/// factor in front of diffusion coefficient: \f$ f = \beta \cdot K_0 \cdot K_P(P,t) / 3 \f$
	double fdiffusion;

	/// Constant \f$ r_e \f$ for scale of radial dependence of diffusion coefficent: distance between Earth and Sun
	double r_e;
};

#endif
