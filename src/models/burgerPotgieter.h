#ifndef BURGERPOTGIETER_H_
#define BURGERPOTGIETER_H_

#include "model.h"
#include "parkerspiral.h"
#include "particle.h"
#include "sde2d.h"

/**
 * Implementation of the model of Burger and Potgieter (1989).
 *
 * Reference:
 * R.A. Burger and M.S. Potgieter, ApJ 339 (1989) 501-511, https://adsabs.harvard.edu/full/1989apj...339..501b
 *
 * To be more precise, this model aims at reproducing Fig. 5 in the Burger&Potgieter paper.
 * The figure was done with the same parameters
 * as used in Kota&Jokipii, ApJ 265 (1983) 573, but with the new implementation
 * of HCS drift by Burger&Potgieter, see also Burger&Hattingh, Astrophys. Space Sci. 230 (1995) 375-382.
 *
 */
class BurgerPotgieter : public Model, RegisteredInFactory<Model, BurgerPotgieter, const ModelInformation&>
{
	public:
	BurgerPotgieter() = delete;
	BurgerPotgieter(const ModelInformation&);
	BurgerPotgieter(const BurgerPotgieter&) = default;
	~BurgerPotgieter() override {}
	std::unique_ptr<Model> clone() const override { return std::make_unique<BurgerPotgieter>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	void calculate(const Particle& particle) override;

	double solarWindSpeed(const Particle& p) const override { return B.solarWindSpeed(p); }
	double getB(const Particle& p) const override { return B.getB(p); }
	double getThetaPrimeHCS(const Particle& p) const override { return B.getThetaPrime(p); }

	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRR() const override { return kappa_rr; }

	double getDriftR(const Particle& p) const override;
	double getDriftTheta(const Particle& p) const override;
	double getDriftSheetR(const Particle& p) const override;

	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override;

	/** Factory function for ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<BurgerPotgieter>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "BP1989"; }

	private:

	/// Magnetic field
	ParkerSpiral B;

	/// Coefficient between the parallel and perpendicular diffusion tensor
	double perPart;

	/// Normalization of the diffusion tensor
	double kappa0;

	/// Parallel component of the diffusion tensor (in field-aligned coordinates)
	double kappa_parallel;

	/// Perpendicular component of the diffusion tensor, in radial and polar direction (in field-aligned coordinates)
	double kappa_perpendicular;

	/// Diffusion coefficient \f$ \kappa_{rr} \f$ (heliospheric coordinates)
	double kappa_rr;

	/// Diffusion coefficient \f$ \kappa_{\theta\theta} \f$ (heliospheric coordinates)
	double kappa_thetatheta;
};

#endif
