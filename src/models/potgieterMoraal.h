#ifndef POTGIETERMORAAL_H_
#define POTGIETERMORAAL_H_

#include "model.h"
#include "parkerspiral.h"
#include "particle.h"
#include "sde2d.h"

class NumericalDerivatives;

/**
 * Implementation for the model of Potgieter and Moraal (1985).
 *
 * Reference: M.S. Potgieter and H. Moraal, ApJ 294 (1985) 425-440, https://adsabs.harvard.edu/abs/1985ApJ...294..425P
 *
 * This model aims at reproducing the main curves in Fig. 6 of the reference paper.
 *
 */
class PotgieterMoraal : public Model, RegisteredInFactory<Model, PotgieterMoraal, const ModelInformation&>
{
	public:
	PotgieterMoraal() = delete;
	PotgieterMoraal(const ModelInformation&);
	PotgieterMoraal(const PotgieterMoraal&) = default;
	~PotgieterMoraal() override {};
	std::unique_ptr<Model> clone() const override { return std::make_unique<PotgieterMoraal>(*this); }
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
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<PotgieterMoraal>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "PM1985"; }

	private:

	/// Magnetic field
	ParkerSpiral B;

	/// smooth function to model HCS drift
	double f;

	/// derivative of smooth function to model HCS drift
	double fprime;

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

	/// Coefficient between the parallel and perpendicular diffusion tensor
	double perPart;

	/// Angle for drift effects in models, based on model BP1989.
	double alphaH;

	/// tan(alphaH)
	double tanAlphaH;

	/// Always use numerical calculation of derivatives? Useful for debugging and testing.
	bool forceNumericalDerivatives = false;
	/// NumericalDerivatives used for numerical calculation of derivatives.
	std::shared_ptr<NumericalDerivatives> nd = nullptr;
};

#endif
