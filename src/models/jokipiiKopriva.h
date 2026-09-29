#ifndef JokipiiKopriva_H_
#define JokipiiKopriva_H_

#include "model.h"
#include "parkerspiral.h"
#include "sde2d.h"

#include <memory>

class NumericalDerivatives;
class Particle;

/**
 * Implementation for the model of Jokipii and Kopriva (1979).
 *
 * Reference:
 * J.R. Jokipii and D.A. Kopriva, ApJ 234 (1979) 384-392, https://adsabs.harvard.edu/pdf/1979ApJ...234..384J
 *
 * This model aims at reproducing the relevant curves in Fig. 2 of the reference paper.
 *
 */
class JokipiiKopriva : public Model, RegisteredInFactory<Model, JokipiiKopriva, const ModelInformation&>
{
	public:
	JokipiiKopriva() = delete;
	JokipiiKopriva(const ModelInformation&);
	JokipiiKopriva(const JokipiiKopriva&) = default;
	~JokipiiKopriva() override {};
	std::unique_ptr<Model> clone() const override { return std::make_unique<JokipiiKopriva>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	void calculate(const Particle& particle) override;

	double solarWindSpeed(const Particle& p) const override { return B.solarWindSpeed(p); }

	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRR() const override { return kappa_rr; }

	double getDriftR(const Particle& p) const override;
	double getDriftTheta(const Particle& p) const override;
	double getDriftSheetR(const Particle& p) const override;

	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override;

	// Functions for numerical drift calculations (only needed if forceNumericalDerivatives option is used).
	double getBOverB2_r(const Particle& p) const override { return B.getBOverB2_r(p); }
	double getBOverB2_theta(const Particle& p) const override { return B.getBOverB2_theta(p); }
	double getBOverB2_phi(const Particle& p) const override { return B.getBOverB2_phi(p); }

	double getB(const Particle& p) const override { return B.getB(p); }
	double getThetaPrimeHCS(const Particle& p) const override { return B.getThetaPrime(p); }
	int getHeaviside(const Particle& p) const override { return B.getHeaviside(p); }

	/** Factory function for ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<JokipiiKopriva>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "JK1979"; }

	private:

	/// Magnetic field
	ParkerSpiral B;

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

	/// Always use numerical calculation of derivatives? Useful for debugging and testing.
	bool forceNumericalDerivatives = false;
	/// NumericalDerivatives used for numerical calculation of derivatives.
	std::shared_ptr<NumericalDerivatives> nd = nullptr;
};

#endif
