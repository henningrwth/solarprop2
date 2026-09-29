#ifndef STANDARD2D_H_
#define STANDARD2D_H_

#include "model.h"
#include "parkerspiral.h"
#include "particle.h"
#include "sde2d.h"

/** Implementation of the standard2D model introduced in vanilla SOLARPROP.
 *
 * \attention The model is retained here for backwards compatibility, but using it is strongly discouraged:
 * The diffusion coefficient is arbitrarily rescaled by a factor of 0.07 for the case of \f$ qA>0 \f$.
 * This mimicks some charge-dependent solar modulation, but has no physics foundation.
 */
class Standard2D : public Model, RegisteredInFactory<Model, Standard2D, const ModelInformation&>
{
	public:
	Standard2D() = delete;
	Standard2D(const ModelInformation&);
	Standard2D(const Standard2D&) = default;
	~Standard2D() override {}
	std::unique_ptr<Model> clone() const override { return std::make_unique<Standard2D>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	void calculate(const Particle& p) override;
	double getB(const Particle& p) const override { return B.getB(p); }
	double getThetaPrimeHCS(const Particle& p) const override { return B.getThetaPrime(p); }

	double solarWindSpeed(const Particle& p) const override { return B.solarWindSpeed(p); }

	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRR() const override { return kappa_rr; }

	double getDriftR(const Particle& p) const override;
	double getDriftTheta(const Particle& p) const override;
	double getDriftSheetR(const Particle& p) const override;

	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override;

	/** Factory function for ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Standard2D>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "standard2D"; }

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

	/// smooth function to model HCS drift
	double f;

	/// derivative of smooth function to model HCS drift
	double fprime;

	/// Angle for drift effects in models based on model BP1989
	double alphaH;

	/// \f$ \tan\alpha_H \f$
	double tanAlphaH;

	bool dt_max_autoset = true;
};

#endif
