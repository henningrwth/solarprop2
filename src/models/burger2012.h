#ifndef BURGER_2012_H_
#define BURGER_2012_H_

#include "model.h"
#include "parkerspiral.h"
#include "particle.h"
#include "sde3d.h"

class NumericalDerivatives;

/**
 * Implementation of the `%Burger2012` model.
 *
 * Reference: R.A. Burger, ApJ 760:60 (2012), https://doi.org/10.1088/0004-637x/760/1/60
 *
 *
 * For the calculation of the diffusion coefficients in heliospheric coordinates and their relevant derivatives,
 * see the \c jupyter notebook \c calculations/parker_diffusion.ipynb. For the calculation of the gradient and
 * curvature drift velocities, see \c calculations/parker_drift.ipynb.
 *
 */
class Burger2012 : public Model, RegisteredInFactory<Model, Burger2012, const ModelInformation&>
{
	public:
	Burger2012() = delete;
	Burger2012(const ModelInformation&);
	Burger2012(const Burger2012&) = default;
	~Burger2012() override {};
	std::unique_ptr<Model> clone() const override { return std::make_unique<Burger2012>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE3d>(this, globalSeed); }

	void calculate(const Particle& particle) override;

	double solarWindSpeed(const Particle& p) const override { return B.solarWindSpeed(p); }
	double getB(const Particle& p) const override { return B.getB(p); }
	double getThetaPrimeHCS(const Particle& p) const override;

	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRR() const override { return kappa_rr; }
	double getKappaRPhi() const override { return kappa_rphi; }
	double getKappaPhiPhi() const override { return kappa_phiphi; }

	double getDriftR(const Particle& p) const override;
	double getDriftTheta(const Particle& p) const override;
	double getDriftPhi(const Particle& p) const override;
	double getDriftSheetR(const Particle& p) const override;
	double getDriftSheetTheta(const Particle& p) const override;
	double getDriftSheetPhi(const Particle& p) const override;

	/** 1/r^2 d/dr (r^2 kappa_rr) */
	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;

	/** 1/(r*sin(theta)) * d/dphi (kappa_r_phi) */
	double getLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const override;

	/** 1/(r^2) d/dcostheta (kappa_thetatheta*sin^2(theta)) */
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override;

	/** 1/(r^2 sin(theta)) * d/dr (r*kappa_r_phi) */
	double getLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const override;

	/** 1/(r^2 sin^2(theta)) * d/dphi (kappa_phi_phi) */
	double getLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const override;

	/** Factory function for ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Burger2012>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "Burger2012"; }

	private:

	/// Magnetic field.
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

	/// Diffusion coefficient \f$ \kappa_{r\phi} \f$ (heliospheric coordinates)
	double kappa_rphi;

	/// Diffusion coefficient \f$ \kappa_{\phi\phi} \f$ (heliospheric coordinates)
	double kappa_phiphi;

	/// Exponent a for scaling of diffusion coefficient with rigidity: \f$ \kappa \f$ is proportional to \f$ R^a \f$, where \f$ R \f$ is rigidity.
	double kappa_exponent;

	/// factor in front of diffusion coefficient: \f$ f = \kappa_0 \cdot K_P(P) \cdot B_\oplus / (B_0  r_0^2) \f$.
	double fdiffusion;

	/// Drift coeffienct \f$ \tanh(k(\theta^\prime-\theta)\cos\nu) \f$.
	double drift_coeff;

	/// Parameter \f$ k \f$.
	double parameter_k;

	/// \f$ \cos\nu \f$, see equation (11) of the reference paper.
	double cos_nu;

	/// Amplitude of current sheet drift velocity vector, see equation (16) of the reference paper.
	double v_dns;

	/// Always use numerical calculation of derivatives? Useful for debugging and testing.
	bool forceNumericalDerivatives = false;
	/// NumericalDerivatives used for numerical calculation of derivatives.
	std::shared_ptr<NumericalDerivatives> nd = nullptr;
};

#endif
