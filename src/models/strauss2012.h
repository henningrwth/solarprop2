#ifndef STRAUSS2012_H_
#define STRAUSS2012_H_

#include "model.h"
#include "parkerspiralhcs.h"
#include "particle.h"

class NumericalDerivatives;

/**
 * Implementation for the `%Strauss2012` model.
 *
 * Reference:
 * R. D. Strauss et al., Astrophys. Space Sci. (2012) 339:223-236, https://doi.org/10.1007/s10509-012-1003-z
 *
 * This model aims at reproducing the results in Figs. 5, 11, and 15 of the reference paper.
 *
 * For the calculation of the diffusion coefficients in heliospheric coordinates and their relevant derivatives,
 * see the \c jupyter notebook \c calculations/parker_diffusion.ipynb. For the calculation of the gradient and
 * curvature drift velocities, see \c calculations/parker_drift.ipynb. For a verification of key expressions used
 * in the reference paper, see \c calculations/strauss2012_hcs_drift.ipynb.
 *
 */
class Strauss2012 : public Model, RegisteredInFactory<Model, Strauss2012, const ModelInformation&>
{
	public:
	Strauss2012() = delete;
	Strauss2012(const ModelInformation&);
	Strauss2012(const Strauss2012&) = default;
	~Strauss2012() override {};

	std::unique_ptr<Model> clone() const override { return std::make_unique<Strauss2012>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override;

	void calculate(const Particle& p) override;

	double solarWindSpeed(const Particle& p) const override { return B.solarWindSpeed(p); }
	double getB(const Particle& p) const override { return B.getB(p); }
	double getThetaPrimeHCS(const Particle& p) const override { return B.getThetaPrime(p); }

	double getKappaRR() const override { return kappa_rr; }
	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRPhi() const override { return kappa_rphi; }
	double getKappaPhiPhi() const override { return kappa_phiphi; }

	double getDriftR(const Particle& p) const override { return B.getDriftVelocity_r(p); }
	double getDriftTheta(const Particle& p) const override { return B.getDriftVelocity_theta(p); }
	double getDriftPhi(const Particle& p) const override { return B.getDriftVelocity_phi(p); }

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

	/** %Model name used by ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Strauss2012>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "Strauss2012"; }

	private:

	/// Magnetic field
	ParkerSpiralHCS B;

	/// Reduce model to two-dimensional case, for testing?
	bool force2D = false;

	/// Mean free path scale parameter
	double lambda0;

	/// Parallel component of the diffusion tensor (in field-aligned coordinates)
	double kappa_parallel;

	/// Perpendicular component of the diffusion tensor (in field-aligned coordinates)
	double kappa_perpendicular;

	/// Coefficient between the parallel and perpendicular diffusion tensor
	double perPart;

	/// Diffusion coefficient \f$ \kappa_{rr} \f$ (heliospheric coordinates)
	double kappa_rr;

	/// Diffusion coefficient \f$ \kappa_{\theta\theta} \f$ (heliospheric coordinates)
	double kappa_thetatheta;

	/// Diffusion coefficient \f$ \kappa_{r\phi} \f$ (heliospheric coordinates)
	double kappa_rphi;

	/// Diffusion coefficient \f$ \kappa_{\phi\phi} \f$ (heliospheric coordinates)
	double kappa_phiphi;

	/// Exponent a for scaling of diffusion coefficient with rigidity: \f$ \kappa \f$ is proportional to \f$ R^a \f$, where \f$ R \f$ is rigidity
	double kappa_exponent;

	/// factor in front of diffusion coefficient: \f$ f = \lambda_0 \cdot v \cdot K_P(P) \cdot B_\oplus / (3B_0 r_0^2) \f$.
	double fdiffusion;

	/// radial component of HCS drift velocity vector
	double hcs_drift_speed_r;
	/// polar component of HCS drift velocity vector
	double hcs_drift_speed_theta;
	/// azimuthal component of HCS drift velocity vector
	double hcs_drift_speed_phi;

	/// use simplified HCS drift scheme: average velocity instead of tracing along HCS?
	bool simple_hcs_drift = false;

	/// fudge factor for HCS drift velocity scale: use 0.6 for match with original results in reference paper
	double hcs_drift_velocity_factor;

	double simpleDriftSheetBP(const Particle& p) const;
	double simpleDriftSheetWNS(const Particle& p) const;

	/// Always use numerical calculation of derivatives? Useful for debugging and testing.
	bool forceNumericalDerivatives = false;
	/// NumericalDerivatives used for numerical calculation of derivatives.
	std::shared_ptr<NumericalDerivatives> nd = nullptr;
};

#endif
