#ifndef CUSTOM_H_
#define CUSTOM_H_

#include "model.h"
#include "particle.h"

// pick the relevant header here:
#include "sde2d.h"
//#include "sde3d.h"

/**
 * Template for custom model implementation.
 *
 * Add class members for the relevant options and parameters, and implement the
 * the needed functions starting with `get...`, as well as `solarWindSpeed` in \c custom.cc.
 * Set the model name in modelName(). Pick the correct set of SDEs in createSDE().
 * Values that are needed several times in a particle step should be cached in class members
 * and calculated in calculate().
 */
class Custom : public Model, RegisteredInFactory<Model, Custom, const ModelInformation&>
{
	public:
	Custom() = delete;
	Custom(const ModelInformation&);
	Custom(const Custom&) = default;
	~Custom() override {}

	/** Clone the model using the copy constructor, thus avoiding repeated extraction of options. */
	std::unique_ptr<Model> clone() const override { return std::make_unique<Custom>(*this); }

	/** Pick the correct set of SDEs (2D or 3D) here. */
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	/** Calculate cached values that are needed several times in a particle step in this function. */
	void calculate(const Particle& particle) override;

	/** Solar wind speed at particle position. */
	double solarWindSpeed(const Particle& p) const override;

	//
	// The following functions may have to be implemented to calculate the
	// components of the diffusion coefficient in heliospheric spherical coordinates:
	//

	/** Value of \f$ \kappa_{rr} \f$. */
	double getKappaRR() const override { return 0.0; }
	/** Value of \f$ \kappa_{\theta\theta} \f$. */
	double getKappaThetaTheta() const override { return 0.0; }
	/** Value of \f$ \kappa_{r\theta} \f$. */
	double getKappaRTheta() const override { return 0.0; }
	/** Value of \f$ \kappa_{\phi\phi} \f$. */
	double getKappaPhiPhi() const override { return 0.0; }
	/** Value of \f$ \kappa_{r\phi} \f$. */
	double getKappaRPhi() const override { return 0.0; }

	//
	// The following functions have to be implemented for models that provide analytic
	// expressions for the drift terms:
	//

	/** Radial component of the gradient and curvature drift velocity: \f$ V_{dr} \f$. */
	double getDriftR(const Particle& p) const override { return 0.0; }
	/** Polar component of the gradient and curvature drift velocity: \f$ V_{d\theta} \f$. */
	double getDriftTheta(const Particle& p) const override { return 0.0; }
	/** Azimuthal component of the gradient and curvature drift velocity: \f$ V_{d\phi} \f$. */
	double getDriftPhi(const Particle& p) const override { return 0.0; }

	/** Radial component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	double getDriftSheetR(const Particle& p) const override { return 0.0; }
	/** Polar component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	double getDriftSheetTheta(const Particle& p) const override { return 0.0; }
	/** Azimuthal component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	double getDriftSheetPhi(const Particle& p) const override { return 0.0; }

	/** Amplitude of magnetic field at particle position. */
	double getB(const Particle& p) const override { return 0.0; }
	/** Location \f$ \theta^\prime(r, \phi) \f$ of the heliospheric current sheet (HCS). */
	double getThetaPrimeHCS(const Particle& p) const override { return Units::pihalf; }

	//
	// The following functions have to be implemented for models that rely on numerical
	// calculation of drifts:
	//

	/** Radial component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	double getBOverB2_r(const Particle& p) const override { return 0.0; }
	/** Polar component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	double getBOverB2_theta(const Particle& p) const override { return 0.0; }
	/** Azimuthal component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	double getBOverB2_phi(const Particle& p) const override { return 0.0; }

	/** Value of Heaviside function H for given particle (H=+1 above HCS, H=-1 below HCS), for numerical drift calculation. */
	int getHeaviside(const Particle& p) const override;

	//
	// The following functions may have to be implemented for models that provide analytic
	// expressions for the derivatives involving the diffusion coefficients:
	//

	/** 1/r^2 d/dr (r^2 kappa_rr)
	 *
	 *  Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial r}(r^2\kappa_{rr}) \f$.
	 */
	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override { return 0.0; }

	/** 1/(r^2) d/dcostheta (kappa_thetatheta*sin^2(theta))
	 *
	 *  Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial\cos\theta}(\kappa_{\theta\theta}\sin^2\theta) \f$.
	 */
	double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const override { return 0.0; }

	/** -d/dcostheta (kappa_r_theta*sin(theta)/r)
	 *
	 *  Calculate \f$ -\frac{\partial}{\partial\cos\theta}\left(\kappa_{r\theta}\frac{\sin\theta}{r}\right) \f$.
	 */
	double getLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const override { return 0.0; }

	/** 1/(r*sin(theta)) * d/dphi (kappa_r_phi)
	 *
	 *  Calculate \f$ \frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}(\kappa_{r\phi}) \f$.
	 */
	double getLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const override { return 0.0; }

	/** -1/r^2 d/dr (r*kappa_r_theta*sin(theta))
	 *
	 *  Calculate \f$ -\frac{1}{r^2}\frac{\partial}{\partial r}(r\kappa_{r\theta}\sin\theta) \f$.
	 */
	double getLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const override { return 0.0; }

	/** 1/(r^2 sin(theta)) * d/dr (r*kappa_r_phi)
	 *
	 *  Calculate \f$ \frac{1}{r^2\sin\theta}\frac{\partial}{\partial r}(r\kappa_{r\phi}) \f$.
	 */
	double getLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const override { return 0.0; }

	/** 1/(r^2 sin^2(theta)) * d/dphi (kappa_phi_phi)
	 *
	 *  Calculate \f$ \frac{1}{r^2\sin^2\theta}\frac{\partial}{\partial\phi}(\kappa_{\phi\phi}) \f$.
	 */
	double getLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const override { return 0.0; }

	/** Override any automatic adjustments to the time step by returning a non-zero value here. */
	double overrideTimeStep(const Particle& p) const override { return 0.0; }

	/** Set the model name, to be used for the \c model parameter. */
	static std::string modelName() { return "custom"; }

	/** This function is called by the model factory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Custom>(modelInfo); }
};

#endif
