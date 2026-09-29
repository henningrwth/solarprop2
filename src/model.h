#ifndef Model_H_
#define Model_H_

#include "modelfactory.h"

#include "units.h"

class ModelInformation;
class Particle;
class SDE;

/** Abstract base class for model implementations.
 *
 * See the Custom model for a template for adding a new
 * model implementation.
 */
class Model
{
	public:
	Model() = delete;
	Model(const ModelInformation& modelInfo);
	Model(const Model&) = default;
	virtual ~Model() {}

	/**
	 *  Clone the model. Inherited classes should call their copy constructor, like this:
	 *
	 *      std::unique_ptr<Model> clone() const override { return std::make_unique<InheritedModel>(*this); }
	 *
	 *  This avoids repeated option parsing.
	 */
	virtual std::unique_ptr<Model> clone() const = 0;

	/** Create the corresponding SDE object, depending on the simplifications assumed for a given model. */
	virtual std::unique_ptr<SDE> createSDE(unsigned int globalSeed) = 0;

	/** Calculate cached values at the beginning of each pseudo-particle step. */
	virtual void calculate(const Particle& p) {}

	/** Solar wind speed at particle position. */
	virtual double solarWindSpeed(const Particle& p) const = 0;

	/** Radial extent of the heliosphere. */
	virtual double getHeliosphereBoundary() const { return heliosphereBoundary; }

	/** Value of \f$ \kappa_{rr} \f$. */
	virtual double getKappaRR() const = 0;
	/** Value of \f$ \kappa_{\theta\theta} \f$. */
	virtual double getKappaThetaTheta() const = 0;
	/** Value of \f$ \kappa_{r\theta} \f$. */
	virtual double getKappaRTheta() const { return 0.0; }
	/** Value of \f$ \kappa_{\phi\phi} \f$. */
	virtual double getKappaPhiPhi() const { return 0.0; }
	/** Value of \f$ \kappa_{r\phi} \f$. */
	virtual double getKappaRPhi() const { return 0.0; }

	/** Radial component of the gradient and curvature drift velocity: \f$ V_{dr} \f$. */
	virtual double getDriftR(const Particle& p) const { return 0.0; }
	/** Polar component of the gradient and curvature drift velocity: \f$ V_{d\theta} \f$. */
	virtual double getDriftTheta(const Particle& p) const { return 0.0; }
	/** Azimuthal component of the gradient and curvature drift velocity: \f$ V_{d\phi} \f$. */
	virtual double getDriftPhi(const Particle& p) const { return 0.0; }

	/** Radial component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	virtual double getDriftSheetR(const Particle& p) const { return 0.0; }
	/** Polar component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	virtual double getDriftSheetTheta(const Particle& p) const { return 0.0; }
	/** Azimuthal component of HCS drift velocity \f$ \vec{v}_\mathrm{hcs} \f$. */
	virtual double getDriftSheetPhi(const Particle& p) const { return 0.0; }

	/** Amplitude of magnetic field at particle position. */
	virtual double getB(const Particle& p) const { return 0.0; }
	/** Location \f$ \theta^\prime(r, \phi) \f$ of the heliospheric current sheet (HCS). */
	virtual double getThetaPrimeHCS(const Particle& p) const { return Units::pihalf; }

	/** Radial component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	virtual double getBOverB2_r(const Particle& p) const { return 0.0; }
	/** Polar component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	virtual double getBOverB2_theta(const Particle& p) const { return 0.0; }
	/** Azimuthal component of \f$ \vec{B}/B^2 \f$, for numerical drift calculation. */
	virtual double getBOverB2_phi(const Particle& p) const { return 0.0; }

	/** Value of Heaviside function H for given particle (H=+1 above HCS, H=-1 below HCS), for numerical drift calculation. */
	virtual int getHeaviside(const Particle& p) const { return 0; }

	/** 1/r^2 d/dr (r^2 kappa_rr)
	 *
	 *  Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial r}(r^2\kappa_{rr}) \f$.
	 */
	virtual double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const { return 0.0; }

	/** 1/(r^2) d/dcostheta (kappa_thetatheta*sin^2(theta))
	 *
	 *  Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial\cos\theta}(\kappa_{\theta\theta}\sin^2\theta) \f$.
	 */
	virtual double getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const { return 0.0; }

	/** -d/dcostheta (kappa_r_theta*sin(theta)/r)
	 *
	 *  Calculate \f$ -\frac{\partial}{\partial\cos\theta}\left(\kappa_{r\theta}\frac{\sin\theta}{r}\right) \f$.
	 */
	virtual double getLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const { return 0.0; }

	/** 1/(r*sin(theta)) * d/dphi (kappa_r_phi)
	 *
	 *  Calculate \f$ \frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}(\kappa_{r\phi}) \f$.
	 */
	virtual double getLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const { return 0.0; }

	/** -1/r^2 d/dr (r*kappa_r_theta*sin(theta))
	 *
	 *  Calculate \f$ -\frac{1}{r^2}\frac{\partial}{\partial r}(r\kappa_{r\theta}\sin\theta) \f$.
	 */
	virtual double getLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const { return 0.0; }

	/** 1/(r^2 sin(theta)) * d/dr (r*kappa_r_phi)
	 *
	 *  Calculate \f$ \frac{1}{r^2\sin\theta}\frac{\partial}{\partial r}(r\kappa_{r\phi}) \f$.
	 */
	virtual double getLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const { return 0.0; }

	/** 1/(r^2 sin^2(theta)) * d/dphi (kappa_phi_phi)
	 *
	 *  Calculate \f$ \frac{1}{r^2\sin^2\theta}\frac{\partial}{\partial\phi}(\kappa_{\phi\phi}) \f$.
	 */
	virtual double getLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const { return 0.0; }

	/** 1/r^2 d/dr (r^2 V_sw)
	 *
	 *  Calculate \f$ \frac{1}{r^2} \frac{\partial}{\partial r} (r^2 V_{sw}) \f$.
	 */
	virtual double getEnergyGainTerm(const Particle& p) const;

	/** Getter for the time step specified by the user or a model default. */
	virtual double getMaximumDeltaTime() const { return dt_max; }

	/** Override any automatic adjustments to the time step by returning a non-zero value. */
	virtual double overrideTimeStep(const Particle& p) const { return 0.0; }

	/** Getter for \c dynamicStep option. */
	virtual bool getDynamicStep() const { return dynamicStep; }

	protected:

	/// Time step specified by the user or a model default.
	double dt_max;

	/// Use dynamic step size?
	bool dynamicStep;

	/// Radial extent of the heliosphere.
	double heliosphereBoundary;
};

#endif
