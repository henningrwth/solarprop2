#ifndef NUMERICALDERIVATIVES_H_
#define NUMERICALDERIVATIVES_H_

#include <memory>

class Model;
class Particle;

/**
 * Helper class to perform numerical differentiation for diffusion and drift terms.
 *
 * This provides independent error checking for terms involving the derivatives of
 * diffusion coefficients and magnetic fields, though it will increase the runtime.
 *
 * See the JokipiiKopriva model for an example how to use this class. Note that the
 * model class contains a \c shared_ptr to its NumericalDerivatives object. This is
 * necessary so that we can clone the model using the default copy constructor, but
 * it only works because the NumericalDerivatives object is only created in the
 * \c calculate() function of the model (if needed), and at that point, no further
 * cloning will take place.
 *
 */
class NumericalDerivatives
{
	public:
	NumericalDerivatives() = delete;
	explicit NumericalDerivatives(const Model*);
	~NumericalDerivatives() {}

	/** Calculate \f$ \frac{1}{r^2} \frac{\partial}{\partial r} (r^2 V_{sw}) \f$ numerically. */
	double calculateEnergyGainTerm(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial r}(r^2\kappa_{rr}) \f$ numerically. */
	double calculateLinearRadialDiffusionTermFromKappaRR(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r\sin\theta}\frac{\partial}{\partial\phi}(\kappa_{r\phi}) \f$ numerically. */
	double calculateLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r^2\sin\theta} \frac{\partial}{\partial\theta} (\sin\theta \kappa_{\theta\theta}) \f$ numerically. */
	double calculateLinearThetaDiffusionTermFromKappaThetaTheta(const Particle& p) const;

	/** Calculate \f$ -\frac{\partial}{\partial\cos\theta}\left(\kappa_{r\theta}\frac{\sin\theta}{r}\right) \f$ numerically. */
	double calculateLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const;

	/** Calculate \f$ -\frac{1}{r^2}\frac{\partial}{\partial r}(r\kappa_{r\theta}\sin\theta) \f$ numerically. */
	double calculateLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r^2}\frac{\partial}{\partial\cos\theta}(\kappa_{\theta\theta}\sin^2\theta) \f$ numerically. */
	double calculateLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r^2\sin\theta}\frac{\partial}{\partial r}(r\kappa_{r\phi}) \f$ numerically. */
	double calculateLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const;

	/** Calculate \f$ \frac{1}{r^2\sin^2\theta}\frac{\partial}{\partial\phi}(\kappa_{\phi\phi}) \f$ numerically. */
	double calculateLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const;


	/** Calculate radial component of \f$ (Rv/3)\cdot\vec{\nabla}\times\vec{B}/B^2 \f$ numerically, excluding delta function from Heaviside function (HCS drift). */
	double calculateRadialGradientCurvatureDriftVelocity(const Particle& p, bool includePolarComponent = false) const;

	/** Calculate polar component of \f$ (Rv/3)\cdot\vec{\nabla}\times\vec{B}/B^2 \f$ numerically, excluding delta function from Heaviside function (HCS drift). */
	double calculatePolarGradientCurvatureDriftVelocity(const Particle& p) const;

	/** Calculate azimuthal component of \f$ (Rv/3)\cdot\vec{\nabla}\times\vec{B}/B^2 \f$ numerically, excluding delta function from Heaviside function (HCS drift). */
	double calculateAzimuthalGradientCurvatureDriftVelocity(const Particle& p, bool includePolarComponent = false) const;

	private:
	static double epsilon_r(const Particle& p);
	static double epsilon_theta(const Particle& p);
	static double epsilon_costheta(const Particle& p);
	static double epsilon_phi(const Particle& p);

	private:
	/// Pointer to the Model for which derivatives are calculated.
	const Model* model;

	/// Clone of the Model instance, which will be used to calculate the infinitesimal shifts, without changing the original model instance.
	std::unique_ptr<Model> modelForDerivatives;
};

#endif
