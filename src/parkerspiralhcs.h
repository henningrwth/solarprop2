#ifndef ParkerspiralHCS_H_
#define ParkerspiralHCS_H_

#include "parkerspiral.h"

#include <gsl/gsl_multimin.h>

class Particle;

/**
 * Simple Parker spiral as magnetic field, but with detailed treatment of HCS drift.
 *
 * HCS drift treated as in Strauss et al., Astrophys. Space Sci. (2012) 339:223-236,
 * https://doi.org/10.1007/s10509-012-1003-z
 *
 * \sa Strauss2012 model
 *
 * Throughout the model, the angle \f$ \beta^\prime \f$ is the angle between the radial
 * direction and a tangent to the HCS.
 *
 * The main difficulty in the model implementaion is finding the point of closest distance
 * to a particle position on the HCS. This problem is solved here by the Simplex minimization
 * algorithm implemented in the GSL.
 *
 * \note In this specialization of the %ParkerSpiral, the full formula for the location of the
 * HCS is used instead of the simplified one.
 *
 */
class ParkerSpiralHCS : public ParkerSpiral
{

	public:
	ParkerSpiralHCS() = delete;
	ParkerSpiralHCS(const ModelInformation& modelInfo);
	ParkerSpiralHCS(const ParkerSpiralHCS& s);
	~ParkerSpiralHCS() override;

	void calculate(const Particle& p) override;

	/** Location of the HCS (full formula). */
	double getThetaPrime(const Particle& p) const override;

	/** Location of the HCS for arbitrary coordinates. */
	double thetaPrimeHcs(double r, double phi) const;

	/** Calculate \f$ \tan\psi \f$, where \f$ \psi \f$ is the spiral angle, for an arbitrary point. */
	double tanPsi(double r, double theta) const;

	/** Calculate \f$ \sin\psi \f$, where \f$ \psi \f$ is the spiral angle, for an arbitrary point. */
	double sinPsi(double r, double theta) const;

	/** Calculate \f$ \cos\psi \f$, where \f$ \psi \f$ is the spiral angle, for an arbitrary point. */
	double cosPsi(double r, double theta) const;

	/** Calculate \f$ \sin^2\psi \f$, where \f$ \psi \f$ is the spiral angle, for an arbitrary point. */
	double sinPsiSqr(double r, double theta) const;

	/** Calculate the sign of the \f$ \beta^\prime \f$ angle for a point on the HCS. */
	double betaSign(double r, double phi) const;

	/** Calculate \f$ \tan^2\beta^\prime \f$ for a point on the HCS. */
	double tanBetaSqr(double r, double theta) const;

	/** Calculate \f$ \cos\beta^\prime \f$ for a point on the HCS. */
	double cosBeta(double r, double theta) const;

	/** Calculate \f$ \sin\beta^\prime \f$ for a point on the HCS. */
	double sinBeta(double r, double theta, double phi) const;

	/** Calculate \f$ \sin\beta^\prime \f$ and \f$ \cos\beta^\prime \f$ for a point on the HCS in one go. */
	std::tuple<double, double> sincosBeta(double r, double theta, double phi) const;

	/** Calculate \f$ \beta^\prime \f$ for a point on the HCS. */
	double beta(double r, double theta, double phi) const;

	/** Calculate the squared distance of a particle to the HCS point defined by the coordinates \f$ (r, \phi) \f$. */
	double distanceToHcsPointSqr(double r, double phi, const Particle& p) const;

	/** Calculate the closest distance of a particle to the HCS.
	 *
	 * @returns tupel of (distance, \f$ r_\mathrm{hcs} \f$, \f$ \theta_\mathrm{hcs} \f$, \f$ \phi_\mathrm{hcs} \f$):
	 * the distance of the particle to the HCS and the spherical coordinates of the closest point on the HCS.
	 */
	std::tuple<double, double, double, double> closestDistanceToHcs(const Particle& p);

	private:
	/// global \f$ \phi \f$ offset for the HCS
	const double phi_0 = 0.;

	/// store particle position during minimization of distance to HCS
	const Particle* particlePosition;

	/// \f$ \Omega / V_{sw} \f$ (assumed to be constant throughout the heliosphere!)
	double omegaOverVsw = 0.;

	/// start point for minimization to find closest point on HCS
	gsl_vector* x0_min;
	/// step sizes for minimization to find closest point on HCS
	gsl_vector* stepsizes_min;
	/// minimizer for finding closest point on HCS
	gsl_multimin_fminimizer* minimizer = NULL;
};

#endif
