#ifndef PARKERSPIRAL_H_
#define PARKERSPIRAL_H_

#include "properties.h"

class Particle;

/**
 *  Magnetic field configuration for models using a simple Parker spiral (with flat or wavy current sheet).
 */
class ParkerSpiral
{
	public:
	ParkerSpiral() = delete;
	ParkerSpiral(const ModelInformation& modelInfo, bool _wavyCurrentSheet);
	ParkerSpiral(const ParkerSpiral&) = default;
	virtual ~ParkerSpiral() {}

	virtual void calculate(const Particle& p);

	/** Location of heliospheric current sheet (simplified formula). */
	virtual double getThetaPrime(const Particle& p) const;

	/** Larmor radius \f$ R_L = |R/B| \f$. */
	double getLarmorRadius(const Particle& p) const;

	/** \f$ |B_\oplus/B| \f$. */
	double getBearthOverB(const Particle& p) const;

	/** Absolute value \f$ B \f$ of magnetic field. */
	double getB(const Particle& p) const;

	/** Heaviside function \f$ H(\theta - \theta^\prime) \f$ where \f$ \theta^\prime \f$ is the location of the HCS. */
	int getHeaviside(const Particle& p) const;

	/** Calculate \f$ \tan\psi \f$, where \f$ psi \f$ is the spiral angle. */
	double getTanPsi(const Particle& p) const;

	/** Calculate \f$ \sin\psi \f$, where \f$ psi \f$ is the spiral angle. */
	double getSinPsi(const Particle& p) const;

	/** Calculate \f$ \cos\psi \f$, where \f$ psi \f$ is the spiral angle. */
	double getCosPsi(const Particle& p) const;

	/** Solar wind speed. */
	double solarWindSpeed(const Particle& p) const;

	/** Component in radial direction of \f$ \vec{B}/B^2 \f$. */
	double getBOverB2_r(const Particle& p) const;

	/** Component in polar direction of \f$ \vec{B}/B^2 \f$. */
	double getBOverB2_theta(const Particle& p) const { return 0.0; }

	/** Component in azimuthal direction of \f$ \vec{B}/B^2 \f$. */
	double getBOverB2_phi(const Particle& p) const;

	/** Radial component of standard gradient and curvature drift velocity field outside HCS. */
	double getDriftVelocity_r(const Particle& p) const;

	/** Polar component of standard gradient and curvature drift velocity field outside HCS. */
	double getDriftVelocity_theta(const Particle& p) const;

	/** Azimuthal component of standard gradient and curvature drift velocity field outside HCS. */
	double getDriftVelocity_phi(const Particle& p) const;

	public:

	/// \f$ \Gamma = \tan\psi \f$.
	double Gamma;

	/// \f$ 1+\Gamma^2 \f$.
	double onePlusGamma2;

	/// Factor needed for calculation of drift velocity components: \f$ 2Ap\beta{}cr\Gamma/\left(3B_0 r_0^2 q(1+\Gamma^2)^2\right) \f$.
	double driftConstant;

	/// Drift factor times Heaviside factor
	double driftConstantTimesHeaviside;

	/// Normalization of the solar magnetic field
	double BfieldNormalization;

	/// Magnitude of solar magnetic field at Earth
	double Bearth;

	/// Square of reference position for magnetic field
	double r0_sqr;

	/// Rotation speed at solar equator
	double OmegaSun;

	/// Solar magnetic field polarity
	int polarity;

	/// HCS tilt angle
	double tiltAngle;

	/// Current sheet flat or wavy for calculation of Heaviside function?
	const bool wavyCurrentSheet;
};

#endif
