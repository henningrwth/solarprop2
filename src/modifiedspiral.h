#ifndef Modifiedspiral_H_
#define Modifiedspiral_H_

#include "properties.h"

class Particle;

/**
 * Properties of a Parker spiral magnetic field with polar modification according to Jokipii and Kóta.
 *
 * The magnetic field is given by
 * \f[
 * \vec{B}=AB_0\frac{r_0^2}{r^2}\left[\,\vec{e}_r - \Gamma\,\vec{e}_\phi + \frac{\delta(\theta)}{r_\odot}r\,\vec{e}_\theta\right] \times \left[1-2H(\theta-\theta^\prime)\right]
 * \f]
 *
 * This form is used by the Helmod2012 and Helmod2018 models.
 *
 * For the derivation of the gradient and curvature drift velocity terms, see the \c jupyter
 * notebook \c calculations/helmod_drift.ipynb.
 */
class ModifiedSpiral
{

	public:
	ModifiedSpiral() = delete;
	ModifiedSpiral(const ModelInformation& modelInfo, bool _isHighSolarActivity, int _polarity, int _phase, double _tiltAngle, double _Bearth,
		       double _vswMin, double _delta0);
	ModifiedSpiral(const ModifiedSpiral&) = default;
	~ModifiedSpiral() {}

	void calculate(const Particle& p);

	double getLarmorRadius(const Particle& p) const;
	double getBearthOverB(const Particle& p) const;
	double getB(const Particle& p) const;
	double getThetaPrime(const Particle& p) const;

	double solarWindSpeed(const Particle& p) const;

	double getBOverB2_r(const Particle& p) const;
	double getBOverB2_theta(const Particle& p) const;
	double getBOverB2_phi(const Particle& p) const;

	double getDriftVelocity_r(const Particle& p) const;
	double getDriftVelocity_theta(const Particle& p) const;

	double deltaFunction(const Particle& p) const;
	double deltaPrimeFunction(const Particle& p) const;

	double dGamma_dcostheta(const Particle& p) const;

	double dvsw_dcostheta(const Particle& p) const;

	public:

	/// cache solar wind speed
	double vsw;

	/// \f$ \Gamma \f$
	double Gamma;

	/// \f$ 1+\Gamma^2 \f$
	double onePlusGamma2;

	/// \f$ \delta(\theta) \f$: correction factor for the magnetic field in polar direction
	double delta;

	/// extension function of the magnetic field: \f$ t(r,\mu) = \delta(\mu) \cdot r/r_\odot \f$
	double t;

	/// \f$ \mathrm{d}t/\mathrm{d}\cos\theta \f$
	double dtdmu;

	/// \f$ \mathrm{d}\Gamma/\mathrm{d}\cos\theta \f$
	double Gammaprime;

	/// \f$ 1+\tan^2\psi + \delta^2\cdot(r/r_\odot)^2 \f$
	double N;

	/// \f$ 1+t^2 \f$
	double tau;

	/// polarity (of magnetic field) times value of Heaviside function (above or below HCS)
	int sign;

	/// Value of Heaviside function \f$ H \f$ for given particle (\f$ H=+1 \f$ above HCS, \f$ H=-1 \f$ below HCS)
	int heavi;

	/// Factor needed for calculation of drift velocity components: \f$ 2Ap\beta cr\Gamma/(3B_0 r_0^2 q (1+\Gamma^2)^2 \f$
	double driftConstant;

	/// Normalization of the solar magnetic field
	double BfieldNormalization;

	/// Magnitude of solar magnetic field at Earth
	double Bearth;

	/// Square of reference position for magnetic field
	const double r0_sqr;

	/// Solar radius
	const double rSun;

	/// Solar wind speed according to OmniWEB data
	double solarWindSpeedMin;

	/// Maximum solar wind speed
	const double solarWindSpeedMax;

	/// When solar activity is high, a constant solar wind speed throughout heliosphere is assumed.
	const bool isHighSolarActivity;

	/// Rotation speed at solar equator
	const double OmegaSun;

	/// Solar magnetic field polarity
	const int polarity;

	/// \brief Solar activity phase.
	///
	/// Possible values:
	/// - \c +1 ascending
	/// - \c -1 descending
	/// - \c 0 unknown
	const int phase;

	/// Smoothed sunspot number
	double ssn;

	/// HCS tilt angle
	const double tiltAngle;

	/// amplitude for magnetic field correction function
	const double delta0;

	/// polar angle at which magnetic field correction saturates
	const double thetaThresholdForCorrection;

	private:
	int getHeaviside(const Particle& p) const;
};

#endif
