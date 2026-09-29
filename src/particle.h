#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "properties.h"

#include <cmath>

/**
 * A single pseudo-particle.
 */
class Particle
{
	public:
	Particle();
	Particle(const Particle&) = default;
	Particle(const ParticleProperties&, double T, int _particleNumber);
	virtual ~Particle();

	/** %Particle number. */
	int getParticleNumber() const { return particleNumber; }

	/** Radial coordinate. */
	double getR() const { return r; }

	/** Kinetic energy. */
	double getT() const { return T; }

	/** Polar coordinate. */
	double getTheta() const { return std::acos(cosTheta); }

	/** Azimuthal coordinate. */
	double getPhi() const { return phi; }

	/** Time coordinate. */
	double getTime() const { return time; }

	/** %Particle momentum. */
	double getP() const { return momentum; }

	/** Is the particle still inside the heliosphere? */
	bool isInsideHeliosphere() const { return insideHeliosphere; }

	/** Mass of the particle. */
	double getMass() const { return mass; }

	/** Closest distance to Sun reached during propagation of the particle. */
	double getRmin() const { return rmin; }

	/** Calculate \f$ \sin\theta \f$, where \f$ \theta \f$ is the polar coordinate. */
	double getSinTheta() const { return sinTheta; }

	/** Calculate \f$ \cos\theta \f$, where \f$ \theta \f$ is the polar coordinate. */
	double getCosTheta() const { return cosTheta; }

	/** Calculate \f$ \cot\theta = \cos\theta / \sin\theta \f$, where \f$ \theta \f$ is the polar coordinate. */
	double getCotTheta() const { return cosTheta / sinTheta; }

	/** Calculate the \f$ x \f$ coordinate of the particle. */
	double getX() const { return r * sinTheta * std::cos(phi); }

	/** Calculate the \f$ y \f$ coordinate of the particle. */
	double getY() const { return r * sinTheta * std::sin(phi); }

	/** Calculate the \f$ z \f$ coordinate of the particle. */
	double getZ() const { return r * cosTheta; }

	/** Calculate the \f$ \rho=\sqrt{x^2+y^2} \f$ coordinate of the particle. */
	double getRho() const { return r * sinTheta; }

	/** Rigidity of the particle. */
	double getRigidity() const { return rigidity; }

	/** Velocity of the particle. */
	double getVelocity() const { return velocity; }

	void setR(double _r);
	void setCosTheta(double _costheta);
	void setPhi(double _phi) { phi = _phi; }
	void setKineticEnergy(double _T);
	void setMomentum(double _p);
	void setTime(double _time);
	void setInsideHeliosphere(bool _val) { insideHeliosphere = _val; }

	std::string asString() const;
	void dump() const;

	private:

	friend class Model;
	friend class ParkerSpiral;
	friend class ModifiedSpiral;

	friend class BurgerPotgieter;
	friend class JokipiiKopriva;
	friend class PotgieterMoraal;
	friend class Standard2D;
	friend class Simple2D;
	friend class Helmod;
	friend class Helmod2012;
	friend class Helmod2018;
	friend class Burger2012;
	friend class Strauss2012;
	friend class Yamada;
	friend class Custom;

	protected:
	/// radial coordinate
	double r;
	/// polar coordinate: \f$ \cos\theta \f$.
	double cosTheta;
	/// azimuthal coordinate
	double phi;
	/// time coordinate
	double time;
	/// still inside heliosphere?
	bool insideHeliosphere;

	/// Charge of particle, in units of the proton charge.
	int charge;
	/// Mass of particle
	double mass;
	/// Kinetic energy of particle
	double T;

	/// Unique number for particle (useful for unique random seed, and when storing individual steps in SDE)
	int particleNumber;

	/// Closest distance to Sun during propagation.
	double rmin;

	/// \f$ \sin\theta \f$
	double sinTheta;

	/// %Particle velocity, in units of speed of light
	double beta;

	/// %Particle velocity
	double velocity;

	/// Absolute value of particle rigidity
	double absRigidity;

	/// %Particle rigidity
	double rigidity;

	/// %Particle momentum
	double momentum;
};

#endif
