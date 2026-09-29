#include "particle.h"
#include "units.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_particle(py::module_& m)
{
	py::class_<Particle>(m, "Particle")
			.def("__repr__", &Particle::asString)
			.def("particleNumber", &Particle::getParticleNumber)
			.def("r", &Particle::getR)
			.def("T", &Particle::getT)
			.def("theta", &Particle::getTheta)
			.def("phi", &Particle::getPhi)
			.def("time", &Particle::getTime)
			.def("p", &Particle::getP)
			.def("x", &Particle::getX)
			.def("y", &Particle::getY)
			.def("z", &Particle::getZ)
			.def("rho", &Particle::getRho)
			.def("rmin", &Particle::getRmin)
			.def("rigidity", &Particle::getRigidity)
			.def("velocity", &Particle::getVelocity)
			.def("isInsideHeliosphere", &Particle::isInsideHeliosphere)
			;
}
#endif

using namespace Units;

/**
 * Default constructor.
 */
Particle::Particle() :
	Particle({0, 0.}, 0., 0)
{
}

/**
 * Constructor with mass, charge, kinetic energy and particle number.
 *
 * The start position is at the location of the Earth.
 *
 * @param sp particle properties (mass and charge)
 * @param T kinetic energy
 * @param _particleNumber particle number
 */
Particle::Particle(const ParticleProperties& sp, double T, int _particleNumber)
{

	// radial position of the earth in the heliosphere
	r = 1.0 * AU;
	rmin = r;

	setCosTheta(0.0);
	phi = 0.0;
	time = 0;
	insideHeliosphere = true;

	mass = sp.mass;
	charge = sp.charge;
	setKineticEnergy(T);

	particleNumber = _particleNumber;
}

/**
 * @brief Destructor
 */
Particle::~Particle() {}

/** Set radial coordinate.
 *
 * This will also keep `rmin` updated.
 */
void Particle::setR(double _r)
{
	r = _r;

	if (r < rmin)
		rmin = r;
}

/** Set \f$ \cos\theta \f$.
 *
 * This will recalculate \f$ \sin\theta \f$ accordingly.
 */
void Particle::setCosTheta(double _costheta)
{
	cosTheta = _costheta;
	sinTheta = sqrt(1.0 - cosTheta*cosTheta);
}

/**
 * Set the kinetic energy.
 *
 * The other kinematic variables will be recalculated accordingly.
 *
 * @param _T Value for the kinetic energy.
 */
void Particle::setKineticEnergy(double _T)
{
	T = _T;

	momentum = sqrt(pow(T/c,2) + 2.*T*mass);
	rigidity = momentum/charge;
	absRigidity = std::abs(rigidity);
	beta = momentum*c/(T+mass*c2);
	velocity = beta * c;
}

/**
 * Set the momentum.
 *
 * The other kinematic variables will be recalculated accordingly.
 *
 * @param _p Value for the momentum.
 */
void Particle::setMomentum(double _p)
{
	momentum = _p;

	double mc2 = mass*c2;
	double E = sqrt(pow(_p*c, 2) + pow(mc2, 2));
	T = E - mc2;

	rigidity = momentum/charge;
	absRigidity = std::abs(rigidity);
	beta = momentum*c/E;
	velocity = beta * c;
}

/** Set the time. */
void Particle::setTime(double _time)
{
	time = _time;
}

/**
 * Get particle state as string.
 *
 * @return A string representation of the particle's state.
 */
std::string Particle::asString() const
{
	std::stringstream sstr;
	// dump i, t, r, x, z, T
	sstr << "# " << std::setw(5) << getParticleNumber() << " "
	     << "t= " << std::setw(12) << getTime() / s << " s, "
	     << "r= " << std::setw(12) << getR() / AU << " "
	     << "x= " << std::setw(12) << getX() / AU << " "
	     << "y= " << std::setw(12) << getY() / AU << " "
	     << "rho= " << std::setw(12) << getRho() / AU << " "
	     << "z= " << std::setw(12) << getZ() / AU << " "
	     << "T= " << std::setw(12) << getT() / GeV << " "
	     << "rmin= " << std::setw(12) << getRmin() / AU << " AU "
	     << "theta= " << std::setw(12) << getTheta() / degree << " deg "
	     << "phi= " << std::setw(12) << getPhi() / degree << " deg";

	return sstr.str();
}

/**
 * Print the particle state.
 */
void Particle::dump() const
{
	std::cout << asString() << std::endl;
}
