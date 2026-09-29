#ifndef SOLARPROP_H
#define SOLARPROP_H

#include "arrows.h"
#include "density.h"
#include "green.h"
#include "particle.h"
#include "properties.h"
#include "spectrum.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#ifdef PYTHON_BINDINGS
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

void bind_arrows(py::module_&);
void bind_modelinformation(py::module_& m);
void bind_particle(py::module_&);
void bind_units(py::module_&);
void bind_solarprop(py::module_&);
#endif

/**
 * \mainpage
 *
 * This documentation covers the C++ part of \c SOLARPROP. Please see the file \c README.md for a general introduction to \c SOLARPROP.
 *
 * #### Quick overview of the structure of the C++ side
 *
 * Typical python code used to calculate a solution of the transport equation for cosmic rays in the heliosphere will first
 * choose the modulation model, set its parameters, define a list of energies and the local interstellar spectrum (LIS)
 * at these energies. Then an instance of the Solarprop C++ interfrace is created like this:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~{.py}
 * s = solarprop.Solarprop(parameters)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Finally, the solution of the transport equation is calculated by, e.g.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~{.py}
 * flux_mod = s.modulate(ekin, LIS(ekin, parameters['mass']))
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * which invokes Solarprop::modulate(). This function will either calculate the Green function matrix or read it from a
 * file stored in a previous run (if the \c green option is used), and then call Solarprop::modulateSpectrum(), which applies
 * the probabilities stored in the Green function matrix to the local interstellar spectrum to obtain the modulated flux at
 * the top of the atmosphere. The calculation of the Green function matrix is done in Solarprop::performModulationRun() which
 * uses a Density object as the main engine: In the constructor Density::Density(), the pseudo-particles are created. In
 * Density::simulate(), the SDE objects containing the stochastic differential equations are set up, one for each thread, and
 * an OpenMP-parallelized loop over the pseudo-particles is executed, tracking each one from its initial position at Earth
 * to the edge of the heliosphere. From the Density object, which keeps track of the final states of the pseudo-particles, the
 * Green function matrix is filled in GreenFunctionMatrix::fillFromDensity().
 *
 * #### Standalone executable
 *
 * When the solarprop code is compiled outside of \c pip, a standalone executable \c sprop_test is created, which can be used
 * for debugging and profiling. See the main() function.
 */

/**
 * Main interface class to the Python side.
 *
 * This class provides the python bindings and the steering code to perform the calculations for solar modulation.
 *
 * The python bindings are generated with [pybind11](https://pybind11.readthedocs.io).
 *
 * The main result of a modulation run is the probability (Green function) matrix that contains the probability
 * for a particle in kinetic energy bin \f$ T_i \f$ at Earth to originate from the kinetic energy bin \f$ T_j \f$
 * at the heliopause. Once this matrix has been calculated, it can be used to %modulate different local
 * interstellar spectra quickly.
 *
 * \attention Be aware that the energy smearing of the initial pseudo-particle states is based on the first LIS
 * that is seen. This means that the results obtained when the modulate() function is called for a new LIS will
 * be slightly different from those that would be obtained if the new LIS had been used for the first call to
 * modulate().
 *
 */
class Solarprop
{
	public:
	Solarprop() = delete;
	Solarprop(const py::dict& parameters);
	~Solarprop();

	void updateParameters(const py::dict& parameters);

	static std::string versionString();
	void printConfig();

#ifdef PYTHON_BINDINGS
	py::array modulate(const std::vector<double>& ekin_lis_in_GeV, const std::vector<double>& flux_lis);
#else
	std::vector<double> modulate(const std::vector<double>& ekin_lis_in_GeV, const std::vector<double>& flux_lis);
#endif

	ParticleHistories trajectories(unsigned int n, double ekin_in_GeV);
	std::vector<Particle> endpoints(unsigned int n, double ekin_in_GeV);
	std::vector<Arrows> arrows(double range_in_AU, double ekin_in_GeV, int grid_points_per_direction);

	double magneticField(double r, double theta, double phi);
	double thetaHCS(double r, double phi);

	private:
	/// The ModelInformation object provides access to all user settings and the relevant lookups for time-dependent heliospheric parameters.
	ModelInformation modelInfo;

	/// Cache model for magneticField() and thetaHCS() functions.
	std::unique_ptr<Model> modelForMagneticField = nullptr;

	/// Cache Green matrix for repeated modulation of different LIS.
	GreenFunctionMatrix greenMatrix;

	/// Read/write Green matrix data in binary mode to save disk space?
	const bool greenFuncBinaryMode = true;

	ParticleHistories trajectories_or_endpoints(unsigned int n, double ekin, int historyLevel);

	GreenFunctionMatrix performModulationRun(const std::vector<double>& ekinTOA, const Spectrum& LIS, const ModelInformation&) const;
	std::vector<double> modulateSpectrum(const std::vector<double>& ekin_toa, const Spectrum& lis, double mass) const;
	double squaredMomentum(double T, double mass) const;
};

#ifdef PYTHON_BINDINGS
PYBIND11_MODULE(solarprop, m)
{

	bind_units(m);
	bind_arrows(m);
	bind_particle(m);
	bind_solarprop(m);
	bind_modelinformation(m);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
#endif

#endif
