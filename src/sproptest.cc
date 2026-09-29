#include "solarprop.h"

#include <pybind11/embed.h>
namespace py = pybind11;

#include <iomanip>
#include <iostream>

/// \file sproptest.cc
/// \brief Main function for a standalone C++ executable that can be used for debugging and profiling.

///
/// Main function for a standalone executable that can be used
/// for debugging and profiling.
///
/// The standalone executable is built when the build is run manually
/// (outside of pip).
///
int main()
{
	// start the python interpreter and keep it alive
	py::scoped_interpreter guard{};

	// set the model parameters
	py::dict parameters;
	parameters["model"] = "JK1979";
	parameters["mass"] = 0.938272;
	parameters["charge"] = 1;
	parameters["polarity"] = 1;
	parameters["forceNumericalDerivatives"] = true;
	parameters["save"] = "sprop_test.green";

	// create the Solarprop instance
	Solarprop s(parameters);

	// define the energy binning
	std::vector<double> ekin(32);
	double E0 = 0.01;
	double E1 = 50.;
	for (size_t i = 0; i < ekin.size(); ++i)
	{
		ekin[i] = E0 * std::pow(E1 / E0, double(i) / (ekin.size() - 1));
	}

	// vector for input LIS
	std::vector<double> flux(ekin.size(), 0.0);

	// set bogus LIS values
	auto lis = [](double T) { return std::pow(T, -3.6); };
	for (size_t i = 0; i < ekin.size(); ++i)
		flux[i] = lis(ekin[i]);

	// perform the modulation run
	auto res = s.modulate(ekin, flux);

	// print the results
	std::cout << "Ekin (GeV)   modulated flux (a.u.)   LIS flux (arbitrary units)" << std::endl;
	for (size_t i = 0; i < ekin.size(); ++i)
	{
		std::stringstream s;
		s << std::setw(10) << ekin[i] << " " << std::setw(16) << res[i] << " " << std::setw(17) << flux[i];
		std::cout << s.str() << std::endl;
	}

	return 0;
}
