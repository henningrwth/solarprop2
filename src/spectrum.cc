#include "spectrum.h"
#include "binning.h"
#include "particle.h"
#include "units.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

using namespace Units;

/**
 * Constructor with energy and flux data.
 *
 * \param _ekin Vector of kinetic energies.
 * \param _flux Vector of corresponding interstellar flux values.
 */
Spectrum::Spectrum(const std::vector<double>& _ekin, const std::vector<double>& _flux)
{
	ekin = _ekin;
	flux = _flux;
}

/**
 * Print energy and flux values.
 */
void Spectrum::print(const std::string& prefix) const
{
	dump(std::cout, prefix);
}

/**
 * Stream energy and flux values.
 *
 * \param out Output stream.
 * \param prefix A prefix to be added at the beginning of each output line.
 */
void Spectrum::dump(std::ostream& out, const std::string& prefix) const
{
	for (unsigned int i = 0; i < ekin.size(); ++i)
	{
		if (!prefix.empty())
			out << prefix << " ";

		out << ekin[i] / GeV << "    " << flux[i] << std::endl;
	}
}

/**
 * Split the energy binning to a finer version and interpolate the flux values at the new bin centers.
 *
 * \attention The Spectrum object is modified in place.
 *
 * @param extra Number of new bins created in place of each original bin.
 */
void Spectrum::quality(unsigned int extra)
{
	if (extra < 1)
		return;

	Binning ekinBinning(ekin);
	ekinBinning.subdivide(extra);

	std::vector<double> newEkin = ekinBinning.calculateBinCenters();
	std::vector<double> newFlux(newEkin.size(), 0.0);

	for (unsigned int i = 0; i < newEkin.size(); ++i)
	{
		newFlux[i] = logEvalFlux(newEkin[i]);
	}

	ekin = newEkin;
	flux = newFlux;
}

/**
 * Evaluate the spectrum at a given kinetic energy \p T using logarithmic interpolation.
 *
 * This assumes that the spectrum is well approximate by a power law locally.
 */
double Spectrum::logEvalFlux(double T) const
{
	const double* fX = ekin.data();
	const double* fY = flux.data();

	// find neighbours simply looping all points
	int low = -1;
	int up = -1;

	if (T < fX[0])
	{
		low = 0;
		up = 1;
	}
	else if (T > ekin.back())
	{
		int nBins = ekin.size();
		low = nBins - 2;
		up = nBins - 1;
	}
	else
	{
		for (unsigned int i = 0; i < ekin.size(); ++i)
		{
			if (std::abs(T - fX[i]) < eV)
			{
				return fY[i]; // no interpolation needed
			}
			else if (fX[i] < T)
			{
				if (low == -1 || fX[i] > fX[low])
					low = i;
			}
			else // (fX[i] > T)
			{
				if (up == -1 || fX[i] < fX[up])
					up = i;
			}
		}
	}

	// use log interpolation
	double logx = log(T);
	double logx1 = log(fX[low]);
	double logx2 = log(fX[up]);
	double logy1 = log(fY[low]);
	double logy2 = log(fY[up]);

	double ylog = logy2 + (logx - logx2) * (logy1 - logy2) / (logx1 - logx2);

	return std::exp(ylog);
}

/** Calculate the approximate spectral index (logarithmic slope) of the spectrum at the given kinetic energy \p T. */
double Spectrum::approximateSpectralIndex(double T) const
{
	assert(ekin.size() > 1);
	int nBins = ekin.size();

	// find neighbours simply looping all points
	int low = -1;
	int up = -1;

	const double* fX = ekin.data();
	const double* fY = flux.data();
	for (unsigned int i = 0; i < ekin.size(); ++i)
	{

		if (fX[i] < T)
		{
			if (low == -1 || fX[i] > fX[low])
				low = i;
		}
		else if (fX[i] > T)
		{
			if (up == -1 || fX[i] < fX[up])
				up = i;
		}
		else // case x == fX[i]
		{
			low = i - 1;
			up = i + 1;
			break;
		}
	}

	if (low < 0)
		low = 0;
	if (up >= nBins)
		up = nBins - 1;

	// use log interpolation
	double logx1 = log(fX[low]);
	double logx2 = log(fX[up]);
	double logy1 = log(fY[low]);
	double logy2 = log(fY[up]);

	return (logy2 - logy1) / (logx2 - logx1);
}
