#include "binning.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

/**
 * Constructor: take a vector of energy values and construct the bin edges.
 *
 * It is assumed that the \p points are more-or-less evenly spaced on a logarithmic scale.
 *
 * @param points Energy values that define the binning. If the values are evenly space on a logarithmic scale, the \p points will be the bin centers on the log
 * scale.
 */
Binning::Binning(const std::vector<double>& points)
{
	// make first bin
	double x0 = log(points[0]);
	double x1 = log(points[1]);
	binEdges.push_back(exp(x0 - 0.5 * (x1 - x0)));

	// make bins in the middle
	auto np = points.size();
	for (unsigned int ip = 0; ip < np - 1; ++ip)
	{
		double xn = log(points[ip]);
		double xm = log(points[ip + 1]);
		binEdges.push_back(exp(0.5 * (xm + xn)));
	}

	// make last bin
	double x2 = log(points[np - 2]);
	double x3 = log(points[np - 1]);
	binEdges.push_back(exp(x3 + 0.5 * (x3 - x2)));
}

/**
 * Modify the binning in-place to obtain a finer binned version.
 *
 * @param nSubdivisionsPerPoint Number of new bins created in place of each original bin.
 */
void Binning::subdivide(unsigned int nSubdivisionsPerPoint)
{
	std::vector<double> newBinEdges;
	for (unsigned int i = 0; i < binEdges.size() - 1; ++i)
	{
		double logBinStart = log(binEdges[i]);
		double logBinWidth = log(binEdges[i + 1]) - logBinStart;

		for (unsigned int j = 0; j < nSubdivisionsPerPoint; ++j)
			newBinEdges.push_back(exp(logBinStart + double(j) / double(nSubdivisionsPerPoint) * logBinWidth));
	}

	newBinEdges.push_back(binEdges[binEdges.size() - 1]);
	binEdges = newBinEdges;
}

/** Obtain a vector of the bin centers. */
std::vector<double> Binning::calculateBinCenters() const
{
	std::vector<double> binCenters;
	binCenters.reserve(binEdges.size() - 1);
	for (unsigned int i = 0; i < binEdges.size() - 1; ++i)
	{
		binCenters.push_back(exp(0.5 * (log(binEdges[i + 1]) + log(binEdges[i]))));
	}

	return binCenters;
}

/**
 * Find the bin number for a given \p value.
 *
 * \param value Value along the binning axis.
 *
 * \returns Bin number. See the \ref bin-numbering convention.
 */
unsigned int Binning::findBin(double value) const
{
	return std::distance(binEdges.cbegin(), std::lower_bound(binEdges.cbegin(), binEdges.cend(), value, [](double a, double b) { return a <= b; }));
}

/** Number of bins. */
unsigned int Binning::nBins() const
{
	if (binEdges.empty())
		return 0;

	return binEdges.size() - 1;
}

/** Print bin edges and centers. */
void Binning::dump() const
{
	auto binCenters = calculateBinCenters();
	for (unsigned int i = 0; i < binEdges.size() - 1; ++i)
	{
		std::cout << "|" << binEdges[i] << "..(" << binCenters[i] << ")..";
	}
	std::cout << "|" << binEdges.back() << std::endl;
}

/** Get a string representation of the given bin number. See the \ref bin-numbering convention. */
std::string Binning::binAsString(unsigned int binNumber) const
{
	if (binNumber == 0)
		return std::string("underflow");
	if (binNumber > nBins())
		return std::string("overflow");
	std::stringstream s;
	s << "[" << binEdges[binNumber - 1] << ".." << binEdges[binNumber] << "]";
	return s.str();
}

/**
 * Get the lower bin edge for a given bin.
 *
 * @param binNumber Bin number, see the \ref bin-numbering convention.
 *
 * @return Lower bin edge.
 */
double Binning::lowerBinEdge(unsigned int binNumber) const
{
	assert(binNumber > 0 && binNumber <= nBins());
	return binEdges[binNumber - 1];
}

/**
 * Get the upper bin edge for a given bin.
 *
 * @param binNumber Bin number, see the \ref bin-numbering convention.
 *
 * @return Upper bin edge.
 */
double Binning::upperBinEdge(unsigned int binNumber) const
{
	assert(binNumber > 0 && binNumber <= nBins());
	return binEdges[binNumber];
}
