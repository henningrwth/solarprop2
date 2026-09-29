#ifndef BINNING_H_
#define BINNING_H_

#include <string>
#include <vector>

/**
 * Helper class for turning a list of energy points into a near-logarithmic binning.
 *
 * The input spectrum to \c SOLARPROP will usually be defined for a list of energies
 * more-or-less evenly spaced on a logarithmic scale. To avoid numerical artifacts in
 * the results, a fine energy binning is needed in the calculation of the Green transfer
 * matrix. Therefore we first use this class to turn a list of points into a binning,
 * conserving the near-logarithmic spacing, and then subdividing the binning to get
 * a sufficiently fine one using subdivide().
 *
 * \anchor bin-numbering
 * #### Bin numbering convention used throughout this class:
 *
 * Let \c N be the number of bins. Then bin numbers are in the range from \c 1 to \c N, or \c 0 for underflow or \p N+1 for overflow.
 */
class Binning
{
	public:

	Binning() = delete;
	Binning(const std::vector<double>& points);

	void subdivide(unsigned int nSubdivisionsPerPoint);

	std::vector<double> calculateBinCenters() const;

	unsigned int findBin(double value) const;

	unsigned int nBins() const;

	void dump() const;
	std::string binAsString(unsigned int binNumber) const;
	double lowerBinEdge(unsigned int binNumber) const;
	double upperBinEdge(unsigned int binNumber) const;

	private:

	/// the list of bin edges
	std::vector<double> binEdges;
};

#endif
