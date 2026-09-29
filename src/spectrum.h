#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <string>
#include <vector>

#include "properties.h"

/**
 * An energy spectrum of a cosmic-ray species.
 *
 * The spectrum is characterized by values of kinetic energies and the
 * cosmic-ray flux at these energies.
 */
class Spectrum
{
	public:
	Spectrum(const Spectrum&) = default;
	Spectrum(const std::vector<double>&, const std::vector<double>&);
	~Spectrum() {}

	/** Values of the kinetic energy. */
	const std::vector<double>& getEkin() const { return ekin; }
	/** Corresponding flux values. */
	const std::vector<double>& getFlux() const { return flux; }

	void print(const std::string& prefix) const;

	void quality(unsigned int);

	double approximateSpectralIndex(double T) const;

	private:
	/// Vector of kinetic energies where the flux is defined.
	std::vector<double> ekin;
	/// Vector of the flux values corresponding to the values in ekin.
	std::vector<double> flux;

	void dump(std::ostream& out, const std::string& prefix) const;

	double logEvalFlux(double T) const;
};

#endif
