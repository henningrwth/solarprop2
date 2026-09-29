#include "green.h"
#include "binning.h"
#include "model.h"
#include "particle.h"
#include "stringtools.h"
#include "units.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std::string_literals;

using namespace Units;

/**
 * Constructor defining the energy values for the TOA and LIS spectra.
 *
 * All probabilities will initially be set to zero.
 *
 * @param ekinTOA Kinetic energies at the top of the atmosphere.
 * @param ekinLIS Kinetic energies for the local interstellar spectrum.
 *
 */
GreenFunctionMatrix::GreenFunctionMatrix(const std::vector<double>& ekinTOA, const std::vector<double>& ekinLIS) :
	fEkinTOA(ekinTOA),
	fEkinLIS(ekinLIS)
{
	// initialize to zero
	for (unsigned int i = 0; i < ekinTOA.size(); ++i)
	{
		fProbabilityMatrix.push_back(std::vector<double>(ekinLIS.size(), 0.0));
	}
}

/**
 * Constructor for reading a previously stored matrix from disk.
 *
 * @param infilename Input file name.
 * @param binaryMode If true, values are stored in binary format instead of plain text. Binary format takes less space, but plain text is human-readable.
 */
GreenFunctionMatrix::GreenFunctionMatrix(const std::string& infilename, bool binaryMode)
{
	std::ios::openmode mode = std::ios::in;
	if (binaryMode)
		mode |= std::ios::binary;

	std::ifstream file(infilename, mode);

	if (!file.is_open())
		throw std::runtime_error("Cannot read Green function matrix from file \""s + infilename + "\"."s);

	if (binaryMode)
	{
		size_t ekin_size_toa;
		file.read(reinterpret_cast<char*>(&ekin_size_toa), sizeof(ekin_size_toa));
		size_t ekin_size_lis;
		file.read(reinterpret_cast<char*>(&ekin_size_lis), sizeof(ekin_size_lis));

		for (size_t i = 0; i < ekin_size_toa; ++i)
		{
			float ek;
			file.read(reinterpret_cast<char*>(&ek), sizeof(ek));
			fEkinTOA.push_back(ek * GeV);
		}

		for (size_t j = 0; j < ekin_size_lis; ++j)
		{
			float ek;
			file.read(reinterpret_cast<char*>(&ek), sizeof(ek));
			fEkinLIS.push_back(ek * GeV);
		}

		for (size_t i = 0; i < ekin_size_toa; ++i)
		{
			std::vector<double> m;

			for (size_t j = 0; j < ekin_size_lis; ++j)
			{
				float p;
				file.read(reinterpret_cast<char*>(&p), sizeof(p));
				m.push_back(double(p));
			}

			fProbabilityMatrix.push_back(m);
		}
	}
	else
	{
		std::string line;

		// read ekin values
		std::getline(file, line);
		const auto& tokenstoa = split(line, " ");
		for (const auto& s : tokenstoa)
			fEkinTOA.push_back(stod(s) * GeV);

		std::getline(file, line);
		const auto& tokenslis = split(line, " ");
		for (const auto& s : tokenslis)
			fEkinLIS.push_back(stod(s) * GeV);

		// read matrix
		while (std::getline(file, line))
		{
			std::vector<double> m;
			const auto& tokens = split(line, " ");
			for (const auto& s : tokens)
				m.push_back(stod(s));
			fProbabilityMatrix.push_back(m);
		}
	}

	file.close();
}

/** Clear the energy vectors and probability matrix. */
void GreenFunctionMatrix::reset()
{
	fEkinTOA.clear();
	fEkinLIS.clear();
	fProbabilityMatrix.clear();
}

/**
 * Fill the probability matrix from a Density object.
 *
 * Count the number of particles with matching energy for each LIS bin.
 * Particles that have not reached the outer edge of the heliosphere
 * are discarded. Normalize by the number of particles simulated per
 * energy bin.
 *
 * @param den Input Density object.
 */
void GreenFunctionMatrix::fillFromDensity(const Density& den)
{
	Binning ekinBinning(fEkinLIS);

	for (unsigned int ekinIndexTOA = 0; ekinIndexTOA < fEkinTOA.size(); ++ekinIndexTOA)
	{
		std::vector<double>& flux = fProbabilityMatrix.at(ekinIndexTOA);

		int red = 0;

		for (unsigned int i = 0; i < den.numberOfParticlesPerEnergyBin(); i++)
		{
			const Particle& p = den.getParticle(ekinIndexTOA, i);

			if (!p.isInsideHeliosphere())
			{
				unsigned int bin = ekinBinning.findBin(p.getT());
				if (bin == 0 || bin > ekinBinning.nBins())
					continue;

				flux[bin - 1] += 1.0;
			}
			else // Take care of the normalization if a particle did not reach the heliosphere in time
			{
				red++;
			}
		}

		for (unsigned int j = 0; j < flux.size(); j++)
		{
			flux[j] /= (den.numberOfParticlesPerEnergyBin() - red);
		}
	}
}

/**
 * Extract stored probability \f$ G(T_i,T_j) \f$.
 * @param ekinIndexTOA Index of TOA energy \f$ T_i \f$.
 * @param ekinIndexLIS Index of LIS energy \f$ T_j \f$.
 *
 * @return Probability \f$ G(T_i,T_j) \f$.
 */
double GreenFunctionMatrix::probability(unsigned int ekinIndexTOA, unsigned int ekinIndexLIS) const
{
	return fProbabilityMatrix.at(ekinIndexTOA).at(ekinIndexLIS);
}

/**
 * Write probability matrix and associated energies to disk.
 *
 * \param filename Output file name.
 * \param binaryMode If true, values are stored in binary format instead of plain text. Binary format takes less space, but plain text is human-readable.
 *
 */
void GreenFunctionMatrix::write(const std::string& filename, bool binaryMode) const
{
	std::ios::openmode mode = std::ios::out;
	if (binaryMode)
		mode |= std::ios::binary;

	std::ofstream file(filename.c_str(), mode);

	if (file.is_open())
	{
		if (binaryMode)
		{
			std::cout << "Writing binary Green function matrix: " << filename << std::endl;
			size_t ekin_size_toa = fEkinTOA.size();
			size_t ekin_size_lis = fEkinLIS.size();
			file.write(reinterpret_cast<char*>(&ekin_size_toa), sizeof(ekin_size_toa));
			file.write(reinterpret_cast<char*>(&ekin_size_lis), sizeof(ekin_size_lis));

			for (size_t i = 0; i < ekin_size_toa; ++i)
			{
				float ek = fEkinTOA[i] / GeV;
				file.write(reinterpret_cast<char*>(&ek), sizeof(ek));
			}

			for (size_t i = 0; i < ekin_size_lis; ++i)
			{
				float ek = fEkinLIS[i] / GeV;
				file.write(reinterpret_cast<char*>(&ek), sizeof(ek));
			}

			for (size_t i = 0; i < ekin_size_toa; ++i)
			{
				for (size_t j = 0; j < ekin_size_lis; ++j)
				{
					float p = probability(i, j);
					file.write(reinterpret_cast<char*>(&p), sizeof(p));
				}
			}
		}
		else
		{
			std::cout << "Writing Green function matrix: " << filename << std::endl;
			file << *this;
		}
	}

	file.close();
}

/**
 * Stream operator.
 */
std::ostream& operator<<(std::ostream& o, const GreenFunctionMatrix& gr)
{
	const auto& ekinTOA = gr.fEkinTOA;
	const auto& ekinLIS = gr.fEkinLIS;

	for (unsigned int i = 0; i < ekinTOA.size(); i++)
		o << std::setw(12) << std::setprecision(8) << ekinTOA[i] / GeV << " ";
	o << std::endl;

	for (unsigned int j = 0; j < ekinLIS.size(); j++)
		o << std::setw(12) << std::setprecision(8) << ekinLIS[j] / GeV << " ";
	o << std::endl;

	for (unsigned int i = 0; i < ekinTOA.size(); i++)
	{
		for (unsigned int j = 0; j < ekinLIS.size(); j++)
		{
			o << std::setw(12) << std::setprecision(8) << gr.probability(i, j) << " ";
		}
		o << std::endl;
	}

	return o;
}

/**
 * Make sure that the externally provided vectors of energies are consistent with the energy values used when
 * the probability matrix stored in this class was calculated.
 *
 * \param otherEkinTOA The vector of energies to be used for the spectrum at the top of the atmosphere.
 * \param otherEkinLIS The vector of energies to be used for the local interstellar spectrum.
 */
bool GreenFunctionMatrix::assert_ekin_consistency(const std::vector<double>& otherEkinTOA, const std::vector<double>& otherEkinLIS) const
{
	if (fEkinTOA.size() != otherEkinTOA.size())
		throw std::runtime_error("Number of TOA Ekin values does not match!");
	if (fEkinLIS.size() != otherEkinLIS.size())
		throw std::runtime_error("Number of LIS Ekin values does not match!");

	for (unsigned int i = 0; i < fEkinTOA.size(); ++i)
	{
		if (std::abs(fEkinTOA[i] - otherEkinTOA[i]) / (0.5 * (fEkinTOA[i] + otherEkinTOA[i])) > 1.e-6)
		{
			std::stringstream s;
			s << "Binning mismatch in TOA Ekin values at index " << i << ": " << std::setprecision(12) << fEkinTOA[i] / GeV << " vs "
			  << otherEkinTOA[i] / GeV << " GeV!";
			throw std::runtime_error(s.str());
		}
	}

	for (unsigned int i = 0; i < fEkinLIS.size(); ++i)
	{
		if (std::abs(fEkinLIS[i] - otherEkinLIS[i]) / (0.5 * (fEkinLIS[i] + otherEkinLIS[i])) > 1.e-6)
		{
			std::stringstream s;
			s << "Binning mismatch in LIS Ekin values at index " << i << ": " << std::setprecision(12) << fEkinLIS[i] / GeV << " vs "
			  << otherEkinLIS[i] / GeV << " GeV!";
			throw std::runtime_error(s.str());
		}
	}

	return true;
}

/** Check if probability matrix is empty. */
bool GreenFunctionMatrix::empty() const
{
	return fEkinTOA.empty();
}
