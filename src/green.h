#ifndef GreenFunctionMatrix_h_
#define GreenFunctionMatrix_h_

#include <string>
#include <vector>

#include "density.h"

/**
 * Green function matrix that is the main result of a simulation run.
 *
 * The entries \f$ G(T_i,T_j) \f$ of the Green matrix are calculated
 * from the probability for a particle in kinetic energy bin \f$ T_i \f$
 * at Earth to originate from the kinetic energy bin \f$ T_j \f$ at the
 * heliopause.
 *
 * The local interstellar spectrum \f$ \Phi_\mathrm{LIS}(T) \f$,
 * i.e. the differential flux of particles as a function of kinetic
 * energy \f$ T \f$ at the heliopause, provides the boundary condition
 * for the solution of the transport equation. The modulated flux at
 * Earth is \f$ \Phi_\mathrm{mod}(T)=p^2(T)f_\mathrm{mod}(T) \f$, with
 * the phase-space density calculated as
 * \f[
 *   f_\mathrm{mod}(T_i)=\sum\limits_{j}G(T_i,T_j)f_\mathrm{LIS}(T_j),
 * \f]
 * where \f$ f_\mathrm{LIS}(T)=\Phi_\mathrm{LIS}(T)/p^2(T) \f$, and
 * \f$ p^2(T)=T^2+2mT \f$, where \f$ m \f$ is the particle mass.
 */
class GreenFunctionMatrix
{
	public:
	GreenFunctionMatrix(const std::vector<double>& ekinTOA, const std::vector<double>& ekinLIS);
	GreenFunctionMatrix(const std::string& infilename, bool binaryMode);
	GreenFunctionMatrix() = default;
	GreenFunctionMatrix(const GreenFunctionMatrix&) = default;
	~GreenFunctionMatrix() {}

	void reset();

	void fillFromDensity(const Density& den);

	double probability(unsigned int ekinIndexTOA, unsigned int ekinIndexLIS) const;

	bool assert_ekin_consistency(const std::vector<double>& otherEkinTOA, const std::vector<double>& otherEkinLIS) const;
	bool empty() const;

	void write(const std::string& filename, bool binaryMode) const;
	friend std::ostream& operator<<(std::ostream& o, const GreenFunctionMatrix& gr);

	private:

	/// probability matrix \f$ G(T_i,T_j) \f$ containing the probabilities for a particle in kinetic energy bin \f$ T_i \f$
	/// at Earth to originate from the kinetic energy bin \f$ T_j \f$ at the heliopause
	std::vector<std::vector<double>> fProbabilityMatrix;

	/// the local energies corresponding to the rows in the matrix
	std::vector<double> fEkinTOA;

	/// the LIS energies corresponding to the columns in the matrix
	std::vector<double> fEkinLIS;
};

std::ostream& operator<<(std::ostream& o, const GreenFunctionMatrix& gr);

#endif
