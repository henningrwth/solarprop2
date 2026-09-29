#ifndef HELMOD_H_
#define HELMOD_H_

#include <vector>

#include "model.h"
#include "modifiedspiral.h"
#include "sde2d.h"

class NumericalDerivatives;
class Particle;

/**
 * Base class with common functions used by all helmod models.
 *
 * For a visualization and numerical investigation of the equations for heliospheric
 * current sheet drift, see \c calculations/helmod_hcs_drift.ipynb.
 */
class Helmod : public Model
{
	public:
	Helmod() = delete;
	Helmod(const ModelInformation&);
	Helmod(const Helmod&) = default;
	virtual ~Helmod() override {};
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	double solarWindSpeed(const Particle& p) const override;

	void calculate(const Particle& particle) override;

	double getKappaThetaTheta() const override { return kappa_thetatheta; }
	double getKappaRR() const override { return kappa_rr; }
	double getKappaRTheta() const override { return kappa_r_theta; }

	double getDriftR(const Particle& p) const override;
	double getDriftTheta(const Particle& p) const override;

	double getDriftSheetR(const Particle& p) const override;

	//
	// Functions for numerical drift calculations (only needed if forceNumericalDerivatives option is used):
	//

	double getBOverB2_r(const Particle& p) const override;
	double getBOverB2_theta(const Particle& p) const override;
	double getBOverB2_phi(const Particle& p) const override;

	double getB(const Particle& p) const override;

	int getHeaviside(const Particle&) const override;

	double overrideTimeStep(const Particle& p) const override;

	protected:

	/** Access to heliospheric sector depending on particle position. */
	unsigned int getSector(const Particle& p) const;

	/** Calculate transition function \f$ f(\theta) \f$ for drift effects and its derivative with respect to \f$ \theta \f$. */
	std::pair<double, double> calculateDriftFactors(const Particle& p) const;

	/// Magnetic field (one for each heliosphere sector)
	std::vector<ModifiedSpiral> v_B;

	/// current heliosphere sector
	unsigned int currentSector;

	/// Normalization of the diffusion tensor (one for each heliosphere sector)
	std::vector<double> v_kappa0;

	/// Diffusion coefficient \f$ \kappa_{rr} \f$ (heliospheric coordinates)
	double kappa_rr;

	/// Diffusion coefficient \f$ \kappa_{\theta\theta} \f$ (heliospheric coordinates)
	double kappa_thetatheta;

	/// Diffusion coefficient \f$ \kappa_{r\theta} \f$ (heliospheric coordinates)
	double kappa_r_theta;

	/// enhancement factor for polar perpendicular diffusion coefficient \f$ K_{\perp\theta} \f$ in polar region (see eq. (8) in Bobik et al. (2012))
	double iotaPolarRegion;

	/// current factor between polar and radial perpendicular diffusion coefficients
	double iota;

	/// \f$ \rho \f$: factor between radial perpendicular and parallel diffusion coefficients (one for each heliosphere sector)
	std::vector<double> v_rho;

	/// value of \f$ \rho \f$ for current sector
	double rho;

	/// Softening term \f$ g_\mathrm{low} \f$ for diffusion at low rigidities (one for each heliosphere sector)
	std::vector<double> v_glow;

	/// driftFactor: drift suppression factor (unity for normal unsuppressed drifts)
	std::vector<double> v_driftFactor;

	/// value of `driftFactor` for current sector
	double driftFactor;

	/// rigidity dependent term of diffusion coefficient
	double KP;

	/// transition function \f$ \tilde{f}(\theta) \f$ for drift effects
	double fdrift;

	/// \f$ \mathrm{d}\tilde{f}(theta)/\mathrm{d}\theta \f$
	double fdriftprime;

	/// use simplified HCS drift scheme: average velocity according to BP model?
	bool simple_hcs_drift = false;

	/// use model-specific time step?
	bool useModelSpecificTimestep;

	/// Always use numerical calculation of derivatives? Useful for debugging and testing.
	bool forceNumericalDerivatives = false;

	/// NumericalDerivatives used for numerical calculation of derivatives.
	std::shared_ptr<NumericalDerivatives> nd = nullptr;
};

#endif
