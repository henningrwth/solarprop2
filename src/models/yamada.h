#ifndef YAMADA_H_
#define YAMADA_H_

#include <memory>

#include "model.h"
#include "particle.h"
#include "sde2d.h"

/**
 * Implementation of the model of %Yamada et al (1998).
 *
 * Reference: Y. %Yamada et al., Geophys. Res. Lett. Vol. 25 No. 13 (1998) 2353-2356, https://doi.org/10.1029/98GL51869
 *
 * This model aims at reproducing Fig. 4 in the reference paper.
 */
class Yamada : public Model, RegisteredInFactory<Model, Yamada, const ModelInformation&>
{
	public:
	Yamada() = delete;
	Yamada(const ModelInformation&);
	Yamada(const Yamada&) = default;
	~Yamada() {}
	std::unique_ptr<Model> clone() const override { return std::make_unique<Yamada>(*this); }
	std::unique_ptr<SDE> createSDE(unsigned int globalSeed) override { return std::make_unique<SDE2d>(this, globalSeed); }

	void calculate(const Particle& p) override;

	double solarWindSpeed(const Particle& p) const override;

	double getKappaThetaTheta() const override { return 0.0; }
	double getKappaRR() const override { return kappa_rr; }

	double getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const override;

	/** %Model name used by ModelFactory. */
	static std::unique_ptr<Model> createModel(const ModelInformation& modelInfo) { return std::make_unique<Yamada>(modelInfo); }
	/** %Model name used by ModelFactory. */
	static std::string modelName() { return "Yamada1998"; }

	private:

	/// Normalization of the diffusion tensor
	double kappa0;

	/// Diffusion coefficient \f$ \kappa_{rr} \f$ (heliospheric coordinates)
	double kappa_rr;
};
#endif
