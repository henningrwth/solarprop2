#include "sde2d.h"

#include "model.h"
#include "random.h"
#include "units.h"

#include <cassert>
#include <iomanip>
#include <iostream>

using namespace Units;

/**
 * Constructor.
 *
 * @param m The model defining the terms in the %SDE.
 * @param _globalSeed Global seed for random numbers. The unique particle number will be added to this to obtain a unique random generator for each particle.
 *
 */
SDE2d::SDE2d(Model* m, unsigned int _globalSeed) :
	SDE(m, _globalSeed)
{
}

/**
 * Define rules for adjusting time step if necessary.
 *
 * If the `dynamicStep` option is used, the time step is calculated such that \f$ r \f$ changes roughly by 5%, or by 0.05 AU if particle is close to the Sun.
 * In any case, the time step will be reduced close to the Sun to avoid large steps and thus the \f$ 1/r \f$ divergence in the calculation of \f$ \mathrm{d}p \f$,
 * and close to the poles to avoid numerical instability.
 */
double SDE2d::adjusted_dt(const Particle& p, double v, double driftR, double driftsheetR, double linearDiffusionTermRR, double kappaRR,
			  double linearCoefficient) const
{
	const double dt_max = model->getMaximumDeltaTime();
	double new_dt = dt_max;

	double r = p.getR();

	// aim for dr given by relative change, or absolute change if close to the Sun
	if (model->getDynamicStep())
	{
		double dr_goal = r < 1.0 * AU ? 0.05 * AU : 0.05 * r;
		double dt_v = dr_goal / std::abs(v);
		double dt_drift = dr_goal / std::abs(driftR);
		double dt_diffusion = dr_goal / std::abs(linearDiffusionTermRR);
		double dt_diffusion2 = dr_goal * dr_goal / 2.0 / std::abs(kappaRR);
		new_dt = dt_v;
		if (dt_drift < new_dt)
			new_dt = dt_drift;
		if (dt_diffusion < new_dt)
			new_dt = dt_diffusion;
		if (dt_diffusion2 < new_dt)
			new_dt = dt_diffusion2;
		if (new_dt < 300.0 * s)
			new_dt = 300.0 * s;
	}

	//
	// limit dt if necessary to avoid numerical instabilities
	//

	// adjust time step close to Sun to avoid large steps and thus 1/r divergence in deltaP
	if (r < 1.0 * AU)
	{
		static constexpr double sMaxRelChange = 0.05;
		double dr1 = linearCoefficient * new_dt;
		if (dr1 / r > sMaxRelChange)
		{
			new_dt = sMaxRelChange * r / linearCoefficient;
		}

		double dr2 = std::sqrt(2. * kappaRR * new_dt);

		if (dr2 / r > sMaxRelChange)
		{
			new_dt = sMaxRelChange * sMaxRelChange * r * r / 2.0 / kappaRR;
		}
	}

	// avoid numerical instability near the poles
	if (std::abs(p.getCosTheta()) > 0.997)
	{
		static constexpr double sMaxChange = 1.0 * AU;
		double dr1 = linearCoefficient * new_dt;
		if (dr1 > sMaxChange)
		{
			new_dt = sMaxChange / linearCoefficient;
		}
	}

	// make sure dt never exceeds predefined maximum value
	if (new_dt > dt_max)
		new_dt = dt_max;

	return new_dt;
}

/**
 * Calculate \f$ \mathrm{d}r \f$, the step in \f$ r \f$ direction, and update \f$ \mathrm{d}t \f$.
 *
 * \attention This function will set the time increment \f$ \mathrm{d}t \f$. So make sure it is called before all other \c getDelta functions.
 */
double SDE2d::getDeltaR(const Particle& p)
{
	double kappaRR = model->getKappaRR();

	//
	// get a proposal for dt from fastest propagation process
	//
	double v = model->solarWindSpeed(p);
	double driftR = model->getDriftR(p);
	double driftsheetR = model->getDriftSheetR(p);
	double linearDiffusionTermRR = model->getLinearRadialDiffusionTermFromKappaRR(p);
	double linearDiffusionTermRTheta = model->getLinearRadialDiffusionTermFromKappaRTheta(p);
	double linearCoefficient = -v - driftR - driftsheetR + linearDiffusionTermRR + linearDiffusionTermRTheta;

	dt = model->overrideTimeStep(p);
	if (dt == 0.0)
		dt = adjusted_dt(p, v, driftR, driftsheetR, linearDiffusionTermRR, kappaRR, linearCoefficient);

	double kappaRTheta = model->getKappaRTheta();
	double kappaThetaTheta = model->getKappaThetaTheta();

	double dr = linearCoefficient * dt;

	if (kappaRTheta != 0.0)
	{
		dr += std::sqrt((kappaRR * kappaThetaTheta - kappaRTheta * kappaRTheta) / (0.5 * kappaThetaTheta) * dt) * wr;
		dr -= kappaRTheta * std::sqrt(2.0 / kappaThetaTheta * dt) * wx;
	}
	else
	{
		dr += std::sqrt(2. * kappaRR * dt) * wr;
	}

	return dr;
}

/** Calculate \f$ \mathrm{d}\cos\theta \f$ for the current step. */
double SDE2d::getDeltaCosTheta(const Particle& p) const
{

	double r = p.getR();
	double kappaThetaTheta = model->getKappaThetaTheta();

	double linearPolarDiffusionTermRTheta = model->getLinearPolarDiffusionTermFromKappaRTheta(p);
	double linearPolarDiffusionTermThetaTheta = model->getLinearPolarDiffusionTermFromKappaThetaTheta(p);

	double driftTheta = model->getDriftTheta(p) + model->getDriftSheetTheta(p);

	double sinTheta = p.getSinTheta();
	double dmu = (driftTheta / r * sinTheta + linearPolarDiffusionTermRTheta + linearPolarDiffusionTermThetaTheta) * dt;
	dmu += std::sqrt(2. * kappaThetaTheta * dt) * std::abs(sinTheta) * wx / r;

	return dmu;
}

/**
 * Calculate the individual contributions of the transport processes for a pseudo-particle at a given position.
 *
 * @param p The pseudo-particle state (position and energy).
 *
 * @return Arrows object.
 */
Arrows SDE2d::calculateArrows(const Particle& p) const
{
	model->calculate(p);

	double r = p.getR();
	double x = p.getX();
	double z = p.getZ();

	double v = model->solarWindSpeed(p);
	double driftR = model->getDriftR(p);
	double driftTheta = model->getDriftTheta(p);
	double driftsheetR = model->getDriftSheetR(p);
	double driftsheetTheta = model->getDriftSheetTheta(p);
	double linearDiffusionTermRR = model->getLinearRadialDiffusionTermFromKappaRR(p);
	double linearDiffusionTermRTheta = model->getLinearRadialDiffusionTermFromKappaRTheta(p);

	double dt = model->getMaximumDeltaTime();

	double kappaRR = model->getKappaRR();
	double kappaRTheta = model->getKappaRTheta();
	double kappaThetaTheta = model->getKappaThetaTheta();

	double linearPolarDiffusionTermRTheta = model->getLinearPolarDiffusionTermFromKappaRTheta(p);
	double linearPolarDiffusionTermThetaTheta = model->getLinearPolarDiffusionTermFromKappaThetaTheta(p);
	double sinTheta = p.getSinTheta();

	Arrows a(p);
	a.dt = dt;

	Particle p_convection(p);
	double dr_convection = -v * dt;
	p_convection.setR(p.getR() + dr_convection);
	a.deltaConvection = {p_convection.getX() - x, p_convection.getZ() - z};

	Particle p_drift(p);
	double dr_drift = -driftR * dt;
	double dmu_drift = driftTheta / r * sinTheta * dt;
	p_drift.setR(p.getR() + dr_drift);
	p_drift.setCosTheta(p.getCosTheta() + dmu_drift);
	a.deltaDrift = {p_drift.getX() - x, p_drift.getZ() - z};

	Particle p_hcs(p);
	double dr_hcsdrift = -driftsheetR * dt;
	double dmu_hcsdrift = driftsheetTheta / r * sinTheta * dt;
	p_hcs.setR(p.getR() + dr_hcsdrift);
	p_hcs.setCosTheta(p.getCosTheta() + dmu_hcsdrift);
	a.deltaHcsDrift = {p_hcs.getX() - x, p_hcs.getZ() - z};

	Random rnd(globalSeed);

	// diffusion
	for (unsigned int i = 0; i < 12; ++i)
	{
		double wr = rnd.getGaussianRandomNumber();
		double wx = rnd.getGaussianRandomNumber();

		Particle p_diff(p);

		double dr_diff = (linearDiffusionTermRR + linearDiffusionTermRTheta) * dt;
		if (kappaRTheta != 0.0)
		{
			dr_diff += std::sqrt((kappaRR * kappaThetaTheta - kappaRTheta * kappaRTheta) / (0.5 * kappaThetaTheta) * dt) * wr;
			dr_diff -= kappaRTheta * std::sqrt(2.0 / kappaThetaTheta * dt) * wx;
		}
		else
		{
			dr_diff += std::sqrt(2. * kappaRR * dt) * wr;
		}

		double dmu_diff = (linearPolarDiffusionTermRTheta + linearPolarDiffusionTermThetaTheta) * dt;
		dmu_diff += std::sqrt(2. * kappaThetaTheta * dt) * std::abs(sinTheta) * wx / r;

		p_diff.setR(p.getR() + dr_diff);
		p_diff.setCosTheta(p.getCosTheta() + dmu_diff);

		a.deltasDiffusion.emplace_back(p_diff.getX() - x, p_diff.getZ() - z);
	}

	return a;
}

void SDE2d::generateRandomNumbers(Random& rndm)
{
	wr = rndm.getGaussianRandomNumber();
	wx = rndm.getGaussianRandomNumber();
}
