#include "sde3d.h"

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
SDE3d::SDE3d(Model* m, unsigned int _globalSeed) :
	SDE(m, _globalSeed)
{
}

/** Define rules for adjusting time step if necessary. */
double SDE3d::adjusted_dt(const Particle&) const
{
	// so far, this is only a placeholder,
	// to make sure particles do not get stuck
	return model->getMaximumDeltaTime();
}

/** Generate the shift in r direction and set the time step.
 */
double SDE3d::getDeltaR(const Particle& p)
{
	double kappaRR = model->getKappaRR();

	double v = model->solarWindSpeed(p);
	double driftR = model->getDriftR(p);
	double driftsheetR = model->getDriftSheetR(p);
	double linearDiffusionTermRR = model->getLinearRadialDiffusionTermFromKappaRR(p);
	double linearDiffusionTermRPhi = model->getLinearRadialDiffusionTermFromKappaRPhi(p);
	double linearCoefficient = -v - driftR - driftsheetR + linearDiffusionTermRR + linearDiffusionTermRPhi;

	dt = model->overrideTimeStep(p);

	// FIXME adjust time step??
	if (dt == 0.0)
		dt = adjusted_dt(p);

	double kappaPhiPhi = model->getKappaPhiPhi();
	double kappaRPhi = model->getKappaRPhi();

	double dr = linearCoefficient * dt;

	dr += std::sqrt(2. * (kappaRR - kappaRPhi * kappaRPhi / kappaPhiPhi) * dt) * w_r;

	dr += kappaRPhi * std::sqrt(2. / kappaPhiPhi * dt) * w_phi;

	return dr;
}

/** Generate step in cos(theta). */
double SDE3d::getDeltaCosTheta(const Particle& p) const
{

	double r = p.getR();
	double kappaThetaTheta = model->getKappaThetaTheta();

	double linearPolarDiffusionTermThetaTheta = model->getLinearPolarDiffusionTermFromKappaThetaTheta(p);

	double driftTheta = model->getDriftTheta(p) + model->getDriftSheetTheta(p);

	double sinTheta = p.getSinTheta();
	double dmu = (driftTheta / r * sinTheta + linearPolarDiffusionTermThetaTheta) * dt;
	dmu += std::sqrt(2. * kappaThetaTheta * dt) * std::abs(sinTheta) * w_mu / r;

	return dmu;
}

/** Generate step in phi. */
double SDE3d::getDeltaPhi(const Particle& p) const
{
	double r = p.getR();
	double sinTheta = p.getSinTheta();

	double kappaPhiPhi = model->getKappaPhiPhi();

	double linearAzimuthalDiffusionTermRPhi = model->getLinearAzimuthalDiffusionTermFromKappaRPhi(p);
	double linearAzimuthalDiffusionTermPhiPhi = model->getLinearAzimuthalDiffusionTermFromKappaPhiPhi(p);

	double driftPhi = model->getDriftPhi(p) + model->getDriftSheetPhi(p);

	double dphi = (-driftPhi / r / sinTheta + linearAzimuthalDiffusionTermRPhi + linearAzimuthalDiffusionTermPhiPhi) * dt;

	dphi += std::sqrt(2. * kappaPhiPhi * dt) * w_phi / r / sinTheta;

	return dphi;
}

/**
 * Calculate the individual contributions of the transport processes for a pseudo-particle at a given position.
 *
 * @param p The pseudo-particle state (position and energy).
 *
 * @return Arrows object.
 */
Arrows SDE3d::calculateArrows(const Particle& p) const
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

	double linearRadialDiffusionTermRR = model->getLinearRadialDiffusionTermFromKappaRR(p);
	double linearRadialDiffusionTermRPhi = model->getLinearRadialDiffusionTermFromKappaRPhi(p);

	double dt = model->getMaximumDeltaTime();

	double kappaRR = model->getKappaRR();
	double kappaPhiPhi = model->getKappaPhiPhi();
	double kappaRPhi = model->getKappaRPhi();
	double kappaThetaTheta = model->getKappaThetaTheta();

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
		double w_r = rnd.getGaussianRandomNumber();
		double w_mu = rnd.getGaussianRandomNumber();
		double w_phi = rnd.getGaussianRandomNumber();

		Particle p_diff(p);

		double dr_diff = (linearRadialDiffusionTermRR + linearRadialDiffusionTermRPhi) * dt;
		dr_diff += std::sqrt(2. * (kappaRR - kappaRPhi * kappaRPhi / kappaPhiPhi) * dt) * w_r;
		dr_diff += kappaRPhi * std::sqrt(2. / kappaPhiPhi * dt) * w_phi;

		double dmu_diff = linearPolarDiffusionTermThetaTheta * dt;
		dmu_diff += std::sqrt(2. * kappaThetaTheta * dt) * std::abs(sinTheta) * w_mu / r;

		p_diff.setR(p.getR() + dr_diff);
		p_diff.setCosTheta(p.getCosTheta() + dmu_diff);

		a.deltasDiffusion.emplace_back(p_diff.getX() - x, p_diff.getZ() - z);
	}

	return a;
}

void SDE3d::generateRandomNumbers(Random& rndm)
{
	w_r = rndm.getGaussianRandomNumber();
	w_mu = rndm.getGaussianRandomNumber();
	w_phi = rndm.getGaussianRandomNumber();
}
