#include "helmod.h"

#include "numdiff.h"
#include "particle.h"
#include "units.h"

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace Units;

/**
 * Constructor
 *
 * @param modelInfo %Model information
 */
Helmod::Helmod(const ModelInformation& modelInfo) :
	Model(modelInfo)
{
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;
	kappa_r_theta = 0.0;

	if (modelInfo.hasOption("dt"))
	{
		dt_max = modelInfo.getFloatOption("dt") * s;
	}
	else
	{
		dt_max = 1000.*s;
	}

	useModelSpecificTimestep = modelInfo.getBooleanOptionWithDefault("modelSpecificTimestep", false, false);
	simple_hcs_drift = modelInfo.getBooleanOptionWithDefault("simpleHcsDrift", false, false);

	forceNumericalDerivatives = modelInfo.getBooleanOptionWithDefault("forceNumericalDerivatives", false, false);

	if (forceNumericalDerivatives)
	{
		std::cout << "[Helmod] Forcing numerical calculation of derivatives in diffusion and drift terms..." << std::endl;
	}
}

/**
 * To account for the propagation of changes in heliospheric conditions with the solar wind,
 * the heliosphere is divided into 15 equidistant sectors. Each consecutive sector goes back one
 * Carrington rotation in time, and the relevant model parameters are extracted for the
 * point in time associated with each sector accordingly.
 *
 * @param p Pseudo-particle.
 *
 * @return Sector index.
 */
unsigned int Helmod::getSector(const Particle& p) const
{
	if (p.r < 0.)
		return 0;

	auto nSectors = v_B.size();
	unsigned int iSector = std::floor(p.r / heliosphereBoundary * nSectors);
	if (iSector >= nSectors)
		iSector = nSectors - 1;

	return iSector;
}

double Helmod::solarWindSpeed(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(getSector(p));
	return B.solarWindSpeed(p);
}


/**
 * Calculate the quantities needed several times in one particle step.
 */
void Helmod::calculate(const Particle& p)
{
	currentSector = getSector(p);
	ModifiedSpiral& B = v_B.at(currentSector);

	rho = v_rho.at(currentSector);
	driftFactor = v_driftFactor.at(currentSector);

	if (!nd && forceNumericalDerivatives)
	{
		nd = std::make_shared<NumericalDerivatives>(this);
	}

	B.calculate(p);

	if (simple_hcs_drift)
	{
		fdrift = B.heavi;
		fdriftprime = 0.;
	}
	else
	{
		std::tie(fdrift, fdriftprime) = calculateDriftFactors(p);
	}
}

double Helmod::getDriftR(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);

	if (forceNumericalDerivatives)
		return fdrift * nd->calculateRadialGradientCurvatureDriftVelocity(p, false) / B.heavi;

	// f(theta) smears out heaviside function and has values from 1 to -1, so we have to cancel the heaviside factor used in the vanilla modified Parker spiral
	return fdrift * B.getDriftVelocity_r(p) / B.heavi;
}

double Helmod::getDriftTheta(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);

	if (forceNumericalDerivatives)
		return fdrift * nd->calculatePolarGradientCurvatureDriftVelocity(p) / B.heavi;

	// see getDriftR() for explanation
	return fdrift * B.getDriftVelocity_theta(p) / B.heavi;
}

double Helmod::getDriftSheetR(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);

	if (simple_hcs_drift)
	{
		double alpha = B.tiltAngle/rad;
		double deltaTheta = 2.*p.absRigidity*B.solarWindSpeed(p) / (B.BfieldNormalization * B.r0_sqr * B.OmegaSun * cos(alpha));

		// Look if the particle is inside the "cone" of the wavy heliospheric current sheet
		if (std::abs(p.getTheta() - pi/2.) < alpha+deltaTheta)
		{
			return B.polarity*p.charge * p.velocity/6.0 * cos(alpha) * deltaTheta/sin(alpha+deltaTheta);
		}
		return 0.;
	}

	return fdriftprime * p.rigidity * p.velocity / 3.0 / p.r * B.getBOverB2_phi(p) / B.heavi;
}

double Helmod::getBOverB2_r(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);
	return B.getBOverB2_r(p);
}

double Helmod::getBOverB2_theta(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);
	return B.getBOverB2_theta(p);
}

double Helmod::getBOverB2_phi(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);
	return B.getBOverB2_phi(p);
}

double Helmod::getB(const Particle& p) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);
	return B.getB(p);
}

int Helmod::getHeaviside(const Particle&) const
{
	const ModifiedSpiral& B = v_B.at(currentSector);
	return B.heavi;
}

/** If the \c useModelSpecificTimestep option is used, set \f$ \mathrm{d}t=r^2/\kappa_{rr} \f$. */
double Helmod::overrideTimeStep(const Particle& p) const
{
	if (useModelSpecificTimestep)
	{
		double r = p.r;
		return r*r / kappa_rr;
	}

	return Model::overrideTimeStep(p);
}

std::pair<double, double> Helmod::calculateDriftFactors(const Particle& p) const
{
	// this is an approximation, we could also use deltaThetaHCS = 2*r_Larmor/r, rLarmor = |R|/|B|, with |B| evaluated at (r, theta=pi/2-tiltAngle),
	// but it was checked that this leads to the same results
	const ModifiedSpiral& B = v_B.at(currentSector);
	Particle pB(p);
	pB.setCosTheta(std::cos(pihalf - B.tiltAngle));
	double Rscale = 0.5 * B.BfieldNormalization * B.OmegaSun * B.r0_sqr / B.solarWindSpeed(pB);
	double deltaThetaHCS = p.absRigidity / Rscale / std::cos(B.tiltAngle);

	double charg = B.tiltAngle + deltaThetaHCS;
	if (charg >= pihalf)
		charg = pihalf;

	double ch = pihalf - 0.5*std::sin(charg);
	if (ch <= pihalf / 2.)
		ch = pihalf / 2. + 0.00001;

	double ah = std::acos(pihalf / ch - 1.0);

	double xi = (1.0 - (2.0 * p.getTheta() / pi)) * std::tan(ah);

	double f = std::atan(xi) / ah;
	double fprime = -2.0 * std::tan(ah) / (pi * ah * (std::pow(xi, 2) + 1.0));

	f *= driftFactor;
	fprime *= driftFactor;

	return std::make_pair(f, fprime);
}
