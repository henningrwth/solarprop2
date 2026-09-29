#include "helmod2018.h"

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
Helmod2018::Helmod2018(const ModelInformation& modelInfo) :
	Helmod(modelInfo)
{
	FACTORY_INIT;

	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;
	kappa_r_theta = 0.0;

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();

	heliosphereBoundary = 100.*AU;

	auto nominalTimestamp = modelInfo.getTimestamp();
	auto nominalCR = modelInfo.getCarringtonRotation();

	double highSolarActivityThreshold = modelInfo.getFloatOptionWithDefault("highSolarActivityThreshold", 150.0);

	double BfieldScaling = modelInfo.getFloatOptionWithDefault("BfieldScaling", 1.0, false);

	iotaPolarRegion = modelInfo.getFloatOptionWithDefault("iotaPolarRegion", 10.0, false);

	r_e = 1.*AU;

	double delta0 = modelInfo.getFloatOptionWithDefault("delta0", 2.e-5, false);

	int verbosity = modelInfo.getIntegerOptionWithDefault("verbosity", 3, false);

	// divide heliosphere into 15 sectors
	const unsigned int nSectors = 15;
	for (unsigned int iSector = 0; iSector < nSectors; ++iSector)
	{
		// each sector goes back 1 Carrington rotation in time
		auto timestamp = nominalTimestamp - iSector * CarringtonRotationDuration;
		auto CR = nominalCR - iSector;

		double tiltAngle = 0.;
		if (modelInfo.hasOption("angle"))
			tiltAngle = modelInfo.getFloatOption("angle") * degree;
		else
			tiltAngle = modelInfo.getTiltAngleFor(CR, false);

		int polarity = 0;
		int phase = 0;
		if (!modelInfo.hasOption("polarity") || !modelInfo.hasOption("phase"))
			std::tie(polarity, phase) = modelInfo.getPolarityAndPhaseFor(CR, false);

		if (modelInfo.hasOption("polarity"))
			polarity = modelInfo.getIntegerOption("polarity");
		if (modelInfo.hasOption("phase"))
			phase = modelInfo.getIntegerOption("phase");

		double ssn = modelInfo.getSunspotNumberFor(timestamp);
		bool isHighSolarActivity = (ssn > highSolarActivityThreshold);

		double Bearth = 0.;
		double vswMin = 0.;
		if (!modelInfo.hasOption("vswMin") || !modelInfo.hasOption("Bearth"))
			std::tie(Bearth, vswMin) = modelInfo.getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebDataFor(timestamp);

		if (modelInfo.hasOption("vswMin"))
			vswMin = modelInfo.getFloatOption("vswMin") * km/s;
		if (modelInfo.hasOption("Bearth"))
			Bearth = modelInfo.getFloatOption("Bearth") * nanotesla;

		Bearth *= BfieldScaling;

		double kappaScaling = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("kappaScaling", 1.0, timestamp, false);

		double kappa0 = kappaScaling * modelInfoKappa0;
		if (kappa0 == 0.0)
		{
			kappa0 = kappaScaling * defaultKappaZeroFromSolarPolarityAndPhase(polarity, phase, ssn);
		}

		v_B.emplace_back(modelInfo, isHighSolarActivity, polarity, phase, tiltAngle, Bearth, vswMin, delta0);
		v_kappa0.push_back(kappa0);

		rho = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("rho", 0.06, timestamp, false);
		v_rho.push_back(rho);

		driftFactor = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("driftFactor", 1.0, timestamp, false);
		v_driftFactor.push_back(driftFactor);

		double g_low = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("glow", 0.3, timestamp, false);
		v_glow.push_back(g_low);

		if (verbosity >= 3 || iSector == 0)
		{
			time_t ttt(timestamp);
			std::cout << "[Helmod2018] sector " << iSector << "/" << nSectors << " for CR " << CR << " "
				  << "(time: " << std::put_time(std::gmtime(&ttt), "%d %b %Y") << "): "
				  << "A=" << polarity << ", phase=" << phase << ", "
				  << "tilt angle: " << tiltAngle/deg << " deg, "
				  << "kappa0= " << kappa0/(1.e-4*AU*AU/s) << " x10^(-4)AU^2/s = " << kappa0/(cm2/s) << " cm^2/s, "
				  << "SSN: " << ssn << ", "
				  << "g_low: " << g_low << ", "
				  << "B_earth = " << Bearth / BfieldScaling / nanotesla << " nT (x " << BfieldScaling << "), "
				  << "vsw = " << vswMin / (km/s) << " km/s, "
				  << "high activity: " << isHighSolarActivity << ", "
				  << "rho: " << rho << ", "
				  << "iota: " << iotaPolarRegion << ", "
				  << "drift factor: " << driftFactor
				  << std::endl;
		}
	}
}

/**
 * Calculate the diffusion coefficient scale \f$ \kappa_0 \f$ from equation (6) and Table 1 of the reference paper.
 *
 * @param polarity Polarity of heliospheric magnetic field.
 * @param phase Phase (`+1` for ascending, `-1` for descending)
 * @param ssn Smoothed sunspot number.
 *
 * @return Diffusion coefficient scale \f$ \kappa_0 \f$.
 */
double Helmod2018::defaultKappaZeroFromSolarPolarityAndPhase(int polarity, int phase, double ssn)
{
	double c0 = 0., c1 = 0., c2 = 0., c3 = 0.;

	if (polarity > 0)
	{
		if (phase > 0)
		{
			c0 = 0.0002262;
			c1 = -5.058e-7;
			c2 = 0.;
			c3 = 0.;
		}
		else
		{
			c0 = 0.0002267;
			c1 = -7.118e-7;
			c2 = 0.;
			c3 = 0.;
		}
	}
	else
	{
		if (phase > 0)
		{
			c0 = 0.0003059;
			c1 = -2.51e-6;
			c2 = 1.284e-8;
			c3 = -2.838e-11;
		}
		else
		{
			c0 = 0.0002876;
			c1 = -3.715e-6;
			c2 = 2.534e-8;
			c3 = -5.689e-11;
		}
	}

	double kappa0_scale = c0 + c1*ssn + c2*ssn*ssn + c3*std::pow(ssn, 3);
	double unit = AU*AU/s;
	return kappa0_scale * unit;
}


/**
 * Calculate the quantities needed several times in one particle step.
 */
void Helmod2018::calculate(const Particle& p)
{
	Helmod::calculate(p);

	const ModifiedSpiral& B = v_B.at(currentSector);

	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);

	iota = 1.0;
	// cos(30deg) = sqrt(3)/2 = 0.8660254
	if (std::abs(p.cosTheta) > 0.8660254) // polar region: theta < 30deg or theta > 150 deg
		iota = iotaPolarRegion;

	double g_low = v_glow.at(currentSector);
	KP = p.absRigidity / (GV/c) + g_low;

	double kappa0 = v_kappa0.at(currentSector);
	fdiffusion = p.beta/3.0 * kappa0 * KP;

	double denominator = B.N*r_e*B.tau;
	kappa_rr = fdiffusion*(p.r + r_e)*(B.N*iota*rho*t2 + Gamma2*rho + B.tau) / denominator;
	kappa_r_theta = fdiffusion*B.heavi*B.t*(p.r + r_e)*(-B.N*iota*rho + Gamma2*rho + B.tau) / denominator;
	kappa_thetatheta = fdiffusion*(p.r + r_e)*(B.N*iota*rho + t2*(Gamma2*rho + B.tau)) / denominator;
}

double Helmod2018::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);
	double N2 = std::pow(B.N, 2);
	double tau2 = std::pow(B.tau, 2);

	return -fdiffusion*(2*B.N*t2*(p.r + r_e)*(B.N*iota*rho*t2 + Gamma2*rho + B.tau)
			    - B.N*B.tau*(p.r*(B.N*iota*rho*t2 + Gamma2*rho + B.tau) + 2*(p.r + r_e)*(B.N*iota*rho*t2 + Gamma2*rho + B.tau) + 2*(p.r + r_e)*(B.N*iota*rho*t2 + Gamma2*rho + iota*rho*t2*(Gamma2 + t2) + t2))
			    + 2*B.tau*(Gamma2 + t2)*(p.r + r_e)*(B.N*iota*rho*t2 + Gamma2*rho + B.tau)) / (N2*p.r*r_e*tau2);
}

double Helmod2018::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double r2 = std::pow(p.r, 2);
	double Gamma2 = std::pow(B.Gamma, 2);
	double mu = p.cosTheta;
	double sintheta2 = std::pow(p.sinTheta, 2);
	double tprime = B.dtdmu;
	double ttprime = B.t*tprime;
	double t2 = std::pow(B.t, 2);
	double N2 = std::pow(B.N, 2);
	double tau2 = std::pow(B.tau, 2);
	double GGp = B.Gammaprime*B.Gamma;

	return -2.*fdiffusion*(p.r + r_e)*(B.N*ttprime*(B.N*iota*rho + t2*(Gamma2*rho + B.tau))*sintheta2
					   + B.N*B.tau*(mu*(B.N*iota*rho + t2*(Gamma2*rho + B.tau)) - (iota*rho*(GGp + ttprime) + t2*(GGp*rho + ttprime) + ttprime*(Gamma2*rho + B.tau))*sintheta2)
					   + B.tau*(GGp + ttprime)*(B.N*iota*rho + t2*(Gamma2*rho + B.tau))*sintheta2) / (N2*r2*r_e*tau2);
}

double Helmod2018::getLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double mu = p.cosTheta;
	double sintheta2 = std::pow(p.sinTheta, 2);
	double t2 = std::pow(B.t, 2);
	double N2 = std::pow(B.N, 2);
	double tau2 = std::pow(B.tau, 2);
	double tprime = B.dtdmu;
	double ttprime = B.t*tprime;
	double absSinTheta = std::abs(p.sinTheta);
	double GGp = B.Gammaprime*B.Gamma;

	return fdiffusion*B.heavi*(p.r + r_e)*(B.N*mu*B.t*B.tau*(-B.N*iota*rho + Gamma2*rho + B.tau)
					       + 2*B.N*t2*tprime*(-B.N*iota*rho + Gamma2*rho + B.tau)*sintheta2 - B.N*B.tau*(2*B.t*(GGp*rho - iota*rho*(GGp + ttprime) + ttprime) + tprime*(-B.N*iota*rho + Gamma2*rho + B.tau))*sintheta2
					       + 2*B.t*B.tau*(GGp + ttprime)*(-B.N*iota*rho + Gamma2*rho + B.tau)*sintheta2) / (N2*p.r*r_e*tau2*absSinTheta);
}

double Helmod2018::getLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaRTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double r2 = std::pow(p.r, 2);
	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);
	double N2 = std::pow(B.N, 2);
	double tau2 = std::pow(B.tau, 2);
	double absSinTheta = std::abs(p.sinTheta);

	return fdiffusion*B.heavi*B.t*absSinTheta*(2*B.N*t2*(p.r + r_e)*(-B.N*iota*rho + Gamma2*rho + B.tau)
						   - B.N*B.tau*(p.r*(-B.N*iota*rho + Gamma2*rho + B.tau) + 2*(p.r + r_e)*(Gamma2*rho - iota*rho*(Gamma2 + t2) + t2) + 2*(p.r + r_e)*(-B.N*iota*rho + Gamma2*rho + B.tau))
						   + 2*B.tau*(Gamma2 + t2)*(p.r + r_e)*(-B.N*iota*rho + Gamma2*rho + B.tau)) / (N2*r2*r_e*tau2);
}
