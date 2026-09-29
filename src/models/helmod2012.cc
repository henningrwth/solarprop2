#include "helmod2012.h"

#include "numdiff.h"
#include "particle.h"
#include "units.h"

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

using namespace Units;

/**
 * @brief Constructor
 *
 * @param modelInfo %Model information
 */
Helmod2012::Helmod2012(const ModelInformation& modelInfo) :
	Helmod(modelInfo)
{
	FACTORY_INIT;

	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;
	kappa_r_theta = 0.0;

	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0);

	// Override diffusion coefficient?
	auto modelInfoKappa0 = modelInfo.getDiffusionCoefficientScale();

	heliosphereBoundary = 100.*AU;

	auto nominalTimestamp = modelInfo.getTimestamp();
	auto nominalCR = modelInfo.getCarringtonRotation();

	double highSolarActivityThreshold = modelInfo.getFloatOptionWithDefault("highSolarActivityThreshold", 150.0);

	double BfieldScaling = modelInfo.getFloatOptionWithDefault("BfieldScaling", 1.0, false);

	iotaPolarRegion = modelInfo.getFloatOptionWithDefault("iotaPolarRegion", 10.0, false);

	double delta0 = modelInfo.getFloatOptionWithDefault("delta0", 8.7e-5, false);

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

		double kappa0 = kappaScaling * modelInfoKappa0;
		if (kappa0 == 0.0)
		{
			kappa0 = kappaScaling * defaultKappaZeroFromSolarPolarityAndPhase(polarity, phase, ssn);
		}

		v_B.emplace_back(modelInfo, isHighSolarActivity, polarity, phase, tiltAngle, Bearth, vswMin, delta0);
		v_kappa0.push_back(kappa0);

		rho = modelInfo.getFloatOptionWithDefault("rho", 0.05, false);
		v_rho.push_back(rho);

		driftFactor = modelInfo.getFloatOptionWithDefault("driftFactor", 1.0, false);
		v_driftFactor.push_back(driftFactor);

		double g_low = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("glow", 0., timestamp, false);
		v_glow.push_back(g_low);

		if (verbosity >= 3 || iSector == 0)
		{
			time_t ttt(timestamp);
			std::cout << "[Helmod2012] sector " << iSector << "/" << nSectors << " for CR " << CR << " "
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
				  << "drift factor: " << driftFactor << ", "
				  << "delta0: " << delta0
				  << std::endl;
		}
	}
}

/**
 * Calculate the diffusion coefficient scale \f$ \kappa_0 \f$ from equation (13) and Table 1 of the reference paper.
 *
 * @param polarity Polarity of heliospheric magnetic field.
 * @param phase Phase (`+1` for ascending, `-1` for descending)
 * @param ssn Smoothed sunspot number.
 *
 * @return Diffusion coefficient scale \f$ \kappa_0 \f$.
 */
double Helmod2012::defaultKappaZeroFromSolarPolarityAndPhase(int polarity, int phase, double ssn)
{
	double c1 = 0., c2 = 0., c3 = 0., c4 = 0.;

	if (polarity > 0)
	{
		if (phase > 0)
		{
			c1 = 2.39708e-4;
			c2 = 0.;
			c3 = -8.28987e-7;
			c4 = 0.;
		}
		else
		{
			c1 = 2.28037e-4;
			c2 = 0.;
			c3 = -1.00984e-6;
			c4 = 0.;
		}
	}
	else
	{
		if (phase > 0)
		{
			c1 = 0.0001686;
			c2 = 0.001488;
			c3 = 0.;
			c4 = -3.164e-9;
		}
		else
		{
			c1 = 8.872e-5;
			c2 = 0.001874;
			c3 = 0.;
			c4 = 0.;
		}
	}

	double kappa0_scale = c1 + c2/ssn + c3*ssn + c4*ssn*ssn;
	double unit = AU*AU/s;
	return kappa0_scale * unit;
}


void Helmod2012::calculate(const Particle& p) {

	Helmod::calculate(p);

	const ModifiedSpiral& B = v_B.at(currentSector);

	double r2 = std::pow(p.r, 2);
	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);

	iota = 1.0;
	// cos(30deg) = sqrt(3)/2 = 0.8660254
	if (std::abs(p.cosTheta) > 0.8660254) // polar region: theta < 30deg or theta > 150 deg
		iota = iotaPolarRegion;

	double g_low = v_glow.at(currentSector);
	KP = p.absRigidity / (GV/c) + g_low;

	// threshold for rigidity-dependent KP according to section 2 of the Helmod paper
	// Ekin = 444 MeV for protons corresponds to R=1.015 GV/c.
	if (KP < 1.015)
		KP = 1.015;

	double kappa0 = v_kappa0.at(currentSector);
	fdiffusion = p.beta * kappa0 * KP / 3.0 * B.Bearth;

	double denominator = std::pow(B.N, 1.5) * B.BfieldNormalization * B.r0_sqr * B.tau;
	kappa_rr = fdiffusion * r2 * (B.N*iota*rho*t2 + Gamma2*rho + B.tau) / denominator;
	kappa_thetatheta = fdiffusion * r2 * (B.N*iota*rho + t2*(Gamma2*rho + B.tau)) / denominator;
	kappa_r_theta = B.heavi * fdiffusion * r2 * B.t * (-B.N*iota*rho + Gamma2*rho + B.tau) / denominator;

}

double Helmod2012::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);

	return -fdiffusion*p.r*(2*B.N*t2*(B.N*iota*rho*t2 + Gamma2*rho + B.tau) - 2*B.N*B.tau*(3*B.N*iota*rho*t2 + 3*Gamma2*rho + iota*rho*t2*(Gamma2 + t2) + 3*t2 + 2)
				+ 3*B.tau*(Gamma2 + t2)*(B.N*iota*rho*t2 + Gamma2*rho + B.tau))/(pow(B.N, 5.0/2.0)*B.BfieldNormalization*B.r0_sqr*std::pow(B.tau, 2));
}

double Helmod2012::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double mu = p.cosTheta;
	double sintheta2 = std::pow(p.sinTheta, 2);
	double tprime = B.dtdmu;
	double t2 = std::pow(B.t, 2);
	double tau2 = std::pow(B.tau, 2);
	double GGp = B.Gammaprime*B.Gamma;

	return -fdiffusion*(2*B.N*B.t*tprime*(B.N*iota*rho + t2*(Gamma2*rho + B.tau))*sintheta2
			    + 2*B.N*B.tau*(mu*(B.N*iota*rho + t2*(Gamma2*rho + B.tau))
					   - (iota*rho*(GGp + B.t*tprime) + t2*(GGp*rho + B.t*tprime)
					      + B.t*tprime*(Gamma2*rho + B.tau))*sintheta2)
			    + 3*B.tau*(GGp + B.t*tprime)*(B.N*iota*rho + t2*(Gamma2*rho + B.tau))*sintheta2) /
			(B.BfieldNormalization*pow(B.N, 5.0/2.0)*B.r0_sqr*tau2);
}

double Helmod2012::getLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double mu = p.cosTheta;
	double sintheta2 = std::pow(p.sinTheta, 2);
	double t2 = std::pow(B.t, 2);
	double tau2 = std::pow(B.tau, 2);
	double tprime = B.dtdmu;
	double absSinTheta = std::abs(p.sinTheta);
	double GGp = B.Gammaprime*B.Gamma;

	return fdiffusion*p.r*(B.N*mu*B.t*B.tau*(-B.N*iota*rho + Gamma2*rho + B.tau) + 2*B.N*t2*tprime*(-B.N*iota*rho + Gamma2*rho + B.tau)*sintheta2
			       - B.N*B.tau*(2*B.t*(GGp*rho - iota*rho*(GGp + B.t*tprime) + B.t*tprime) + tprime*(-B.N*iota*rho + Gamma2*rho + B.tau))*sintheta2
			       + 3*B.t*B.tau*(GGp + B.t*tprime)*(-B.N*iota*rho + Gamma2*rho + B.tau)*sintheta2) /
			(B.BfieldNormalization*pow(B.N, 5.0/2.0)*B.heavi*B.r0_sqr*tau2*absSinTheta);
}

double Helmod2012::getLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaRTheta(p);

	const ModifiedSpiral& B = v_B.at(currentSector);
	double Gamma2 = std::pow(B.Gamma, 2);
	double t2 = std::pow(B.t, 2);
	double absSinTheta = std::abs(p.sinTheta);

	return fdiffusion*B.t*absSinTheta*(2*B.N*t2*(-B.N*iota*rho + Gamma2*rho + B.tau)
					   - 2*B.N*B.tau*(-2*B.N*iota*rho + 3*Gamma2*rho - iota*rho*(Gamma2 + t2) + t2 + 2*B.tau)
					   + 3*B.tau*(Gamma2 + t2)*(-B.N*iota*rho + Gamma2*rho + B.tau))/(pow(B.N, 5.0/2.0)*B.BfieldNormalization*B.heavi*B.r0_sqr*std::pow(B.tau, 2));
}


