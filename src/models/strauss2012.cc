#include "strauss2012.h"
#include "numdiff.h"
#include "sde2d.h"
#include "sde3d.h"
#include "units.h"

#include <cmath>
#include <iostream>

using namespace Units;

/**
 * Constructor.
 *
 * @param modelInfo %Model information.
 */
Strauss2012::Strauss2012(const ModelInformation& modelInfo) :
	Model(modelInfo),
	B(modelInfo)
{
	FACTORY_INIT;

	kappa_parallel = 0.0;
	kappa_perpendicular = 0.0;
	kappa_rr = 0.0;
	kappa_thetatheta = 0.0;
	kappa_phiphi = 0.0;
	fdiffusion = 0.;
	hcs_drift_speed_r = 0.;
	hcs_drift_speed_theta = 0.;
	hcs_drift_speed_phi = 0.;

	dt_max = modelInfo.getFloatOptionWithDefault("dt", 0.004 * days/s) * s;
	heliosphereBoundary = modelInfo.getFloatOptionWithDefault("heliosphereBoundary", 100., false) * AU;

	// we use the solarprop default value for OmegaSun, even though Strauss et al. probably used a slightly different value

	lambda0 = 0.05*AU;

	// Scale diffusion coefficient?
	double kappaScaling = modelInfo.getFloatOptionWithDefault("kappaScaling", 1.0, false);
	lambda0 *= kappaScaling;

	perPart = modelInfo.getFloatOptionWithDefault("rho", 0.02, false);
	kappa_exponent = modelInfo.getFloatOptionWithDefault("kappaExponent", 1.0, false);

	hcs_drift_velocity_factor = modelInfo.getFloatOptionWithDefault("hcsFactor", 1.0, false);
	simple_hcs_drift = modelInfo.getBooleanOptionWithDefault("simpleHcsDrift", false, false);
	forceNumericalDerivatives = modelInfo.getBooleanOptionWithDefault("forceNumericalDerivatives", false, false);

	force2D = modelInfo.getBooleanOptionWithDefault("force2D", false, false);

	double kappa0 = lambda0 * Units::c / 3.0;
	std::cout << "[Strauss2012] lambda0= " << lambda0/AU << " AU (kappa0= " << kappa0/(cm2/s) << " cm^2/s), "
		  << "a= " << kappa_exponent << ", "
		  << "lambda_perp/lambda_par= " << perPart << ", R= " << heliosphereBoundary/AU << " AU, "
		  << "HCS drift factor: " << hcs_drift_velocity_factor << ", full HCS drift: " << !simple_hcs_drift << std::endl;
}

std::unique_ptr<SDE> Strauss2012::createSDE(unsigned int globalSeed)
{
	if (force2D)
		return std::make_unique<SDE2d>(this, globalSeed);

	return std::make_unique<SDE3d>(this, globalSeed);
}

/**
 * Calculate the quantities needed several times in one particle step.
 */
void Strauss2012::calculate(const Particle& p)
{
	if (!nd && forceNumericalDerivatives)
	{
		nd = std::make_shared<NumericalDerivatives>(this);
	}

	B.calculate(p);

	double rigidityValue = p.absRigidity / (GV/c);
	double kappa_R = std::pow(rigidityValue, kappa_exponent);

	double lambda_parallel = lambda0 * B.getBearthOverB(p) * kappa_R;
	kappa_parallel = lambda_parallel * p.getVelocity() / 3.0;
	kappa_perpendicular = perPart * kappa_parallel;

	double Gamma2 = B.onePlusGamma2 - 1.0;
	kappa_rr = (kappa_parallel + kappa_perpendicular*Gamma2) / B.onePlusGamma2;
	kappa_thetatheta = kappa_perpendicular;
	kappa_rphi = (kappa_perpendicular - kappa_parallel) * B.Gamma / B.onePlusGamma2;
	kappa_phiphi = (kappa_parallel * Gamma2 + kappa_perpendicular) / B.onePlusGamma2;

	fdiffusion = lambda0 * p.getVelocity() * kappa_R * B.Bearth / (3.0 * B.BfieldNormalization * B.r0_sqr);

	// HCS drift
	hcs_drift_speed_r = 0.;
	hcs_drift_speed_theta = 0.;
	hcs_drift_speed_phi = 0.;

	if (!simple_hcs_drift)
	{
		try
		{
			auto [d_hcs, r_hcs, theta_hcs, phi_hcs] = B.closestDistanceToHcs(p);

			double r_larmor = B.getLarmorRadius(p);
			double x = d_hcs / r_larmor;

			if (x < 0.)
			{
				x = 0.;
			}

			if (x >= 2.)
			{
				hcs_drift_speed_r = 0.;
				hcs_drift_speed_theta = 0.;
				hcs_drift_speed_phi = 0.;
			}
			else
			{
				double velocity = (0.457 - 0.412*x + 0.0915*x*x) * p.getVelocity();
				velocity *= hcs_drift_velocity_factor;

				double sinbeta = 0.;
				double cosbeta = 1.;
				if (B.tiltAngle > 0.)
				{
					auto sincosbeta = B.sincosBeta(r_hcs, theta_hcs, phi_hcs);
					sinbeta = std::get<0>(sincosbeta);
					cosbeta = std::get<1>(sincosbeta);
				}

				double sinpsi = B.sinPsi(r_hcs, theta_hcs);
				double cospsi = B.cosPsi(r_hcs, theta_hcs);

				// speed vector at the HCS
				double v_r_prime = cosbeta * sinpsi * velocity;
				double v_theta_prime = sinbeta * velocity;
				double v_phi_prime = cosbeta * cospsi * velocity;

				// transform vector field to particle position
				//
				// from the Strauss+ paper, it is not entirely clear if this is needed/done in the reference implementation,
				// but Figure 3 of the paper implies that the velocity vector is evaluated at the HCS and then used at the
				// particle's position
				double dphi = p.getPhi() - phi_hcs;
				double sindphi = std::sin(dphi);
				double cosdphi = std::cos(dphi);
				double sintheta = p.getSinTheta();
				double costheta = p.getCosTheta();
				double sinthetahcs = std::sin(theta_hcs);
				double costhetahcs = std::cos(theta_hcs);
				// transformations: see calculations/vector_fields.ipynb
				hcs_drift_speed_r = v_phi_prime*sintheta*sindphi + v_r_prime * (sintheta*sinthetahcs*cosdphi + costheta*costhetahcs) + v_theta_prime * (sintheta*costhetahcs*cosdphi - sinthetahcs*costheta);
				hcs_drift_speed_theta = v_phi_prime*sindphi*costheta + v_r_prime * (-sintheta*costhetahcs + sinthetahcs*costheta*cosdphi) + v_theta_prime * (sintheta*sinthetahcs + costheta*costhetahcs*cosdphi);
				hcs_drift_speed_phi = v_phi_prime*cosdphi - v_r_prime*sinthetahcs*sindphi - v_theta_prime*sindphi*costhetahcs;

			}
		}
		catch (const std::exception& e)
		{
			std::cout << "[Strauss2012] Exception: " << e.what() << " for particle " << p.asString() << std::endl;
			throw(e);
			hcs_drift_speed_r = 0.;
			hcs_drift_speed_theta = 0.;
			hcs_drift_speed_phi = 0.;
		}
	}
}

/**
 * Simplified treatment of current sheet drifts, according to the "BP" model.
 *
 * BP model (see section 4.2 of Burger and Hattingh, Astrophys Space Sci 230 (1995) 375-382).
 */
double Strauss2012::simpleDriftSheetBP(const Particle& p) const
{
	double alpha = B.tiltAngle/rad;

	if (alpha > 0.)
	{
		double deltaTheta = 2. * p.absRigidity * B.solarWindSpeed(p) / (B.BfieldNormalization * B.r0_sqr * B.OmegaSun * std::cos(alpha));

		// Look if the particle is inside the "cone" of the wavy heliospheric current sheet
		if (std::abs(p.getTheta() - pi/2.) < alpha + deltaTheta)
		{
			double vhcs = B.polarity * p.charge * p.velocity / 6.0 * std::cos(alpha) * deltaTheta / std::sin(alpha + deltaTheta);
			return hcs_drift_velocity_factor * vhcs;
		}
		return 0.;
	}
	else
	{
		// Distance to the current sheet
		double d = fabs(p.r * p.cosTheta);
		double rL = B.getLarmorRadius(p);

		// There is only a drift effect close to the current sheet
		if (d < 2.*rL)
		{
			double vhcs = B.polarity * p.charge * p.velocity / 6.0;
			return hcs_drift_velocity_factor * vhcs;
		}

		return 0.0;
	}
	return 0.0;
}

/**
 * Simplified treatment of current sheet drifts, according to the "WNS" model.
 *
 * WNS model (see section 4.3 of Burger and Hattingh, Astrophys Space Sci 230 (1995) 375-382.
 */
double Strauss2012::simpleDriftSheetWNS(const Particle& p) const
{
	double alpha = B.tiltAngle/rad;
	double deltaTheta = 2. * p.absRigidity * B.solarWindSpeed(p) / (B.BfieldNormalization * B.r0_sqr * B.OmegaSun * std::cos(alpha));
	double GammaPrime = B.OmegaSun / B.solarWindSpeed(p) * p.r * std::cos(alpha / std::sqrt(2.));
	double vhcs = B.polarity * p.charge * p.velocity / 6.0 * 2. * p.absRigidity * p.r * GammaPrime / (B.BfieldNormalization * B.r0_sqr * (1. + GammaPrime*GammaPrime) * (alpha + deltaTheta));
	return hcs_drift_velocity_factor * vhcs;
}


double Strauss2012::getDriftSheetR(const Particle& p) const
{
	if (simple_hcs_drift)
	{
		return simpleDriftSheetBP(p);
	}

	return B.polarity * p.charge * hcs_drift_speed_r;
}

double Strauss2012::getDriftSheetTheta(const Particle& p) const
{
	if (simple_hcs_drift)
		return 0.0;

	return B.polarity * p.charge * hcs_drift_speed_theta;
}

double Strauss2012::getDriftSheetPhi(const Particle& p) const
{
	if (simple_hcs_drift)
		return 0.0;

	return B.polarity * p.charge * hcs_drift_speed_phi;
}


double Strauss2012::getLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRR(p);

	double Gamma2 = B.onePlusGamma2 - 1.0;

	return fdiffusion*p.r*(6*B.onePlusGamma2*Gamma2*perPart + 4*B.onePlusGamma2 - 3*Gamma2*Gamma2*perPart - 3*Gamma2) / pow(B.onePlusGamma2, 2.5);
}

double Strauss2012::getLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearRadialDiffusionTermFromKappaRPhi(p);

	return 0.0;
}

double Strauss2012::getLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearPolarDiffusionTermFromKappaThetaTheta(p);

	return -fdiffusion*p.cosTheta*perPart*(B.onePlusGamma2 + 1) / pow(B.onePlusGamma2, 1.5);
}

double Strauss2012::getLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearAzimuthalDiffusionTermFromKappaRPhi(p);

	return fdiffusion*B.Gamma*(B.onePlusGamma2 + 3)*(perPart - 1) / (pow(B.onePlusGamma2, 2.5)*p.sinTheta);
}

double Strauss2012::getLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const
{
	if (forceNumericalDerivatives)
		return nd->calculateLinearAzimuthalDiffusionTermFromKappaPhiPhi(p);

	return 0.0;
}
