#include "numdiff.h"

#include "model.h"
#include "particle.h"
#include "units.h"

using namespace Units;

/**
 * Constructor.
 *
 * \param m Model for which derivatives are calculated.
 */
NumericalDerivatives::NumericalDerivatives(const Model* m) :
	model(m)
{
	modelForDerivatives = model->clone();
}

/** Infinitesimal shift in \f$ r \f$ for radial derivatives. */
double NumericalDerivatives::epsilon_r(const Particle& p)
{
	return 0.0001*p.getR();
}

/** Infinitesimal shift in \f$ \theta \f$ for polar derivatives. */
double NumericalDerivatives::epsilon_theta(const Particle& p)
{
	double dtheta = 1.e-7*rad;
	if (p.getCosTheta() < 0.0)
		dtheta = -dtheta;
	return dtheta;
}

/** Infinitesimal shift in \f$ \cos\theta \f$ for polar derivatives. */
double NumericalDerivatives::epsilon_costheta(const Particle& p)
{
	double dmu = -1.e-7;
	if (p.getCosTheta() < 0.0)
		dmu = -dmu;
	return dmu;
}

/** Infinitesimal shift in \f$ \phi \f$ for azimuthal derivatives. */
double NumericalDerivatives::epsilon_phi(const Particle& p)
{
	double dphi = 1.e-7*rad;
	if (p.getPhi() > pi)
		dphi = -dphi;
	return dphi;
}

double NumericalDerivatives::calculateEnergyGainTerm(const Particle& p) const
{
	double r = p.getR();
	double V_sw = model->solarWindSpeed(p);

	Particle p2(p);
	double dr = epsilon_r(p);
	p2.setR(r+dr);
	modelForDerivatives->calculate(p2);
	double dVsw_dr = (modelForDerivatives->solarWindSpeed(p2) - V_sw) / dr;

	return 2.0*V_sw/r + dVsw_dr;
}

double NumericalDerivatives::calculateLinearRadialDiffusionTermFromKappaRR(const Particle& p) const
{
	double r = p.getR();
	double kappaRR = model->getKappaRR();

	Particle p2(p);
	double dr = epsilon_r(p);
	p2.setR(r+dr);
	modelForDerivatives->calculate(p2);
	double kappaRR_2 = modelForDerivatives->getKappaRR();
	double d_kappa_rr_d_r = (kappaRR_2 - kappaRR) / dr;
	return 2.*kappaRR / r + d_kappa_rr_d_r;
}

double NumericalDerivatives::calculateLinearRadialDiffusionTermFromKappaRPhi(const Particle& p) const
{
	// 1/(r*sin(theta)) * d/dphi (kappa_r_phi)

	double r = p.getR();
	double sinTheta = p.getSinTheta();
	double phi = p.getPhi();

	double pdphi = 0.5*epsilon_phi(p);

	Particle p1(p);
	p1.setPhi(phi - pdphi);
	modelForDerivatives->calculate(p1);
	double kappa_r_phi_1 = modelForDerivatives->getKappaRPhi();

	Particle p2(p);
	p2.setPhi(phi + pdphi);
	modelForDerivatives->calculate(p2);
	double kappa_r_phi_2 = modelForDerivatives->getKappaRPhi();

	return 1./(r*sinTheta) * (kappa_r_phi_2 - kappa_r_phi_1) / (2*pdphi);
}

double NumericalDerivatives::calculateLinearThetaDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	double r = p.getR();
	double cosTheta = p.getCosTheta();
	double sinTheta = p.getSinTheta();
	double kappaThetaTheta = model->getKappaThetaTheta();

	Particle p2(p);
	double pdtheta = epsilon_theta(p);
	p2.setCosTheta(cosTheta - sinTheta*pdtheta);
	modelForDerivatives->calculate(p2);
	double kappaThetaTheta_2 = modelForDerivatives->getKappaThetaTheta();
	double d_kappa_thetatheta_d_theta = (kappaThetaTheta_2 - kappaThetaTheta) / pdtheta;
	return 1./(r*r) * (p.getCotTheta() * kappaThetaTheta + d_kappa_thetatheta_d_theta);
}

double NumericalDerivatives::calculateLinearRadialDiffusionTermFromKappaRTheta(const Particle& p) const
{
	double r = p.getR();
	double cosTheta = p.getCosTheta();
	double pdcostheta = 0.5*epsilon_costheta(p);

	Particle p1(p);
	p1.setCosTheta(cosTheta - pdcostheta);
	modelForDerivatives->calculate(p1);
	double kappaRTheta_x_sinTheta_1 = modelForDerivatives->getKappaRTheta() * p1.getSinTheta();

	Particle p2(p);
	p2.setCosTheta(cosTheta + pdcostheta);
	modelForDerivatives->calculate(p2);
	double kappaRTheta_x_sinTheta_2 = modelForDerivatives->getKappaRTheta() * p2.getSinTheta();

	return -1./r * (kappaRTheta_x_sinTheta_2 - kappaRTheta_x_sinTheta_1) / (2*pdcostheta);
}

double NumericalDerivatives::calculateLinearPolarDiffusionTermFromKappaRTheta(const Particle& p) const
{
	double r = p.getR();
	double sinTheta = p.getSinTheta();
	double kappaRTheta = model->getKappaRTheta();

	Particle p2(p);
	double dr = epsilon_r(p);
	p2.setR(r+dr);
	modelForDerivatives->calculate(p2);
	double kappaRTheta_2 = modelForDerivatives->getKappaRTheta();
	double d_kappa_rtheta_d_r = (kappaRTheta_2 - kappaRTheta) / dr;
	return -(sinTheta/r) * (kappaRTheta / r + d_kappa_rtheta_d_r);
}

double NumericalDerivatives::calculateLinearPolarDiffusionTermFromKappaThetaTheta(const Particle& p) const
{
	double r = p.getR();
	double cosTheta = p.getCosTheta();
	double pdcostheta = 0.5*epsilon_costheta(p);

	Particle p1(p);
	p1.setCosTheta(cosTheta - pdcostheta);
	modelForDerivatives->calculate(p1);
	double kappaThetaTheta_x_sin2theta_1 = modelForDerivatives->getKappaThetaTheta() * std::pow(p1.getSinTheta(), 2);

	Particle p2(p);
	p2.setCosTheta(cosTheta + pdcostheta);
	modelForDerivatives->calculate(p2);
	double kappaThetaTheta_x_sin2theta_2 = modelForDerivatives->getKappaThetaTheta() * std::pow(p2.getSinTheta(), 2);

	return 1./(r*r) * (kappaThetaTheta_x_sin2theta_2 - kappaThetaTheta_x_sin2theta_1) / (2*pdcostheta);
}

double NumericalDerivatives::calculateLinearAzimuthalDiffusionTermFromKappaRPhi(const Particle& p) const
{
	// 1/(r^2 sin(theta)) * d/dr (r*kappa_r_phi)

	double r = p.getR();
	double sinTheta = p.getSinTheta();

	double r_kappa_r_phi = r * model->getKappaRPhi();

	Particle p2(p);
	double dr = epsilon_r(p);
	double r_2 = r + dr;
	p2.setR(r_2);
	modelForDerivatives->calculate(p2);
	double r_kappa_r_phi_2 = r_2 * modelForDerivatives->getKappaRPhi();

	return 1./(r*r*sinTheta) * (r_kappa_r_phi_2 - r_kappa_r_phi) / dr;
}

double NumericalDerivatives::calculateLinearAzimuthalDiffusionTermFromKappaPhiPhi(const Particle& p) const
{
	// 1/(r^2 sin^2(theta)) * d/dphi (kappa_phi_phi)

	double r = p.getR();
	double sinTheta = p.getSinTheta();
	double phi = p.getPhi();

	double pdphi = 0.5*epsilon_phi(p);

	Particle p1(p);
	p1.setPhi(phi - pdphi);
	modelForDerivatives->calculate(p1);
	double kappa_phi_phi_1 = modelForDerivatives->getKappaPhiPhi();

	Particle p2(p);
	p2.setPhi(phi + pdphi);
	modelForDerivatives->calculate(p2);
	double kappa_phi_phi_2 = modelForDerivatives->getKappaPhiPhi();

	return 1./(r*r*sinTheta*sinTheta) * (kappa_phi_phi_2 - kappa_phi_phi_1) / (2*pdphi);
}

double NumericalDerivatives::calculateRadialGradientCurvatureDriftVelocity(const Particle& p, bool includePolarComponent) const
{
	double r = p.getR();

	double cosTheta = p.getCosTheta();
	double sinTheta = p.getSinTheta();
	double F_phi = model->getBOverB2_phi(p);
	double F_phi_sintheta = F_phi * sinTheta;

	Particle p2theta(p);
	double dtheta = epsilon_theta(p);
	p2theta.setCosTheta(cosTheta - sinTheta*dtheta);
	modelForDerivatives->calculate(p2theta);

	double F_phi_sintheta_2 = modelForDerivatives->getBOverB2_phi(p2theta) * p2theta.getSinTheta();
	if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2theta) < 0)
		F_phi_sintheta_2 = -F_phi_sintheta_2;
	double dF_phi_sintheta_dtheta = (F_phi_sintheta_2 - F_phi_sintheta) / dtheta;

	double dF_theta_dphi = 0.0;
	if (includePolarComponent)
	{
		double phi = p.getPhi();
		double F_theta = model->getBOverB2_theta(p);

		Particle p2phi(p);
		double dphi = epsilon_phi(p);
		p2phi.setPhi(phi + dphi);
		modelForDerivatives->calculate(p2phi);
		double F_theta_2 = modelForDerivatives->getBOverB2_theta(p2phi);
		if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2phi) < 0)
			F_theta_2 = -F_theta_2;
		dF_theta_dphi = (F_theta_2 - F_theta) / dphi;
	}

	double factor = p.getRigidity() * p.getVelocity() / 3.;

	double vDr = (factor / r / sinTheta) * (dF_phi_sintheta_dtheta - dF_theta_dphi);

	return vDr;
}

double NumericalDerivatives::calculatePolarGradientCurvatureDriftVelocity(const Particle& p) const
{
	double r = p.getR();
	double sinTheta = p.getSinTheta();
	double F_phi = model->getBOverB2_phi(p);

	Particle p2r(p);
	double dr = epsilon_r(p);
	p2r.setR(r + dr);
	modelForDerivatives->calculate(p2r);
	double F_phi_2 = modelForDerivatives->getBOverB2_phi(p2r);
	if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2r) < 0)
		F_phi_2 = -F_phi_2;
	double dF_phi_dr = (F_phi_2 - F_phi) / dr;

	// disable phi dependence
	//
	//double phi = p.getPhi();
	//double F_r = model->getBOverB2_r(p);

	//Particle p2phi(p);
	//double dphi = epsilon_phi(p);
	//p2phi.setPhi(phi + dphi);
	//modelForDerivatives->calculate(p2phi);
	//double F_r_2 = modelForDerivatives->getBOverB2_r(p2phi);
	//if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2phi) < 0)
	//	F_r_2 = -F_r_2;
	//double dF_r_dphi = (F_r_2 - F_r) / dphi;

	double dF_r_dphi = 0.0;

	double factor = p.getRigidity() * p.getVelocity() / 3.;
	double vDtheta = factor * (dF_r_dphi / r / sinTheta - dF_phi_dr - F_phi / r);
	return vDtheta;
}

double NumericalDerivatives::calculateAzimuthalGradientCurvatureDriftVelocity(const Particle& p, bool includePolarComponent) const
{
	double r = p.getR();
	double cosTheta = p.getCosTheta();
	double sinTheta = p.getSinTheta();

	double F_r = model->getBOverB2_r(p);
	double F_theta = model->getBOverB2_theta(p);

	Particle p2theta(p);
	double dtheta = epsilon_theta(p);
	p2theta.setCosTheta(cosTheta - sinTheta*dtheta);
	modelForDerivatives->calculate(p2theta);

	double F_r_2 = modelForDerivatives->getBOverB2_r(p2theta);
	if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2theta) < 0)
		F_r_2 = -F_r_2;
	double dF_r_dtheta = (F_r_2 - F_r) / dtheta;


	double dF_theta_dr = 0.0;
	if (includePolarComponent)
	{
		Particle p2r(p);
		double dr = epsilon_r(p);
		p2r.setR(r + dr);
		modelForDerivatives->calculate(p2r);
		double F_theta_2 = modelForDerivatives->getBOverB2_theta(p2r);
		if (model->getHeaviside(p) * modelForDerivatives->getHeaviside(p2r) < 0)
			F_theta_2 = -F_theta_2;
		dF_theta_dr = (F_theta_2 - F_theta) / dr;
	}

	double factor = p.getRigidity() * p.getVelocity() / 3.;
	double vDphi = factor * (dF_theta_dr + F_theta / r  - dF_r_dtheta / r);
	return vDphi;
}
