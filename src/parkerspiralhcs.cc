#include "parkerspiralhcs.h"

#include "particle.h"
#include "units.h"

#include <iostream>

using namespace Units;

/**
 * Constructor
 * @param modelInfo %Model information.
 */
ParkerSpiralHCS::ParkerSpiralHCS(const ModelInformation& modelInfo) :
	ParkerSpiral(modelInfo, true)
{
	// error handling for GSL
	gsl_set_error_handler(
		[](const char* reason, const char* file, int line, int gsl_errno)
		{
			std::cout << file << " at line " << line << ": " << reason << " errno " << gsl_errno << " (" << gsl_strerror(gsl_errno) << ")"
				  << std::endl;
			throw std::runtime_error(reason);
		});

	const gsl_multimin_fminimizer_type* minimizer_type = gsl_multimin_fminimizer_nmsimplex2;

	x0_min = gsl_vector_alloc(2);
	stepsizes_min = gsl_vector_alloc(2);
	minimizer = gsl_multimin_fminimizer_alloc(minimizer_type, 2);
}

/**
 * Copy constructor.
 *
 * Makes a deep copy of the GSL minimizer.
 */
ParkerSpiralHCS::ParkerSpiralHCS(const ParkerSpiralHCS& s) :
	ParkerSpiral(s)
{
	const gsl_multimin_fminimizer_type* minimizer_type = gsl_multimin_fminimizer_nmsimplex2;

	x0_min = gsl_vector_alloc(2);
	stepsizes_min = gsl_vector_alloc(2);
	minimizer = gsl_multimin_fminimizer_alloc(minimizer_type, 2);
}

ParkerSpiralHCS::~ParkerSpiralHCS()
{
	gsl_vector_free(x0_min);
	gsl_vector_free(stepsizes_min);
	gsl_multimin_fminimizer_free(minimizer);
}

/** Calculate cached values at the beginning of each pseudo-particle step.
 *
 * \attention Make sure this function is called before any of the getters in this class.
 */
void ParkerSpiralHCS::calculate(const Particle& p)
{
	ParkerSpiral::calculate(p);

	// we assume that the solar wind speed is constant and radially directed throughout the heliosphere
	omegaOverVsw = OmegaSun / solarWindSpeed(p);
}

double ParkerSpiralHCS::getThetaPrime(const Particle& p) const
{
	return thetaPrimeHcs(p.getR(), p.getPhi());
}

double ParkerSpiralHCS::thetaPrimeHcs(double r, double phi) const
{
	// see eq. (11) of the Strauss+ (2012) paper
	return pihalf + std::asin(std::sin(tiltAngle) * std::sin(phi - phi_0 + omegaOverVsw * r));
}

double ParkerSpiralHCS::tanPsi(double r, double theta) const
{
	return omegaOverVsw * r * std::sin(theta);
}

double ParkerSpiralHCS::sinPsi(double r, double theta) const
{
	double x = tanPsi(r, theta);
	return x / std::sqrt(x*x + 1.);
}

double ParkerSpiralHCS::cosPsi(double r, double theta) const
{
	double x = tanPsi(r, theta);
	return 1. / std::sqrt(x*x + 1.);
}

double ParkerSpiralHCS::sinPsiSqr(double r, double theta) const
{
	double x = tanPsi(r, theta);
	double x2 = x*x;
	return x2 / (x2 + 1.);
}

double ParkerSpiralHCS::betaSign(double r, double phi) const
{
	// see eq. (15) of the Strauss+ (2012) paper
	return std::copysign(1.0, std::cos(phi + phi_0 + omegaOverVsw * r));
}

double ParkerSpiralHCS::tanBetaSqr(double r, double theta) const
{
	// see eq. (14) of the Strauss+ (2012) paper

	double sinThetaSqr = std::pow(std::sin(theta), 2);
	double cosThetaSqr = 1. - sinThetaSqr;
	double numer = std::pow(std::sin(tiltAngle), 2) - cosThetaSqr;
	if (numer < 0.)
		numer = 0.;
	return std::pow(omegaOverVsw * r, 2) / sinPsiSqr(r, theta) * numer / sinThetaSqr;
}

double ParkerSpiralHCS::cosBeta(double r, double theta) const
{
	double x2 = tanBetaSqr(r, theta);
	return 1. / std::sqrt(x2+1.);
}

double ParkerSpiralHCS::sinBeta(double r, double theta, double phi) const
{
	double sign = betaSign(r, phi);
	double x2 = tanBetaSqr(r, theta);
	return sign * std::sqrt(x2 / (x2 + 1.));
}

std::tuple<double, double> ParkerSpiralHCS::sincosBeta(double r, double theta, double phi) const
{
	double sign = betaSign(r, phi);
	double x2 = tanBetaSqr(r, theta);
	double cosbeta = 1. / std::sqrt(x2 + 1.);
	double sinbeta = sign * std::sqrt(x2) * cosbeta;
	return {sinbeta, cosbeta};
}

double ParkerSpiralHCS::beta(double r, double theta, double phi) const
{
	double sign = betaSign(r, phi);
	return sign * std::atan(std::sqrt(tanBetaSqr(r, theta)));
}

double ParkerSpiralHCS::distanceToHcsPointSqr(double r, double phi, const Particle& p) const
{
	double theta = thetaPrimeHcs(r, phi);
	double rp = p.getR();
	double sinthetap = p.getSinTheta();
	double costhetap = p.getCosTheta();
	double phip = p.getPhi();

	return r*r + rp*rp - 2*r*rp*(std::sin(theta)*sinthetap*std::cos(phi-phip) + std::cos(theta)*costhetap);
}

std::tuple<double, double, double, double> ParkerSpiralHCS::closestDistanceToHcs(const Particle& p)
{
	const bool verbose = false;
	gsl_multimin_function minex_func;

	int iter = 0;
	int status;

	particlePosition = &p;

	// Starting point
	gsl_vector_set(x0_min, 0, p.getR() / AU);
	gsl_vector_set(x0_min, 1, std::fmod(p.getPhi(), twopi));

	// Set initial step sizes.
	gsl_vector_set(stepsizes_min, 0, 0.05);
	gsl_vector_set(stepsizes_min, 1, 0.01);

	// Initialize method and iterate.
	minex_func.n = 2;
	auto myf = [](const gsl_vector* v, void* params) -> double
	{
		const ParkerSpiralHCS* thisSpiral = static_cast<ParkerSpiralHCS*>(params);

		double r = gsl_vector_get(v, 0) * AU;
		double phi = gsl_vector_get(v, 1);

		if (r < 0.)
			return GSL_NAN;

		return thisSpiral->distanceToHcsPointSqr(r, phi, *(thisSpiral->particlePosition));
	};
	minex_func.f = myf;
	minex_func.params = this;

	gsl_multimin_fminimizer_set(minimizer, &minex_func, x0_min, stepsizes_min);

	do
	{
		++iter;
		status = gsl_multimin_fminimizer_iterate(minimizer);
		if (status)
			break;

		double size = gsl_multimin_fminimizer_size(minimizer);
		status = gsl_multimin_test_size(size, 2e-3);

		if (verbose)
		{
			if (status == GSL_SUCCESS)
			{
				printf("converged to minimum at\n");
			}

			printf("%5d %10.3f AU %10.3f deg: D = %7.4f AU, size = %.4f\n",
			       iter,
			       gsl_vector_get(minimizer->x, 0),
			       gsl_vector_get(minimizer->x, 1) / deg,
			       std::sqrt(minimizer->fval) / AU, size);
		}
	} while (status == GSL_CONTINUE && iter < 80);

	double mindist = std::sqrt(minimizer->fval);
	double r_hcs = gsl_vector_get(minimizer->x, 0) * AU;
	double phi_hcs = gsl_vector_get(minimizer->x, 1);
	double theta_hcs = thetaPrimeHcs(r_hcs, phi_hcs);

	return {mindist, r_hcs, theta_hcs, phi_hcs};
}
