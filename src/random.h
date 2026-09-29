#ifndef RANDOM_H_
#define RANDOM_H_

#include <random>

/**
 * Wrapper class for random number generation.
 */
class Random
{
	public:

	Random();
	Random(unsigned int seed);
	~Random() {}

	double getGaussianRandomNumber();
	double getUniformRandomNumber();

	private:

	/// Random number engine
	std::mt19937 generator;

	/// Normal (gaussian) distribution of random numbers
	std::normal_distribution<double> normal_dist;

	/// Uniform distribution of random numbers
	std::uniform_real_distribution<double> uniform_dist;
};

#endif
