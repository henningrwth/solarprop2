#include "random.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

/**
 * Constructor with seed calculated from current time.
 */
Random::Random() :
	generator(),
	normal_dist(0., 1.),
	uniform_dist(0., 1.)
{
	unsigned int seed = std::chrono::system_clock().now().time_since_epoch() / std::chrono::microseconds(1);

	//std::cout << "[Random] Setting random seed: " << seed << std::endl;
	generator.seed(seed);
}

/**
 * Constructor with fixed seed.
 *
 * @param seed Random seed.
 */
Random::Random(unsigned int seed) :
	generator(),
	normal_dist(0., 1.),
	uniform_dist(0., 1.)
{
	//std::cout << "[Random] Setting random seed: " << seed << std::endl;
	generator.seed(seed);
}

/**
 * Generate a Gaussian-distributed random number.
 *
 * @return Gaussian-distributed random number.
 */
double Random::getGaussianRandomNumber()
{
	return normal_dist(generator);
}

/**
 * Generate a uniformly distributed random number between 0 and 1.
 *
 * @return Uniformly distributed random number between 0 and 1.
 */
double Random::getUniformRandomNumber()
{
	return uniform_dist(generator);
}
