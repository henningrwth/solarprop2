#ifndef ARROWS_H
#define ARROWS_H

#include "particle.h"

#include <vector>

/// \file arrows.h

/// \typedef std::pair<double, double> PosXZ
/// \brief A pair of coordinates indicating a 2D position.
using PosXZ = std::pair<double, double>;

/// \typedef std::pair<double, double> DeltaXZ
/// \brief A pair of components of a 2D vector.
using DeltaXZ = std::pair<double, double>;

/** Helper class used to store the contributions of the individual transport processes to the overall step of a pseudo-particle.
 *
 * This class is used to store the results of Solarprop::arrows() for \c solarprop_arrows.
 *
 * \sa SDE2d::calculateArrows
 * \sa SDE3d::calculateArrows
 *
 */
class Arrows
{
	public:
	Arrows(const Particle&);

	const Particle& getStartParticle() const { return start; }
	PosXZ getStartPosition() const;

	DeltaXZ getDeltaConvection() const { return deltaConvection; }
	DeltaXZ getDeltaDrift() const { return deltaDrift; }
	DeltaXZ getDeltaHcsDrift() const { return deltaHcsDrift; }
	std::vector<DeltaXZ> getDeltasDiffusion() const { return deltasDiffusion; }
	double getDt() const { return dt; }

	private:

	friend class SDE2d;
	friend class SDE3d;

	/// Particle used for the start position (base of arrows)
	const Particle start;

	/// displacement due to (gradient and curvature) drift, excluding HCS drift
	DeltaXZ deltaDrift;

	/// displacement due to HCS drift
	DeltaXZ deltaHcsDrift;

	/// displacement due to convection
	DeltaXZ deltaConvection;

	/// random sample of displacements due to diffusive transport
	std::vector<DeltaXZ> deltasDiffusion;

	/// time step
	double dt;
};

#endif
