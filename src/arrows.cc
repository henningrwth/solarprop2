#include "arrows.h"

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void bind_arrows(py::module_& m)
{
	py::class_<Arrows>(m, "Arrows")
			.def("startPosition", &Arrows::getStartPosition)
			.def("deltaConvection", &Arrows::getDeltaConvection)
			.def("deltaDrift", &Arrows::getDeltaDrift)
			.def("deltaHcsDrift", &Arrows::getDeltaHcsDrift)
			.def("deltasDiffusion", &Arrows::getDeltasDiffusion)
			.def("dt", &Arrows::getDt)
			;
}
#endif

/**
 * @brief Constructor.
 * @param _start Start position: base point of all step vectors.
 */
Arrows::Arrows(const Particle& _start) :
	start(_start)
{
}

/** Get start position (for base of arrows). */
PosXZ Arrows::getStartPosition() const
{
	return {start.getX(), start.getZ()};
}
