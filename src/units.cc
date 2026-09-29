#include "units.h"

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
namespace py = pybind11;

void bind_units(py::module_& m)
{
	auto m_units = m.def_submodule("units");
	m_units.attr("meter") = Units::meter;
	m_units.attr("AU") = Units::AU;
	m_units.attr("meter2") = Units::meter2;
	m_units.attr("millimeter") = Units::millimeter;
	m_units.attr("millimeter2") = Units::millimeter2;
	m_units.attr("centimeter") = Units::centimeter;
	m_units.attr("centimeter2") = Units::centimeter2;
	m_units.attr("kilometer") = Units::kilometer;
	m_units.attr("kilometer2") = Units::kilometer2;
	m_units.attr("parsec") = Units::parsec;
	m_units.attr("mm") = Units::mm;
	m_units.attr("mm2") = Units::mm2;
	m_units.attr("cm") = Units::cm;
	m_units.attr("cm2") = Units::cm2;
	m_units.attr("m2") = Units::m2;
	m_units.attr("km") = Units::km;
	m_units.attr("km2") = Units::km2;
	m_units.attr("pc") = Units::pc;
	m_units.attr("pi") = Units::pi;
	m_units.attr("twopi") = Units::twopi;
	m_units.attr("pihalf") = Units::pihalf;
	m_units.attr("radian") = Units::radian;
	m_units.attr("milliradian") = Units::milliradian;
	m_units.attr("degree") = Units::degree;
	m_units.attr("steradian") = Units::steradian;
	m_units.attr("rad") = Units::rad;
	m_units.attr("mrad") = Units::mrad;
	m_units.attr("sr") = Units::sr;
	m_units.attr("deg") = Units::deg;
	m_units.attr("second") = Units::second;
	m_units.attr("days") = Units::days;
	m_units.attr("CarringtonRotationDuration") = Units::CarringtonRotationDuration;
	m_units.attr("s") = Units::s;
	m_units.attr("c") = Units::c;
	m_units.attr("c2") = Units::c2;
	m_units.attr("eplus") = Units::eplus;
	m_units.attr("gigaelectronvolt") = Units::gigaelectronvolt;
	m_units.attr("megaelectronvolt") = Units::megaelectronvolt;
	m_units.attr("electronvolt") = Units::electronvolt;
	m_units.attr("eV") = Units::eV;
	m_units.attr("MeV") = Units::MeV;
	m_units.attr("GeV") = Units::GeV;
	m_units.attr("gigavolt") = Units::gigavolt;
	m_units.attr("megavolt") = Units::megavolt;
	m_units.attr("volt") = Units::volt;
	m_units.attr("GV") = Units::GV;
	m_units.attr("MV") = Units::MV;
	m_units.attr("tesla") = Units::tesla;
	m_units.attr("nanotesla") = Units::nanotesla;
	m_units.attr("gauss") = Units::gauss;
	m_units.attr("kilogauss") = Units::kilogauss;
}
#endif
