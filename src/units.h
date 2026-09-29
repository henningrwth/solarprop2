#ifndef UNITS_H_
#define UNITS_H_

/// \file units.h
/// System of units used in \c SOLARPROP.
///
/// A HEP coherent system of Units (adapted from CLHEP/Geant4)
/// is used throught the C++ side of SOLARPROP. Function
/// arguments used to pass information from the Python side
/// are supposed to have a name explicitly stating the required
/// unit. For example:
///
///     #include "units.h"
///     using namespace Units;
///
///     double Solarprop::some_position(double ekin_in_GeV)
///     {
///         double ekin = ekin_in_GeV * GeV;
///         // (...)
///     }
///
/// Return values are given in units defined by this unit system
/// and have to be converted manually, though this will often not
/// be necessary due to the choices for the basic units (see below).
/// Example:
///
///     from solarprop.units import AU
///
///     s = solarprop.Solarprop(parameters)
///
///     value = s.some_position(5.0)
///     print(value/AU)
///
///
/// The basic units are :
///  - meter                   (`meter`)
///  - second                  (`second`)
///  - giga electron volt      (`GeV`)
///  - positron charge         (`eplus`)
///  - radian                  (`radian`)
///  - steradian               (`steradian`)
///
/// The relevant derived units are defined in this file.
///

namespace Units
{

//
// Length [L]
//

static constexpr double meter = 1.0;

/// astronomical unit
static constexpr double AU = 1.495978707e+11 * meter;

static constexpr double meter2 = meter * meter;
static constexpr double millimeter = 1.e-3 * meter;
static constexpr double millimeter2 = millimeter * millimeter;

static constexpr double centimeter = 1.e-2 * meter;
static constexpr double centimeter2 = centimeter * centimeter;

static constexpr double kilometer = 1000. * meter;
static constexpr double kilometer2 = kilometer * kilometer;

static constexpr double parsec = 3.0856775807e+16 * meter;

// abbreviations
static constexpr double mm = millimeter;
static constexpr double mm2 = millimeter2;

static constexpr double cm = centimeter;
static constexpr double cm2 = centimeter2;

static constexpr double m2 = meter2;

static constexpr double km = kilometer;
static constexpr double km2 = kilometer2;

static constexpr double pc = parsec;

//
// Angle
//
static constexpr double pi = 3.14159265358979323846;
static constexpr double twopi = 2.0 * pi;
static constexpr double pihalf = pi / 2.0;
static constexpr double radian = 1.;
static constexpr double milliradian = 1.e-3 * radian;
static constexpr double degree = (pi / 180.0) * radian;

static constexpr double steradian = 1.;

// abbreviations
static constexpr double rad = radian;
static constexpr double mrad = milliradian;
static constexpr double sr = steradian;
static constexpr double deg = degree;

//
// Time [T]
//
static constexpr double second = 1.0;
static constexpr double hour = 3600.0 * second;
static constexpr double days = 24.0 * hour;

static constexpr double CarringtonRotationDuration = 27.2752316 * days;

// abbreviation
static constexpr double s = second;

/// speed of light
static constexpr double c = 2.99792458e+8 * meter / s;

/// speed of light squared
static constexpr double c2 = c * c;

//
// Electric charge [Q]
//
static constexpr double eplus = 1.; // positron charge

//
// Energy [E]
//
static constexpr double gigaelectronvolt = 1.;
static constexpr double megaelectronvolt = 1.e-3 * gigaelectronvolt;
static constexpr double electronvolt = 1.e-6 * megaelectronvolt;

// abbreviations
static constexpr double eV = electronvolt;
static constexpr double MeV = megaelectronvolt;
static constexpr double GeV = gigaelectronvolt;

//
// Electric potential [E][Q^-1]
//
static constexpr double gigavolt = gigaelectronvolt / eplus;
static constexpr double megavolt = 1.e-3 * gigavolt;
static constexpr double volt = 1.e-6 * megavolt;

static constexpr double GV = gigavolt;
static constexpr double MV = megavolt;

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static constexpr double tesla = volt * second / meter2;
static constexpr double nanotesla = 1.e-9 * tesla;

static constexpr double gauss = 1.e-4 * tesla;
static constexpr double kilogauss = 1.e-1 * tesla;

}

#endif
