#include "properties.h"

#include "model.h"
#include "modelfactory.h"
#include "units.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace std::string_literals;

using namespace Units;

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void bind_modelinformation(py::module_& m)
{
	py::class_<ModelInformation>(m, "ModelInformation")
			.def(py::init<>())
			.def("calculateCR", &ModelInformation::calculateCR)
			.def("getPolarityAndPhaseFor", &ModelInformation::getPolarityAndPhaseFor)
			;
}
#endif

/**
 * Fill the lookups for time-dependent data from the input files specified by the respective options, or clear them.
 *
 * The user has to specify the location of the relevant input file for each lookup. This allows the user to either use
 * the default data files shipped with the solarprop code or provide their own version. The following lookups are filled
 * from data files:
 *
 * | Option for input file name       | Lookup|
 * |--------------------------------- | ----------------------|
 * | ModelInformation::ssnFileOption  | Smoothed sunspot number|
 * | ModelInformation::omniFileOption | OMNIweb data on magnetic field amplitude and solar wind speed|
 * | `nmFile`                         | Modulation potential inferred from neutron monitor data|
 * | `angleFile`                      | HCS tilt angle|
 * | `kappaFile`                      | A user-provided lookup for the diffusion coeffient scale \f$ \kappa_0 \f$|
 * | `parameterLookup`                | User-provided lookup for arbitrary model parameters, see getFloatParameterFromOptionOrLookupOrDefault().|
 *
 * If an option for the input file name is not given, the respective lookup is cleared.
 *
 * Data files will only be re-read if the source file name has changed.
 */
void ModelInformation::readDatafiles()
{
	if (hasOption(ssnFileOption))
	{
		auto fname = getStringOption(ssnFileOption);
		if (fname != ssnLookupFilename)
		{
			ssnLookup = FileHelper::readSunspotNumberData(fname);
			if (ssnLookup.size() < 2)
				throw std::runtime_error("ERROR reading SSN lookup: "s + fname);
			ssnLookupFilename = fname;
		}
	}
	else
	{
		ssnLookup.clear();
		ssnLookupFilename.clear();
	}

	if (hasOption(omniFileOption))
	{
		auto fname = getStringOption(omniFileOption);
		if (fname != omniFilename)
		{
			std::tie(bmagLookup, vswLookup) = FileHelper::readOmniWebData(fname);
			if (bmagLookup.size() < 2)
				throw std::runtime_error("ERROR reading B field lookup: "s + fname);
			if (vswLookup.size() < 2)
				throw std::runtime_error("ERROR reading solar wind lookup: "s + fname);
			omniFilename = fname;
		}
	}
	else
	{
		bmagLookup.clear();
		vswLookup.clear();
		omniFilename.clear();
	}

	if (hasOption("nmFile"))
	{
		auto fname = getStringOption("nmFile");
		if (fname != nmFilename)
		{
			nmData = FileHelper::readNeutronMonitorData(fname);
			if (nmData.size() < 2)
				throw std::runtime_error("ERROR reading NM data lookup: "s + fname);
			nmFilename = fname;
		}
	}
	else
	{
		nmData.clear();
		nmFilename.clear();
	}

	if (hasOption("kappaFile"))
	{
		auto fname = getStringOption("kappaFile");
		if (fname != kappaFilename)
		{
			kappaLookup = FileHelper::readKappaFile(fname);
			if (kappaLookup.size() < 2)
				throw std::runtime_error("ERROR reading kappa lookup: "s + fname);
			kappaFilename = fname;
		}
	}
	else
	{
		kappaLookup.clear();
		kappaFilename.clear();
	}

	if (hasOption("angleFile"))
	{
		std::string fname = getStringOption("angleFile");
		if (fname != angleFilename)
		{

			// Get tilt angle model (L or R model)
			std::string tiltModel = getStringOptionWithDefault("tiltModel", "R");

			std::cout << "Reading tilt angle data (model \"" << tiltModel << "\") from file: " << fname << std::endl;
			if (tiltModel == "R")
			{
				tiltAngleData = FileHelper::readFile<int, double>(fname, 2, 5);
			}
			else if (tiltModel == "L")
			{
				tiltAngleData = FileHelper::readFile<int, double>(fname, 2, 8);
			}
			else
				throw std::runtime_error("Invalid tilt model: "s + tiltModel);

			// set unit
			std::for_each(tiltAngleData.begin(), tiltAngleData.end(), [](std::pair<const int, double>& kv) { kv.second *= degree; });

			angleFilename = fname;
		}
	}
	else
	{
		tiltAngleData.clear();
		angleFilename.clear();
	}

	if (hasOption("parameterLookup"))
	{
		auto fname = getStringOption("parameterLookup");
		if (fname != parameterLookupFilename)
		{
			parameterLookup = FileHelper::readParameterLookupFile(fname);
			parameterLookupFilename = fname;
		}
	}
	else
	{
		parameterLookup.clear();
		parameterLookupFilename.clear();
	}
}

/** Use the \c ModelFactory to create an instance of a model implementation, based on the value of the \c model option provided by the user. */
std::unique_ptr<Model> ModelInformation::createModel() const
{
	auto model = getStringOption("model");

	std::unique_ptr<Model> p = Factory<Model, const ModelInformation&>::Create(model, *this);
	if (!p)
		throw std::runtime_error("Unknown model: "s + model);
	return p;
}

/** Check if option \p key was specified by the user. */
bool ModelInformation::hasOption(const std::string& key) const
{
	return hasStringOption(key) || hasIntegerOption(key) || hasFloatOption(key);
}

/** Internal helper function: Check if the option \p key was provided with an integer value. */
bool ModelInformation::hasIntegerOption(const std::string& key) const
{
	return integerOptions.find(key) != integerOptions.end();
}

/** Internal helper function: Check if the option \p key was provided with a floating-point value. */
bool ModelInformation::hasFloatOption(const std::string& key) const
{
	return floatOptions.find(key) != floatOptions.end();
}

/** Internal helper function: Check if the option \p key was provided with a string value. */
bool ModelInformation::hasStringOption(const std::string& key) const
{
	return stringOptions.find(key) != stringOptions.end();
}

/** Store option \p key that was provided with an integer value of \p value. */
void ModelInformation::storeIntegerOption(const std::string& key, int value)
{
	//std::cout << "[DEBUG] New integer option \"" << key << "\" = " << value << std::endl;
	usedOptionKeys.erase(key);
	integerOptions[key] = value;
}

/** Store option \p key that was provided with a floating-point value of \p value. */
void ModelInformation::storeFloatOption(const std::string& key, double value)
{
	//std::cout << "[DEBUG] New float option \"" << key << "\" = " << value << std::endl;
	usedOptionKeys.erase(key);
	floatOptions[key] = value;
}

/** Store option \p key that was provided with a string value of \p value. */
void ModelInformation::storeStringOption(const std::string& key, const std::string& value)
{
	// std::cout << "[DEBUG] New string option \"" << key << "\" = \"" << value << "\"" << std::endl;
	usedOptionKeys.erase(key);
	stringOptions[key] = value;
}

/** Query the value of option \p key and return it as a boolean. */
bool ModelInformation::getBooleanOption(const std::string& key) const
{
	return bool(getIntegerOption(key));
}

/**
 * Query the value of option \p key and return it as a boolean, or a default value if the option was not set.
 *
 * \param key Option name.
 * \param defaultValue Value to be returned in case the option was not used.
 * \param printMessageIfDefault Print a message alerting the user of the fact that the default value was used?
 *
 * \return Option value as a boolean, or \p defaultValue if the option was not set.
 */
bool ModelInformation::getBooleanOptionWithDefault(const std::string& key, bool defaultValue, bool printMessageIfDefault) const
{
	return bool(getIntegerOptionWithDefault(key, defaultValue, printMessageIfDefault));
}

/** Query the value of option \p key and return it as an integer. */
int ModelInformation::getIntegerOption(const std::string& key) const
{
	int opt;

	if (hasIntegerOption(key))
	{
		opt = integerOptions.at(key);
	}
	else if (hasFloatOption(key))
	{
		opt = std::lround(floatOptions.at(key));
	}
	else if (hasStringOption(key))
	{
		try
		{
			opt = std::stoi(stringOptions.at(key));
		}
		catch (const std::invalid_argument& e)
		{
			std::stringstream err;
			err << e.what() << ": Invalid argument for integer option \"" << key << "\": \"" << stringOptions.at(key) << "\".";
			throw std::runtime_error(err.str());
		}
		catch (const std::out_of_range& e)
		{
			std::stringstream err;
			err << e.what() << ": Out-of-range error for integer option \"" << key << "\": \"" << stringOptions.at(key) << "\".";
			throw std::runtime_error(err.str());
		}
	}
	else
		throw std::runtime_error(std::string("Cannot determine integer value for option: ") + key);

	//std::cout << "[DEBUG] Get integer option \"" << key << "\" = " << opt << std::endl;
	usedOptionKeys.insert(key);
	return opt;
}

/** Query the value of option \p key and return it as a floating-point number. */
double ModelInformation::getFloatOption(const std::string& key) const
{
	double opt;
	if (hasFloatOption(key))
	{
		opt = floatOptions.at(key);
	}
	else if (hasIntegerOption(key))
	{
		opt = float(integerOptions.at(key));
	}
	else if (hasStringOption(key))
	{
		try
		{
			opt = std::stod(stringOptions.at(key));
		}
		catch (const std::invalid_argument& e)
		{
			std::stringstream err;
			err << e.what() << ": Invalid argument for float option \"" << key << "\": \"" << stringOptions.at(key) << "\".";
			throw std::runtime_error(err.str());
		}
		catch (const std::out_of_range& e)
		{
			std::stringstream err;
			err << e.what() << ": Out-of-range error for float option \"" << key << "\": \"" << stringOptions.at(key) << "\".";
			throw std::runtime_error(err.str());
		}
	}
	else
		throw std::runtime_error("Cannot determine float value for option: "s + key);

	//std::cout << "[DEBUG] Get float option \"" << key << "\" = " << opt << std::endl;
	usedOptionKeys.insert(key);
	return opt;
}

/** Query the value of option \p key and return it as a string. */
std::string ModelInformation::getStringOption(const std::string& key) const
{
	std::string opt;
	if (hasStringOption(key))
	{
		opt = stringOptions.at(key);
	}
	else if (hasIntegerOption(key))
	{
		opt = std::to_string(integerOptions.at(key));
	}
	else if (hasFloatOption(key))
	{
		opt = std::to_string(floatOptions.at(key));
	}
	else
		throw std::runtime_error(std::string("Cannot determine string value for option: ") + key);

	//std::cout << "[DEBUG] Get string option \"" << key << "\" = " << opt << std::endl;
	usedOptionKeys.insert(key);
	return opt;
}

/**
 * Query the value of option \p key and return it as an integer, or a default value if the option was not set.
 *
 * \param key Option name.
 * \param defaultValue Value to be returned in case the option was not used.
 * \param printMessageIfDefault Print a message alerting the user of the fact that the default value was used?
 *
 * \return Option value as an integer, or \p defaultValue if the option was not set.
 */
int ModelInformation::getIntegerOptionWithDefault(const std::string& key, int defaultValue, bool printMessageIfDefault) const
{
	if (!hasOption(key))
	{
		if (printMessageIfDefault)
		{
			std::cout << "No value given for integer option \"" << key << "\". Using default value: " << defaultValue << std::endl;
		}
		return defaultValue;
	}

	return getIntegerOption(key);
}

/**
 * Query the value of option \p key and return it as a floating-point number, or a default value if the option was not set.
 *
 * \param key Option name.
 * \param defaultValue Value to be returned in case the option was not used.
 * \param printMessageIfDefault Print a message alerting the user of the fact that the default value was used?
 *
 * \return Option value as a floating-point number, or \p defaultValue if the option was not set.
 */
double ModelInformation::getFloatOptionWithDefault(const std::string& key, double defaultValue, bool printMessageIfDefault) const
{
	if (!hasOption(key))
	{
		if (printMessageIfDefault)
		{
			std::cout << "No value given for float option \"" << key << "\". Using default value: " << defaultValue << std::endl;
		}
		return defaultValue;
	}

	return getFloatOption(key);
}

/**
 * Query the value of option \p key and return it as a string, or a default value if the option was not set.
 *
 * \param key Option name.
 * \param defaultValue Value to be returned in case the option was not used.
 * \param printMessageIfDefault Print a message alerting the user of the fact that the default value was used?
 *
 * \return Option value as a string, or \p defaultValue if the option was not set.
 */
std::string ModelInformation::getStringOptionWithDefault(const std::string& key, const std::string& defaultValue, bool printMessageIfDefault) const
{
	if (!hasOption(key))
	{
		if (printMessageIfDefault)
		{
			std::cout << "No value given for string option \"" << key << "\". Using default value: " << defaultValue << std::endl;
		}
		return defaultValue;
	}

	return getStringOption(key);
}

/** Check if the user-provided parameter lookup contains values for option \p key. */
bool ModelInformation::hasLookupFor(const std::string& key) const
{
	return parameterLookup.contains(key);
}

/**
 * Interpolate the user-provided parameter lookup to get the value of a given option at given point in time.
 *
 * \param key Option name.
 * \param timestamp Unix time stamp.
 *
 * \returns Interpolated value for option \p key at the time defined by the given \p timestamp.
 *
 */
double ModelInformation::getLookupValueFor(const std::string& key, time_t timestamp) const
{
	if (!parameterLookup.contains(key))
		throw std::runtime_error("No lookup available for parameter \""s + key + "\"!"s);

	return FileHelper::interpolateTimestampLookup(parameterLookup.at(key), timestamp, false);
}

/**
 * Get an option value either from a specified value, or a time-dependent lookup, or a default value.
 *
 * Source for parameter value (with decreasing priority):
 *   1. numeric value declared for option, e.g. in parameter string or dictionary passed to Solarprop instance.
 *   2. lookup value interpolated from user-provided file specified by `parameterLookup` option, interpolated to given \p timestamp, if option is set to "lookup".
 *   3. the value passed as \p defaultValue.
 *
 * For example, during model setup, we want to extract the value of the `kappaScaling` parameter that is multiplied
 * to the value of the diffusion coefficient by calling
 *
 *     auto timestamp = modelInfo.getTimestamp();
 *     double ks = modelInfo.getFloatParameterFromOptionOrLookupOrDefault("kappaScaling", 1.0, timestamp);
 *
 * Then at run time, if the user does not set a value for the `kappaScaling` option, `ks` will be 1. The user can instead
 * set the `kappaScaling` option to a certain value, e.g.:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~{.py}
 * parameters = { 'kappaScaling' = 0.9 }
 * s = solarprop.Solarprop(parameters)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Finally, the user can provide a time-dependent lookup file and direct this function to interpolate the value of `kappaScaling`
 * based on this file:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~{.py}
 * parameters = { 'kappaScaling' = 'lookup',
 *                'parameterLookup' = '/path/to/kappa_lookup.txt' }
 * s = solarprop.Solarprop(parameters)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * where the user-provided lookup file, here called `kappa_lookup.txt`, must have the following format:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * timestamp  kappaScaling
 * 1306800000 0.95
 * 1307145600 0.97
 * 1307491200 0.995
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * and so on.
 *
 * \sa FileHelper::readParameterLookupFile()
 *
 */
double ModelInformation::getFloatParameterFromOptionOrLookupOrDefault(const std::string& key, double defaultValue, time_t timestamp,
								      bool printMessageIfDefault) const
{
	// case 3.
	if (!hasOption(key))
		return getFloatOptionWithDefault(key, defaultValue, printMessageIfDefault);

	// case 2.
	if (getStringOption(key) == "lookup")
	{
		return getLookupValueFor(key, timestamp);
	}

	// case 1.
	return getFloatOption(key);
}

/**
 * Check for any unused options and print warnings for them.
 *
 * \returns \p true if no unused options are found.
 */
bool ModelInformation::checkForUnusedOptions() const
{
	int nUnusedOptions = 0;

	for (const auto& [key, value] : integerOptions)
	{
		if (!usedOptionKeys.contains(key))
		{
			++nUnusedOptions;
			std::cout << "WARNING: Option \"" << key << "\" with value " << value << " was never used!" << std::endl;
		}
	}

	for (const auto& [key, value] : floatOptions)
	{
		if (!usedOptionKeys.contains(key))
		{
			++nUnusedOptions;
			std::cout << "WARNING: Option \"" << key << "\" with value " << value << " was never used!" << std::endl;
		}
	}

	for (const auto& [key, value] : stringOptions)
	{
		if (!usedOptionKeys.contains(key))
		{
			++nUnusedOptions;
			std::cout << "WARNING: Option \"" << key << "\" with value \"" << value << "\" was never used!" << std::endl;
		}
	}

	return (nUnusedOptions == 0);
}

/**
 * Extract option value for the scale \f$ \kappa_0 \f$ of the diffusion tensor.
 *
 * If the option `kappa0` is set, its value (to be given in units of \f$ (\mathrm{cm}^2\,\mathrm{s}^{-1}\,\mathrm{GV}^{-1}) \f$) will be returned.
 * If the option `kappaFile` is set to a lookup file containing time-dependent values of \f$ \kappa_0 \f$ instead, the value is interpolated
 * based on the current time stamp.
 *
 * If neither of these options is used, a value of zero will be returned, signalling the model to use its default value
 * of the diffusion coefficient.
 *
 * \sa FileHelper::readKappaFile()
 *
 * \returns Value of \f$ \kappa_0 \f$, or zero if none of the relevant options is set.
 */
double ModelInformation::getDiffusionCoefficientScale() const
{
	// Set default (to be set by models later)
	double kappaZero = 0.;

	// Get scale of diffusion coefficient
	if (hasOption("kappa0"))
	{
		if (hasOption("kappaFile"))
			std::cout << "WARNING: Mutually exclusive options \"kappa0\" and \"kappaFile\" given. Diffusion coefficient scale will be taken from \"kappa0\" option." << std::endl;

		// 1. either set it directly, or ...
		kappaZero = getFloatOption("kappa0") * cm2/s/GV;
		std::cout << "Overriding diffusion constant: " << kappaZero / (cm2/s/GV) << " cm^2 s^-1 GV^-1" << std::endl;
	}
	else if (hasOption("kappaFile"))
	{
		// 2. ... evaluate it from given time series, using linear interpolation
		auto timestamp = getTimestamp();

		double kappaZero = FileHelper::interpolateTimestampLookup(kappaLookup, timestamp);
		auto prec = std::cout.precision();
		std::cout << "Interpolated kappa0 value at time " << std::setprecision(12) << timestamp << std::setprecision(prec)
			  << ": " << kappaZero / (cm2/s/GV) << " cm^2/s/GV" << std::endl;
	}

	return kappaZero;
}

/**
 * Get magnetic field amplitude and solar wind speed for current time stamp from OMNIweb data.
 *
 * @return Tuple of magnetic field amplitude and solar wind speed.
 */
std::tuple<double, double> ModelInformation::getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebData() const
{
	return getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebDataFor(time_t(getTimestamp()));
}

/**
 * Get magnetic field amplitude and solar wind speed for given \p timestamp from OMNIweb data.
 *
 * @return Tuple of magnetic field amplitude and solar wind speed.
 */
std::tuple<double, double> ModelInformation::getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebDataFor(time_t timestamp) const
{
	if (!bmagLookup.size() || !vswLookup.size())
		throw std::runtime_error("OMNIweb data lookup is empty! Make sure \""s + omniFileOption + "\" option is given.");

	double bmag = FileHelper::interpolateTimestampLookup(bmagLookup, double(timestamp));
	double vsw = FileHelper::interpolateTimestampLookup(vswLookup, double(timestamp));

	return {bmag, vsw};
}

/** Get lookup value of the smoothed sunspot number for the current time stamp. */
double ModelInformation::getSunspotNumber() const
{
	return getSunspotNumberFor(time_t(getTimestamp()));
}

/** Get lookup value of the smoothed sunspot number for the given \p timestamp. */
double ModelInformation::getSunspotNumberFor(time_t timestamp) const
{
	if (!ssnLookup.size())
		throw std::runtime_error("Sunspot number lookup is empty! Make sure \""s + ssnFileOption + "\" option is given.");

	return FileHelper::interpolateTimestampLookup(ssnLookup, double(timestamp));
}

/** Extract cosmic-ray particle properties (charge and mass) from user options \c mass and \c charge. */
ParticleProperties ModelInformation::getParticleProperties() const
{
	// Read particle properties
	ParticleProperties spec;
	spec.mass = getFloatOption("mass") * GeV / c2;
	spec.charge = getIntegerOptionWithDefault("charge", 1) * eplus;

	return spec;
}

/**
 * Calculate Carrington rotation number from date.
 *
 * @param year The year.
 * @param month The month.
 * @param day The day.
 *
 * @return Carrington rotation number.
 */
int ModelInformation::calculateCR(int year, int month, int day)
{
	// Calculate Julian date
	int A = floor(year / 100.);
	int B = 2 - A + floor(A / 4.);
	double jd = 0;
	if (month > 2)
	{
		jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + B + day - 1524.5;
	}
	else
	{
		jd = floor(365.25 * ((year - 1.) + 4716)) + floor(30.6001 * ((month + 12) + 1)) + B + day - 1524.5;
	}
	// Calculate Carrington rotation
	return (int)(jd - 2398140.227) / (CarringtonRotationDuration / days);
}

/**
 * Calculate average tilt angle from Carrington rotations between two dates.
 *
 * @return Averaged tilt angle.
 */
double ModelInformation::averageAngle(int startYear, int startMonth, int startDay, int endYear, int endMonth, int endDay) const
{
	double average = 0.0;
	// Convert Julian date to Carrington rotation
	int carringtonStart = calculateCR(startYear, startMonth, startDay);
	int carringtonEnd = calculateCR(endYear, endMonth, endDay);
	// Sum and calculate average
	int counter = 0;
	for (int i = carringtonStart; i <= carringtonEnd; i++)
	{
		average += tiltAngleData.at(i);
		counter++;
	}
	average /= counter;

	return average;
}

/**
 * Average modulation potential from neutron monitor data over given timerange.
 *
 * @return Averaged modulation potential.
 */
double ModelInformation::averageNeutronMonitorData(int startYear, int startMonth, int endYear, int endMonth) const
{
	double average = 0;
	int counter = 0;
	for (int i = startYear; i <= endYear; i++)
	{
		if (i == startYear && i == endYear)
		{
			for (int k = startMonth; k <= endMonth; k++)
			{
				average += nmData.at({i, k});
				counter++;
			}
		}
		else if (i == startYear)
		{
			for (int k = startMonth; k <= 12; k++)
			{
				average += nmData.at({i, k});
				counter++;
			}
		}
		else if (i == endYear)
		{
			for (int k = 1; k <= endMonth; k++)
			{
				average += nmData.at({i, k});
				counter++;
			}
		}
		else
		{
			for (int k = 1; k <= 12; k++)
			{
				average += nmData.at({i, k});
				counter++;
			}
		}
	}
	average /= counter;

	return average;
}

/** Calculate the Carrington Rotation number for the current time stamp. */
int ModelInformation::getCarringtonRotation() const
{
	const auto [year, month, day] = getYMD();
	return calculateCR(year, month, day);
}

/**
 * Extract heliospheric polarity and phase from the user options.
 *
 * The user can either set the `polarity` and optionally, the `phase` options to
 * set the values manually, or they will be extracted based on the current
 * time stamp.
 *
 * Possible values for the
 * - polarity: \c +1 for \f$ A>0 \f$, \c -1 for \f$ A<0 \f$,
 * - heliospheric phase: \c +1 ascending, \c -1 descending, \c 0 unknown.
 *
 * @return Pair of polarity and phase.
 */
std::pair<int, int> ModelInformation::getPolarityAndPhase() const
{
	// Determine polarity from Carrington rotation if not given
	int polarity = 0;
	int phase = 0;
	if (hasOption("polarity"))
	{
		polarity = getIntegerOption("polarity");
		phase = getIntegerOptionWithDefault("phase", 0, false);
	}
	else
	{
		std::tie(polarity, phase) = getPolarityAndPhaseFor(getCarringtonRotation());
	}

	if (polarity != 1 && polarity != -1)
	{
		std::stringstream err;
		err << "Invalid polarity: " << polarity;
		throw std::runtime_error(err.str());
	}

	return {polarity, phase};
}

/**
 * Calculate heliospheric polarity and phase for a given Carrington Rotation.
 *
 * Possible values for the
 * - polarity: \c +1 for \f$ A>0 \f$, \c -1 for \f$ A<0 \f$,
 * - heliospheric phase: \c +1 ascending, \c -1 descending, \c 0 unknown.
 *
 * The time intervals for these values have been roughly read off a plot of the
 * smoothed sunspot number vs time. They are valid for years 1970 to 2026.
 *
 * @param carringtonRotation Carrington Rotation number.
 * @param verbose Also print the results?
 *
 * @return Pair of polarity and phase.
 */
std::pair<int, int> ModelInformation::getPolarityAndPhaseFor(int carringtonRotation, bool verbose)
{
	int polarity = 0;
	int phase = 0;

	// start with solar maximum ca 1959 (CR 1409)
	if (carringtonRotation < 1409)
	{
		throw std::runtime_error("Cannot determine polarity for CR "s + std::to_string(carringtonRotation) + "!"s);
	}
	else if (carringtonRotation >= 1409 && carringtonRotation < 1536) // A<0 until ca mid 1968 (CR 1536)
	{
		polarity = -1;

		// descending activity up to solar minimum ca April 1964
		phase = 1 - 2 * (carringtonRotation < 1481);
	}
	else if (carringtonRotation >= 1536 && carringtonRotation < 1689) // A>0 until 1979-11-30
	{
		polarity = 1;

		// descending activity up to solar minimum ca 1976.5
		phase = 1 - 2 * (carringtonRotation < 1643);
	}
	else if (carringtonRotation >= 1689 && carringtonRotation < 1837) // A<0 until 1990-12-19
	{
		polarity = -1;

		// descending activity up to solar minimum ca 1986.5
		phase = 1 - 2 * (carringtonRotation < 1777);
	}
	else if (carringtonRotation >= 1837 && carringtonRotation < 1955) // A>0 until 1999-10-11
	{
		polarity = 1;

		// descending activity up to solar minimum ca 1996.5
		phase = 1 - 2 * (carringtonRotation < 1911);
	}
	else if (carringtonRotation >= 1955 && carringtonRotation < 2143) // A<0 until 2013-10-26 (CR 2143)
	{
		// estimated reversal in October 2013 according to
		// X. Sun et al., https://dx.doi.org/10.1088/0004-637X/798/2/114

		polarity = -1;

		// descending activity up to solar minimum ca 2009
		phase = 1 - 2 * (carringtonRotation < 2079);
	}
	else if (carringtonRotation >= 2143 && carringtonRotation < 2295) // A>0 until estimated reversal during solar cycle 25: March 2025
	{
		polarity = 1;

		// descending activity up to solar minimum ca 2020
		phase = 1 - 2 * (carringtonRotation < 2226);
	}
	else
	{
		polarity = -1; // next: A<0
		phase = -1;    // first: descending phase
	}

	if (verbose)
	{
		std::string polarityWord = "invalid";
		if (polarity == +1)
			polarityWord = "A>0";
		else if (polarity == -1)
			polarityWord = "A<0";

		std::string phaseWord = "unknown";
		if (phase == +1)
			phaseWord = "ascending";
		else if (phase == -1)
			phaseWord = "descending";
		std::cout << "Computed polarity for CR " << carringtonRotation << ": " << polarityWord << " (phase: " << phaseWord << ")" << std::endl;
	}

	return {polarity, phase};
}

/**
 * Extract option value for the tilt angle of the heliospheric current sheet.
 *
 * If the option `angle` is set, its value (to be given in units of degrees) will be returned.
 * If the option `angleFile` is set to a lookup file containing time-dependent values of the tilt angle instead,
 * the value is calculated based on the current time stamp: If the time stamp is defined using only the `year` option, the
 * values per Carrington Rotation are averaged over the given year. If the time stamp is defined using the `timestamp` options or using both
 * `year` and `month` options, the value valid for the current Carrington Rotation is returned.
 *
 * If neither of these options is used, a default value of zero will be returned, indicating a flat current sheet.
 *
 * @return Value of the current sheet tilt angle.
 */
double ModelInformation::getTiltAngle() const
{
	// Set default
	double angle = 0.;

	// Check if file with tilt angle information is provided. If yes, take tilt angle from file
	if (hasOption("angle"))
	{
		if (hasOption("angleFile"))
			std::cout << "WARNING: Mutually exclusive options \"angle\" and \"angleFile\" given. Tilt angle will be taken from \"angle\" option."
				  << std::endl;

		angle = getFloatOption("angle") * degree;
		std::cout << "Override tilt angle: " << angle / deg << " deg" << std::endl;
	}
	else if (hasOption("angleFile"))
	{
		// If year but no month is provided, average over the whole year.
		if (hasOption("year") && !hasOption("month") && !hasOption("timestamp"))
		{
			int year = getIntegerOption("year");
			angle = averageAngle(year, 1, 1, year, 12, 31);
			std::cout << "Computed average tilt angle for year " << year << ": " << angle / deg << " deg" << std::endl;
		}
		else
		{
			angle = getTiltAngleFor(getCarringtonRotation());
		}
	}

	return angle;
}

/**
 * Get lookup value for the tilt angle of the heliospheric current sheet for a given Carrington Rotation.
 *
 * \param carringtonRotation Carrington Rotation number.
 * \param verbose Also print the result?
 *
 * @return Value of the current sheet tilt angle.
 */
double ModelInformation::getTiltAngleFor(int carringtonRotation, bool verbose) const
{
	if (tiltAngleData.contains(carringtonRotation))
	{
		double angle = tiltAngleData.at(carringtonRotation);
		if (verbose)
			std::cout << "Tilt angle for CR " << carringtonRotation << ": " << angle / deg << " deg" << std::endl;
		return angle;
	}

	std::cout << "WARNING: CR " << carringtonRotation << " outside time range for tilt angle data. Assuming flat HCS." << std::endl;
	return 0. * degree;
}

/**
 * Extract option value for the modulation potential.
 *
 * If the option `nmValue` is set, its value (to be given in units of \f$ (\mathrm{MV}) \f$) will be returned.
 * If the option `nmFile` is set to a lookup file containing time-dependent values of the modulation potential instead,
 * the value is calculated based on the current time stamp: If the time stamp is defined using only the `year` option, the
 * monthly values are averaged over the given year. If the time stamp is defined using the `timestamp` options or using both
 * `year` and `month` options, the value valid for the given month is returned.
 *
 * If neither of these options is used, a default value of \f$ 500\,\mathrm{MV} \f$ will be returned.
 *
 * @return Value of the modulation potential.
 */
double ModelInformation::getModulationPotential() const
{
	// Set default
	double nmModulationPotential = 500.0 * MV;

	// Get neutron monitor data
	if (hasOption("nmValue"))
	{
		if (hasOption("nmFile"))
			std::cout << "WARNING: Mutually exclusive options \"nmValue\" and \"nmFile\" given. Modulation potential will be taken from "
				     "\"nmValue\" option."
				  << std::endl;

		nmModulationPotential = getFloatOption("nmValue") * MV;
		std::cout << "Override neutron monitor value: " << nmModulationPotential / MV << " MV" << std::endl;
		return nmModulationPotential;
	}
	else if (hasOption("nmFile"))
	{
		// If year but no month is provided, average over the whole year.
		if (hasOption("year") && !hasOption("month") && !hasOption("timestamp"))
		{
			int year = getIntegerOption("year");
			nmModulationPotential = averageNeutronMonitorData(year, 1, year, 12);
			std::cout << "Averaged neutron monitor value for year " << year << ": " << nmModulationPotential / MV << " MV" << std::endl;
		}
		else
		{
			auto ymd = getYMD();
			const auto [year, month, day] = ymd;
			nmModulationPotential = averageNeutronMonitorData(year, month, year, month);
			std::cout << "Computed neutron monitor value for date " << ymd << ": " << nmModulationPotential / MV << " MV" << std::endl;
		}

		return nmModulationPotential;
	}

	std::cout << "Using default neutron monitor value: " << nmModulationPotential / MV << " MV" << std::endl;
	return nmModulationPotential;
}

/**
 * Calculate the timestamp (unix time: number of seconds since epoch) that quantifies the point in time specified by user options.
 *
 * There are two ways the user may set the current time for their model:
 *  1. by providing the `timestamp` option to specify the timestamp (unix time) directly,
 *  2. or by setting the `year` option, and maybe the `month` and `day` options.
 *
 *  In the second case, if only the year is given, the timestamp will be 1 July of that year. If `month` and `year` are given, the
 *  timestamp will be the 15th day of the month. If `day`, `month`, and `year` are given, the timestamp will be the given calendar
 *  day at midnight UTC.
 */
double ModelInformation::getTimestamp() const
{
	if (hasOption("timestamp"))
	{
		if (hasOption("year") || hasOption("month") || hasOption("day"))
			std::cout << "WARNING: Mutually exclusive options \"timestamp\" and at least one of \"year\", \"month\", and \"day\" given. Will use "
				     "\"timestamp\" value."
				  << std::endl;

		return getFloatOption("timestamp");
	}

	if (!hasOption("year"))
		throw std::runtime_error(std::string("Cannot determine timestamp. Need at least option \"year\"."));

	// Determine date
	int year = getIntegerOption("year");

	double timestamp = 0.0;

	if (hasOption("month"))
	{
		int month = getIntegerOption("month");
		int day = getIntegerOptionWithDefault("day", 15);
		timestamp = FileHelper::calculateTimestamp(year, month, day);
	}
	else
	{
		timestamp = FileHelper::calculateTimestamp(year, 7, 1);
	}

	return timestamp;
}

/** Get the timestamp as a YearMonthDay tuple. */
YearMonthDay ModelInformation::getYMD() const
{
	time_t ttt(getTimestamp());
	auto* ptm = std::gmtime(&ttt);
	return std::make_tuple(ptm->tm_year + 1900, ptm->tm_mon + 1, ptm->tm_mday);
}

/**
 * Extract heliospheric polarity and phase from the user options.
 *
 * See getPolarityAndPhase().
 */
int ModelInformation::getPolarity() const
{
	auto polarityAndPhase = getPolarityAndPhase();
	return polarityAndPhase.first;
}

/** Streaming operator for a YearMonthDay tuple. */
std::ostream& operator<<(std::ostream& o, const YearMonthDay& ymd)
{
	std::stringstream s;
	s << std::get<0>(ymd) << "-" << std::setw(2) << std::setfill('0') << std::get<1>(ymd) << "-" << std::get<2>(ymd);
	o << s.str();
	return o;
}
