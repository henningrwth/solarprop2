#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>

#include "filehelper.h"
#include "model.h"

/// \file properties.h

/// \typedef std::tuple<int, int, int> YearMonthDay
/// \brief A tuple of `(year, month, day)`.
using YearMonthDay = std::tuple<int, int, int>;

/// Cosmic-ray particle properties (charge and mass).
struct ParticleProperties
{
	/// Charge of the particle
	int charge;

	/// Mass of the particle
	double mass;
};

/// Parameters for smearing of particle energies.
struct EnergySmearingProperties
{

	/// Lower kinetic energy bound for smearing
	double T_low = 0.0;

	/// Upper kinetic energy bound for smearing
	double T_up = 0.0;

	/// Approximate spectral index used for energy smearing
	double spectralIndex = 0.0;

	/// %Random seed to be used for smearing
	unsigned int randomSeed = 999999;
};

/**
 * %Model information.
 *
 * This class serves three purposes:
 * - Extract the values of options (flags or parameters) provided by the user.
 * - Read lookup files of time-dependent data and interpolate or extract the values for the relevant point(s) in time.
 * - Create model implementation classes, based on the value of the \c model option provided by the user.
 *
 * Options provided by the user are stored in three different maps depending on the type of the option value initially
 * received in the python-dictionary passed to the Solarprop interface class (integer, floating-point or string). This
 * is useful because option values may be defined in config files or on the command line and as such be treated as
 * strings.
 * However, when an option value is queried here using one of the `get...Option...()` functions, type conversions are
 * performed to provide the type implied by the function name, e.g., getBooleanOption() will return a \c bool.
 *
 * The class keeps track of the options whose values have been queried so far (usually during model setup). Calling
 * checkForUnusedOptions() towards the end of the program is a good idea because the function will print a warning
 * for each unused option, which usually means that the user misspelled an option name, or used an option not
 * applicable for the chosen model.
 *
 * In most cases, the time-dependent values are interpolated or read from the lookups for the current time stamp
 * (Unix time: number of seconds since the epoch), which can be set in different ways, see getTimestamp().
 */
class ModelInformation
{

	public:
	ModelInformation() = default;
	void readDatafiles();

	std::unique_ptr<Model> createModel() const;

	bool hasOption(const std::string& key) const;

	void storeIntegerOption(const std::string& key, int value);
	void storeFloatOption(const std::string& key, double value);
	void storeStringOption(const std::string& key, const std::string& value);

	bool getBooleanOption(const std::string& key) const;
	bool getBooleanOptionWithDefault(const std::string& key, bool defaultValue, bool printMessageIfDefault = true) const;
	int getIntegerOption(const std::string& key) const;
	int getIntegerOptionWithDefault(const std::string& key, int defaultValue, bool printMessageIfDefault = true) const;
	double getFloatOption(const std::string& key) const;
	double getFloatOptionWithDefault(const std::string& key, double defaultValue, bool printMessageIfDefault = true) const;
	std::string getStringOption(const std::string& key) const;
	std::string getStringOptionWithDefault(const std::string& key, const std::string& defaultValue, bool printMessageIfDefault = true) const;

	bool hasLookupFor(const std::string& key) const;
	double getLookupValueFor(const std::string& key, time_t timestamp) const;
	double getFloatParameterFromOptionOrLookupOrDefault(const std::string& key, double defaultValue, time_t timestamp,
							    bool printMessageIfDefault = true) const;

	bool checkForUnusedOptions() const;

	ParticleProperties getParticleProperties() const;
	double getTimestamp() const;
	YearMonthDay getYMD() const;
	int getCarringtonRotation() const;
	int getPolarity() const;
	std::pair<int, int> getPolarityAndPhase() const;
	double getTiltAngle() const;
	double getTiltAngleFor(int carringtonRotation, bool verbose = true) const;
	double getModulationPotential() const;
	double getDiffusionCoefficientScale() const;
	std::tuple<double, double> getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebData() const;
	std::tuple<double, double> getMagneticFieldAmplitudeAndSolarWindSpeedFromOmniWebDataFor(time_t timestamp) const;
	double getSunspotNumber() const;
	double getSunspotNumberFor(time_t timestamp) const;

	static int calculateCR(int year, int month, int day);
	static std::pair<int, int> getPolarityAndPhaseFor(int carringtonRotation, bool verbose = true);

	private:
	/// Store options (map option name to option value) for values initially provided as integers.
	std::map<std::string, int> integerOptions;
	/// Store options (map option name to option value) for values initially provided as floating-point numbers.
	std::map<std::string, double> floatOptions;
	/// Store options (map option name to option value) for values initially provided as strings.
	std::map<std::string, std::string> stringOptions;
	/// Keep track of options whose value has been requested.
	mutable std::set<std::string> usedOptionKeys;

	/// Lookup mapping unix time to smoothed sunspot number (SSN).
	TimestampValueMap ssnLookup;

	/// Option for name of SSN data file.
	///
	/// The name of the option is defined here because it is needed in several places.
	const std::string ssnFileOption = "ssnFile";

	/// Source (filename) of current sunspot number data.
	std::string ssnLookupFilename;

	/// Lookup mapping unix time to magnetic field amplitude at Earth.
	TimestampValueMap bmagLookup;

	/// Lookup mapping unix time to solar wind speed.
	TimestampValueMap vswLookup;

	/// Option for name of OMNIweb data file.
	///
	/// The name of the option is defined here because it is needed in several places.
	const std::string omniFileOption = "omniFile";

	/// Source (filename) of current OMNIweb data.
	std::string omniFilename;

	/// Lookup mapping Carrington Rotation number to HCS tilt angle.
	std::map<int, double> tiltAngleData;

	/// Source (filename) of current sheet tilt angle data.
	std::string angleFilename;

	/// Lookup mapping year and month to modulation potential inferred from neutron monitor data.
	std::map<std::pair<int, int>, double> nmData;

	/// Source (filename) of current neutron monitor data.
	std::string nmFilename;

	/// Lookup mapping unix time to kappa scale, based on user-provided input.
	TimestampValueMap kappaLookup;

	/// Source (filename) of current kappa scale data.
	std::string kappaFilename;

	/// Lookups (time to parameter value) for a set of parameters.
	TimedependentParameters parameterLookup;

	/// Source (filename) of current parameter lookup.
	std::string parameterLookupFilename;

	double averageNeutronMonitorData(int, int, int, int) const;
	double averageAngle(int, int, int, int, int, int) const;

	bool hasIntegerOption(const std::string& key) const;
	bool hasFloatOption(const std::string& key) const;
	bool hasStringOption(const std::string& key) const;
};

std::ostream& operator<<(std::ostream& o, const YearMonthDay& ymd);

#endif
