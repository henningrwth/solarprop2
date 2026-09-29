#ifndef FILEHELPER_H_
#define FILEHELPER_H_

#include <map>
#include <string>
#include <vector>

/// \file filehelper.h

/// \typedef std::map<std::pair<int, int>, double> YearMonthValueMap
/// \brief Map a year and month to a value.
///
/// The key is a pair of year and month (first entry: year, second entry: month).
using YearMonthValueMap = std::map<std::pair<int, int>, double>;

/// \typedef std::map<double, double> TimestampValueMap
/// \brief Parameter lookup: Map a unix timestamp to a value.
using TimestampValueMap = std::map<double, double>;

/// \typedef std::map<std::string, TimestampValueMap> TimedependentParameters
/// \brief Map parameter name to a parameter lookup (map of timestamp to value).
using TimedependentParameters = std::map<std::string, TimestampValueMap>;

/**
 * Static functions for reading data files in various formats.
 */
class FileHelper
{
	public:
	FileHelper() {}
	~FileHelper() {}

	static double calculateTimestamp(int year, int month, int day, int hour = 0);

	static double interpolateTimestampLookup(const TimestampValueMap& lookup, double timestamp, bool extrapolationIsError = true);

	template <typename T1, typename T2>
	static std::map<T1, T2> readFile(const std::string& input, unsigned int first, unsigned int second);

	template <typename T>
	static std::vector<T> readColumnInFile(const std::string&, unsigned int);

	static std::vector<double> readNumbersFromColumnInFile(const std::string&, unsigned int, double);
	static std::vector<std::vector<double>> readNumbersFromFileColumns(const std::string& input, const std::vector<unsigned int>& columnNumbers,
									   const std::vector<double>& units);
	static std::map<std::pair<int, int>, double> readNeutronMonitorData(const std::string&);
	static std::tuple<TimestampValueMap, TimestampValueMap> readOmniWebData(const std::string&);
	static TimestampValueMap readSunspotNumberData(const std::string&);
	static TimestampValueMap readKappaFile(const std::string&);
	static TimedependentParameters readParameterLookupFile(const std::string&);
};

#endif
