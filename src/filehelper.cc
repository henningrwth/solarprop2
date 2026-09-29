#include "filehelper.h"

#include "stringtools.h"
#include "units.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std::string_literals;
using namespace Units;

/**
 * Read input file and store result in a map.
 *
 * Lines starting with '#' will be skipped.
 *
 * @param input File name of the input file.
 * @param first Number of the first column to be imported (start counting at 1).
 * @param second Number of the second column to be imported (start counting at 1).
 *
 * @return Map of the columns from the input file.
 */
template <typename T1, typename T2>
std::map<T1, T2> FileHelper::readFile(const std::string& input, unsigned int first, unsigned int second)
{
	std::string line;
	std::string dummy;
	std::ifstream file(input);

	std::map<T1, T2> data;

	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			if (line.starts_with('#'))
				continue;

			std::stringstream s(line);
			T1 key;
			T2 value;
			// Skip all columns before first
			for (unsigned int i = 1; i < first; i++)
			{
				s >> dummy;
			}
			// Store this column as key
			s >> key;
			// Skip all columns before second
			for (unsigned int i = first + 1; i < second; i++)
			{
				s >> dummy;
			}
			// Store this column as value
			s >> value;

			// Attach key value pair to data
			data[key] = value;
		}
		//std::cout << "All data read from file " << input << std::endl;
		file.close();
	}
	else
		throw std::runtime_error("Cannot open file: "s + input);

	return data;
}

/**
 * Read given column from input file and store result in a vector.
 *
 * Lines starting with '#' will be skipped.
 *
 * @param input File name of the input file.
 * @param col Number of the column to be imported (start counting at 1).
 *
 * @return Vector of the values from the input file, or empty or partially empty vector if column number is invalid.
 */
template <typename T>
std::vector<T> FileHelper::readColumnInFile(const std::string& input, unsigned int col)
{
	std::string line;
	std::ifstream file(input);

	std::vector<T> data;
	if (col < 1)
		return data;

	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			if (line.starts_with('#'))
				continue;

			const auto& tokens = split(line, " ");
			if (tokens.size() < col)
				continue;

			std::stringstream s(tokens[col - 1]);
			T value;
			s >> value;

			data.push_back(value);
		}
		file.close();
	}
	else
		throw std::runtime_error("Cannot open file: "s + input);

	return data;
}

/**
 * Read given column from input file, multiply the values by the given unit, and store result in a vector.
 *
 * Lines starting with '#' will be skipped.
 *
 * @param input File name of the input file.
 * @param col Number of the column to be imported (start counting at 1).
 * @param unit Unit the values are given in.
 *
 * @return Vector of the values from the input file, or empty or partially empty vector if column number is invalid.
 */
std::vector<double> FileHelper::readNumbersFromColumnInFile(const std::string& input, unsigned int col, double unit)
{
	std::vector<double> data = readColumnInFile<double>(input, col);
	for (double& d : data)
		d *= unit;

	return data;
}

/**
 * Read data from one or more columns of an input file in one go.
 *
 * Lines starting with '#' will be skipped.
 *
 * @param input File name of the input file.
 * @param columnNumbers Vector of column numbers to be imported (start counting at 1).
 * @param units Vector of corresponding units, one for each column.
 *
 * @return Vector of value vectors, one for each column.
 */
std::vector<std::vector<double>> FileHelper::readNumbersFromFileColumns(const std::string& input, const std::vector<unsigned int>& columnNumbers,
									const std::vector<double>& units)
{
	std::string line;
	std::ifstream file(input);

	std::vector<std::vector<double>> alldata;

	assert(columnNumbers.size() == units.size());

	for (unsigned int i = 0; i < columnNumbers.size(); ++i)
	{
		if (columnNumbers[i] < 1)
		{
			throw std::runtime_error("Illegal column number "s + std::to_string(columnNumbers[i]));
		}
		alldata.push_back(std::vector<double>());
	}

	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			// Lines starting with '#' will be skipped.
			if (line.starts_with('#'))
				continue;

			const auto& tokens = split(line, " ");
			//for (unsigned int i = 0; i < tokens.size(); ++i)
			//std::cout << i+1 << ": \"" << tokens[i] << "\" ";
			//std::cout << std::endl;

			for (unsigned int i = 0; i < columnNumbers.size(); ++i)
			{
				auto col = columnNumbers[i];
				if (col > tokens.size())
					continue;
				std::vector<double>& data = alldata[i];
				double v = std::stod(tokens.at(col-1)) * units.at(i);
				//std::cout << col << ": " << "\"" << tokens.at(col-1) << "\" * " << units.at(i) << " = " << v << std::endl;
				data.push_back(v);
			}
		}
		file.close();
	}
	else
		throw std::runtime_error("Cannot open file: "s + input);

	return alldata;
}

// explicit instantiation
template std::map<int, double> FileHelper::readFile<int, double>(const std::string&, unsigned int, unsigned int);

/**
 * Read neutron monitor data from file.
 *
 * Lines starting with '#' will be skipped.
 *
 * For an example of the file format, see the standard neutron monitor data file in
 * \c data/nm.dat .
 *
 * @param input File name for input file that contains neutron monitor data in units of MV.
 *
 * @return Map with modulation potential value for each year and month.
 */
YearMonthValueMap FileHelper::readNeutronMonitorData(const std::string& input)
{
	std::cout << "Reading neutron monitor data from file: " << input << std::endl;

	std::string line;
	std::ifstream file(input);

	std::map<std::pair<int, int>, double> data;

	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			if (line.starts_with('#'))
				continue;

			std::stringstream s(line);
			std::pair<int, int> key;
			// Initialization necessary to avoid bad data in the case of blank lines
			int year = 0;
			double value;
			s >> year;
			// Parse line
			// FIXME first line (year 1936) is not read correctly because it contains '-' entries
			for (unsigned int i = 0; i < 12; i++)
			{
				s >> value;
				key = {year, i + 1};
				data[key] = value * MV;
				//std::cout << "key: " << key.first << "/" << key.second << " value (MV): " << data[key] / MV << std::endl;
			}
		}
		file.close();
	}
	else
		throw std::runtime_error("Cannot open NM file: "s + input);

	return data;
}

/**
 * Convert date to unix time stamp (seconds since epoch).
 *
 * @param year Year (1900 onwards)
 * @param month Month (1 to 12)
 * @param day Day of the month (1 to 31), or Day in Year (1 to 365 or 366) if \c month=1
 * @param hour (optional) Hour (0 to 23)
 *
 * @return Unix time stamp.
 */
double FileHelper::calculateTimestamp(int year, int month, int day, int hour)
{

	tm t;
	t.tm_year = year - 1900;
	t.tm_mon = month - 1;
	t.tm_mday = day; // timegm will wrap months correctly if day is day-of-year
	t.tm_hour = hour;
	t.tm_min = 0;
	t.tm_sec = 0;

	double timestamp = double(timegm(&t));
	// std::cout << year << "-" << month << "-" << day << ":" << hour << " -> " << int(timestamp) << std::endl;
	return timestamp;
}

/**
 * Interpolate a lookup at a given point in time.
 *
 * @param lookup Map of Unix time stamps to values.
 * @param timestamp Unix time stamp at which to interpolate lookup.
 * @param extrapolationIsError Raise exception in case timestamp is outside lookup range?
 *
 * @return Interpolated value. First or last value, respectively, if timestamp is out of range and \p extrapolationIsError is false.
 */
double FileHelper::interpolateTimestampLookup(const TimestampValueMap& lookup, double timestamp, bool extrapolationIsError)
{
	if (lookup.size() < 2)
		throw std::runtime_error("Too few entries in lookup!");

	static const char* fmt = "%d %b %Y";

	auto it = lookup.lower_bound(timestamp);

	if (it == lookup.end())
	{
		double t2 = lookup.rbegin()->first;
		double y2 = lookup.rbegin()->second;

		if (extrapolationIsError)
		{
			std::stringstream msg;
			time_t t2t(t2);
			msg << "Timestamp " << timestamp << " too late for lookup ending at " << std::put_time(std::gmtime(&t2t), fmt) << "!";
			throw std::runtime_error(msg.str());
		}

		// debug message
		//time_t ttt(timestamp);
		//std::cout << "t = " << std::put_time(std::gmtime(&ttt), fmt) << " -> " << y2 << std::endl;
		return y2;
	}


	double t2 = it->first;
	double y2 = it->second;

	if (it == lookup.begin())
	{
		if (timestamp < t2 && extrapolationIsError)
		{
			std::stringstream msg;
			time_t t2t(t2);
			msg << "Timestamp " << timestamp << " too early for lookup starting at " << std::put_time(std::gmtime(&t2t), fmt) << "!";
			throw std::runtime_error(msg.str());
		}

		// debug message
		//time_t ttt(timestamp);
		//std::cout << "t = " << std::put_time(std::gmtime(&ttt), fmt) << " -> " << y2 << std::endl;
		return y2;
	}

	auto pr = std::prev(it);
	double t1 = pr->first;
	double y1 = pr->second;

	double y = (t1 == t2 ? y1 : (y2 * (timestamp - t1) - y1 * (timestamp - t2)) / (t2 - t1));

	// debug message
	// time_t t1t(t1);
	// time_t t2t(t2);
	// time_t ttt(timestamp);
	// std::cout << "t: " << std::put_time(std::gmtime(&t1t), fmt) << " <= " << std::put_time(std::gmtime(&ttt), fmt) << " <= " <<
	// std::put_time(std::gmtime(&t2t), fmt) << " -> " << y1 << " .. " << y2 << " -> " << y << std::endl;

	return y;
}

/**
 * Read data on magnitude of magnetic field and solar wind plasma speed extracted from OmniWEB as a function of time.
 *
 * To extract the data file from OMNIweb, visit https://omniweb.gsfc.nasa.gov/form/dx1.html
 *
 * Then:
 *   - Select "Create file".
 *   - Select resolution: 27-day averaged.
 *   - Enter start date: 19640101.
 *   - Select: Magnetic field -> "IMF Magnitude Avg, nT" and Plasma -> "Flow Speed, km/sec".
 *   - Submit.
 *
 * @return A tuple whose first entry is a map of Unix time stamp to the magnitude of magnetic field, and whose second entry is a map of Unix time stamp to the
 * solar wind plasma speed.
 */
std::tuple<TimestampValueMap, TimestampValueMap> FileHelper::readOmniWebData(const std::string& input)
{
	TimestampValueMap magField;
	TimestampValueMap swSpeed;

	std::string line;
	std::ifstream file(input);

	std::cout << "Reading OMNI web data from file: " << input << std::endl;
	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			if (line.starts_with('#'))
				continue;

			std::stringstream ss(line);
			int year, doy, hour;
			std::string bstr, vstr;
			double b, v;
			ss >> year >> doy >> hour >> bstr >> vstr;

			double ts = calculateTimestamp(year, 1, doy, hour);
			//time_t ttt(ts); // debug
			if (bstr != "999.9"s)
			{
				b = std::stod(bstr) * nanotesla;
				magField[ts] = b;
				//std::cout << "t: " << ts << " " << std::put_time(std::gmtime(&ttt), "%F") << " B = " << b / nanotesla << " nT" << std::endl;
			}

			if (vstr != "9999."s)
			{
				v = std::stod(vstr) * km/s;
				swSpeed[ts] = v;
				//std::cout << "t: " << ts << " " << std::put_time(std::gmtime(&ttt), "%F") << " v = " << v / (km/s) << " km/s" << std::endl;
			}
		}
	}

	return {magField, swSpeed};
}

/**
 * Read time-dependent sunspot number data from file.
 *
 * Lines starting with '#' will be skipped.
 *
 * For an example of the file format, see the standard sunspot number data file in
 * \c data/ssn.dat .
 *
 * @param input File name for input file that contains sunspot number data in SIDC format.
 *
 * @return Map of unix timestamps to sunspot number.
 */
TimestampValueMap FileHelper::readSunspotNumberData(const std::string& input)
{
	std::cout << "Reading SSN data from file: " << input << std::endl;
	auto data = readNumbersFromFileColumns(input, {3, 4}, {1.0, 1.0});
	const std::vector<double>& fracYears = data[0];
	const std::vector<double>& ssnValues = data[1];

	TimestampValueMap ssnLookup;
	for (size_t i = 0; i < fracYears.size(); ++i)
	{
		// 3rd column of SSN data file has fractional year (e.g., 2020.54), we split this into
		// the full year (2020) and the remainder (0.54) is converted to day-of-year.
		double fracYear = fracYears[i];
		double fullYear;
		double yearFraction = std::modf(fracYear, &fullYear);
		long doy = std::lround(yearFraction * 365.);
		double t = calculateTimestamp(std::lround(fullYear), 1, doy);
		//std::cout << i << " " << fracYear << " " << std::lround(fullYear) << "+" << yearFraction << " " << doy << " " << t << " -> " << ssnValues[i] << std::endl;
		ssnLookup[t] = ssnValues[i];
	}

	return ssnLookup;
}

/**
 * Read time-dependent data on diffusion parameter scale from a file.
 *
 * Format of input file: First column contains unix time stamps, second
 * column contains \f$ \kappa_0 \f$ in units of \f$ \mathrm{cm}^2/\mathrm{s}/\mathrm{GV} \f$.
 *
 * @param input File name of input file.
 *
 * @return Map of unix timestamps to \f$ \kappa_0 \f$.
 */
TimestampValueMap FileHelper::readKappaFile(const std::string& input)
{
	std::cout << "Reading diffusion coefficient data from file: " << input << std::endl;
	std::map<double, double> kappaData = FileHelper::readFile<double, double>(input, 1, 2);

	TimestampValueMap kappaLookup;

	for (const auto& x : kappaData)
	{
		double t = x.first;
		double kappa = x.second * cm2/s/GV;
		kappaLookup[t] = kappa;
	}

	return kappaLookup;
}

/**
 * Read time-dependent data on arbitrary parameters from a file.
 *
 * Format of lookup file:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * timestamp  parname1  parname2  parname3
 * 1306800000 value value value
 * 1307145600 value value value
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * and so on. The word \c timestamp has to be present in the first column on the first line,
 * while the parameter names can be chose at will.
 *
 * @param input File name of input file.
 *
 * @return Map of parameter names to (time-dependent) parameter lookup.
 */
TimedependentParameters FileHelper::readParameterLookupFile(const std::string& input)
{
	std::cout << "Reading parameter lookup file: " << input << std::endl;
	TimedependentParameters lookupMap;

	std::string line;
	std::ifstream file(input);

	if (file.is_open())
	{
		// get parameter names from first line
		std::getline(file, line);
		auto parNames = split(line, " ");
		if (parNames.at(0) != "timestamp")
			throw std::runtime_error("Cannot parse parameter file. Expected token \"timestamp\", but received \""s + parNames.at(0) + "\""s);
		for (size_t i = 1; i < parNames.size(); ++i)
			lookupMap[parNames[i]] = TimestampValueMap();

		// parse remaining lines
		int lineCounter = 1;
		while (std::getline(file, line))
		{
			++lineCounter;
			if (line.starts_with('#'))
				continue;
			auto tokens = split(line, " ");
			if (tokens.empty())
				continue;
			if (tokens.size() != lookupMap.size() + 1)
			{
				std::stringstream message;
				message << "Invalid line " << lineCounter << " in parameter file: Found " << tokens.size() << " tokens, expected "
					<< lookupMap.size() + 1 << "!" << std::endl;
				throw std::runtime_error(message.str());
			}

			double timestamp = std::stod(tokens[0]);
			for (size_t ipar = 1; ipar < parNames.size(); ++ipar)
			{
				TimestampValueMap& lookup = lookupMap[parNames[ipar]];
				lookup[timestamp] = std::stod(tokens[ipar]);
			}
		}
		file.close();
	}
	else
		throw std::runtime_error("Cannot open file: "s + input);

	// dump
	for (const auto& elem : lookupMap)
	{
		std::cout << elem.first << ": " << std::endl;
		const auto& lookup = elem.second;
		size_t counter = 0;
		for (const auto& p : lookup)
		{
			if (counter < 3 || counter > lookup.size() - 4)
				std::cout << (long long int)(p.first) << ": " << p.second << std::endl;
			if (counter == 3)
				std::cout << " ... " << std::endl;
			++counter;
		}
	}

	return lookupMap;
}
