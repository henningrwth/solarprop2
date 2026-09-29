#include "stringtools.h"

/**
 * Tokenize a string.
 *
 * The input string will be split at each character contained in the \p delimiters string.
 *
 * @param str Input string.
 * @param delimiters A string of delimiters.
 *
 * @return Vector of tokens.
 */
std::vector<std::string> split(const std::string& str, const std::string& delimiters)
{
	std::vector<std::string> tokens;

	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}
