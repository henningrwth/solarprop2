#ifndef STRINGTOOLS_H_
#define STRINGTOOLS_H_

/// \file stringtools.h
/// \brief Define a helper function for tokenizing a string.

#include <string>
#include <vector>

std::vector<std::string> split(const std::string& str, const std::string& delimiters);

#endif
