#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <vector>
#include <string>

namespace StringUtils {

const std::string RANGE_DELIMS = "{:}"; // must have 3 chars

int stoi(const std::string &s);
double stod(const std::string &s);
float stof(const std::string &s);
std::string itos(int i);
std::string findDelimiters(const std::string &s, const std::string &c);

// will not return blanks
std::vector <std::string> tokenizeMultipleDelimiters(const std::string &s, const std::string &c);

// basic range template: expand "{start:end}" to vector <string> with one entry per range element
// if end==start-1, will return empty
std::vector <std::string> expandRangeTemplate(const std::string &str);
std::vector <std::string> expandRangeTemplates(const std::vector <std::string> &rangeTemplates);
}

#endif
