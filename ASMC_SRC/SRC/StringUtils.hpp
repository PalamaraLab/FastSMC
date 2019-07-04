//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.


#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <vector>
#include <string>

namespace StringUtils {

const std::string RANGE_DELIMS = "{:}";

double stod(const std::string &s);
float stof(const std::string &s);
std::string itos(int i);
std::string findDelimiters(const std::string &s, const std::string &c);

std::vector <std::string> tokenizeMultipleDelimiters(const std::string &s, const std::string &c);
std::vector <std::string> expandRangeTemplate(const std::string &str);
std::vector <std::string> expandRangeTemplates(const std::vector <std::string> &rangeTemplates);
}

#endif
