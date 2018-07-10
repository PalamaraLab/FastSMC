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


#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <string>
#include <vector>

using namespace std;

class Individual {

public:
  string famId;
  string IId;
  string name;
  vector <bool> genotype1;
  vector <bool> genotype2;

public:
  Individual(string _famId = "", string _IId = "", int numOfSites = 0);
  void setGenotype(int hap, int pos, bool val);

};

#endif
