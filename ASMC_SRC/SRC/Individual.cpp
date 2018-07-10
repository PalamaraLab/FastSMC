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


#include "Individual.hpp"

Individual::Individual(string _famId, string _IId, int numOfSites) {
	famId = _famId;
	IId = _IId;
	name = famId + "\t" + IId;
	genotype1 = vector <bool> (numOfSites);
	genotype2 = vector <bool> (numOfSites);
}

void Individual::setGenotype(int hap, int pos, bool val) {
	if (hap == 1) genotype1[pos] = val;
	else genotype2[pos] = val;
}
