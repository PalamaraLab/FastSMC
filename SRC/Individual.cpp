/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */

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
