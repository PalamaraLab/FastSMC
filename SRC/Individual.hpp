#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */

#include <string>
#include <vector>

using namespace std;

class Individual {

  /* **************************** */
  /* **************************** */
  // contains individual names and data
  /* **************************** */
  /* **************************** */
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
