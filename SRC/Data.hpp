#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <string>

#include "Individual.hpp"

using namespace std;

class Data {

  /* **************************** */
  /* **************************** */
  // loads haps/samples
  /* **************************** */
  /* **************************** */

public:
  /* **************************** */
  // structures for sample names and individuals
  vector <string> FamIDList;
  vector <string> IIDList;
  vector <string> famAndIndNameList;
  vector <Individual> individuals;

  // some variables about the data
  int sampleSize;
  int haploidSampleSize;
  int sites;
  int totalSamplesBound;
  bool decodingUsesCSFS = false;
  vector<float> geneticPositions;
  vector<int> physicalPositions;
  vector<bool> siteWasFlippedDuringFolding;
  vector<float> recRateAtMarker;
  vector< vector<int> > undistinguishedCounts;

  Data(string hapsFileRoot, int numOfSites, int totalSamplesBound, bool foldToMinorAlleles, bool decodingUsesCSFS);
  static int countHapLines(string hapsFileRoot);

private:
  void readSamplesList(string hapsFileRoot);
  void readHaps(string hapsFileRoot, bool foldToMinorAlleles);
  int readMap(string hapsFileRoot);
  void makeUndistinguished(bool foldToMinorAlleles);
  vector<int> totalSamplesCount;
  vector<int> derivedAlleleCounts;
  vector<string> SNP_IDs;
  int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);

};

#endif
