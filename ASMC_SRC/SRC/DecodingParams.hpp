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

#ifndef DECODINGPARAMS_HPP
#define DECODINGPARAMS_HPP

#include <string>
#include <vector>

using namespace std;

enum class DecodingMode { sequenceFolded, arrayFolded, sequence, array };

enum class DecodingModeOverall {
    sequence, array
};

class DecodingParams {

public:


  public:
  string hapsFileRoot;
  string decodingQuantFile;
  string outFileRoot;
  int jobs;
  int jobInd;
  string decodingModeString;
  DecodingModeOverall decodingModeOverall;
  DecodingMode decodingMode;
  bool decodingSequence = false;
  bool foldData;
  bool usingCSFS;
  bool compress;
  bool useAncestral;
  float skipCSFSdistance;
  bool noBatches;

  // main tasks
  bool doPosteriorSums;
  bool doPerPairMAP;           // output MAP for each pair
  bool doPerPairPosteriorMean; // output posterior mean for each pair
  string expectedCoalTimesFile; // expected coalescence times within each interval
  bool withinOnly;      // only compute decoding within individuals
  bool doMajorMinorPosteriorSums = false;

  bool processOptions();
  bool processCommandLineArgs(int argc, char* argv[]);

  /// constructor with default parameters set
  DecodingParams();
      
};

#endif
