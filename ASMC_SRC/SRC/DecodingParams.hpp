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

enum class DecodingMode {
    sequenceFolded, arrayFolded, sequence, array
};

// individual ids and XOR of genotypes
struct DecodingReturnValues {
    vector < vector <float> > sumOverPairs; // output for sum over all pairs
    vector < vector <float> > sumOverPairs00; // output for sum over all pairs with genotype 00
    vector < vector <float> > sumOverPairs01; // output for sum over all pairs with genotype 01 or 10
    vector < vector <float> > sumOverPairs11; // output for sum over all pairs with genotype 11
};


class DecodingParams {

public:

    string hapsFileRoot;
    string decodingQuantFile;
    string outFileRoot;
    int jobs;
    int jobInd;
    string decodingModeString;
    DecodingMode decodingMode;
    bool decodingSequence = false;
    bool foldData = false;
    bool usingCSFS = false;
    bool compress = false;
    bool useAncestral = false;
    float skipCSFSdistance;
    bool noBatches = false;

    // main tasks
    bool doPosteriorSums = false;
    bool doPerPairMAP = false; // output MAP for each pair
    bool doPerPairPosteriorMean = false; // output posterior mean for each pair
    string expectedCoalTimesFile; // expected coalescence times within each interval
    bool withinOnly = false; // only compute decoding within individuals
    bool doMajorMinorPosteriorSums = false;

    bool processCommandLineArgs(int argc, char *argv[]);

};

#endif
