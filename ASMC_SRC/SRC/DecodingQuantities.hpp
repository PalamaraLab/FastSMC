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


#ifndef DECODINGQUANTITIES_HPP
#define DECODINGQUANTITIES_HPP

#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

enum class DataType {
    TransitionType, States, CSFSSamples, TimeVector, SizeVector, Discretization, ExpectedTimes, CSFS, FoldedCSFS, ClassicEmission, AscertainedCSFS, FoldedAscertainedCSFS, CompressedAscertainedEmission, initialStateProb, ColumnRatios, RowRatios, Uvectors, Bvectors, Dvectors, HomozygousEmissions, None
};

class DecodingQuantities {

public:

  unsigned int states;
  int CSFSSamples;
  vector <float> initialStateProb;
  vector <float> expectedTimes;
  vector <float> discretization;
  vector <float> timeVector;
  vector <float> columnRatios;
  vector < vector <float> > classicEmissionTable;
  vector < vector <float> > compressedEmissionTable;
  unordered_map<float, vector<float> > Dvectors;
  unordered_map<float, vector<float> > Bvectors;
  unordered_map<float, vector<float> > Uvectors;
  unordered_map<float, vector<float> > rowRatioVectors;
  unordered_map<int, vector<float> > homozygousEmissionMap;
  vector < vector < vector <float> > > CSFSmap;
  vector < vector < vector <float> > > foldedCSFSmap;
  vector < vector < vector <float> > > ascertainedCSFSmap;
  vector < vector < vector <float> > > foldedAscertainedCSFSmap;

  DecodingQuantities(const char *fileName);

private:
  // implemented, but need to update other code
  // void createFromBinary(const char *fileName);
  void createFromGzippedText(const char *fileName);
};

#endif
