  #ifndef DECODINGQUANTITIES_HPP
#define DECODINGQUANTITIES_HPP

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */

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
  vector < vector < vector <float> > > CSFSmap; // index is number of undistinguished
  vector < vector < vector <float> > > foldedCSFSmap; // index is number of undistinguished -- folded
  vector < vector < vector <float> > > ascertainedCSFSmap; // index is number of undistinguished
  vector < vector < vector <float> > > foldedAscertainedCSFSmap; // index is number of undistinguished -- folded

  DecodingQuantities(const char *fileName);

private:
  // TODO: restore this
  // void createFromBinary(const char *fileName);
  void createFromGzippedText(const char *fileName);
};

#endif
