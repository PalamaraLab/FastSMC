#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "StringUtils.hpp"
#include "FileUtils.hpp"
#include "Types.hpp"

#include "DecodingQuantities.hpp"

using namespace std;

DecodingQuantities::DecodingQuantities(const char *fileName) {
  cout << "Using text decoding quantities from " << fileName << endl;
  createFromGzippedText(fileName);
}

void DecodingQuantities::createFromGzippedText(const char *fileName) {
  FileUtils::AutoGzIfstream br; br.openOrExit(fileName);
  string line;
  DataType currentType = DataType::None;
  bool parsedStates = false;
  bool parsedCSFSSamples = false;
  while (getline(br, line)) {
    vector <string> splitString; istringstream iss(line); string buf; while (iss >> buf) splitString.push_back(buf);
    if (splitString.size() == 0 || (splitString.size() == 1 && splitString[0] == "")) {
      continue;
    }
    string firstToken = splitString[0]; boost::algorithm::to_lower(firstToken);
    if (firstToken == string("states")) {
      currentType = DataType::States;
      getline(br, line);
      states = StringUtils::stoi(line);
      parsedStates = true;
      cout << "Parsed states " << states << endl;
      continue;
    } else if (firstToken == string("transitiontype")) {
      currentType = DataType::TransitionType;
      getline(br, line);
      cout << "Ignored TransitionType\n";
      continue;
    } else if (firstToken == string("csfssamples")) {
      currentType = DataType::CSFSSamples;
      parsedCSFSSamples = true;
      getline(br, line);
      CSFSSamples = StringUtils::stoi(line);
      cout << "Parsed CSFSSamples " << CSFSSamples << endl;
      CSFSmap = vector < vector < vector <float> > > (CSFSSamples - 1);
      foldedCSFSmap = vector < vector < vector <float> > > (CSFSSamples - 1);
      ascertainedCSFSmap = vector < vector < vector <float> > > (CSFSSamples - 1);
      foldedAscertainedCSFSmap = vector < vector < vector <float> > > (CSFSSamples - 1);
      continue;
    } else if (firstToken == string("timevector")) {
      currentType = DataType::TimeVector;
      getline(br, line);
      vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
      timeVector = vector <float> (inSplitString.size());
      for (unsigned int i = 0; i < inSplitString.size(); i++) {
        // cout << inSplitString[i] << endl;
        timeVector[i] = StringUtils::stof(inSplitString[i]);
      }
      cout << "Parsed TimeVector " << timeVector.size() << endl;
      continue;
    } else if (firstToken == string("sizevector")) {
      currentType = DataType::SizeVector;
      getline(br, line);
      cout << "Ignored SizeVector\n";
      continue;
    } else if (firstToken == string("expectedtimes")) {
      currentType = DataType::ExpectedTimes;
      if (!parsedStates) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states." << endl;
        exit(1);
      }
      getline(br, line);
      vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
      expectedTimes = vector <float> (inSplitString.size());
      if (inSplitString.size() != states) {
        cerr << "ERROR. Parsed " << inSplitString.size() + " ExpectedTimes entries for " << states << " states." << endl;
        exit(1);
      }
      for (unsigned int i = 0; i < inSplitString.size(); i++) {
        // cout << inSplitString[i] << endl;
        expectedTimes[i] = StringUtils::stof(inSplitString[i]);
      }
      cout << "Parsed ExpectedTimes " << expectedTimes.size() << endl;
      continue;
    } else if (firstToken == string("discretization")) {
      currentType = DataType::Discretization;
      if (!parsedStates) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states." << endl;
        exit(1);
      }
      getline(br, line);
      discretization = vector <float> (states + 1);
      vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
      if (inSplitString.size() != states + 1) {
        cerr << "ERROR. Parsed " << inSplitString.size() + " Discretization entries for " << states << " states." << endl;
        exit(1);
      }
      for (unsigned int i = 0; i < inSplitString.size(); i++) {
        discretization[i] = StringUtils::stof(inSplitString[i]);
      }
      cout << "Parsed Discretization " << discretization.size() << endl;
      continue;
    } else if (firstToken == string("classicemission")) {
      currentType = DataType::ClassicEmission;
      if (!parsedStates) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states." << endl;
        exit(1);
      }
      classicEmissionTable = vector < vector <float> > (2);
      for (int k = 0; k < 2; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " ClassicEmission entries for " << states << " states." << endl;
          exit(1);
        }
        vector <float> thisEmission = vector <float> (states);
        for (unsigned int i = 0; i < states; i++) {
          thisEmission[i] = StringUtils::stof(inSplitString[i]);
        }
        classicEmissionTable[k] = thisEmission;
      }
      cout << "Parsed ClassicEmission " << classicEmissionTable.size() << endl;
      continue;
    } else if (firstToken == string("compressedascertainedemission")) {
      currentType = DataType::CompressedAscertainedEmission;
      if (!parsedStates) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states." << endl;
        exit(1);
      }
      compressedEmissionTable = vector < vector <float> > (2);
      for (int k = 0; k < 2; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " CompressedAscertainedEmission entries for " << states << " states." << endl;
          exit(1);
        }
        vector <float> thisEmission = vector <float> (states);
        for (unsigned int i = 0; i < states; i++) {
          thisEmission[i] = StringUtils::stof(inSplitString[i]);
        }
        compressedEmissionTable[k] = thisEmission;
      }
      cout << "Parsed CompressedAscertainedEmission " << compressedEmissionTable.size() << endl;
      continue;
    } else if (firstToken == string("csfs")) {
      currentType = DataType::CSFS;
      if (!parsedStates || !parsedCSFSSamples) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states and number of CSFS samples." << endl;
        exit(1);
      }
      int undistinguishedIndex = StringUtils::stoi(splitString[1]);
      vector < vector <float> > thisCSFS = vector < vector <float> > ();
      for (int k = 0; k < 3; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " CSFS entries for " << states << " states: " << line << endl;
          exit(1);
        }
        vector <float> CSFSline = vector <float> (states);
        thisCSFS.push_back(CSFSline);
        for (unsigned int i = 0; i < states; i++) {
          thisCSFS[k][i] = StringUtils::stof(inSplitString[i]);
        }
      }
      CSFSmap[undistinguishedIndex] = thisCSFS;
      // cout << "Parsed CSFS " << undistinguishedIndex << " " << CSFSmap[undistinguishedIndex].size() << " " << CSFSmap[undistinguishedIndex][0].size() << endl;
      continue;
    } else if (firstToken == string("foldedcsfs")) {
      currentType = DataType::FoldedCSFS;
      if (!parsedStates || !parsedCSFSSamples) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states and number of CSFS samples." << endl;
        exit(1);
      }
      int undistinguishedIndex = StringUtils::stoi(splitString[1]);
      vector < vector <float> > thisCSFS = vector < vector <float> > ();
      for (int k = 0; k < 2; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " FoldedCSFS entries for " << states << " states: " << line << endl;
          exit(1);
        }
        vector <float> CSFSline = vector <float> (states);
        thisCSFS.push_back(CSFSline);
        for (unsigned int i = 0; i < states; i++) {
          thisCSFS[k][i] = StringUtils::stof(inSplitString[i]);
        }
      }
      foldedCSFSmap[undistinguishedIndex] = thisCSFS;
      // cout << "Parsed FoldedCSFS " << undistinguishedIndex << " " << foldedCSFSmap[undistinguishedIndex].size() << " " << foldedCSFSmap[undistinguishedIndex][0].size() << endl;
      continue;
    } else if (firstToken == string("ascertainedcsfs")) {
      currentType = DataType::AscertainedCSFS;
      if (!parsedStates || !parsedCSFSSamples) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states and number of CSFS samples." << endl;
        exit(1);
      }
      int undistinguishedIndex = StringUtils::stoi(splitString[1]);
      vector < vector <float> > thisCSFS = vector < vector <float> > ();
      for (int k = 0; k < 3; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " AscertainedCSFS entries for " << states << " states: " << line << endl;
          exit(1);
        }
        vector <float> CSFSline = vector <float> (states);
        thisCSFS.push_back(CSFSline);
        for (unsigned int i = 0; i < states; i++) {
          thisCSFS[k][i] = StringUtils::stof(inSplitString[i]);
        }
      }
      ascertainedCSFSmap[undistinguishedIndex] = thisCSFS;
      // cout << "Parsed AscertainedCSFS " << undistinguishedIndex << " " << ascertainedCSFSmap[undistinguishedIndex].size() << " " << ascertainedCSFSmap[undistinguishedIndex][0].size() << endl;
      continue;
    } else if (firstToken == string("foldedascertainedcsfs")) {
      currentType = DataType::FoldedAscertainedCSFS;
      if (!parsedStates || !parsedCSFSSamples) {
        cerr << "ERROR. Parsed " << int(currentType) << " before parsing states and number of CSFS samples." << endl;
        exit(1);
      }
      int undistinguishedIndex = StringUtils::stoi(splitString[1]);
      vector < vector <float> > thisCSFS = vector < vector <float> > ();
      for (int k = 0; k < 2; k++) {
        getline(br, line);
        vector <string> inSplitString; istringstream iss(line); string buf; while (iss >> buf) inSplitString.push_back(buf);
        if (inSplitString.size() != states) {
          cerr << "ERROR. Parsed " << inSplitString.size() + " FoldedAscertainedCSFS entries for " << states << " states: " << line << endl;
          exit(1);
        }
        vector <float> CSFSline = vector <float> (states);
        thisCSFS.push_back(CSFSline);
        for (unsigned int i = 0; i < states; i++) {
          thisCSFS[k][i] = StringUtils::stof(inSplitString[i]);
        }
      }
      foldedAscertainedCSFSmap[undistinguishedIndex] = thisCSFS;
      // cout << "Parsed FoldedAscertainedCSFS " << undistinguishedIndex << " " << foldedAscertainedCSFSmap[undistinguishedIndex].size() << " " << foldedAscertainedCSFSmap[undistinguishedIndex][0].size() << endl;
      continue;
    } else if (firstToken == string("homozygousemissions")) {
      currentType = DataType::HomozygousEmissions;
    } else if (firstToken == string("initialstateprob")) {
      currentType = DataType::initialStateProb;
    } else if (firstToken == string("columnratios")) {
      currentType = DataType::ColumnRatios;
    } else if (firstToken == string("rowratios")) {
      currentType = DataType::RowRatios;
    } else if (firstToken == string("uvectors")) {
      currentType = DataType::Uvectors;
    } else if (firstToken == string("bvectors")) {
      currentType = DataType::Bvectors;
    } else if (firstToken == string("dvectors")) {
      currentType = DataType::Dvectors;
    } else {
      // reading a content line
      if (currentType == DataType::ColumnRatios) {
        columnRatios = vector <float> (states, 0.0);
        for (unsigned int i = 0; i < splitString.size(); i++) {
          columnRatios[i] = StringUtils::stof(splitString[i]);
        }
        // cout << "Parsed ColumnRatios" << endl;
      } else if (currentType == DataType::initialStateProb) {
        initialStateProb = vector <float> (states, 0.0);
        for (unsigned int i = 0; i < splitString.size(); i++) {
          initialStateProb[i] = StringUtils::stof(splitString[i]);
        }
        // cout << "Parsed initialStateProb" << endl;
      } else if (currentType == DataType::RowRatios) {
        vector <float> thisRow(states, 0.0);
        for (unsigned int i = 1; i < splitString.size(); i++) {
          thisRow[i-1] = StringUtils::stof(splitString[i]);
        }
        float index = StringUtils::stof(splitString[0]);
        rowRatioVectors[index] = thisRow;
        // cout << "Parsed RowRatios" << endl;
      } else if (currentType == DataType::Uvectors) {
        vector <float> thisRow(states, 0.0);
        for (unsigned int i = 1; i < splitString.size(); i++) {
          thisRow[i-1] = StringUtils::stof(splitString[i]);
        }
        float index = StringUtils::stof(splitString[0]);
        Uvectors[index] = thisRow;
        // cout << "Parsed Uvectors" << endl;
      } else if (currentType == DataType::Bvectors) {
        vector <float> thisRow(states, 0.0);
        for (unsigned int i = 1; i < splitString.size(); i++) {
          thisRow[i-1] = StringUtils::stof(splitString[i]);
        }
        float index = StringUtils::stof(splitString[0]);
        Bvectors[index] = thisRow;
        // cout << "Parsed Bvectors" << endl;
      } else if (currentType == DataType::Dvectors) {
        vector <float> thisRow(states, 0.0);
        for (unsigned int i = 1; i < splitString.size(); i++) {
          thisRow[i-1] = StringUtils::stof(splitString[i]);
        }
        float index = StringUtils::stof(splitString[0]);
        Dvectors[index] = thisRow;
        // cout << "Parsed Dvectors" << endl;
      } else if (currentType == DataType::HomozygousEmissions) {
        vector <float> thisRow(states, 0.0);
        for (unsigned int i = 1; i < splitString.size(); i++) {
          thisRow[i-1] = StringUtils::stof(splitString[i]);
        }
        int index = StringUtils::stoi(splitString[0]);
        homozygousEmissionMap[index] = thisRow;
        // cout << "Parsed HomozygousEmissions" << endl;
      }
    }
  }
}
