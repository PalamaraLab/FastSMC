#include <iostream>
#include <string>
#include <zlib.h>
#include <vector>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

  //make sure parameters are ok
  if ( argc != 2 ) {
    cout << "Number of parameters is wrong." << endl;
    cout << "Only one parameter (name of binary file) is required." << endl;
    exit(1);
  }

  // FILE
  gzFile gzinBIBD;
  gzinBIBD = gzopen( argv[1] , "rb" );

  bool outputIbdSegmentLength;
  gzread( gzinBIBD , (char*) &outputIbdSegmentLength, sizeof(bool) );
  bool doPerPairPosteriorMean;
  gzread( gzinBIBD , (char*) &doPerPairPosteriorMean, sizeof(bool) );
  bool doPerPairMAP;
  gzread( gzinBIBD , (char*) &doPerPairMAP, sizeof(bool) );
  int chrNumber;
  gzread( gzinBIBD , (char*) &chrNumber , sizeof(int) );
  unsigned int nbInd;
  unsigned int lengthFamid;
  unsigned int lengthIid;
  gzread( gzinBIBD , (char*) &nbInd , sizeof(unsigned int) );
  vector <string> idsFam;
  vector <string> idsId;
  for (unsigned int i = 0; i < nbInd; i++) {
    char famID[50];
    char iID[50];
    gzread( gzinBIBD , (char*) &lengthFamid , sizeof(unsigned int) );
    gzread( gzinBIBD , &famID[0] , lengthFamid );
    gzread( gzinBIBD , (char*) &lengthIid , sizeof(unsigned int) );
    gzread( gzinBIBD , &iID[0] , lengthIid );
    std::string idsFam_elmt(famID, lengthFamid);
    std::string idsId_elmt(iID, lengthIid);
    idsFam.push_back(idsFam_elmt);
    idsId.push_back(idsId_elmt);
  }

  int ind[2];
  int hap[2];
  unsigned int pos[2];
  float score;
  float length_cM;

  while ( gzread( gzinBIBD , (char*) &ind[0] , sizeof(int) ) != 0 ) {
    gzread( gzinBIBD , &hap[0] , sizeof( char ) );
    gzread( gzinBIBD , (char*) &ind[1] , sizeof(int) );
    gzread( gzinBIBD , &hap[1] , sizeof( char ) );
    gzread( gzinBIBD , (char*) &pos[0] , sizeof( unsigned int ) );
    gzread( gzinBIBD , (char*) &pos[1] , sizeof( unsigned int ) );
    cout << idsFam[ind[0]] << "\t" << idsId[ind[0]] << "\t" << hap[0] << "\t";
    cout << idsFam[ind[1]] << "\t" << idsId[ind[1]] << "\t" << hap[1] << "\t";
    cout << chrNumber << "\t" << pos[0] << "\t" << pos[1] << "\t";

    if (outputIbdSegmentLength) {
      gzread( gzinBIBD , (char*) &length_cM , sizeof( float ) );
      cout << length_cM << "\t";
    }

    gzread( gzinBIBD , (char*) &score , sizeof( float ) );
    cout << score;

    if (doPerPairPosteriorMean) {
      float postMean;
      gzread( gzinBIBD , (char*) &postMean , sizeof( float ) );
      cout << "\t" << postMean;
    }

    if (doPerPairMAP) {
      float map;
      gzread( gzinBIBD , (char*) &map , sizeof( float ) );
      cout << "\t" << map;
    }

    cout << "\n";
  }

  gzclose(gzinBIBD);

  return 0;
}
